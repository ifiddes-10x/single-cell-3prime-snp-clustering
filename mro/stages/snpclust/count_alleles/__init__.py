#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#
import collections
import itertools
import pyfasta
import shutil
import martian
import numpy as np
import vcf
import tenkit.bam as tk_bam
import tenkit.bio_io as tk_io
import cellranger.matrix as cr_matrix
import cellranger.utils as cr_utils

__MRO__ = '''
stage COUNT_ALLELES(
    in  path      reference_path,
    in  bam[]     input_bams,
    in  vcf[]     variants,
    in  tsv       cell_barcodes,
    in  int       min_bcs_per_snp,
    in  int       min_snp_obs,
    in  int       min_snp_base_qual,
    in  float     base_error_rate,
    out h5        raw_allele_bc_matrices_h5,
    out path      raw_allele_bc_matrices_mex,
    src py        "stages/snpclust/count_alleles",
) split using (
    in  bam       chunk_bam,
    in  vcf       chunk_variants,
    in  fasta     fasta,
    in  int       num_alleles,
    out  h5       allele_matrix,
)
'''

def split(args):

    # prepare a local pyfasta-capable reference
    genome_fasta_path = cr_utils.get_reference_genome_fasta(args.reference_path)
    local_fasta = martian.make_path('genome.fa')
    shutil.copy(genome_fasta_path, local_fasta)
    _ = pyfasta.Fasta(local_fasta)  # will flatten the file

    # hacky way to find the largest number of alleles we will expect to see
    max_alt_alleles = 0
    for v in args.vcf:
        for rec in vcf.Reader(open(v)):
            max_alt_alleles = max(max_alt_alleles, len(rec.ALT))
    num_alleles = max_alt_alleles + 1

    chunks = []
    for chunk_bam, chunk_variants in zip(args.input_bams, args.variants):
        chunks.append({'chunk_variants': chunk_variants, 'chunk_bam': chunk_bam, 'fasta': local_fasta,
                       'num_alleles': num_alleles, '__mem_gb': 16})

    return {'chunks': chunks, 'join': {'__mem_gb': 128}}


def main(args, outs):
    in_bam = tk_bam.create_bam_infile(args.chunk_bam)
    fasta = pyfasta.Fasta(args.fasta)
    # hack to strip fasta comments
    fasta = {name.split()[0]: rec for name, rec in fasta.iteritems()}

    # barcodes
    bcs = cr_utils.load_barcode_tsv(args.cell_barcodes)
    # efficient barcode index
    bc_index = {x: i for i, x in enumerate(bcs)}

    # load all records so we know the matrix shape
    records = list(vcf.Reader(open(args.chunk_variants)))

    # hacky way to create fake 'genes'
    class Pos(object):
        def __init__(self, rec):
            self.id = '{}_{}'.format(rec.CHROM, rec.POS)

    positions = [Pos(rec) for rec in records]

    # initialize matrix -- replacing 'genomes' with num_alleles and 'genes' with snps
    mat = cr_matrix.GeneBCMatrices(range(4), itertools.repeat(positions, 4), bcs)

    # begin iterating over records and assigning allele calls
    for i, record in enumerate(records):
        alleles = tk_io.get_record_alt_alleles(record)
        ref = tk_io.get_record_ref(record)
        all_alleles = [ref] + alleles
        chrom = tk_io.get_record_chrom(record)
        pos = tk_io.get_record_pos(record)
        ref = tk_io.get_record_ref(record)

        # iterating over single column means this loop only actually happens once
        # this is required or pysam gets very upset when you try to look at the records
        for col in in_bam.pileup(chrom, pos - 1, pos, truncate=True):
            # keep track of mapping of CB -> allele calls in order to populate matrix
            cb_map = collections.defaultdict(lambda: np.zeros(args.num_alleles))

            # strip N-cigar operations (introns) to avoid needless computation
            # also check for having a cell barcode
            reads = [x.alignment for x in col.pileups if not x.is_refskip and x.alignment.has_tag('CB')
                     and x.alignment.get_tag('CB') in bc_index]

            for read in reads:
                allele = tk_bam.read_contains_allele_sw(ref, all_alleles, pos, read, fasta[chrom],
                                                        match=6, mismatch=-4, gap_open=-6, gap_extend=-1)
                if allele != -1:
                    cb = read.get_tag('CB')
                    cb_map[cb][allele] += 1

            # transform CB map into matrix values -- each cell contains the number of reads of that allele seen
            for cb, cb_m in cb_map.iteritems():
                j = bc_index[cb]
                for k, val in enumerate(cb_m):
                    mat.get_matrix(k).m[i, j] = val

    mat.save_h5(outs.allele_matrix)


def join(args, outs, chunk_defs, chunk_outs):
    outs.coerce_strings()

    raw_chunk_h5s = [chunk_out.raw_allele_bc_matrices_h5 for chunk_out in chunk_outs]
    raw_allele_bc_matrices = cr_matrix.merge_matrices(raw_chunk_h5s)

    raw_allele_bc_matrices.save_h5(outs.raw_allele_bc_matrices_h5)
    raw_allele_bc_matrices.save_mex(outs.raw_allele_bc_matrices_mex)
