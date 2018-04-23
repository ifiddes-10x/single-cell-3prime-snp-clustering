#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#
import martian
import subprocess
import shutil
import itertools
import tenkit.bam as tk_bam
import cellranger.utils as cr_utils
import snpclust.constants as snp_constants
import os

__MRO__ = '''
stage FILTER_MODIFY_BAM(
    in  bam       input,
    in  path      reference_path,
    in  path      bed_file,
    in  tsv       cell_barcodes,
    out bam[]     output_bams, 
    out string[]  loci,
    src py     "stages/snpclust/filter_modify_bam",
) split (
    in  string locus,
    in path    genome_fasta,
    out bam    output,
)
'''

def split(args):
    # bring in genome fasta and index it -- cellranger references have no fasta index or dict file
    genome_fasta_path = cr_utils.get_reference_genome_fasta(args.reference_path)
    local_path = martian.make_path('genome.fa')
    try:
        os.symlink(genome_fasta_path, local_path)
    except OSError:
        shutil.copy(genome_fasta_path, local_path)
    subprocess.check_call(['samtools', 'faidx', local_path])
    with open(local_path.replace('.fa', '.dict'), 'w') as outf:
        subprocess.check_call(['samtools', 'dict', local_path], stdout=outf)

    if args.bed_file is not None:
        loci = open(args.bed_file).readlines()
    else:
        in_bam = tk_bam.create_bam_infile(args.input)
        loci = build_loci(in_bam, snp_constants.REGION_SPLIT_SIZE)
    chunks = [{'locus': locus, 'genome_fasta': local_path, '__mem_gb': 16} for locus in loci]
    return {'chunks': chunks}


def main(args, outs):

    in_bam = tk_bam.create_bam_infile(args.input)
    tmp_bam = martian.make_path('tmp.bam')
    out_bam, _ = tk_bam.create_bam_outfile(tmp_bam, None, None, template=in_bam)

    cell_bcs = set(cr_utils.load_barcode_tsv(args.cell_barcodes))
    loci = [x.split() for x in args.locus.split('\n')]
    for chrom, start, stop in loci:
        bam_iter = in_bam.fetch(chrom, int(start), int(stop), multiple_iterators=True)
        for (tid, pos), reads_iter in itertools.groupby(bam_iter, key=cr_utils.pos_sort_key):
            dupe_keys = set()
            for read in reads_iter:
                if cr_utils.get_read_barcode(read) not in cell_bcs:
                    continue

                if cr_utils.is_read_dupe_candidate(read, cr_utils.get_high_conf_mapq({'high_conf_mapq': 255})):
                    dupe_key = (cr_utils.si_pcr_dupe_func(read), cr_utils.get_read_umi(read))
                    if dupe_key in dupe_keys:
                        continue

                    dupe_keys.add(dupe_key)
                    read.is_duplicate = False
                    out_bam.write(read)

    out_bam.close()

    # Correct the STAR mapping from 255 to 60 and take care of split reads
    output_bam = martian.make_path(outs.output)
    star_args = ['gatk-launch', 'SplitNCigarReads',
                 '-R', args.genome_fasta,
                 '-I', tmp_bam,
                 '-O', output_bam,
                 '--skip-mapping-quality-transform', 'false',
                 '--create-output-bam-index', 'false',
                 '--TMP_DIR', os.getcwd()]

    subprocess.check_call(star_args)
    os.remove(tmp_bam)
    tk_bam.index(output_bam)


def join(args, outs, chunk_defs, chunk_outs):
    outs.coerce_strings()
    # pass along every BAM produced
    outs.output_bams = [chunk.output for chunk in chunk_outs]
    outs.loci = [chunk_def.locus for chunk_def in chunk_defs]


def build_loci(in_bam, chunk_size):
    """
    Given a BAM path, builds a set of loci
    :param in_bam: Open BAM handle
    :param chunk_size: Number of bases to pack into one bin
    :return: BED-like string
    """
    # greedy implementation of the bin packing problem
    loci = []
    this_locus = []
    bin_size = 0
    for chrom, chrom_length in zip(in_bam.header.references, in_bam.header.lengths):
        region_start = 0
        while region_start < chrom_length:
            start = region_start
            end = min(region_start + chunk_size, chrom_length)
            this_locus.append([chrom, start, end])
            bin_size += end - start
            if bin_size >= chunk_size:
                loci.append('\n'.join(['\t'.join(map(str, x)) for x in this_locus]))
                this_locus = []
                bin_size = 0
            region_start = end
    return loci