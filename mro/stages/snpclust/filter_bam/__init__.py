#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#
import martian
import subprocess
import collections
import itertools
import tenkit.bam as tk_bam
import cellranger.utils as cr_utils
import snpclust.constants as snp_constants
import fidlib.fidlib.intervals as fl_intervals

__MRO__ = '''
stage FILTER_BAM(
    in  bam       input,
    in  path      reference_path,
    in  path      bed_file,
    in  tsv       cell_barcodes,
    out bam[]     output_bams, 
    out string[]  loci,
    src py     "stages/snpclust/filter_bam",
) split (
    in  string locus,
    out bam    output,
)
'''

def split(args):
    if args.bed_file is not None:
        loci = open(args.bed_file).readlines()
    else:
        # split up the genome, but only into exonic chunks.
        ref_gtf = cr_utils.get_reference_genes_gtf(args.reference_path)
        exons = find_exon_loci(ref_gtf)
        loci = build_loci(exons, snp_constants.REGION_SPLIT_SIZE)
    chunks = [{'locus': locus, '__mem_gb': 16} for locus in loci]
    return {'chunks': chunks}


def main(args, outs):

    in_bam = tk_bam.create_bam_infile(args.input)
    tmp_bam = martian.make_path('tmp.bam')
    out_bam, _ = tk_bam.create_bam_outfile(tmp_bam, None, None, template=in_bam)
    #out_bam, _ = tk_bam.create_bam_outfile(outs.output, None, None, template=in_bam)

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
    # not sure why I have to do this, it should be sorted already. Maybe due to nearby exonic ranges?
    subprocess.check_call(['sambamba', 'sort', '-o', outs.output, tmp_bam])
    tk_bam.index(outs.output)


def join(args, outs, chunk_defs, chunk_outs):
    outs.coerce_strings()
    # pass along every BAM produced
    outs.output_bams = [chunk.output for chunk in chunk_outs]
    outs.loci = [chunk_def.locus for chunk_def in chunk_defs]


def find_exon_loci(ref_gtf):
    """
    Extracts exonic regions of a GTF into merged, sorted intervals
    """
    regions = collections.defaultdict(list)
    for l in open(ref_gtf):
        if '\texon\t' in l:
            l = l.split()
            i = fl_intervals.ChromosomeInterval(l[0], int(l[3]), int(l[4]), '.')
            if len(i) > 0:
                regions[l[0]].append(i)
    merged_regions = {}
    for chrom, region_list in regions.iteritems():
        merged_regions[chrom] = fl_intervals.gap_merge_intervals(region_list, 0)
    # report a flat list
    return [item for sublist in merged_regions.itervalues() for item in sublist]


def build_loci(exons, chunk_size):
    """
    Given a list of intervals, pack bins with as many until chunk_size is exceeded
    :param chunk_size: Number of bases to pack into one bin
    :return: BED-like string
    """
    loci = []
    this_locus = []
    bin_size = 0
    for e in exons:
        this_locus.append(e)
        bin_size += len(e)
        if bin_size >= chunk_size:
            loci.append('\n'.join(['\t'.join(map(str, [e.chromosome, e.start, e.stop])) for e in this_locus]))
            this_locus = []
            bin_size = 0
    return loci