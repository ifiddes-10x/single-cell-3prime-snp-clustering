#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#
import martian
import subprocess
import tenkit.bam as tk_bam
import tenkit.bio_io as tk_io
import cellranger.utils as cr_utils
import os

__MRO__ = '''
stage CALL_SNPS(
    in  path     reference_path,
    in  bam[]    input_bams,
    in  string[] loci,
    out vcf[]    output,
    out vcf      raw_variants,
    src py       "stages/snpclust/call_snps",
) split using (
    in  bam    input,
    in  string locus,
    out vcf    raw_variant_chunk,
)
'''
def split(args):
    chunks = []
    for bam, locus in zip(args.input_bams, args.loci):
        chunks.append({'locus': locus, 'input_bam': bam, '__mem_gb': 32, '__threads': 8})
    return {'chunks': chunks, 'join': {'__mem_gb': 16}}


def main(args, outs):
    genome_fasta_path = cr_utils.get_reference_genome_fasta(args.reference_path)

    bed_path = martian.make_path('region.bed')
    with open(bed_path, 'w') as f:
        f.write(args.locus)

    # get the sample. Making the assumption that this is a single-sample BAM
    in_bam = tk_bam.create_bam_infile(args.input_bam)
    sample = in_bam.header['RG'][0]['SM']
    raw_variant_chunk = martian.make_path(outs.raw_variant_chunk)
    subprocess.check_call(['gatk-launch', 'Mutect2', '-R', genome_fasta_path, '--intervals',
                           bed_path, '-I', args.input_bam, '-tumor', sample,
                           '-O', raw_variant_chunk, '--TMP_DIR', os.getcwd(),
                           '--native-pair-hmm-threads', str(args.__threads)])

        
def join(args, outs, chunk_defs, chunk_outs):
    outs.coerce_strings()
    outs.raw_variants = [chunk.raw_variant_chunk for chunk in chunk_outs]
    
    tk_io.combine_vcfs(outs.raw_variants, outs.output)
