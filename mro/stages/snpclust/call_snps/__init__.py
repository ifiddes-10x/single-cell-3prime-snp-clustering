#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#
import martian
import subprocess
import shutil
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
    out vcf.gz   raw_variants,
    src py       "stages/snpclust/call_snps",
) split using (
    in  bam    input,
    in  string locus,
    in  path   genome_fasta,
    out vcf    raw_variant_chunk,
)
'''
def split(args):
    # bring in genome fasta and index it -- cellranger references have no fasta index
    genome_fasta_path = cr_utils.get_reference_genome_fasta(args.reference_path)
    local_path = martian.make_path('genome.fa')
    try:
        os.symlink(genome_fasta_path, local_path)
    except OSError:
        shutil.copy(genome_fasta_path, local_path)
    subprocess.check_call(['samtools', 'faidx', local_path])
    #with open(local_path.replace('.fa', '.dict'), 'w') as outf:
    #    subprocess.check_call(['samtools', 'dict', local_path], stdout=outf)

    chunks = []
    for bam, locus in zip(args.input_bams, args.loci):
        chunks.append({'locus': locus, 'input_bam': bam, 'genome_fasta': local_path, '__mem_gb': 32, '__threads': 1})
    return {'chunks': chunks, 'join': {'__mem_gb': 64}}


def main(args, outs):
    bed_path = martian.make_path('region.bed')
    with open(bed_path, 'w') as f:
        f.write(args.locus)

    raw_variant_chunk = martian.make_path(outs.raw_variant_chunk)
    with open(raw_variant_chunk, 'w') as outf:
        subprocess.check_call(['freebayes',
                               '-t', bed_path,
                               '-f', args.genome_fasta,
                               '--haplotype-length', '0',
                               '--min-alternate-count', '1',
                               '--min-alternate-fraction', '0',
                               '--pooled-continuous',
                               '--use-best-n-alleles', '2',
                               args.input_bam], stdout=outf)

    # get the sample. Making the assumption that this is a single-sample BAM
    #in_bam = tk_bam.create_bam_infile(args.input_bam)
    #sample = in_bam.header['RG'][0]['SM']

    # TODO: fix hardcoded path
    # had to do this to deal with memory issues
    #subprocess.check_call(['java', '-Xmx{}g'.format(args.__mem_gb),
    #                       '-jar', '/mnt/home/stephen/miniconda2/share/gatk4-4.0.1.2-0/gatk-package-4.0.1.2-local.jar',
    #                       'Mutect2', '-R', args.genome_fasta, '--intervals',
    #                       bed_path, '-I', args.input_bam, '-tumor', sample,
    #                       '-O', raw_variant_chunk, '--TMP_DIR', os.getcwd(),
    #                       '--native-pair-hmm-threads', str(args.__threads)])
    #subprocess.check_call(['gatk-launch', 'Mutect2', '-R', args.genome_fasta, '--intervals',
    #                       bed_path, '-I', args.input_bam, '-tumor', sample,
    #                       '-O', raw_variant_chunk, '--TMP_DIR', os.getcwd(),
    #                       '--native-pair-hmm-threads', str(args.__threads)])

        
def join(args, outs, chunk_defs, chunk_outs):
    outs.coerce_strings()
    outs.output = [chunk.raw_variant_chunk for chunk in chunk_outs]

    raw_variants = martian.make_path('raw_variants.vcf')
    tk_io.combine_vcfs(raw_variants, outs.output)
    outs.raw_variants = raw_variants + '.gz'