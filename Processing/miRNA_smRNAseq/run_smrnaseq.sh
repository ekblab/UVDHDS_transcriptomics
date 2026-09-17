#!/bin/bash

# Script to execute the nf-core/smrnaseq pipeline for the UV-DHDS project
# using NEXTFLEX Small RNA-Seq Kit v4-compatible trimming.
#
# IMPORTANT:
# Do NOT use the nf-core "nextflex" profile for NEXTFLEX v4 libraries.
# The nf-core profile applies the older NEXTFLEX v3 4N + 4N trimming strategy.
#
# NEXTFLEX v4 requires:
#   - removal of the 3' adapter
#   - no additional 5' clipping
#   - no additional 3' clipping
#   - minimum read length of 16 nt
#
# ~/RNAseq/smRNAseq/UV_DHDS should be replaced by the respective
# working directory if necessary.

~/bin/nextflow-25.10.6 run ~/nextflow/smrnaseq_v2.4.1 \
  -profile docker \
  -work-dir ~/RNAseq/smRNAseq/UV_DHDS/work_nextflex_v4 \
  --input ~/RNAseq/smRNAseq/UV_DHDS/samplesheet.csv \
  --skip_mirdeep \
  --genome GRCh38 \
  --mirtrace_species hsa \
  --three_prime_adapter TGGAATTCTCGGGTGCCAAGG \
  --clip_r1 0 \
  --three_prime_clip_r1 0 \
  --fastp_min_length 16 \
  --outdir ~/RNAseq/smRNAseq/UV_DHDS/output_nextflex_v4
