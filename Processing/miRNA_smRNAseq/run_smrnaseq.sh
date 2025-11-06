#!/bin/bash

# script to execute the nf-core smRnaseq pipeline with custom options for the UV-DHDS project
# "~RNAseq/smRNAseq/UV_DHDS should be replaced by your working directory

~/bin/nextflow-24.10.6-dist run ~/nextflow/smrnaseq \
   -profile docker,nextflex \
  --input ~/RNAseq/smRNAseq/UV_DHDS/samplesheet.csv \
  --skip_mirdeep \
  --genome 'GRCh38' \
  --mirtrace_species 'hsa' \
  --outdir ~/RNAseq/smRNAseq/UV_DHDS/output
  
