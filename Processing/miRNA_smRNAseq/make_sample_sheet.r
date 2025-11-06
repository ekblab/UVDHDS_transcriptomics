# script to prepare sample sheet based on fastq files in the fastq folder

library(tidyverse)

## set working directory
setwd("~/RNAseq/smRNAseq/UV_DHDS")

## list all fastq files in the input folder
fastqs <- list.files("~/RNAseq/smRNAseq/UV_DHDS/fastq")
path <- "/home/zellbio/RNAseq/smRNAseq/UV_DHDS/fastq/"

samplesheet <- data.frame(
  sample = "",
  fastq_1 = paste0(path, fastqs[seq(1, length(fastqs), 2)])
) %>%
  mutate(sample = str_extract(fastq_1, "[dm]\\d+"))

## order sample sheet
samplesheet <- samplesheet[stringr::str_order(samplesheet$sample, numeric = T),]

## output
write.csv(samplesheet, "samplesheet.csv")

