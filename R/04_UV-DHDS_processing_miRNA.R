## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<HEAD>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

##*********************************************************************************************************

## load packages
library(tidyverse)
library(DESeq2)
library(readxl)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ekbSeq) ## install with devtools::install_github("https://github.com/MBender1992/ekbSeq")

##*********************************************************************************************************
## Load data

## load functions for the UV-DHDS project
source("R/UV-DHDS_functions.R")

## load data in the same order as the metadata rows
counts <- read.table(file = 'Data/UV-DHDS_mirna_counts_mirtop.tsv', sep = '\t', header = TRUE)
# counts_hairpin_edgeR <- read_edgeR_counts('Data/hairpin_counts_edgeR.csv')

## convert counts into matrix
counts <- column_to_rownames(counts, "miRNA") 

## reorder counts based on numerical order
counts <- counts[,stringr::str_order(colnames(counts), numeric = T)]

## assign coldata
colData <- read_csv("Data/metadata_miRNA.csv") %>% 
  filter(!is.na(run)) %>% 
  filter(irradiation == "control") %>%
  mutate(sample = str_remove(sample, "_\\d+h")) %>%
  dplyr::select(-c(condition, irradiation)) 

## remove unused counts
counts_filtered <- counts[,colnames(counts) %in% tolower(colData$ID)]

## reorder colData based on numerical values (only the )
colData <- colData[stringr::str_order(colData$ID, numeric = T),]

## check whether counts names and colData are ordered in the same way
if(all(toupper(colnames(counts_filtered)) == toupper(colData$ID))) {
  message("Counts matrix and colData object share the same identifiers in the correct order")
} else {
  stop("There are discrepancies between the counts matrix and the colData object.")
}

## construct summarized experiment object
se <- SummarizedExperiment(assays= as.matrix(counts_filtered), colData=colData)
seColl <- collapseReplicates(se, se$sample, se$time)

##*********************************************************************************************************
## Set thresholds and variable names and perform DE analysis

## always specify this filter criterion to take the desired difference into account in significance testing
pThres <- 0.05
lfcThres <- round(log2(1), 4)

## define contrast
dds <- DESeqDataSet(seColl, design = ~ cell)

## run DESeq
dds <- DESeq(dds, parallel = TRUE)

## add metadata to DESeq object with attributes
attr(dds, "thresholds") <- list(
  pvalue = pThres, 
  lfc = lfcThres
)

##*********************************************************************************************************
## Data transformation

## Transform data
vsd_raw <- vst(dds, blind = FALSE, nsub = sum( rowMeans( counts(dds, normalized=TRUE)) > 5 ))

## remove batch effect (not necessary due to low batch)
# vsd_adjusted <- vsd_raw
# mat <- assay(vsd_adjusted)
# mat <- limma::removeBatchEffect(mat, vsd_adjusted[["run"]])
# assay(vsd_adjusted) <- mat
# mat <- NULL

## store processed objects as list
processed_data <- list(dds = dds, 
                       se = seColl, 
                       vsd_unadjusted = vsd_raw
                       )

## processsed data
saveRDS(processed_data, "Data/processed_data_miRNA.rds")
