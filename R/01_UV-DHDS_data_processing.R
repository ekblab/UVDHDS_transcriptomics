## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<HEAD>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

##*********************************************************************************************************

## load packages
library(tidyverse)
library(DESeq2)
library(readxl)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ekbSeq) ## install with devtools::install_github("https://github.com/MBender1992/ekbSeq")
library(biomaRt)

## load functions for the UV-DHDS project
source("R/UV-DHDS_functions.R")

##*********************************************************************************************************
## Load data

## load marker genes
markerGenes <- read_csv2("Data/marker_genes.csv")

## extract marker genes from csv
McMarkerGenes    <- markerGenes[markerGenes$TYPE == "Melanocytes",]$SYMBOL
DscMarkerGenes   <- markerGenes[markerGenes$TYPE == "DSCs",]$SYMBOL
FibroMarkerGenes <- markerGenes[markerGenes$TYPE == "Fibroblasts",]$SYMBOL
SkinMarkerGenes  <- c(McMarkerGenes, DscMarkerGenes, FibroMarkerGenes)

## load data in the same order as the metadata rows
counts_20240924 <- read.table(file = 'Data/UV-DHDS_count_matrix_EKB_20240924_lexogen_20241016.tsv', sep = '\t', header = TRUE)[,-1]
counts_20250410 <- read.table(file = 'Data/UV-DHDS_count_matrix_EKB_20250410_lexogen_20250428.tsv', sep = '\t', header = TRUE)[,-1]

## combine data into merged counts
counts_merged <- merge(counts_20240924, counts_20250410, by = "gene_id")

## convert counts into matrix
counts_merged <- column_to_rownames(counts_merged, "gene_id") 

## assign coldata and filter colData to only include controls samples
colData <- read_csv("Data/metadata_mRNA.csv") %>% 
  filter(!run %in% "TUD" & !is.na(ID) & !is.na(run) & run != "library failed") %>% 
  mutate(ID = ifelse(run %in% c("EKB_20250410", "EKB_20250417"), tolower(ID), ID)) %>%  
  filter(irradiation == "control") %>%
  mutate(sample = str_remove(sample, "_\\d+h")) %>%
  dplyr::select(-c(condition, irradiation)) 

## reformat colnames and coldata to share the same naming convention
colnames(counts_merged) <- str_replace_all(colnames(counts_merged), "\\.|-", "_")
colData$ID <- str_replace_all(colData$ID, "\\.|-", "_")

## remove unused counts
counts_filtered <- counts_merged[,colnames(counts_merged) %in% colData$ID]
  
## reorder colData based on numerical values (only the )
colData[colData$run %in% c("EKB_20250410", "EKB_20250417"),] <- colData[colData$run %in% c("EKB_20250410", "EKB_20250417"),][stringr::str_order(colData[colData$run %in% c("EKB_20250410", "EKB_20250417"),]$ID, numeric = T),]

## assign collapsing string for collapseReplicates
colData$collapse <- paste(colData$run, colData$time, sep = "_")

## check whether counts names and colData are ordered in the same way
if(all(toupper(colnames(counts_filtered)) == toupper(colData$ID))) {
  message("Counts matrix and colData object share the same identifiers in the correct order")
} else {
  stop("There are discrepancies between the counts matrix and the colData object.")
}

## construct summarized experiment object
se <- SummarizedExperiment(assays= as.matrix(counts_filtered), colData=colData)
seColl <- collapseReplicates(se, se$sample, se$collapse)

##*********************************************************************************************************
## Set thresholds and variable names and perform DE analysis

## always specify this filter criterion to take the desired difference into account in significance testing
pThres <- 0.01
lfcThres <- round(log2(1), 4)

## define contrast
dds <- DESeqDataSet(seColl, design = ~ run + cell)

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
vsd_raw <- vst(dds, blind = FALSE)

## remove batch effect
vsd_adjusted <- vsd_raw
mat <- assay(vsd_adjusted)
mat <- limma::removeBatchEffect(mat, vsd_adjusted[["run"]])
assay(vsd_adjusted) <- mat
mat <- NULL

## Construct annotation object
anno <- AnnotationDbi::select(org.Hs.eg.db, rownames(dds), 
                              columns=c("ENSEMBL", "ENTREZID", "SYMBOL", "GENENAME", "GENETYPE"), 
                              keytype="ENSEMBL")

## remove duplicated ENSEMBL ID entries
anno <- anno %>% filter(!duplicated(ENSEMBL))
colnames(anno)[colnames(anno) == "GENETYPE"] <- "GENETYPE_AnnoDBI"

## annotate gene IDs
mart <- useMart('ENSEMBL_MART_ENSEMBL')
mart <- useDataset('hsapiens_gene_ensembl', mart)
anno_bm <- getBM(mart = mart, attributes = c('ensembl_gene_id',  'hgnc_symbol', 'entrezgene_id', 'gene_biotype'), uniqueRows = TRUE) %>% 
  dplyr::select(ensembl_gene_id, gene_biotype) %>%
  dplyr::filter(!duplicated(ensembl_gene_id)) %>%
  setNames(c("ENSEMBL", "GENETYPE_biomaRt")) 

## combine annotations from annoDBI and biomaRt. Add RP4 as a name as it was detected in IPA workflow
anno <- anno %>% left_join(anno_bm)
anno[anno$ENSEMBL == "ENSG00000277287",]$SYMBOL <- "RP4_794I64"

## store processed objects as list
processed_data <- list(dds = dds, 
                       seColl = seColl, 
                       vsd_unadjusted = vsd_raw, 
                       vsd_adj = vsd_adjusted, 
                       anno = anno,
                       SkinMarkerGenes = SkinMarkerGenes)

## processsed data
saveRDS(processed_data, "Data/processed_data_mRNA.rds")
