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

counts <- read.table(file = 'Data/UV-DHDS_mRNA_counts.tsv', sep = '\t', header = TRUE)

## convert counts into matrix
counts <- as.matrix(column_to_rownames(as.data.frame(counts), "gene_id"))

## load metadata
colData <- read_csv("Data/metadata_mRNA_GEO.csv")

## assign collapsing string for collapseReplicates
colData$collapse <- paste(colData$run, colData$replicate, sep = "_")

## check whether counts names and colData are ordered in the same way
if(all(colnames(counts) == colData$ID)) {
  message("Counts matrix and colData object share the same identifiers in the correct order")
} else {
  stop("There are discrepancies between the counts matrix and the colData object.")
}

## construct summarized experiment object
se <- SummarizedExperiment(assays = as.matrix(counts), colData = colData)
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
design_viz <- model.matrix( ~ cell,  data = as.data.frame(colData(vsd_adjusted)))
mat <- limma::removeBatchEffect(mat, vsd_adjusted[["run"]], design = design_viz)
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

# colData_geo <- colData %>% mutate(time = ifelse(time == "6h", "rep1", "rep2"))
# colnames(colData_geo)[5] <- "replicate"
# colData_geo <- colData_geo %>% dplyr::select(-c(`Irradiation date`,cell_passage, entity))
# colData_geo <- colData_geo
# colData_geo <- colData_geo %>% mutate(sample = str_replace_all(.$sample, "_control", "")) %>%
#   mutate(ID_new = paste(.$sample, .$run, sep = "_"))
# colData_geo$ID_old <- colData$ID
# colData_geo$ID <- colData_geo$ID_new
# colData_geo$ID_new <- NULL
# colData_geo$ID <- str_remove(paste(colData_geo$ID, colData_geo$replicate, sep = "_"), "rep")
# 
# colData_geo %>% dplyr::select(-c(paired_end, pf_reads_sample_percent, seed_date, harvest_date, rna_extract_date, collapse)) %>% write.csv("metadata_mRNA_GEO.csv")
# 
# colnames(counts_filtered) <- colData_geo$ID
# 
# counts_filtered %>% rownames_to_column("gene_id") %>% write_tsv("UV-DHDS_counts.tsv")

