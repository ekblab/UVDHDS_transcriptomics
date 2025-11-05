## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<HEAD>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
##*********************************************************************************************************

## clear workspace
rm(list = ls())
gc()

## load packages
library(tidyverse)
library(DESeq2)
library(ekbSeq)
library(ashr)
library(ComplexHeatmap)
library(EnhancedVolcano)
library(circlize)
library(RColorBrewer)
library(ggprism)
library(ggvenn)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

## load functions for the UV-DHDS project
source("R/UV-DHDS_functions.R")

##*********************************************************************************************************
## Differential expression analysis for differences in cell lines

## load data
processed_data <- readRDS("Data/processed_data_mRNA.rds")
dds <- processed_data$dds
se  <- processed_data$se
anno <- processed_data$anno
vsd_adjusted <- processed_data$vsd_adj
pThres <- dds@thresholds$pvalue
lfcThres <- dds@thresholds$lfc
SkinMarkerGenes <- processed_data$SkinMarkerGenes
rm(processed_data)

## extract results
res <- results(dds, lfcThreshold = lfcThres, alpha = pThres, 
               contrast = contraster(dds,
                                     group1 = list(c("cell", "Melanocytes")),
                                     group2 = list(c("cell", "DSC"))))

##*********************************************************************************************************
## Ma plots

## shrink LFCs and annotate results object
resShrunk <- lfcShrink(dds, res = res, contrast = c("cell", "Melanocytes", "DSC"), type="ashr")
allRes <- cbind(ENSEMBL = rownames(resShrunk), resShrunk)
allRes <- if (dim(anno)[1] == dim(allRes)[1]) {
  left_join(as.data.frame(allRes), anno)
} else {
  stop("Dimensions of annotation object and result object are different.")
}
allRes <- allRes[!str_detect(allRes$ENSEMBL, "\\."), ]
sigRes <- allRes[!is.na(allRes$padj) & allRes$padj < pThres,]

## save results
openxlsx::write.xlsx(allRes, "Results/mRNA/Melanocytes_vs_DSCs_shrunken_LFC.xlsx")

# plot with original lfc
export_plot_dual("Results/mRNA/MA_plot_original", quote({
  plotMA(res, ylim = c(-6, 6), alpha = pThres, colSig = "red2", colNonSig = "grey30", main = "MA plot with original log fold changes")
}), width = 9, height = 6, dpi = 600)

# shrink estimates and exlpore differences in MA plots
export_plot_dual("Results/mRNA/MA_plot_shrunken", quote({
  plotMA(resShrunk, ylim = c(-6, 6), alpha = pThres, colSig = "red2", colNonSig = "grey30", main = "MA plot with shrunken log fold changes")
}), width = 9, height = 6, dpi = 600)

##*********************************************************************************************************
## Transcript distribution

## Type of transcripts
p_AnnoDBI <- plot_transcript_dist(data = sigRes, anno.col = "GENETYPE_AnnoDBI")
p_bm      <- plot_transcript_dist(data = sigRes, anno.col = "GENETYPE_biomaRt")

## save plot
export_plot_dual("Results/mRNA/transcript_type", ggpubr::ggarrange(p_AnnoDBI, p_bm, nrow = 2), width = 10, height = 15)

##*********************************************************************************************************
## Volcano plots

##******
##* Volcano plot with labelling of most strongly DE genes

## get results without preliminary filter
res_volcano <- results(dds)
 
# adding ENSEMBL gene ID as a column in significant differentially expression gene table.
dat_volcano <- cbind(ENSEMBL = rownames(res_volcano), res_volcano)
dat_volcano <- left_join(as.data.frame(dat_volcano), anno)
res_volcano <- NULL

## parameters for Volcano plot
pointSize <- 2
labSize   <- 5

# plot Volcano
# this only serves as diagnostic and representation tool. the actual fold change thresholds and pvalues change later due to the way the results function works
p <- plot_volcano(dat_volcano, title = 'Melanocytes vs. DSCs', show.genes = c("GDF15", "DCT", "CDH1", "CDH2", "NES", "SOX2", "PAX3", "RAB38"), pThres = pThres, labSize = labSize, pointsize = pointSize, raster = TRUE)
export_plot_dual("Results/mRNA/volcano_topSignif", p, width = 10, height = 8)

##******
##* Volcano plot with labelling of skin marker genes

## put marker genes at the end of the dataframe for better plotting
dat_volcano_aux <- dat_volcano[!dat_volcano$SYMBOL %in% SkinMarkerGenes, ]
dat_volcano_markers <- rbind(dat_volcano_aux, dat_volcano[dat_volcano$SYMBOL %in% SkinMarkerGenes, ])

## plot and savevolcano
p <- plot_volcano(dat_volcano_markers, title = 'Melanocytes vs. DSCs', show.genes = SkinMarkerGenes, pThres = pThres, labSize = labSize, pointsize = pointSize)
export_plot_dual("Results/mRNA/volcano_markers", p, width = 10, height = 8)

##*********************************************************************************************************
## Venn diagram of exclusively expressed genes in DSCs and Melanocytes

## plo and save venn diagram
p_venn <-  extract_and_plot_venn(dds, results_object = allRes, pThres = pThres)
export_plot_dual("Results/mRNA/venn_diagram_dsc_vs_melanocytes", p_venn, width = 9, height = 6)

##*********************************************************************************************************
## Heatmap
Ht <- plot_biotype_heatmap(dds, vsd_adjusted, results_object = sigRes)
export_plot_dual("Results/mRNA/Heatmap", quote({
  draw(Ht, merge_legend = TRUE)
}), width = 12, height = 15)

##*********************************************************************************************************
## GO analysis

##******
## store upregulated genes in a list
genes_upregulated <- allRes[allRes$padj < pThres & allRes$log2FoldChange > 0 & !is.na(allRes$padj) & !is.na(allRes$ENTREZID), ]$ENTREZID
## Do GO analysis
go_upregulated <- enrichGO(gene = genes_upregulated, OrgDb = org.Hs.eg.db,  ont = "BP",  keyType = "ENTREZID",  pAdjustMethod = "BH", readable = TRUE)
go_upregulated_reduced <- reduce_go_rrvgo(go_upregulated)

## plot and save bubble plots
p <- bubble_plot_clusterprofiler_style(go_upregulated_reduced, name_col = "term", score_col = "GeneRatio", pval_col = "p.adjust", genes_col = "geneID", 
                                  name_label = "enriched GO terms in melanocytes", rotate_x = TRUE)
export_plot_dual("Results/mRNA/GOBP/upregulated_GOBP", p, width = 5.6, height = 4)

## save results lists
write.csv(go_upregulated, "Results/mRNA/GOBP/upregulated_GOBP_all_terms.csv")
write.csv(go_upregulated_reduced, "Results/mRNA/GOBP/upregulated_GOBP_reduced_terms.csv")

##******
## store downregulated genes in a list
genes_downregulated <- allRes[allRes$padj < pThres & allRes$log2FoldChange < 0 & !is.na(allRes$padj) & !is.na(allRes$ENTREZID), ]$ENTREZID
## plot enriched pathways
go_downregulated <- enrichGO(gene = genes_downregulated, OrgDb = org.Hs.eg.db,  ont = "BP",  keyType = "ENTREZID",  pAdjustMethod = "BH", readable = TRUE)
go_downregulated_reduced <- reduce_go_rrvgo(go_downregulated)

## plot and save bubble plots
p <- bubble_plot_clusterprofiler_style(go_downregulated_reduced, name_col = "term", score_col = "GeneRatio", pval_col = "p.adjust", genes_col = "geneID", 
                                       name_label = "enriched GO terms in DSCs", rotate_x = TRUE)
export_plot_dual("Results/mRNA/GOBP/downregulated_GOBP", p, width = 6.2, height = 4)

## save results lists
write.csv(go_downregulated, "Results/mRNA/GOBP/downregulated_GOBP_all_terms.csv")
write.csv(go_downregulated_reduced, "Results/mRNA/GOBP/downregulated_GOBP_reduced_terms.csv")

##*********************************************************************************************************
## lncRNA analysis

##******
##* Volcano plot of lncRNA
dat_volcano_lncRNA <- dat_volcano[dat_volcano$GENETYPE_biomaRt == "lncRNA",]
p <- plot_volcano(results = dat_volcano_lncRNA, title = "Melanocytes vs. DSCs", show.genes = "GAS5",
                  pThres = pThres, pointsize = pointSize, labSize = labSize, raster = TRUE)
export_plot_dual("Results/mRNA/volcano_topSignif_lncRNA", p, width = 10, height = 8)

##******
##* Venn Diagramm of overlapping genes
p_venn <- extract_and_plot_venn(dds, results_object = allRes, pThres = pThres, biotype_filter = "lncRNA")
export_plot_dual("Results/mRNA/venn_diagram_dsc_vs_melanocytes_lncRNA", p_venn, width = 9, height = 6)

##******
##* Heatmap
Ht <- plot_biotype_heatmap(dds, vsd_adjusted, results_object = sigRes, biotype_filter = "lncRNA")
export_plot_dual("Results/mRNA/Heatmap_lncRNA", quote({
  draw(Ht, merge_legend = TRUE)
}), width = 12, height = 15)

##*********************************************************************************************************
## Gene expression of genes involved in the TGF-beta, miRNA, lncRNA network 

p <- plot_control_expression_comparison(vsd_adjusted, allRes, merge.plots = TRUE, genes = c("TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "TGFBR3",
                                                                                       "SMAD2", "SMAD3", "SMAD4", "MITF", "SOX2")) 
export_plot_dual("Results/mRNA/TGF_beta_network_boxplots_genes", p, width = 5.2, height = 2.5)

p <- plot_control_expression_comparison(vsd_adjusted, allRes, merge.plots = TRUE, genes = c("MEG3", "LINC00520")) 
export_plot_dual("Results/mRNA/TGF_beta_network_boxplots_lncRNA", p, width = 1.8, height = 2.6)

##*********************************************************************************************************
## TGFbeta analysis

## KEGG pathway information for TGF beta retrieved from: https://www.gsea-msigdb.org/gsea/msigdb/cards/KEGG_TGF_BETA_SIGNALING_PATHWAY
# Read KEGG pathway file
KEGG_TGF_beta <- read.table(file = 'Data/KEGG_TGF_BETA_SIGNALING_PATHWAY.v2025.1.Hs.tsv', sep = '\t', header = TRUE)
REACTOME_TGF_beta <- read.table(file = 'Data/REACTOME_SIGNALING_BY_TGF_BETA_RECEPTOR_COMPLEX.v2025.1.Hs.tsv', sep = '\t', header = TRUE)

# Extract and split the gene symbols
KEGG_TGF_beta_genes <- KEGG_TGF_beta[KEGG_TGF_beta$STANDARD_NAME == "GENE_SYMBOLS", ]
REACTOME_TGF_beta_genes <- REACTOME_TGF_beta[REACTOME_TGF_beta$STANDARD_NAME == "GENE_SYMBOLS", ]
KEGG_gene_vector <- unlist(strsplit(KEGG_TGF_beta_genes[, 2], split = ","))
REACTOME_gene_vector <- unlist(strsplit(REACTOME_TGF_beta_genes[, 2], split = ","))
gene_vector <- unique(c(KEGG_gene_vector, REACTOME_gene_vector))

## plot KEGG genes
Ht <- plot_biotype_heatmap(dds, vsd_adjusted, results_object = allRes %>% filter(SYMBOL %in% KEGG_gene_vector), 
                           show_row_dend = TRUE, show_row_names = TRUE, split = 3, plot_title = "Expression of genes involved in KEGG TGF-beta signalling")
export_plot_dual("Results/mRNA/Heatmap_TGF_beta_KEGG", quote({
  draw(Ht, merge_legend = TRUE)
}), width = 12, height = 15)

## plot Reactome genes
Ht <- plot_biotype_heatmap(dds, vsd_adjusted, results_object = allRes %>% filter(SYMBOL %in% REACTOME_gene_vector), 
                           show_row_dend = TRUE, show_row_names = TRUE, split = 3, plot_title = "Expression of genes involved in Reactome TGF-beta signalling")
export_plot_dual("Results/mRNA/Heatmap_TGF_beta_Reactome", quote({
  draw(Ht, merge_legend = TRUE)
}), width = 12, height = 15)

## plot TGFB1/2/3 and TGFBR1/2/3
p <- plot_control_expression_comparison(vsd_adjusted, allRes, merge.plots = FALSE, genes = c("TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "TGFBR3", "BMP1", "BMP2", "BMP4",  
                                                                        "SMAD2", "SMAD3", "SMAD4", "SMAD1", "SMAD5", "SMAD6"), nrow = 3) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) 
export_plot_dual("Results/mRNA/TGF_beta_boxplots", p, width = 8, height = 4)





