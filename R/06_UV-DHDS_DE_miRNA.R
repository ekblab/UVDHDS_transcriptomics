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

## load functions for the UV-DHDS project
source("R/UV-DHDS_functions.R")

##*********************************************************************************************************
## Differential expression analysis for differences in cell lines

## load data
processed_data <- readRDS("Data/processed_data_miRNA.rds")
dds <- processed_data$dds
se  <- processed_data$se
vsd_raw <- processed_data$vsd_unadjusted
pThres <- dds@thresholds$pvalue
lfcThres <- dds@thresholds$lfc
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
resShrunk$SYMBOL <- rownames(resShrunk)
resShrunk$ENSEMBL <- resShrunk$SYMBOL
allRes <- data.frame(resShrunk)
sigRes <- allRes[!is.na(allRes$padj) & allRes$padj < pThres,]

## save results
openxlsx::write.xlsx(rownames_to_column(allRes, "miRNA"), "Results/miRNA/Melanocytes_vs_DSCs_shrunken_LFC.xlsx")

# plot with original lfc
svg("Results/miRNA/MA_plot_original.svg", width=9, height=6)
plotMA(res, ylim = c(-6, 6), alpha = pThres, colSig = "red2", colNonSig = "grey30", main = "MA plot with original log fold changes")
dev.off()

# shrink estimates and exlpore differences in MA plots
svg("Results/miRNA/MA_plot_shrunken.svg", width=9, height=6)
plotMA(resShrunk, ylim = c(-6, 6), alpha = pThres, colSig = "red2", colNonSig = "grey30", main = "MA plot with shrunken log fold changes")
dev.off()

##*********************************************************************************************************
## Volcano plots

##******
##* Volcano plot with labelling of most strongly DE genes

## get results without preliminary filter
dat_volcano <- results(dds)
dat_volcano$SYMBOL <- rownames(dat_volcano)

## parameters for Volcano plot
labSize   <- 5

# plot Volcano
dat_volcano <- as.data.frame(dat_volcano) %>%
  mutate(SYMBOL = str_replace_all(.$SYMBOL, "hsa-", ""))
# this only serves as diagnostic and representation tool. the actual fold change thresholds and pvalues change later due to the way the results function works
p <- plot_volcano(dat_volcano, title = "Melanocytes vs. DSCs", xlim = c(-12, 12), labSize = labSize, raster = TRUE)
export_plot_dual("Results/miRNA/volcano_topSignif", p, width = 10, height = 8)

## plot targets of lncRNA
p <- EnhancedVolcano(as.data.frame(dat_volcano), lab = rownames(dat_volcano), selectLab = c("hsa-miR-125b-5p", "hsa-miR-125a-3p", "hsa-miR-125a-5p", "hsa-miR-27b-3p", "hsa-miR-34a-3p", "hsa-miR-34a-5p"),
                     x = 'log2FoldChange', y = 'padj', gridlines.major = FALSE, gridlines.minor = FALSE, raster = TRUE
                     xlim = c(-12, 12), ylim = c(-5, 60), title = 'Melanocytes vs. DSCs', drawConnectors = TRUE, arrowheads = FALSE,
                     pCutoff = pThres, FCcutoff = 1, pointSize = 2, labSize = labSize) 

## save plot
export_plot_dual("Results/miRNA/volcano_lncTargets", p, width = 10, height = 8)


##*********************************************************************************************************
## Heatmap

## specify fontsize
heatmapFontsize <- 18

## specify colors for heatmap
colors  <-  colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
col_fun <- colorRamp2(c(-2, 0, 2), c(c(colors[1], colors[51], colors[100])))

## Heatmap of all DE genes
matAll <- assay(vsd_raw)[rownames(assay(vsd_raw)) %in% rownames(sigRes), ]

## scale matrix
datScaled <- matAll %>%  t() %>% scale() %>%  t()
datScaled <- datScaled[complete.cases(datScaled) & rowSums(is.infinite(datScaled)) == 0, ]

## specify annotation bar
annoHeatmap <- as.data.frame(colData(vsd_raw)[, c("cell", "donor")])
colnames(annoHeatmap) <- c("Cell_type", "Donor")
donors <- colorRampPalette(rev(brewer.pal(n = 7, name = "Set1")))(6)
names(donors) <- unique(colData(vsd_raw)$donor)

## specify annotation bar colors
annColors <- HeatmapAnnotation(df =annoHeatmap,
                               col = list(Cell_type = c("DSC" = ggsci::pal_npg("nrc")(4)[1], "Melanocytes" = ggsci::pal_npg("nrc")(4)[2]),
                                          Donor     = donors),
                               annotation_legend_param = list(Cell_type = list(nrow=1), 
                                                              title_gp = gpar(fontsize = heatmapFontsize), 
                                                              labels_gp = gpar(fontsize = heatmapFontsize)),
                               annotation_name_gp = gpar(fontsize = heatmapFontsize, fontface = "bold"))

## plot Heatmap
Ht <- Heatmap(datScaled,col= col_fun,
              top_annotation = annColors,
              clustering_method_row = "average",
              clustering_method_columns = "average",
              clustering_distance_row = "pearson",
              clustering_distance_column = "euclidean",
              show_row_dend = FALSE,
              show_row_names =  FALSE,
              show_column_names = TRUE,
              column_names_gp = gpar(fontsize = 10),
              column_title = paste("Differentially Expressed genes (", nrow(sigRes), ")", sep =""),
              heatmap_legend_param = list(
                title = "row Z-score",
                at = seq(-2,2,by=1),
                color_bar="continuous",
                title_position ="topcenter",
                legend_direction = "horizontal",
                legend_width = unit(4, "cm")
              ))

export_plot_dual("Results/miRNA/Heatmap", quote({
  draw(Ht, merge_legend = TRUE)
}), width = 12, height = 15)

##*********************************************************************************************************
## Venn diagram of exclusively expressed genes in DSCs and Melanocytes
p_venn <-  extract_and_plot_venn(dds, results_object = allRes, pThres = pThres)
export_plot_dual("Results/miRNA/venn_diagram_dsc_vs_melanocytes", p_venn, width = 9, height = 6)

## plot targets of lncRNA
p <- plot_control_expression_comparison(vsd_raw, allRes, genes = c("hsa-miR-125a-3p","hsa-miR-125b-5p", "hsa-miR-125a-5p", "hsa-miR-27b-3p", "hsa-miR-34a-3p", "hsa-miR-34a-5p"))
export_plot_dual("Results/miRNA/boxplots_lnc_targets", p, width = 6, height = 5)

##*********************************************************************************************************
## Gene expression of genes involved in the TGF-beta, miRNA, lncRNA network 

rownames(vsd_raw) <- str_remove(rownames(vsd_raw), "hsa-")
allRes$SYMBOL <- str_remove(allRes$SYMBOL, "hsa-")
allRes$ENSEMBL <- str_remove(allRes$ENSEMBL, "hsa-")
rownames(allRes) <- str_remove(rownames(allRes), "hsa-")

p <- plot_control_expression_comparison(vsd_raw, allRes, merge.plots = TRUE, genes = c("miR-204-5p", "miR-1268a", "miR-338-3p", "miR-27b-3p", 
                                                                                       "miR-18a-5p",  "miR-31-5p", "miR-146a-5p",  "miR-146b-5p", 
                                                                                       "miR-125a-3p", "miR-125b-5p", "miR-125b-2-3p",  "miR-137",
                                                                                       "miR-210-3p", "miR-340", "miR-381-3p","miR-211-5p", "miR-127-3p")) 
export_plot_dual("Results/miRNA/TGF_beta_network_boxplots_miRNA", p, width = 7, height = 2.8)

