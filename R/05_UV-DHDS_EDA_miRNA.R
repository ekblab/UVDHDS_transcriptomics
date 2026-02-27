## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<HEAD>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
##*********************************************************************************************************

## clear workspace
rm(list = ls())
gc()

## load packages
library(tidyverse)
library(DESeq2)
library(ekbSeq) ## install with devtools::install_github("https://github.com/MBender1992/ekbSeq")
library(pheatmap)
library(ggpubr)
library(PoiClaClu)

##*********************************************************************************************************
## Exploratory data analysis

## load data
processed_data <- readRDS("Data/processed_data_miRNA.rds")
dds <- processed_data$dds
se  <- processed_data$se
vsd_raw <- processed_data$vsd_unadjusted
# vsd_adjusted <- processed_data$vsd_adj
rm(processed_data)

##******
## plot correlation matrix

## calculate Poisson distances and define parameters for graphical representation
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix(poisd$dd)
rownames(samplePoisDistMatrix) <- colData(se)$sample
colnames(samplePoisDistMatrix) <- NULL
colorsPoisd <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(255)

## plot heatmap
export_plot_dual("Results/miRNA/heatmap_poisson_distance", quote({
  p <- pheatmap(samplePoisDistMatrix,
           clustering_distance_rows=poisd$dd,
           clustering_distance_cols=poisd$dd,
           col=colorsPoisd)
  print(p)
}), width = 10, height = 8, dpi = 300)

##******
## check for outliers
export_plot_dual("Results/miRNA/cooks_distance", quote({
  par(mar = c(16, 5, 2, 2))
  boxplot(log10(assays(dds)[["cooks"]]), range = 0, las = 2)
  title("Cooks distance across samples")
}), width = 7, height = 5, dpi = 300)

##******
## plot PCA for potential batch effects
p <- pca_plot(data = vsd_raw, pcsToUse = 1:2, intgroup = c("donor", "cell"), 
         title = "PCA with unadjusted VST data", subtitle = "PC1 vs. PC2")
export_plot_dual("Results/miRNA/PCA",  p,  width = 4, height = 3)

##******
##* Plot dispersion
export_plot_dual("Results/miRNA/dispersion", quote({
  plotDispEsts(dds)
}), width = 6, height = 4, dpi = 300)

