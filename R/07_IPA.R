## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<HEAD>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
##*********************************************************************************************************

## clear workspace
rm(list = ls())
gc()

## load packages
library(tidyverse)
library(DESeq2)
library(ekbSeq)
library(readxl)

## load functions for the UV-DHDS project
source("R/UV-DHDS_functions.R")

##*********************************************************************************************************
## Load data
dat_upstream <- read_xlsx("Data/IPA_upstream_regulators.xlsx")
dat_pw <- read_xlsx("Data/IPA_canonical_pathways.xlsx") %>%
  mutate(padj = 10^-`-log(B-H p-value)`) %>%
  select(-`-log(B-H p-value)`)


## function to replace space and dash in column names
replace_dash_space <- function(str){
  str <- str_replace_all(str, " ", "_")
  str <- str_replace_all(str, "-", "_")
  return(str)
}

## replace column names
names(dat_upstream) <- replace_dash_space(names(dat_upstream))
names(dat_pw) <- replace_dash_space(names(dat_pw))

##*********************************************************************************************************
## Pathways

dat_pw_up <- dat_pw %>% filter(z_score > 0 & padj < 0.05)
dat_pw_down <- dat_pw %>% filter(z_score < 0 & padj < 0.05)

## 
p1 <- ipa_bubble_plot(dat_pw_up, name_col = "Ingenuity_Canonical_Pathways", name_label = "upregulated IPA pathways")
p2 <- ipa_bubble_plot(dat_pw_down, name_col = "Ingenuity_Canonical_Pathways", name_label = "downregulated IPA pathways")

export_plot_dual("Results/mRNA/IPA/upregulated_pathways", p1, width = 6.5, height = 4)
export_plot_dual("Results/mRNA/IPA/downregulated_pathways", p2, width = 6.2, height = 4)

##*********************************************************************************************************
## Upstream regulators

dat_upstream_up <- dat_upstream %>% filter(Activation_z_score >= 2 & Expr_Log_Ratio > 1 & p_value_of_overlap < 0.05)
dat_upstream_down <- dat_upstream %>% filter(Activation_z_score <=-2 & Expr_Log_Ratio < -1 & p_value_of_overlap < 0.05)
## 
p1 <- ipa_bubble_plot(dat_upstream_up, name_col = "Upstream_Regulator", score_col = "Activation_z_score", pval_col = "p_value_of_overlap", 
                                  genes_col = "Target_Molecules_in_Dataset", name_label = "upregulated upstream regulators")
p2 <- ipa_bubble_plot(dat_upstream_down, name_col = "Upstream_Regulator", score_col = "Activation_z_score", pval_col = "p_value_of_overlap",
                                  genes_col = "Target_Molecules_in_Dataset", name_label = "downregulated upstream regulators")

export_plot_dual("Results/mRNA/IPA/upregulated_upstream_regulators", p1, width = 6.5, height = 4)
export_plot_dual("Results/mRNA/IPA/downregulated_upstream_regulators", p2, width = 4.1, height = 4)
