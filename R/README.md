# R Scripts for UVDHDS_transcriptomics

This folder contains all R scripts used for the analysis and visualization within the *UVDHDS_transcriptomics* project, which investigates transcriptomic differences between human dermal stem cells (DSCs) and their melanocyte progeny using bulk RNA-seq profiling.

---

## Folder Structure & Script Overview

- **01_UV-DHDS_data_processing.R**  
  Importing, pre-processing, and merging raw/processed RNA-seq data (counts, metadata, etc.) for mRNA and lncRNA 

- **02_UV-DHDS_EDA.R**  
  Exploratory Data Analysis (EDA) of mRNA/lncRNA data, including QC, normalization assessment, batch inspection, and data visualization.

- **03_UV-DHDS_DE.R**  
  Visualisation of the differential expression analysis of mRNA/lncRNA using DESeq2, as well as post-DE annotation and table preparation.

- **04_UV-DHDS_processing_miRNA.R**  
  Importing, preprocessing, and initial normalization for miRNA data.

- **05_UV-DHDS_EDA_miRNA.R**  
  Exploratory analysis, normalization checks, and PCA for miRNA samples.

- **06_UV-DHDS_DE_miRNA.R**  
  Differential expression analysis and visualisation of miRNA expression.

- **07_IPA.R**  
  Ingenuity Pathway Analysis (IPA) workflows and figure generation, including import of IPA results and visualization with custom R scripts.

- **UV-DHDS_functions.R**  
  Collection of custom R functions and wrapper utilities used throughout the analysis pipeline. Many analytical steps call reusable functions from this script for enhanced reproducibility and code modularity.

---

## Usage

Each script is annotated and intended to be run sequentially as the analysis progresses, unless otherwise noted.  
Paths, input, and output files are set relative to the repository root or `/Data` subfolder.  
Custom functions are loaded from `UV-DHDS_functions.R`; if additional utility functions are required, they can be sourced at the top of each script.

All external R package dependencies are documented in the root `sessionInfo.txt` file.

---

## Notes

- Scripts were written for R (≥4.5.1).
- Outputs include summary tables, processed count data, differential expression results, and figures.
- IPA results (Qiagen) must be exported from the IPA GUI as `.xls(x)` and are imported for further plotting/interpretation.

---

## Citation

Please cite the main UVDHDS_transcriptomics publication and this repository if you use these scripts or derivatives in your own work.

---

## Contact

For questions, clarifications, or collaboration inquiries, contact:  
marc.bender@elbekliniken.de
