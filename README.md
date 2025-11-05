# UVDHDS_transcriptomics
Part of the **UVDHDS** project investigating the Effects of **UV** exposition on **D**ifferentiation of **H**uman **D**ermal **S**tem cells in melanomagenesis.

---

## Overview
This repository contains the analysis pipeline, processed data, and figure scripts for the first publication associated with the *UVDHDS* project. 
The study describes for the first time transcriptomic (mRNA, lncRNA, miRNA) differences between dermal stem cells (DSCs) and their melanocyte progeny. This transcriptomic profiling is the foundation for further investigations focussing differentiation process itself as well as UV induced changes. 

This repository specifically covers **bulk RNA-seq transcriptomic profiling** of donor-matched DSCs and melanocytes to establish a transcriptomic "baseline".

---

## Project context
The **UVDHDS** project (*"Effects of UV exposition on differentiation of human dermal stem cells in melanomagenesis"*) is part of the molecular cell biology program at the  
**Laboratory for Molecular Cell Biology (LMZ), Elbe Klinikum Buxtehude, Germany.**

### Planned publications
- **UVDHDS_transcriptomics:** Bulk RNA-seq profiling of dermal stem cells compared to their melanocyte progeny
- **UVDHDS_irradiation:** Cellular and transcriptional dynamics after (repeated) UV irradiation  
- **UVDHDS_scRNAseq:** Single-cell RNA-seq tracing of DSC differentiation trajectories  

---

## Repository Structure
```
UVDHDS_transcriptomics/
│
├── Data/ # Processed data (counts, metadata, etc.)
├── R/ # R analysis scripts
├── Results/ 
└── README.md # You are here
```

The data folder contains counts of 2 runs which were merged as technical replicates. The metadata object contains combined data from the UVDHDS_transcriptomics and UVDHDS_irradiation study as some of the datasets were sequenced together to save resources and increase cost-efficiency. Pathway information for KEGG/Reactome pathways were retrieved from: https://www.gsea-msigdb.org/gsea/msigdb. Results from ingenuity pathway analysis (IPA, Qiagen) were exported as .xls(x). The R folder contains analysis scripts for data processing, exploratory data analysis (EDA) and differential expression (DE) analysis as well as visualization of IPA results with clusterProfiler-style plots. Due to restriction in file size the results folder only includes dummy folders to mimic the structure of the analysis for easier reproducibility. On request original figures will be provided.  

## Main analyses included
- Quality control and normalization of bulk RNA-seq data  
- Differential expression analysis (DESeq2)  
- Gene Ontology and pathway enrichment (clusterProfiler)  
- Transcription factor and upstream regulator inference (IPA)  

---

## Reproducibility
All analyses were performed using R (≥4.5.1).  
Dependencies are listed in `sessionInfo.txt`.  
Most packages are available on CRAN or bioconductor. Custom functions and wrappers for the analysis are either available in R/functions.R or in the custom package ekbSeq which can be installed via *devtools::install_github("mbender1992/ekbSeq")*

---

## Contact
**ekblab** – Laboratory for Molecular Cell Biology.  
Elbe Klinikum Buxtehude, Germany.  
marc.bender@elbekliniken.de 

---

## LICENSE
This project is distributed under the MIT License (see LICENSE).



