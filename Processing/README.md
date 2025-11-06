# Processing Folder for UVDHDS_transcriptomics

This folder contains all pipelines, command-line scripts, and configuration files used for processing raw sequencing data into count matrices for the UVDHDS_transcriptomics project.  
It is subdivided into workflows for mRNA/lncRNA (using Kangaroo-optimized scripts) and miRNA (using nf-core/smRNAseq), fully documenting every reproducible processing step from raw FASTQ files to quantification output.

---

## Folder Structure

```
Processing/
├── mRNA_Kangaroo/
│   ├── 01_trimming_cutadapt.txt
│   ├── 02_alignment_star.txt
│   ├── 03_sorting_samtools.txt
│   ├── 04_qc_multiqc.txt
│   └── 05_quantification_featureCounts.txt
|
└── miRNA_smRNAseq/
    ├── run_smrnaseq.sh
    ├── make_sample_sheet.R
    ├── samplesheet.csv
    └── Logs/ # containing pipeline execution logs and software information

```

---

## mRNA_Kangaroo: Bulk mRNA/lncRNA Processing

The **mRNA_Kangaroo** subfolder contains shell command templates used in sequence to process polyA+ bulk RNA-seq reads:

- **01_trimming_cutadapt.txt**  
  Multi-step adapter and polyA trimming using `cutadapt` to ensure high-quality, adapter-free reads.
  
- **02_alignment_star.txt**  
  Alignment to the GRCh38 genome with `STAR` using parameters optimized for sensitivity and reproducibility.

- **03_sorting_samtools.txt**  
  Sorting of unsorted BAM alignments with `samtools`; original unsorted files are removed to save space.

- **04_qc_multiqc.txt**  
  Aggregation of quality control metrics with `multiQC` (using exported LC_ALL and LANG for robustness).

- **05_quantification_featureCounts.txt**  
  Gene-level quantification using `featureCounts` in three modes:
    - Unique-mapping reads
    - Multimapping reads (all counts)
    - Multimapping reads (fractionally assigned)

Each script is a plain text file containing bash commands or pipelines requiring variable substitution where indicated (`${...}`). See inline comments and variable names for customization to your system/environment.

---

## miRNA_smRNAseq: Small RNA Pipeline

The **miRNA_smRNAseq** subfolder contains:
- **run_smrnaseq.sh**  
  Template bash script to execute the nf-core/smRNAseq pipeline (e.g., Nextflow with Docker), specifying project settings and sample sheet locations.

- **make_sample_sheet.R** and **samplesheet.csv**  
  R script to auto-generate a samplesheet compatible with nf-core/smRNAseq, and an example samplesheet as input.

- **Logs/**  
  Subfolder storing all Nextflow, software, and pipeline output logs for full transparency and reproducibility.

---

## Usage Notes

- *mRNA/lncRNA pipeline*:  
  Execute each step in order (trimming → alignment → sorting → QC → quantification). Edit input/output paths and global variables specific to your data. NOTE: we used the Kangaroo workflow and uploaded our fastq files to their server, where all of these steps were performed in an automatised way. These commands were provided by lexogen on request for reproducibility. 

- *miRNA processing (nf-core/smRNAseq)*:  
  - Customize `run_smrnaseq.sh` by updating the path variables to match your filesystem.
  - Use `make_sample_sheet.R` to generate an up-to-date `samplesheet.csv` from a directory of FASTQ files.
  - Check the `Logs` directory for records of pipeline runs and software versions used.

---

## Software Requirements

- See specific `.txt` files or the manuscript for tool versions and parameterization.
- For miRNA, install [Nextflow](https://www.nextflow.io/), [nf-core/smRNAseq](https://nf-co.re/smrnaseq), and required containers or dependencies as specified in `run_smrnaseq.sh`.

---

## Contact

For questions about running or modifying these pipelines, please contact:  
marc.bender@elbekliniken.de

---