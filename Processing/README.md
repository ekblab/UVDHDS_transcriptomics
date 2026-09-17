# Processing Folder for UVDHDS_transcriptomics

This folder contains all pipelines, command-line scripts, and configuration files used for processing raw sequencing data into count matrices for the UVDHDS_transcriptomics project.  
It is subdivided into workflows for mRNA/lncRNA (using Kangaroo-optimized scripts) and miRNA (using nf-core/smRNAseq), fully documenting every reproducible processing step from raw FASTQ files to quantification output.

---

## Folder Structure

```text
Processing/
├── mRNA_Kangaroo/
│   ├── 01_trimming_cutadapt.txt
│   ├── 02_alignment_star.txt
│   ├── 03_sorting_samtools.txt
│   ├── 04_qc_multiqc.txt
│   └── 05_quantification_featureCounts.txt
│
└── miRNA_smRNAseq/
    ├── run_smrnaseq.sh
    ├── make_sample_sheet.R
    ├── samplesheet.csv
    ├── local_docker.config
    └── Logs/ # pipeline execution logs, software versions, and run information
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
  Aggregation of quality control metrics with `MultiQC` (using exported `LC_ALL` and `LANG` variables for robustness).

- **05_quantification_featureCounts.txt**  
  Gene-level quantification using `featureCounts` in three modes:
    - Unique-mapping reads
    - Multimapping reads (all counts)
    - Multimapping reads (fractionally assigned)

Each script is a plain text file containing bash commands or pipelines requiring variable substitution where indicated (`${...}`). See inline comments and variable names for customization to your system/environment.

---

## miRNA_smRNAseq: Small RNA Pipeline

The **miRNA_smRNAseq** subfolder contains the workflow used to process small RNA-seq data generated with the **NEXTFLEX Small RNA-Seq Kit v4**.

### Files

- **run_smrnaseq.sh**  
  Bash script used to execute `nf-core/smrnaseq` with Nextflow and Docker. The script contains the project-specific input paths and the trimming parameters required for NEXTFLEX v4 libraries.

- **make_sample_sheet.R** and **samplesheet.csv**  
  R script to generate a samplesheet compatible with `nf-core/smrnaseq`, together with an example samplesheet used as pipeline input. Paths specified in the samplesheet should point to the corresponding FASTQ files.
  
- **local_docker.config**  
Local Nextflow configuration used to limit pipeline resource usage to a maximum of 6 CPU cores, 48 GB of memory, and 72 hours runtime per process when executing the workflow with Docker.

- **Logs/**  
  Subfolder storing Nextflow logs, software versions, pipeline reports, and other run information for reproducibility.

### NEXTFLEX v4-specific preprocessing

The libraries in this project were generated using the **NEXTFLEX Small RNA-Seq Kit v4**. These libraries require removal of the 3' adapter but **do not require the additional 4-nt clipping at the 5' and 3' ends that was used for earlier NEXTFLEX library chemistries**.

For this reason, the built-in `nextflex` profile (v.2.4.1) of `nf-core/smrnaseq` is **not used** for this project. Instead, the NEXTFLEX v4 preprocessing parameters are set explicitly:

```text
--three_prime_adapter TGGAATTCTCGGGTGCCAAGG
--clip_r1 0
--three_prime_clip_r1 0
--fastp_min_length 16
```

The relevant pipeline configuration is therefore equivalent to:

```bash
-profile docker --three_prime_adapter TGGAATTCTCGGGTGCCAAGG --clip_r1 0 --three_prime_clip_r1 0 --fastp_min_length 16
```

This ensures that only the NEXTFLEX v4 3' adapter is removed and that authentic small-RNA insert sequence is not additionally clipped.

### Pipeline execution

The analysis was performed using:

- **nf-core/smrnaseq:** v2.4.1
- **Nextflow:** v25.10.6
- **Container execution:** Docker
- **Reference genome:** GRCh38
- **miRTrace species:** `hsa`
- **miRDeep:** skipped

The exact command used for execution is documented in `run_smrnaseq.sh`.

---

## Usage Notes

- *mRNA/lncRNA pipeline:*  
  Execute each step in order (trimming → alignment → sorting → QC → quantification). Edit input/output paths and global variables specific to your data. NOTE: We used the Kangaroo workflow and uploaded our FASTQ files to their server, where all of these steps were performed in an automated manner. These commands were provided by Lexogen on request for reproducibility.

- *miRNA processing (`nf-core/smrnaseq`):*  
  - Customize `run_smrnaseq.sh` by updating filesystem paths if necessary.
  - Use `make_sample_sheet.R` to generate an up-to-date `samplesheet.csv` from the directory containing the FASTQ files.
  - For NEXTFLEX v4 libraries, use the explicit adapter and clipping parameters documented above rather than the built-in `nextflex` profile.
  - Check the `Logs` directory for records of pipeline runs, reports, and software versions.

---

## Software Requirements

- See the individual scripts and the manuscript for tool versions and detailed parameterization.
- For miRNA processing, install [Nextflow](https://www.nextflow.io/), [nf-core/smrnaseq](https://nf-co.re/smrnaseq), and Docker or another supported container runtime.
- The miRNA workflow documented here was validated with `nf-core/smrnaseq` v2.4.1 and Nextflow v25.10.6.

---

## Contact

For questions about running or modifying these pipelines, please contact:  
marc.bender@elbekliniken.de

---
