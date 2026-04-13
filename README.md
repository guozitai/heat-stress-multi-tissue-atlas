# Multi-tissue transcriptome atlas of heat-stressed dairy cows

Code repository for the manuscript:

**Multi-tissue mRNA-seq, miRNA-seq and single-cell RNA-seq data from heat-stressed lactating Holstein cows**

Zitai Guo, Shengtao Gao, Lu Ma, Dengpan Bu*

Submitted to *Nature Scientific Data* (2026)

## Overview

This repository contains bioinformatics analysis code for a multi-tissue, multi-omics transcriptome dataset from lactating Holstein cows under heat stress (HS) versus pair-fed thermoneutral (PFTN) conditions.

## Dataset

- Bulk mRNA-seq: 37 samples across 5 tissues (rumen, liver, mammary, adipose, muscle)
- miRNA-seq: 16 samples (adipose, muscle)
- scRNA-seq: 6 PBMC samples (69,078 cells after QC)
- Arteriovenous blood biochemistry

## Software

- R (v4.x)
- Key packages: DESeq2, HISAT2, StringTie, Seurat (v5.1.0), hdWGCNA (v0.4.11), CellChat (v2.2.0), clusterProfiler, miRDeep2, SingleR, SCINA

## Data availability

Raw sequencing data will be deposited in NCBI SRA. Processed data will be available via figshare. Accession numbers will be updated upon publication.

## Contact

Dengpan Bu - budengpan@caas.cn
