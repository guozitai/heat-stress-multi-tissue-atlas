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
- scRNA-seq: 6 PBMC samples (82,013 cells after QC)
- Arteriovenous blood biochemistry

## Repository contents

| File | Covers |
|---|---|
| `SciData_Master_RCode_Final.R` | Main analysis script: count matrix → DESeq2 differential expression (Fig. 2), GO/KEGG enrichment (Fig. 3), miRNA differential expression and target networks (Fig. 4), scRNA-seq integration and SingleR annotation (Fig. 5), hdWGCNA co-expression modules (Fig. 6a–c), CellChat cell–cell communication (Fig. 6d) |
| `SciData_Figure_RCode_final.R` | Final figure-rendering script used to produce the published panels. Its plotting sections overlap with the master script; where the two differ, this file reflects the version used for the submitted figures. |
| `run_dnbc4tools_all.sh` | scRNA-seq upstream processing (DNBC4tools v2.1.3) for the six PBMC libraries |
| `sessionInfo.txt` | Full R session and package versions for the analyses above |

### Not included

Alignment and read counting for the bulk mRNA-seq libraries (FASTQ → count matrix) were carried out by the sequencing provider and are not reproduced here. The resulting count matrix is the starting point of `SciData_Master_RCode_Final.R` and is archived at figshare (below), so the analysis can be reproduced from that point onward.

## Running the code

Input files are expected under a `data/` directory, or wherever `HS_DATA_DIR` points:

```bash
export HS_DATA_DIR=/path/to/data      # default: ./data
export HS_OUT_DIR=/path/to/figures    # default: ./figures
Rscript SciData_Master_RCode_Final.R
```

The input files are archived in the figshare record listed under Data availability.

## Software

- R 4.4.2 on Ubuntu 22.04
- Key packages: Seurat 5.1.0, SingleR 2.8.0, SCINA 1.2.0, celldex 1.16.0, hdWGCNA 0.4.11, CellChat 2.2.0, clusterProfiler 4.14.3, DESeq2 1.46.0
- Upstream: DNBC4tools 2.1.3 (scRNA-seq), HISAT2 / StringTie (bulk mRNA-seq), miRDeep2 (miRNA-seq)
- Full details in `sessionInfo.txt`

## Data availability

- **Sequence data**: NCBI Gene Expression Omnibus, accession [GSE338587](https://identifiers.org/geo:GSE338587) (samples GSM9877815–GSM9877864)
- **Processed data and supplementary tables**: figshare, <https://doi.org/10.6084/m9.figshare.32993768>

## Contact

Dengpan Bu - budengpan@caas.cn
