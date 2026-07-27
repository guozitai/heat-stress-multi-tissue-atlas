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
- scRNA-seq: 6 PBMC samples (82,013 cells; no cell-level filtering downstream of cell calling)
- Arteriovenous blood biochemistry

## Repository contents

| File | Covers |
|---|---|
| `SciData_Master_RCode_Final.R` | Main analysis script. Starting from the compiled count matrix and the intermediate result objects of the original analysis, it produces the differential-expression summaries (Fig. 2), GO/KEGG enrichment (Fig. 3), miRNA target networks (Fig. 4), cell-type annotation and cell-type-resolved comparisons (Fig. 5), hdWGCNA co-expression modules (Fig. 6a–c) and CellChat cell–cell communication (Fig. 6d). See **What this code does and does not run** below. |
| `figures/` | Scripts that rendered the **final submitted versions** of Fig. 2c, Fig. 3, Fig. 5a–d and Fig. 6b. See `figures/README.md` — it also documents the text-rendering requirement and the excluded unassigned cluster. |
| `SciData_Figure_RCode_final.R` | A later partial collection covering Fig. 2a, the Fig. 2 and Fig. 4 combined legends, Fig. 6b and Fig. 6c, plus the colour update applied to the Fig. 4b network. It is **not** a standalone script for all figures, and for the panels listed under `figures/` it has been superseded. |
| `run_dnbc4tools_all.sh` | scRNA-seq upstream processing (DNBC4tools v2.1.3) for the six PBMC libraries |
| `sessionInfo.txt` | Full R session and package versions for the analyses above |

Panels not covered by `figures/` come from `SciData_Master_RCode_Final.R` unchanged. Fig. 1 is a schematic with no R code.

### What this code does and does not run

`SciData_Master_RCode_Final.R` reproduces the reported results and figure panels from data
that already exist; it is not a pipeline that runs every analysis from raw reads. Four steps of the original workflow were performed earlier: one was carried out by the
sequencing provider and does not appear in the code, and three are **loaded, not
recomputed**:

| Step | How it appears in the script | Status |
|---|---|---|
| FASTQ → count matrix (bulk mRNA-seq) | not present | performed by the sequencing provider; the resulting count matrix is archived at figshare |
| DESeq2 differential expression | `readRDS("DEG_results_all_tissues.rds")` | result object from the original analysis is read in |
| miRNA target prediction (TargetScan, miRWalk) | `read_excel("targetgene_enrich.xlsx")` | prediction table from the original analysis is read in |
| Single-cell integration (Seurat anchors) | `readRDS("integrated_data_with_singleR.rds")` | integrated object from the original analysis is read in |

Everything downstream of those inputs — the enrichment analyses, the figures, the
cell-type-resolved comparisons, the co-expression modules and the cell–cell communication
analysis — is executed by the deposited code. The single-cell integration parameters are
reported in the Methods of the manuscript.

## Running the code

Input files are expected under a `data/` directory, or wherever `HS_DATA_DIR` points:

```bash
export HS_DATA_DIR=/path/to/data      # default: ./data
export HS_OUT_DIR=/path/to/figures    # default: ./figures
Rscript SciData_Master_RCode_Final.R
```

The bulk count matrix is archived at figshare (see Data availability). The other inputs are
the archived intermediate objects of the original analysis listed in the table above; the
single-cell objects are intermediates of the original analysis rather than deposited files,
and the deposited per-sample single-cell matrices at GEO allow the single-cell workflow to be
repeated with the integration parameters given in the Methods of the manuscript.

## Software

- R 4.4.2 on Ubuntu 22.04
- Key packages: Seurat 5.1.0, SingleR 2.8.0, SCINA 1.2.0, celldex 1.16.0, hdWGCNA 0.4.11, CellChat 2.2.0.9001, clusterProfiler 4.14.3, DESeq2 1.46.0
- Upstream: DNBC4tools 2.1.3 (scRNA-seq), HISAT2 / StringTie (bulk mRNA-seq), miRDeep2 (miRNA-seq)
- Full details in `sessionInfo.txt`

## Data availability

- **Sequence data**: NCBI Gene Expression Omnibus, accession [GSE338587](https://identifiers.org/geo:GSE338587) (samples GSM9877815–GSM9877864)
- **Processed data and supplementary tables**: figshare, <https://doi.org/10.6084/m9.figshare.32993768>

## Contact

Dengpan Bu - budengpan@caas.cn
