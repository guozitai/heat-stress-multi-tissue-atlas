#!/bin/bash
#
# scRNA-seq upstream processing — DNBelab C4 (DNBC4tools v2.1.3)
#
# Runs the standard DNBC4tools RNA pipeline on the six PBMC libraries
# (3 heat-stressed, 3 pair-fed thermoneutral). Each library has paired cDNA
# and Oligo FASTQ files. Output feeds the Seurat integration in
# SciData_Master_RCode_Final.R.
#
# Sequence data: NCBI GEO GSE338587
#
# Set the paths below to your local copies before running.
# 运行前请把下面的路径改成本地实际路径。

# 定义输入和输出的基本路径
# Base input/output paths
BASE_DIR="${BASE_DIR:-/path/to/fastq}"            # contains cDNA/ and Oligo/ subdirectories
GENOME_DIR="${GENOME_DIR:-/path/to/genomes/ref}"  # DNBC4tools reference index (ARS-UCD1.2)
OUTPUT_DIR="${OUTPUT_DIR:-./analysis_results}"
DNBC4TOOLS="${DNBC4TOOLS:-dnbc4tools}"            # path to the dnbc4tools executable
THREADS="${THREADS:-8}"

# 定义样本名称列表（根据你的实际情况调整）
# Sample list (3 PFTN + 3 HS)
SAMPLES=("PFTN0305" "PFTN0205" "PFTN0350" "HS1698" "HS1730" "HS1839")

# 循环处理每个样本
# Process each sample in turn
for SAMPLE in "${SAMPLES[@]}"; do
  echo "Processing sample: $SAMPLE"

  # 定义文件路径
  # File paths for this sample
  cDNA_R1="${BASE_DIR}/cDNA/${SAMPLE}/${SAMPLE}_S1_L001_R1_001.fastq.gz"
  cDNA_R2="${BASE_DIR}/cDNA/${SAMPLE}/${SAMPLE}_S1_L001_R2_001.fastq.gz"
  Oligo_R1="${BASE_DIR}/Oligo/${SAMPLE}/${SAMPLE}_S1_L001_R1_001.fastq.gz"
  Oligo_R2="${BASE_DIR}/Oligo/${SAMPLE}/${SAMPLE}_S1_L001_R2_001.fastq.gz"

  # 输出目录
  # Per-sample output directory
  SAMPLE_OUTPUT_DIR="${OUTPUT_DIR}/${SAMPLE}"

  # 运行 dnbc4tools 分析
  # Run the DNBC4tools RNA pipeline
  "${DNBC4TOOLS}" rna run \
    --name "${SAMPLE}" \
    --cDNAfastq1 "${cDNA_R1}" \
    --cDNAfastq2 "${cDNA_R2}" \
    --oligofastq1 "${Oligo_R1}" \
    --oligofastq2 "${Oligo_R2}" \
    --genomeDir "${GENOME_DIR}" \
    --outdir "${SAMPLE_OUTPUT_DIR}" \
    --threads "${THREADS}"

  echo "Sample ${SAMPLE} completed!"
done
