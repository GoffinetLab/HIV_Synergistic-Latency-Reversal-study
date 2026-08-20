#!/usr/bin/env bash

set -euo pipefail

########################################
# Script: 00_QC_and_Trimming.sh
# Description: Quality control of raw reads, adapter and quality trimming,
#              and quality control of the trimmed output.
#
# NOTE ON WHICH READS FEED WHICH ANALYSIS
#   Raw reads are used for the host alignment feeding VIRUSBreakend
#   (scripts 01, 02), since VIRUSBreakend performs its own read processing.
#   Trimmed paired reads (*_tr_1P / *_tr_2P) are used for the viral alignment,
#   coverage and consensus analyses (rest of the scripts).
########################################

########## CONFIG ##########
RAW_DIR="data/rawdata"
TRIM_DIR="data/trimmed"
RAW_QC_DIR="data/QC_raw"
TRIM_QC_DIR="data/QC_trimmed"

ADAPTER="ref/adapters/TruSeq3-PE.fa"        # set to the adapter file used
THREADS=8
MEMORY="16G"

########## SETUP ##########
mkdir -p "$TRIM_DIR" "$RAW_QC_DIR" "$TRIM_QC_DIR"

echo "========================================"
echo "Starting QC and trimming"
date
SECONDS=0

########## STEP 1: QC OF RAW READS ##########

for FQ in "$RAW_DIR"/*.fq.gz; do
    echo "FastQC (raw): $(basename "$FQ")"
    fastqc "$FQ" --outdir "$RAW_QC_DIR"
done

multiqc "$RAW_QC_DIR" -o "$RAW_QC_DIR" -n raw_multiqc_report.html

########## STEP 2: TRIMMOMATIC ##########

for R1 in "$RAW_DIR"/*_R1_001_1.fq.gz; do

    SAMPLE=$(basename "$R1" _R1_001_1.fq.gz)
    R2="$RAW_DIR/${SAMPLE}_R2_001_2.fq.gz"

    echo "----------------------------------------"
    date
    echo "Trimming: $SAMPLE"

    trimmomatic PE \
        -threads "$THREADS" \
        -Xmx"$MEMORY" \
        "$R1" "$R2" \
        "$TRIM_DIR/${SAMPLE}_R1_tr_1P.fq.gz" "$TRIM_DIR/${SAMPLE}_R1_tr_1U.fq.gz" \
        "$TRIM_DIR/${SAMPLE}_R2_tr_2P.fq.gz" "$TRIM_DIR/${SAMPLE}_R2_tr_2U.fq.gz" \
        ILLUMINACLIP:"$ADAPTER":2:30:10 \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:15 \
        MINLEN:36

done

########## STEP 3: QC OF TRIMMED PAIRED READS ##########

for FQ in "$TRIM_DIR"/*_tr_1P.fq.gz "$TRIM_DIR"/*_tr_2P.fq.gz; do
    echo "FastQC (trimmed): $(basename "$FQ")"
    fastqc "$FQ" --outdir "$TRIM_QC_DIR"
done

multiqc "$TRIM_QC_DIR" -o "$TRIM_QC_DIR" -n final_multiqc_report.html

########## RUNTIME ##########
duration=$SECONDS
echo "========================================"
echo "Pipeline completed"
echo "Total time: $((duration / 60)) min $((duration % 60)) sec"
date

