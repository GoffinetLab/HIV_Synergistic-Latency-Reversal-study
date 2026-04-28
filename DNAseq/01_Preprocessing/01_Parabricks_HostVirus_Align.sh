#!/usr/bin/env bash

set -euo pipefail

########################################
# Script: 01_Parabricks_Align.sh
# Description: Align FASTQ files to hg38 using Parabricks fq2bam
########################################

########## CONFIG ##########
INPUT_DIR="data/input"
OUTPUT_DIR="data/output"
REFERENCE="ref/host.fa"
THREADS=8

########## SETUP ##########
mkdir -p "$OUTPUT_DIR"

########## PROCESSING ##########

for R1 in "$INPUT_DIR"/*_R1_001_1.fq.gz; do
    
    # Extract sample name
    SAMPLE=$(basename "$R1" | sed 's/_R1_001_1.fq.gz//')
    R2="$INPUT_DIR/${SAMPLE}_R2_001_2.fq.gz"

    echo "========================================"
    date
    echo "Processing sample: $SAMPLE"

    ########## RUN ALIGNMENT ##########
    
    docker run --rm --gpus all \
        -v "$INPUT_DIR":/workdir \
        -v "$OUTPUT_DIR":/outputdir \
        -w /workdir \
        nvcr.io/nvidia/clara/clara-parabricks:4.0.1-1 \
        pbrun fq2bam \
        --ref /workdir/$REFERENCE \
        --in-fq /workdir/$(basename "$R1") /workdir/$(basename "$R2") \
        --out-bam /outputdir/${SAMPLE}.bam \
        --logfile /outputdir/${SAMPLE}.log

done