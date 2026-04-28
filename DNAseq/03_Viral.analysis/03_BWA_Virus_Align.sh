#!/usr/bin/env bash

set -euo pipefail

########################################
# Script: 03_wgs_align_virus.sh
# Description: Align FASTQ reads to viral reference (HIV),
#              extract mapped reads, and prepare GATK-ready BAMs
########################################

########## CONFIG ##########
INPUT_DIR="data/raw_fastq"
OUTPUT_DIR="data/virus_alignment"
FINAL_OUTPUT_DIR="data/virus_alignment_final"
REFERENCE="ref/hiv_reference.fasta"
THREADS=8

########## SETUP ##########
mkdir -p "$OUTPUT_DIR"
mkdir -p "$FINAL_OUTPUT_DIR"

echo "========================================"
echo "Starting viral alignment pipeline"
date

SECONDS=0

########## PROCESSING ##########

for R1 in "$INPUT_DIR"/*_R1_001_1.fq.gz; do

    SAMPLE=$(basename "$R1" | sed 's/_R1_001_1.fq.gz//')
    R2="$INPUT_DIR/${SAMPLE}_R2_001_2.fq.gz"

    SORTED_BAM="$OUTPUT_DIR/${SAMPLE}.sorted.bam"
    MAPPED_BAM="$OUTPUT_DIR/${SAMPLE}.virus_mapped.bam"
    FINAL_BAM="$FINAL_OUTPUT_DIR/${SAMPLE}.virus_final.bam"

    echo "----------------------------------------"
    date
    echo "Processing sample: $SAMPLE"

    ########## STEP 1: ALIGN TO VIRAL GENOME ##########

    bwa mem -t "$THREADS" "$REFERENCE" "$R1" "$R2" | \
        samtools view -bS - | \
        samtools sort -@ "$THREADS" -o "$SORTED_BAM" -

    samtools index "$SORTED_BAM"

    ########## STEP 2: EXTRACT MAPPED READS ##########

    samtools view -F 4 -b "$SORTED_BAM" > "$MAPPED_BAM"

    ########## STEP 3: ADD READ GROUPS (GATK READY) ##########

    samtools addreplacerg \
        -r "@RG\tID:${SAMPLE}\tSM:${SAMPLE}" \
        "$MAPPED_BAM" \
        -o "$FINAL_BAM"

    samtools index "$FINAL_BAM"

    echo "Final BAM ready: $FINAL_BAM"

done

########## RUNTIME ##########

duration=$SECONDS
echo "========================================"
echo "Pipeline completed"
echo "Total time: $(($duration / 60)) min $(($duration % 60)) sec"
date