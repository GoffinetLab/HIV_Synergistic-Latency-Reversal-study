#!/usr/bin/env bash

set -euo pipefail

########################################
# Script: 03_BWA_Virus_Align.sh
# Description: Align trimmed FASTQ reads to a COMBINED human + HIV reference,
#              mark duplicates, and extract HIV-mapped reads for downstream
#              variant calling and coverage analysis.
########################################

########## CONFIG ##########
INPUT_DIR="data/trimmed_fastq"
OUTPUT_DIR="data/virus_alignment"
FINAL_OUTPUT_DIR="data/virus_alignment_final"

# combined human + HIV (identical-LTR) reference, indexed with bwa index
REFERENCE="ref/combined_human_HIVnew.fasta"
VIRUS_CONTIG="HIV_ENVGFP"

THREADS=32

########## SETUP ##########
mkdir -p "$OUTPUT_DIR" "$FINAL_OUTPUT_DIR"

echo "========================================"
echo "Starting viral alignment pipeline"
date
SECONDS=0

########## PROCESSING ##########

for R1 in "$INPUT_DIR"/*_R1_tr_1P.fq.gz; do

    SAMPLE=$(basename "$R1" _R1_tr_1P.fq.gz)
    R2="$INPUT_DIR/${SAMPLE}_R2_tr_2P.fq.gz"

    MARKDUP_BAM="$OUTPUT_DIR/${SAMPLE}.markdup.bam"
    VIRUS_BAM="$OUTPUT_DIR/${SAMPLE}.virus_mapped.bam"
    FINAL_BAM="$FINAL_OUTPUT_DIR/${SAMPLE}.virus_final.bam"

    echo "----------------------------------------"
    date
    echo "Processing sample: $SAMPLE"

    ########## STEP 1: ALIGN, FIX MATES, MARK DUPLICATES ##########
    # -M   flags split hits as secondary (required by some downstream tools)
    # fixmate -m  adds the mate score tags that markdup needs
    # markdup     flags PCR/optical duplicates so depth reflects real molecules

    bwa mem -t "$THREADS" -M "$REFERENCE" "$R1" "$R2" \
        | samtools fixmate -m -u - - \
        | samtools sort -@ 5 -u - \
        | samtools markdup -@ 5 - "$MARKDUP_BAM"

    samtools index "$MARKDUP_BAM"
    samtools quickcheck -v "$MARKDUP_BAM"

    ########## STEP 2: EXTRACT HIV-MAPPED READS ##########
    # Restrict to the viral contig. Reads whose mate maps to the host are
    # retained by region extraction; no separate mate-chasing step is needed
    # for coverage, and appending BAMs with shell redirection would corrupt
    # the file (a second header and EOF marker mid-file).

    samtools view -b "$MARKDUP_BAM" "$VIRUS_CONTIG" > "$VIRUS_BAM"
    samtools index "$VIRUS_BAM"

    ########## STEP 3: ADD READ GROUPS (GATK READY) ##########

    samtools addreplacerg \
        -r "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${SAMPLE}" \
        "$VIRUS_BAM" \
        -o "$FINAL_BAM"

    samtools index "$FINAL_BAM"

    echo "Final BAM ready: $FINAL_BAM"

done

########## RUNTIME ##########
duration=$SECONDS
echo "========================================"
echo "Pipeline completed"
echo "Total time: $((duration / 60)) min $((duration % 60)) sec"
date
