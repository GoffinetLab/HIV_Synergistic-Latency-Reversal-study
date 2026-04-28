#!/usr/bin/env bash

set -euo pipefail

########################################
# Script: 06_bcftools_Consensus.sh
# Description: Clean filtered VCFs, generate consensus FASTA sequences,
#              and combine into a single multi-FASTA file
########################################

########## CONFIG ##########
INPUT_DIR="data/annotation"
CLEAN_VCF_DIR="data/vcf_clean"
FASTA_DIR="data/consensus_fasta"
FINAL_DIR="data/final"
REFERENCE="ref/hiv_reference.fasta"

########## SETUP ##########
mkdir -p "$CLEAN_VCF_DIR"
mkdir -p "$FASTA_DIR"
mkdir -p "$FINAL_DIR"

echo "========================================"
echo "Starting consensus sequence generation"
date

SECONDS=0

########## STEP 1: REMOVE FILTER FIELD ##########

for VCF in "$INPUT_DIR"/*.ann.vcf; do

    SAMPLE=$(basename "$VCF" .ann.vcf)
    CLEAN_VCF="$CLEAN_VCF_DIR/${SAMPLE}.clean.vcf"

    echo "----------------------------------------"
    date
    echo "Cleaning VCF: $SAMPLE"

    bcftools annotate -x FILTER \
        -o "$CLEAN_VCF" \
        -O v "$VCF"

done

########## STEP 2: INDEX VCF ##########

for VCF in "$CLEAN_VCF_DIR"/*.vcf; do
    echo "Indexing VCF: $VCF"
    bcftools index "$VCF"
done

########## STEP 3: GENERATE CONSENSUS FASTA ##########

for VCF in "$CLEAN_VCF_DIR"/*.vcf; do

    SAMPLE=$(basename "$VCF" .clean.vcf)
    FASTA_OUT="$FASTA_DIR/${SAMPLE}.fasta"

    echo "Generating consensus for: $SAMPLE"

    bcftools consensus \
        -f "$REFERENCE" \
        "$VCF" > "$FASTA_OUT"

done

########## STEP 4: FORMAT FASTA ##########

for FILE in "$FASTA_DIR"/*.fasta; do

    SAMPLE=$(basename "$FILE" .fasta)
    MODIFIED="$FASTA_DIR/${SAMPLE}_modified.fasta"

    echo "Formatting FASTA: $SAMPLE"

    grep -v "^>" "$FILE" | \
        tr -d '\n' | \
        sed "1 s/^/>${SAMPLE}\n/" > "$MODIFIED"

done

########## STEP 5: COMBINE FASTA ##########

COMBINED="$FINAL_DIR/combined.fasta"

cat "$FASTA_DIR"/*_modified.fasta > "$COMBINED"

echo "Combined FASTA created: $COMBINED"

########## RUNTIME ##########

duration=$SECONDS
echo "========================================"
echo "Pipeline completed"
echo "Total time: $(($duration / 60)) min $(($duration % 60)) sec"
date