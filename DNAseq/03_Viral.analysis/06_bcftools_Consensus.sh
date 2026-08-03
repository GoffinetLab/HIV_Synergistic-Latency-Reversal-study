#!/usr/bin/env bash

set -euo pipefail

########################################
# Script: 06_bcftools_Consensus.sh
# Description: Clean filtered VCFs, generate consensus FASTA sequences with
#              low-coverage regions masked, and combine into a multi-FASTA.
########################################

########## CONFIG ##########
INPUT_DIR="data/annotation"
BAM_DIR="data/virus_alignment_final"
CLEAN_VCF_DIR="data/vcf_clean"
MASK_DIR="data/lowcov_mask"
FASTA_DIR="data/consensus_fasta"
FINAL_DIR="data/final"
REFERENCE="ref/hiv_reference.fasta"

MIN_COV=10          # positions below this depth are masked as N

########## SETUP ##########
mkdir -p "$CLEAN_VCF_DIR" "$MASK_DIR" "$FASTA_DIR" "$FINAL_DIR"

echo "========================================"
echo "Starting consensus sequence generation"
date
SECONDS=0

########## STEP 1: CLEAN, COMPRESS AND INDEX VCFs ##########

for VCF in "$INPUT_DIR"/*.ann.vcf; do

    SAMPLE=$(basename "$VCF" .ann.vcf)
    CLEAN_VCF="$CLEAN_VCF_DIR/${SAMPLE}.clean.vcf.gz"

    echo "----------------------------------------"
    date
    echo "Cleaning and indexing VCF: $SAMPLE"

    # -O z writes bgzip-compressed VCF, which bcftools index requires
    bcftools annotate -x FILTER -O z -o "$CLEAN_VCF" "$VCF"
    bcftools index -f "$CLEAN_VCF"

done

########## STEP 2: BUILD LOW-COVERAGE MASKS ##########

for BAM in "$BAM_DIR"/*.bam; do

    SAMPLE=$(basename "$BAM" .virus_final.bam)
    MASK="$MASK_DIR/${SAMPLE}.lowcov.bed"

    echo "Building low-coverage mask (<${MIN_COV}x): $SAMPLE"

    bedtools genomecov -ibam "$BAM" -bga \
        | awk -v m="$MIN_COV" 'BEGIN{OFS="\t"} $4 < m {print $1,$2,$3}' \
        > "$MASK"

done

########## STEP 3: GENERATE MASKED CONSENSUS FASTA ##########

for VCF in "$CLEAN_VCF_DIR"/*.clean.vcf.gz; do

    SAMPLE=$(basename "$VCF" .clean.vcf.gz)
    MASK="$MASK_DIR/${SAMPLE}.lowcov.bed"
    FASTA_OUT="$FASTA_DIR/${SAMPLE}.fasta"

    echo "Generating consensus for: $SAMPLE"

    if [[ -s "$MASK" ]]; then
        bcftools consensus -f "$REFERENCE" -m "$MASK" "$VCF" > "$FASTA_OUT"
    else
        bcftools consensus -f "$REFERENCE" "$VCF" > "$FASTA_OUT"
    fi

done

########## STEP 4: RENAME FASTA HEADERS ##########

for FILE in "$FASTA_DIR"/*.fasta; do

    SAMPLE=$(basename "$FILE" .fasta)
    [[ "$SAMPLE" == *_renamed ]] && continue

    RENAMED="$FASTA_DIR/${SAMPLE}_renamed.fasta"

    echo "Renaming header: $SAMPLE"

    # keep the sequence line-wrapped rather than collapsing to one long line
    awk -v s="$SAMPLE" '/^>/{print ">"s; next}{print}' "$FILE" > "$RENAMED"

done

########## STEP 5: COMBINE ##########

COMBINED="$FINAL_DIR/combined.fasta"
cat "$FASTA_DIR"/*_renamed.fasta > "$COMBINED"
echo "Combined FASTA created: $COMBINED"

########## RUNTIME ##########
duration=$SECONDS
echo "========================================"
echo "Pipeline completed"
echo "Total time: $((duration / 60)) min $((duration % 60)) sec"
date
