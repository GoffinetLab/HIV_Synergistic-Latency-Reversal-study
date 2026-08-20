#!/usr/bin/env bash

set -euo pipefail

########################################
# Script: 05_snpEff_Variant_Annotate.sh
# Description: Annotate SNPs using SnpEff and generate tsv output tables
########################################

########## CONFIG ##########
INPUT_DIR="data/variant_calling"
OUTPUT_DIR="data/annotation"
SNPEFF_DIR="tools/snpEff"
GENOME_ID="HIV_NL43"
MEMORY="4g"

########## SETUP ##########
mkdir -p "$OUTPUT_DIR"

echo "========================================"
echo "Starting SnpEff annotation pipeline"
date

SECONDS=0

########## STEP 0: BUILD DATABASE (run once) ##########

if [[ ! -d "$SNPEFF_DIR/data/$GENOME_ID" ]]; then
    echo "Building SnpEff database..."

    java -jar "$SNPEFF_DIR/snpEff.jar" build \
        -gtf22 \
        -noGenome \
        -noCheckCds \
        -noCheckProtein \
        -v "$GENOME_ID"
fi

########## PROCESSING ##########

for VCF in "$INPUT_DIR"/*_filteredsnps.vcf; do

    SAMPLE=$(basename "$VCF" _filteredsnps.vcf)

    ANN_VCF="$OUTPUT_DIR/${SAMPLE}.ann.vcf"
    TIDY_OUT="$OUTPUT_DIR/${SAMPLE}.ann.tsv"

    echo "----------------------------------------"
    date
    echo "Processing sample: $SAMPLE"

    ########## STEP 1: RUN SNPEFF ##########

    java -Xmx"$MEMORY" -jar "$SNPEFF_DIR/snpEff.jar" \
        -v "$GENOME_ID" \
        "$VCF" > "$ANN_VCF"

    ########## STEP 2: TIDY OUTPUT ##########

    grep -v '^#' "$ANN_VCF" | \
    awk -F'\t' '{
        split($5, alt, ",");
        split($8, ann, "|");
        print $1, $2, $4, alt[1], ann[3], ann[4], ann[11], ann[12], ann[13], ann[14], ann[15], $7
    }' OFS='\t' > "$TIDY_OUT"

    echo "Annotated VCF: $ANN_VCF"
    echo "Tidy table: $TIDY_OUT"

done

########## RUNTIME ##########

duration=$SECONDS
echo "========================================"
echo "Pipeline completed"
echo "Total time: $(($duration / 60)) min $(($duration % 60)) sec"
date