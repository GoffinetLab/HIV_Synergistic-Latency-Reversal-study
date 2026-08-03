#!/usr/bin/env bash

set -euo pipefail

########################################
# Script: 02_Gridss_Intsites.sh
# Description: Detect viral integration sites using GRIDSS VirusBreakend
########################################

########## CONFIG ##########
INPUT_DIR="data/processed_bam"
OUTPUT_DIR="data/output_virusbreakend"
REFERENCE="ref/host.fa"
GRIDSS_DB="ref/virusbreakenddb"
GRIDSS_JAR="/opt/gridss/gridss-2.13.2-gridss-jar-with-dependencies.jar"
THREADS=16

########## SETUP ##########
mkdir -p "$OUTPUT_DIR"

# Activate conda safely
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gridss

echo "========================================"
echo "Starting GRIDSS VirusBreakend pipeline"
date
SECONDS=0

########## PROCESSING ##########

for BAM in "$INPUT_DIR"/*.bam; do

    SAMPLE=$(basename "$BAM" .bam)

    RAW_VCF="$OUTPUT_DIR/${SAMPLE}.vcf"
    PASS_VCF="$OUTPUT_DIR/${SAMPLE}.pass.vcf"
    SUMMARY="$OUTPUT_DIR/${SAMPLE}.integration_sites.tsv"

    echo "----------------------------------------"
    date
    echo "Processing sample: $SAMPLE"

    ########## STEP 1: RUN VIRUSBREAKEND ##########

    docker-gridss virusbreakend \
        --reference "$REFERENCE" \
        --threads "$THREADS" \
        --workingdir "$OUTPUT_DIR" \
        --jar "$GRIDSS_JAR" \
        --output "$RAW_VCF" \
        --db "$GRIDSS_DB" \
        "$BAM"

    ########## STEP 2: KEEP JUNCTIONS PASSING QUALITY FILTERS ##########

    echo "Selecting PASS junctions for: $SAMPLE"

    bcftools view -f PASS "$RAW_VCF" > "$PASS_VCF"

    ########## STEP 3: SUMMARISE INTEGRATION SITES ##########
    # One row per breakend, with the evidence needed to judge each call:
    #   QUAL - GRIDSS variant quality score
    #   BVF  - independent DNA fragments supporting the junction
    #   BSC  - split reads crossing the virus-host boundary
    # A genuine integration in a clonal line typically yields two junctions
    # (one per LTR end) with support orders of magnitude above background.

    {
      printf "sample\tchrom\tpos\tid\talt\tqual\tBVF\tBSC\n"
      bcftools query -f '%CHROM\t%POS\t%ID\t%ALT\t%QUAL\t%INFO/BVF\t%INFO/BSC\n' \
          "$PASS_VCF" \
        | awk -v s="$SAMPLE" 'BEGIN{OFS="\t"} {print s, $0}'
    } > "$SUMMARY"

    echo "PASS VCF:  $PASS_VCF"
    echo "Summary:   $SUMMARY"

done

########## RUNTIME ##########
duration=$SECONDS
echo "========================================"
echo "Pipeline completed"
echo "Total time: $((duration / 60)) min $((duration % 60)) sec"
date
