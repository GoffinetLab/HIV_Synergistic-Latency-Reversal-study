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
    FILTERED_VCF="$OUTPUT_DIR/filtered_${SAMPLE}.vcf"

    echo "----------------------------------------"
    date
    echo "Processing sample: $SAMPLE"

    ########## STEP 1: RUN GRIDSS ##########

    docker-gridss virusbreakend \
        --reference "$REFERENCE" \
        --threads "$THREADS" \
        --workingdir "$OUTPUT_DIR" \
        --jar /opt/gridss/gridss-2.13.2-gridss-jar-with-dependencies.jar \
        --output "$RAW_VCF" \
        --db "$GRIDSS_DB" \
        "$BAM"

    ########## STEP 2: FILTER INTEGRATION SITES ##########

    echo "Filtering integration sites for: $SAMPLE"

    bcftools view \
        -i '((POS>=180 && POS<=182) || (POS>=8630 && POS<=8632) || \
             (ALT~"NC_001802.1:180" || ALT~"NC_001802.1:181" || ALT~"NC_001802.1:182") || \
             (ALT~"NC_001802.1:8630" || ALT~"NC_001802.1:8631" || ALT~"NC_001802.1:8632"))' \
        "$RAW_VCF" > "$FILTERED_VCF"

    echo "Filtered VCF saved: $FILTERED_VCF"

done

########## RUNTIME ##########

duration=$SECONDS
echo "========================================"
echo "Pipeline completed"
echo "Total time: $(($duration / 60)) min $(($duration % 60)) sec"
date