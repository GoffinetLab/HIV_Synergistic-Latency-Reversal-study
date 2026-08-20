#!/usr/bin/env bash

set -euo pipefail

########################################
# Script: 04_GATK_Variant_Call_Filter
# Description: Call SNPs from viral BAMs using GATK HaplotypeCaller
#              and apply SNP filtering
########################################

########## CONFIG ##########
INPUT_DIR="data/virus_alignment_final"
OUTPUT_DIR="data/variant_calling"
REFERENCE="ref/hiv_reference.fasta"
THREADS=8

GATK_IMAGE="broadinstitute/gatk:4.1.3.0"

########## SETUP ##########
mkdir -p "$OUTPUT_DIR"

echo "========================================"
echo "Starting GATK variant calling pipeline"
date

SECONDS=0

########## STEP 0: INDEX REFERENCE ##########

if [[ ! -f "${REFERENCE}.fai" ]]; then
    echo "Indexing reference..."
    samtools faidx "$REFERENCE"
fi

if [[ ! -f "${REFERENCE%.fasta}.dict" ]]; then
    echo "Creating sequence dictionary..."
    gatk CreateSequenceDictionary -R "$REFERENCE"
fi

########## PROCESSING ##########

for BAM in "$INPUT_DIR"/*.bam; do

    SAMPLE=$(basename "$BAM" .bam)

    RAW_VCF="$OUTPUT_DIR/${SAMPLE}.vcf"
    SNP_VCF="$OUTPUT_DIR/${SAMPLE}_rawsnps.vcf"
    FILTERED_VCF="$OUTPUT_DIR/${SAMPLE}_filteredsnps.vcf"

    echo "----------------------------------------"
    date
    echo "Processing sample: $SAMPLE"

    ########## STEP 1: INDEX BAM ##########
    samtools index "$BAM"

    ########## STEP 2: HAPLOTYPECALLER ##########

    docker run --rm \
        -v "$(pwd)":/data \
        -w /data \
        "$GATK_IMAGE" \
        gatk HaplotypeCaller \
        -R "$REFERENCE" \
        -I "$BAM" \
        -O "$RAW_VCF" \
        --pcr-indel-model NONE \
        -ploidy 1 \
        -stand-call-conf 30 \
        -mbq 20 \
        -A QualByDepth

    ########## STEP 3: SELECT SNPs ##########

    docker run --rm \
        -v "$(pwd)":/data \
        -w /data \
        "$GATK_IMAGE" \
        gatk SelectVariants \
        -R "$REFERENCE" \
        -V "$RAW_VCF" \
        --select-type-to-include SNP \
        -O "$SNP_VCF"

    ########## STEP 4: FILTER SNPs ##########

    docker run --rm \
        -v "$(pwd)":/data \
        -w /data \
        "$GATK_IMAGE" \
        gatk VariantFiltration \
        -R "$REFERENCE" \
        -V "$SNP_VCF" \
        --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0" \
        --filter-name "LowConf" \
        -O "$FILTERED_VCF"

    echo "Filtered SNPs: $FILTERED_VCF"

done

########## RUNTIME ##########

duration=$SECONDS
echo "========================================"
echo "Pipeline completed"
echo "Total time: $(($duration / 60)) min $(($duration % 60)) sec"
date