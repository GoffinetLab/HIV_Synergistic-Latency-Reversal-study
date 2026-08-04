#!/usr/bin/env bash

set -euo pipefail

########################################
# Script: 06_bcftools_Consensus.sh
# Description: Call variants against the proviral reference and build
#              consensus sequences for each clone.
########################################

########## CONFIG ##########
INPUT_DIR="data/virus_alignment_final"     # *.hiv.rg.bam from script 03
VCF_DIR="data/variant_calling"
PASS_VCF_DIR="data/vcf_pass"
FASTA_DIR="data/consensus_fasta"
FINAL_DIR="data/final"

REFERENCE="ref/NL43_HIVENVGFP.fasta"       # identical-LTR proviral reference
CONTIG="HIV_ENVGFP"

GATK_IMAGE="broadinstitute/gatk:4.1.3.0"

REFDIR=$(cd "$(dirname "$REFERENCE")" && pwd)
REFBASE=$(basename "$REFERENCE")

########## SETUP ##########
mkdir -p "$VCF_DIR" "$PASS_VCF_DIR" "$FASTA_DIR" "$FINAL_DIR"

echo "========================================"
echo "Starting variant calling and consensus generation"
date
SECONDS=0

########## STEP 0: INDEX REFERENCE ##########

if [[ ! -f "${REFERENCE}.fai" ]]; then
    echo "Indexing reference..."
    samtools faidx "$REFERENCE"
fi

if [[ ! -f "${REFERENCE%.fasta}.dict" ]]; then
    echo "Creating sequence dictionary..."
    docker run --rm -v "$REFDIR":/ref -w /ref "$GATK_IMAGE" \
        gatk CreateSequenceDictionary -R "/ref/$REFBASE"
fi

########## STEP 1: VARIANT CALLING ##########

for BAM in "$INPUT_DIR"/*.hiv.rg.bam; do

    SAMPLE=$(basename "$BAM" .hiv.rg.bam)

    echo "----------------------------------------"
    date
    echo "Calling variants: $SAMPLE"

    docker run --rm \
        -v "$(pwd)":/data -v "$REFDIR":/ref -w /data \
        "$GATK_IMAGE" \
        gatk HaplotypeCaller \
        -R "/ref/$REFBASE" \
        -I "$BAM" \
        -O "$VCF_DIR/${SAMPLE}.vcf" \
        -L "$CONTIG" \
        --pcr-indel-model NONE \
        -ploidy 1 \
        -stand-call-conf 30 \
        -mbq 20 \
        -A QualByDepth

    docker run --rm \
        -v "$(pwd)":/data -v "$REFDIR":/ref -w /data \
        "$GATK_IMAGE" \
        gatk SelectVariants \
        -R "/ref/$REFBASE" \
        -V "$VCF_DIR/${SAMPLE}.vcf" \
        --select-type-to-include SNP \
        -L "$CONTIG" \
        -O "$VCF_DIR/${SAMPLE}_rawsnps.vcf"

    docker run --rm \
        -v "$(pwd)":/data -v "$REFDIR":/ref -w /data \
        "$GATK_IMAGE" \
        gatk VariantFiltration \
        -R "/ref/$REFBASE" \
        -V "$VCF_DIR/${SAMPLE}_rawsnps.vcf" \
        -L "$CONTIG" \
        --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 4.0" \
        --filter-name "LowConf" \
        -O "$VCF_DIR/${SAMPLE}_filteredsnps.vcf"

    echo -n "  raw calls: "
    grep -vc '^#' "$VCF_DIR/${SAMPLE}.vcf" || true

done

########## STEP 2: BUILD CONSENSUS ##########

for VCF in "$VCF_DIR"/*_filteredsnps.vcf; do

    SAMPLE=$(basename "$VCF" _filteredsnps.vcf)

    echo "Generating consensus: $SAMPLE"

    bcftools view -f PASS -O z \
        -o "$PASS_VCF_DIR/${SAMPLE}.pass.vcf.gz" "$VCF"
    bcftools index -f "$PASS_VCF_DIR/${SAMPLE}.pass.vcf.gz"

    echo -n "  variants applied: "
    bcftools view -H "$PASS_VCF_DIR/${SAMPLE}.pass.vcf.gz" | wc -l

    bcftools consensus -f "$REFERENCE" \
        "$PASS_VCF_DIR/${SAMPLE}.pass.vcf.gz" \
        > "$FASTA_DIR/${SAMPLE}.consensus.fasta"

    awk -v s="$SAMPLE" '/^>/{print ">"s; next}{print}' \
        "$FASTA_DIR/${SAMPLE}.consensus.fasta" > tmp && \
        mv tmp "$FASTA_DIR/${SAMPLE}.consensus.fasta"

done

########## STEP 3: CHECKS ##########

echo "----------------------------------------"
echo "Consensus summary:"
for FILE in "$FASTA_DIR"/*.consensus.fasta; do
    SAMPLE=$(basename "$FILE" .consensus.fasta)
    LEN=$(awk '!/^>/{n+=length($0)} END{print n}' "$FILE")
    NS=$(awk '!/^>/{gsub(/[^Nn]/,""); n+=length($0)} END{print n+0}' "$FILE")
    printf "  %-20s %s bp, %s ambiguous bases\n" "$SAMPLE" "$LEN" "$NS"
done

cat "$FASTA_DIR"/*.consensus.fasta > "$FINAL_DIR/all_clones.fasta"

########## RUNTIME ##########
duration=$SECONDS
echo "========================================"
echo "Pipeline completed"
echo "Total time: $((duration / 60)) min $((duration % 60)) sec"
date

# RESULT (this dataset): no high-confidence variants were called in any clone.
# Confirmed independently with bcftools mpileup/call. The consensus sequences
# are therefore identical to the reference construct.
