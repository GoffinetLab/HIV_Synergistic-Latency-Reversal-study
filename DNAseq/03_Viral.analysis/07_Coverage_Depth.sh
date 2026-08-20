#!/usr/bin/env bash

set -euo pipefail

########################################
# Script: 07_Coverage_Depth.sh
# Description: Per-base sequencing depth across the provirus, breadth-of-
#              coverage statistics, and bedgraph tracks for plotting.
#
# Reported statistics use median depth plus breadth rather than a raw minimum:
# the raw minimum is dominated by the read-placement ramp at the very ends of
# a linear reference and understates the data.
########################################

########## CONFIG ##########
BAM_DIR="data/virus_alignment"
OUTPUT_DIR="data/coverage"
VIRUS_CONTIG="HIV_ENVGFP"
MIN_DEPTH=30

########## SETUP ##########
mkdir -p "$OUTPUT_DIR"

echo "========================================"
echo "Starting coverage analysis"
date
SECONDS=0

########## PROCESSING ##########

for BAM in "$BAM_DIR"/*.markdup.bam; do

    SAMPLE=$(basename "$BAM" .markdup.bam)

    DEPTH="$OUTPUT_DIR/${SAMPLE}.depth.txt"
    BEDGRAPH="$OUTPUT_DIR/${SAMPLE}.bedgraph"

    echo "----------------------------------------"
    echo "Processing sample: $SAMPLE"

    ########## STEP 1: PER-BASE DEPTH ##########
    # -a            report every position, including zeros, so breadth is
    #               calculated over the full contig length
    # --excl-flags  count primary, non-duplicate reads only

    samtools depth -a --excl-flags UNMAP,SECONDARY,QCFAIL,DUP \
        -r "$VIRUS_CONTIG" "$BAM" > "$DEPTH"

    ########## STEP 2: BEDGRAPH FOR PLOTTING ##########
    # zero depths are floored to 0.1 so a log-scaled y-axis can render them

    awk -v c="$VIRUS_CONTIG" 'BEGIN{OFS="\t"}
        {d = ($3 == 0) ? 0.1 : $3; print c, $2-1, $2, d}' \
        "$DEPTH" > "$BEDGRAPH"

    ########## STEP 3: DEPTH AND BREADTH STATISTICS ##########

    gawk -v s="$SAMPLE" -v t="$MIN_DEPTH" '
        {d[NR]=$3; sum+=$3; n++; if($3>=1)c1++; if($3>=t)ct++}
        END{
            asort(d);
            med = (n%2) ? d[(n+1)/2] : (d[n/2]+d[n/2+1])/2;
            printf "%s  mean=%.1f  median=%d  breadth>=1x=%.1f%%  breadth>=%dx=%.1f%%\n",
                   s, sum/n, med, 100*c1/n, t, 100*ct/n;
        }' "$DEPTH"

done

########## RUNTIME ##########
duration=$SECONDS
echo "========================================"
echo "Coverage analysis completed"
echo "Total time: $((duration / 60)) min $((duration % 60)) sec"
date
