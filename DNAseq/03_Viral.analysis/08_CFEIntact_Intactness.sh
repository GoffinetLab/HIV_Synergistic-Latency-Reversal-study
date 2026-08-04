#!/usr/bin/env bash

set -euo pipefail

########################################
# Script: 08_CFEIntact_Intactness.sh
# Description: Assess genetic intactness of the consensus proviral sequences
#              with CFEIntact, and run a stop-codon positive control.
#
# NOTES
#
# 1. --subtype NL43 uses the NL4-3 reference added by the CFEIntact developers
#    at our request; more appropriate for this construct than generic subtype B.
#
# 2. The positive control inserts in-frame stop codons into gag and pol
#    (introduce_stop_codons.py) and repeats the analysis, confirming the tool
#    detects genuine coding defects rather than defaulting to an intact call.
########################################

########## CONFIG ##########
FASTA_DIR="data/consensus_fasta"           # output of script 06
RESULT_DIR="data/intactness"
SCRIPT_DIR="scripts"

CFEINTACT_IMAGE="cfelab/cfeintact:1.26.1-206-g3ed0801"
SUBTYPE="NL43"

########## SETUP ##########
mkdir -p "$RESULT_DIR"

echo "========================================"
echo "Starting CFEIntact intactness analysis"
date
SECONDS=0

docker pull "$CFEINTACT_IMAGE"
docker run --rm "$CFEINTACT_IMAGE" version

########## STEP 1: REAL CONSENSUS SEQUENCES ##########

for FILE in "$FASTA_DIR"/*.consensus.fasta; do

    SAMPLE=$(basename "$FILE" .consensus.fasta)

    echo "----------------------------------------"
    echo "Checking: $SAMPLE"

    docker run --rm -v "$(pwd)":/data -w /data \
        "$CFEINTACT_IMAGE" \
        check "$FILE" \
        --subtype "$SUBTYPE" \
        --output-csv \
        --output "$RESULT_DIR/real_${SAMPLE}"

done

########## STEP 2: BUILD THE POSITIVE CONTROL ##########
# introduce_stop_codons.py inserts in-frame TAA at gag 892 and pol 2502

for FILE in "$FASTA_DIR"/*.consensus.fasta; do
    echo "Inserting stop codons: $(basename "$FILE")"
    python3 "$SCRIPT_DIR/introduce_stop_codons.py" "$FILE"
done

########## STEP 3: POSITIVE CONTROL ##########

for FILE in "$FASTA_DIR"/*.consensus_stoptest.fasta; do

    SAMPLE=$(basename "$FILE" .consensus_stoptest.fasta)

    echo "----------------------------------------"
    echo "Checking control: $SAMPLE"

    docker run --rm -v "$(pwd)":/data -w /data \
        "$CFEINTACT_IMAGE" \
        check "$FILE" \
        --subtype "$SUBTYPE" \
        --output-csv \
        --output "$RESULT_DIR/control_${SAMPLE}"

done

########## STEP 4: COLLATE ##########

echo "========================================"
echo "Real sequences:"
cat "$RESULT_DIR"/real_*/defects.csv
echo
echo "Positive control:"
cat "$RESULT_DIR"/control_*/defects.csv

########## RUNTIME ##########
duration=$SECONDS
echo "========================================"
echo "Pipeline completed"
echo "Total time: $((duration / 60)) min $((duration % 60)) sec"
date

# RESULT (this dataset):
#   Real sequences   - the only defect reported in any clone is the env
#                      insertion, i.e. the engineered GFP reporter cassette
#                      (720 insertions against a tolerance of 123). No stop
#                      codons, frameshifts or deletions in gag or pol.
#   Positive control - internal stop codons recovered at gag 892 and pol 2502,
#                      plus a consequent gag frameshift call, in all three
#                      clones. Absent from the real sequences.
