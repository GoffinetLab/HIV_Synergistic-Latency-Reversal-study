#!/bin/bash

## count feature expression from bam files using featurecount from subread package

#run with 8 cores, rev stranded, paired end

## run in directory with sorted bam files

#list=$(echo *_sorted.bam)
list=$(echo `ls -1 | grep ".*_S[0-9]\{1,3\}.bam"`)

featureCounts -T 8 -s 2 -t exon -p --countReadPairs -a "/path/to/reference_gtf_file/Homo_sapiens.GRCh38.109.gtf" \
 	-o "/path/to/countmatrix_output/BryoJQ1_AllSamples.Rmatrix.txt" \
 	$list
