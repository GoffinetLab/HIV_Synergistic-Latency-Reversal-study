#!/bin/bash

## name sort bams from STAR output for counting in featurecounts


for bam in `ls -1 | grep .*.bam`;

do


echo "___________________________________________________________________"
date
echo "Now sorting"
echo $bam

#v=$(echo $bc | sed 's/\(.*\).bam/\1/') ;
v=$(echo $bam | sed 's/\(^P[0-9]\{4\}_[0-9]\{3\}\)_S[0-9]\{1,2\}.bam/\1/') ;
#echo "$v"_sorted.bam

samtools sort -n -o "$v"_sorted.bam -@ 8 $bc ;

done
