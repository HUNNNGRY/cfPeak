#!/bin/bash
INPUT=$1 # use cfpeak_by_sample not cfpeakCNN
OUTPUT=$2
awk '{print $1 "\t" $2 "\t" $3 "\t" "peak_" NR "\t" $5 "\t" "+"}' ${INPUT} > ${OUTPUT}
