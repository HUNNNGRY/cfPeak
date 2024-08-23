#!/bin/bash
INPUT=$1 # use cfpeak_by_sample not cfpeakCNN
awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $8 "\t" $10 "\t" $10}' ${INPUT} > ${INPUT}.tmp
