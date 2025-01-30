#!/bin/bash

# Note
# default input is forward lib: R1 is RNA template
# if take reverse lib as input like SMARTer Pico v2/v3, need switch R1 R2 in parameter: -1 test_2.fq.gz -2 test_1.fq.gz 

# Usage function
usage() {
    echo "Usage: $0 -t THREADS -1 FASTQ1 -2 FASTQ2 -o MERGED_FQ -d DOWNSAMPLE -m MIN_FQ -M MAX_FQ -s SEED -T TMP_FQ"
    exit 1
}

# Parse arguments
while getopts "t:1:2:o:d:m:M:s:T:" opt; do
    case "$opt" in
        t) threads=$OPTARG ;;
        1) fastq1=$OPTARG ;;
        2) fastq2=$OPTARG ;;
        o) fastp_fq=$OPTARG ;;
        d) downsample_fq=$OPTARG ;;
        m) min_fq_num=$OPTARG ;;
        M) max_fq_num=$OPTARG ;;
        s) seed=$OPTARG ;;
        T) tmp_fq=$OPTARG ;;
        *) usage ;;
    esac
done

# Ensure all required arguments are provided
if [ -z "$threads" ] || [ -z "$fastq1" ] || [ -z "$fastq2" ] || [ -z "$fastp_fq" ] || 
   [ -z "$downsample_fq" ] || [ -z "$min_fq_num" ] || [ -z "$max_fq_num" ] || 
   [ -z "$seed" ] || [ -z "$tmp_fq" ]; then
    usage
fi

# Run fastp for merging
fastp -w "$threads" -i "$fastq1" -I "$fastq2" --merge --merged_out "$fastp_fq"
if [ $? -ne 0 ]; then
    echo "Error: fastp failed"
    exit 1
fi

# Get read count
num=$(seqkit stats --basename --tabular -j "$threads" "$fastp_fq" | grep -v "file" | cut -f 4)
if [ -z "$num" ]; then
    echo "Error: Failed to retrieve read count"
    exit 1
fi

# Calculate subsampling number
sub_num=$(( num * downsample_fq / 100 ))
if [ "$num" -lt "$min_fq_num" ]; then
    sub_num=$num
elif [ "$num" -gt "$min_fq_num" ] && [ "$sub_num" -lt "$min_fq_num" ]; then
    sub_num=$min_fq_num
elif [ "$sub_num" -gt "$max_fq_num" ]; then
    sub_num=$max_fq_num
fi

echo "Using fq num: $sub_num"

# Check if subsampling count is zero
if [ "$sub_num" -eq 0 ]; then
    echo "Error: Zero num of fq"
    exit 1
fi

# Perform subsampling and compression
seqtk sample -s "$seed" "$fastp_fq" "$sub_num" | pigz -c -p "$threads" > "$tmp_fq"
if [ $? -ne 0 ]; then
    echo "Error: seqtk sampling failed"
    exit 1
fi

echo "Processing completed successfully!"

