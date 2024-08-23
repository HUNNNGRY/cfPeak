#!/bin/bash

## test
#./getBed12NearestMRNA.sh -m /path/to/mRNA.bed12 -q /path/to/query.bed12 -c /path/to/chrom.sizes -d 10000 -o output.bed


# Default parameters
default_chrSizePath="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/transcriptome_genome_sort_uniq_newTxID"
default_mRNApath="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/bed_by_biotype/mRNA.bed"
nearGeneDist=100000

# Usage function
usage() {
    echo "Usage: $0 -m mRNA_bed12_path -q query_bed12_path [-c chrSizePath] [-d nearGeneDist] [-o outFile]"
    echo "  -m  Path to the mRNA BED12 file (default: $default_mRNApath)"
    echo "  -q  Path to the query BED12 file"
    echo "  -c  Path to the chromosome sizes file (default: $default_chrSizePath)"
    echo "  -d  Distance cutoff (default: $nearGeneDist)"
    echo "  -o  Output file (default: stdout)"
    exit 1
}

# Parse command-line arguments
while getopts "m:q:c:d:o:" opt; do
    case $opt in
        m) mRNA_bed12_Path=$OPTARG ;;
        q) query_Bed12_Path=$OPTARG ;;
        c) chrSizePath=$OPTARG ;;
        d) nearGeneDist=$OPTARG ;;
        o) outFile=$OPTARG ;;
        *) usage ;;
    esac
done

# Check if mandatory parameters are provided
if [ -z "$mRNA_bed12_Path" ] || [ -z "$query_Bed12_Path" ]; then
    usage
fi

# Set default chromosome sizes path if not provided
chrSizePath=${chrSizePath:-$default_chrSizePath}

# Extract TSS positions and sort
awk 'BEGIN {OFS="\t"} {
    if ($6 == "+") {
        $3 = $2 + 1;
    } else if ($6 == "-") {
        $2 = $3 - 1;
    }
    print $0;
}' $mRNA_bed12_Path | LC_COLLATE=C sort -k1,1 -k2,2n | \
#bedtools sort -g $chrSizePath -i stdin
# Find the nearest gene
bedtools closest -D ref -t first -a $query_Bed12_Path -b stdin | \
awk -v dist=$nearGeneDist '{
    if ($16 != "." && $25 != "." && ($25 <= dist && $25 >= -dist)) {
        print $0;
    }
}' > "${outFile:-/dev/stdout}" 

