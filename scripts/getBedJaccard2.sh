#!/bin/bash

# More advanced version

# Calculate overlap
cmb_file=$1  # cmb.txt
jac_dir=$2   # "test/$dst/jaccard"
strand_info=${3:-stranded}  # Strand information, default is stranded if not set
f_value=${4:-}  # -f flag value, default is empty if not set
F_value=${5:-}  # -F flag value, default is empty if not set
mkdir -p $jac_dir

echo "Strand: $strand_info"
echo "-f: ${f_value:-default}"
echo "-F: ${F_value:-default}"
# awk -f scripts/combinations.awk <<< `ls test/${dst}/*_*.bed_optimal_cfpeak.bed | tr "\n" " "` > ${jac_dir}/cmb.txt

rm -f ${jac_dir}/jac.txt
echo -e "from\tto\tfrom_n\tto_n\tlen_intersection\tlen_union\tjaccard\tn_intersections\tfrom_only_n\tto_only_n\toverlap_n_from\toverlap_n_to" >> ${jac_dir}/jac.txt

for tmp in `cat $cmb_file | awk '{print $1 "," $2}'`
do
    i=`echo $tmp | cut -d ',' -f 1`
    j=`echo $tmp | cut -d ',' -f 2`
    echo -e "${i}:${j}"
    i_n=`wc -l $i | cut -d ' ' -f 1`
    j_n=`wc -l $j | cut -d ' ' -f 1`

    # Construct the bedtools jaccard command with optional flags
    jaccard_cmd="bedtools jaccard -split -a ${i} -b ${j}"
    intersect_cmd_from="bedtools intersect -split -v -wa -a ${i} -b ${j}"
    intersect_cmd_to="bedtools intersect -split -v -wa -a ${j} -b ${i}"
    intersect_cmd_overlap_from="bedtools intersect -split -u -wa -a ${i} -b ${j}" # same in py.pkg: intervene
    intersect_cmd_overlap_to="bedtools intersect -split -u -wa -a ${j} -b ${i}"

    [ -n "$f_value" ] && jaccard_cmd+=" -f $f_value" && intersect_cmd_from+=" -f $f_value" && intersect_cmd_to+=" -f $f_value" && intersect_cmd_overlap_from+=" -f $f_value" && intersect_cmd_overlap_to+=" -f $f_value"
    [ -n "$F_value" ] && jaccard_cmd+=" -F $F_value" && intersect_cmd_from+=" -F $F_value" && intersect_cmd_to+=" -F $F_value" && intersect_cmd_overlap_from+=" -F $F_value" && intersect_cmd_overlap_to+=" -F $F_value"
    [ "$strand_info" == "stranded" ] && jaccard_cmd+=" -s" && intersect_cmd_from+=" -s" && intersect_cmd_to+=" -s" && intersect_cmd_overlap_from+=" -s" && intersect_cmd_overlap_to+=" -s"

    # Execute the commands and process the output
    jaccard_output=$($jaccard_cmd | grep -v 'jaccard' | awk -v i=$i -v j=$j -v i_n=$i_n -v j_n=$j_n 'BEGIN {FS=OFS="\t"} {print i,j,i_n,j_n,$0}')
    from_only_n=$(eval "$intersect_cmd_from | wc -l")
    to_only_n=$(eval "$intersect_cmd_to | wc -l")
    overlap_n_from=$(eval "$intersect_cmd_overlap_from | wc -l")
    overlap_n_to=$(eval "$intersect_cmd_overlap_to | wc -l")

    # Append results to jac.txt
    echo -e "${jaccard_output}\t${from_only_n}\t${to_only_n}\t${overlap_n_from}\t${overlap_n_to}" >> ${jac_dir}/jac.txt
done

