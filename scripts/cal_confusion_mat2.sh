
#!/bin/bash

while getopts "i:s:t:m:r:" op
do
	case "$op" in
		i)  sample="$OPTARG";;
		s)  standard="$OPTARG";;
		t)  threads="$OPTARG";;
		m)  method="$OPTARG";;
		r)  replicate="$OPTARG";;
		\?) exit 1;;
	esac
done
# f)  standardNum="$OPTARG";;

# echo $sample
# echo $standard
# echo $threads

# variables and functions ---------------------------------------------------------------

calculate_rate () {
	awk -v term_1=$1 -v term_2=$2 'BEGIN {print (term_1 + 0.0001) / (term_1 + term_2 + 0.0001)}' # + 0.00001 avoid piranha fdr division by zero (cut -b1-5 : truncate first 5 number)
}
# why use cut
calculate_f1_score () {
	awk -v term_1=$1 -v term_2=$2 'BEGIN {print 2 * (term_1 * term_2 + 0.0001) / (term_1 + term_2 + 0.0001)}' # | cut -b1-5
}
# f1 division zero error

threshold() {
	counts=$1
	predicted_truth=$(sort -n -k5 $sample | awk -v s=$1 '$5 >= s')
	predicted_false=$(sort -n -k5 $sample | awk -v s=$1 '$5 < s')
	TP=$(echo "$predicted_truth" | cut -f1-3 | bedtools intersect -u -a - -b $standard | wc -l) # -f 0.1 -F 0.1: require both >10% overlap, when piranha perform best
	FP=$(echo "$predicted_truth" | cut -f1-3 | bedtools intersect -v -a - -b $standard | wc -l)
	# FN=$(echo "$predicted_false" | cut -f1-3 | bedtools intersect -u -f 0.1 -a - -b $standard | wc -l)
	FN=$(($RealTrueNum-$TP))
	TN=$(echo "$predicted_false" | cut -f1-3 | bedtools intersect -v -a - -b $standard | wc -l)
#higher -F better piranha, poorer clipper
#higher -f better clipper (a little effect)
	# RealFalseNum=$FP+$TN

	# calculate precision, recall, F1. export results.
	precision=$(calculate_rate $TP $FP)
	recall=$(calculate_rate $TP $FN)
	fpr=$(calculate_rate $FP $TN)
	f1=$(calculate_f1_score $precision $recall)

	echo -e "$method\t$replicate\t$counts\t$TP\t$FP\t$FN\t$TN\t$precision\t$recall\t$fpr\t$f1"
}
#\t$f1


# file I/O ------------------------------------------------------------------------------

# identify method,condition,replicate,mark
#/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/$dst/call_domain_withRepeats_all/domains_${method}/intersect/${id}.bed.count

method=$method
# method=$(echo $sample | cut -d "/" -f10) # usual
# method=$(basename $sample | sed s/".bed"/""/g) # test params

#condition=$(basename $sample | cut -d_ -f2)
# replicate=$(basename $sample | cut -d "." -f1)
replicate=$replicate
#mark=$(basename $sample | cut -d_ -f4 | sed 's/.bed//')

echo -e "method\treplicate\tcounts\tTP\tFP\tFN\tTN\tprecision\trecall\tfpr\tf1" # header, rm \tf1 

# file prep -----------------------------------------------------------------------------

# identify all counts within a sample.
declare -i max=$(cut -f5 $sample | sort -n | uniq | tail -n 1)
max=max+2

sample_signals=$(cut -f5 $sample | sort -n | uniq)
sample_signals=$(echo -e "${sample_signals}\n${max}")

# export vars and fx to child processes
export -f calculate_rate
export -f calculate_f1_score
export -f threshold
# echo $standardNum

declare -i standardNum=$((`wc -l ${standard} | cut -d " " -f1`))
# echo "truth num: $standardNum"

declare -i RealTrueNum=$standardNum # export standardNum=1000 before run this script to treat it as global var
export RealTrueNum

if [ $method != "" ]; then export method; else echo "method $method is undefined" && exit 1; fi
# if [ $condition != "" ]; then export condition; else echo "condition $condition is undefined" && exit 1; fi
if [ $replicate != "" ]; then export replicate; else echo "replicate $replicate is undefined" && exit 1; fi
# if [ $mark != "" ]; then export mark; else echo "mark $mark is undefined" && exit 1; fi
if [ $sample != "" ]; then export sample; else echo "sample $sample is undefined" && exit 1; fi
if [ $standard != "" ]; then export standard; else echo "standard $standard is undefined" && exit 1; fi



parallel --keep-order --trim lr -j $threads threshold ::: $sample_signals # process peaks in parallel!
#--bibtex: silence ref