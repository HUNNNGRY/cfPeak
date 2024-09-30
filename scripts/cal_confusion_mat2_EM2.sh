
#!/bin/bash

# newest for EM testing 

while getopts "i:s:f:m:r:t:" op
do
	case "$op" in
		i)  sample="$OPTARG";;
		s)  standard="$OPTARG";;
        f)  standard_false="$OPTARG";;
		m)  method="$OPTARG";;
		r)  replicate="$OPTARG";;
		t)  threads="$OPTARG";;
		\?) exit 1;;
	esac
done
# echo $sample
# echo $standard
# echo $threads

# variables and functions ---------------------------------------------------------------

calculate_rate () {
	#func(A,B): A/(A+B)
	awk -v term_1=$1 -v term_2=$2 'BEGIN {print (term_1 + 0.0001) / (term_1 + term_2 + 0.0001)}' # + 0.00001 avoid piranha fdr division by zero (cut -b1-5 : truncate first 5 number)
}

calculate_rate2 () {
	#func(A,B): A/B
	awk -v term_1=$1 -v term_2=$2 'BEGIN {print (term_1 + 0.0001) / (term_2 + 0.0001)}' # + 0.00001 avoid piranha fdr division by zero (cut -b1-5 : truncate first 5 number)
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
	declare -i TP=$(echo "$predicted_truth" | cut -f1-6 | bedtools intersect -split -s -u -a - -b $standard | wc -l) # -f 0.1 -F 0.1: require both >10% overlap, when piranha perform best
	tmp1=$(echo "$predicted_truth" | cut -f1-6 | bedtools intersect -split -s -v -a - -b $standard)
	declare -i FP=$(echo "$tmp1" | wc -l)
	declare -i FP1=$(echo "$tmp1" | bedtools intersect -split -s -v -a - -b $standard_false | wc -l)
	declare -i FP2=$(echo "$tmp1" | bedtools intersect -split -s -u -a - -b $standard_false | wc -l)
	declare -i TN=$(echo "$predicted_false" | cut -f1-6 | bedtools intersect -split -s -u -a - -b $standard_false | wc -l) # -f 0.1 -F 0.1: require both >10% overlap, when piranha perform best
	tmp2=$(echo "$predicted_false" | cut -f1-6 | bedtools intersect -split -s -v -a - -b $standard_false)
	declare -i FN=$(echo "$tmp2" | wc -l)
	declare -i FN1=$(echo "$tmp2" | bedtools intersect -split -s -v -a - -b $standard | wc -l)
	declare -i FN2=$(echo "$tmp2" | bedtools intersect -split -s -u -a - -b $standard | wc -l)

	declare -i RealTrueNum2=$(($RealTrueNum+$FN1))
	declare -i RealFalseNum2=$(($RealFalseNum+$FP1))
	declare -i totalN=$(($RealTrueNum+$RealFalseNum))
	declare -i totalN2=$(($RealTrueNum2+$RealFalseNum2))
	# declare -i FN=$(($totalN-$TP-$FP-$TN))
	# FN=$(echo "$predicted_false" | cut -f1-3 | bedtools intersect -u -f 0.1 -a - -b $standard | wc -l)
	# FN=$(($RealTrueNum-$TP))
	# TN=$(echo "$predicted_false" | cut -f1-3 | bedtools intersect -v -a - -b $standard | wc -l)
#higher -F better piranha, poorer clipper
#higher -f better clipper (a little effect)
	# RealFalseNum=$FP+$TN


	## calculate precision, recall, F1. export results.

	# precision=$(calculate_rate $TP $FP)
	# recall=$(calculate_rate $TP $FN)
	# fpr=$(calculate_rate $FP $TN)
	# f1=$(calculate_f1_score $precision $recall)
	precision=$(calculate_rate $TP $FP)
	recall=$(calculate_rate2 $TP $RealTrueNum2)
	fpr=$(calculate_rate2 $FP $RealFalseNum2)
	f1=$(calculate_f1_score $precision $recall)

	echo -e "$method\t$replicate\t$counts\t$RealTrueNum\t$RealTrueNum2\t$RealFalseNum\t$RealFalseNum2\t$TP\t$FP\t$FP1\t$FP2\t$FN\t$FN1\t$FN2\t$TN\t$precision\t$recall\t$fpr\t$f1"
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

echo -e "method\treplicate\tcounts\tRealTrueNum\tRealTrueNum2\tRealFalseNum\tRealFalseNum2\tTP\tFP\tFP1\tFP2\tFN\tFN1\tFN2\tTN\tprecision\trecall\tfpr\tf1" # header, rm \tf1 

# file prep -----------------------------------------------------------------------------

# identify all counts within a sample.
#declare -i max=$(cut -f5 $sample | sort -n | uniq | tail -n 1) 
max=$(cut -f5 $sample | sort -n | uniq | tail -n 1)
declare -i max=$( printf "%.0f" $max ) # float to integer
max=max+2

sample_signals=$(cut -f5 $sample | sort -n | uniq)
sample_signals=$(echo -e "${sample_signals}\n${max}")

# export vars and fx to child processes
export -f calculate_rate
export -f calculate_rate2
export -f calculate_f1_score
export -f threshold
# echo $standardNum

declare -i RealTrueNum=$((`wc -l ${standard} | cut -d " " -f1`))
export RealTrueNum
declare -i RealFalseNum=$((`wc -l ${standard_false} | cut -d " " -f1`))
export RealFalseNum

# declare -i RealTrueNum=$standardNum # export standardNum=1000 before run this script to treat it as global var
# export RealTrueNum

if [ $method != "" ]; then export method; else echo "method $method is undefined" && exit 1; fi
# if [ $condition != "" ]; then export condition; else echo "condition $condition is undefined" && exit 1; fi
if [ $replicate != "" ]; then export replicate; else echo "replicate $replicate is undefined" && exit 1; fi
# if [ $mark != "" ]; then export mark; else echo "mark $mark is undefined" && exit 1; fi
if [ $sample != "" ]; then export sample; else echo "sample $sample is undefined" && exit 1; fi
if [ $standard != "" ]; then export standard; else echo "standard $standard is undefined" && exit 1; fi
if [ $standard_false != "" ]; then export standard_false; else echo "standard_false $standard_false is undefined" && exit 1; fi



parallel --keep-order --trim lr -j $threads threshold ::: $sample_signals # process peaks in parallel!
#--bibtex: silence ref