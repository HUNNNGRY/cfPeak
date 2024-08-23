#!/bin/bash

#TODO:
# final test in local and Acluster
# need update ref to _submitA.sh

#eg. clusterA sub cmd:
# sbatch -p Acluster -N 1 -n 1 --job-name=test --output=test.outerr --error=test.outerr get_chipr_consensus_peak_submitA.sh \

#eg. local run cmd:
# bash get_chipr_consensus_peak.sh \
# lulab_oscc_tissue_diff \
# /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/lulab_oscc_tissue_diff/sample_table.txt \
# /data2/lulab1/bpf/projects/WCHSU-FTC/output/lulab_oscc_tissue_diff/call_peak_dedup/cfpeak_by_sample/b5_d50_p1


## fixed var
#maxJobs=200 # cluster only
chrSize="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/transcriptome_genome_sort_uniq_newTxID"
cov_num=1
codeDir="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/scripts"

## not fixed var
dst=$1     # lulab_oscc_tissue_diff
smpTab=$2  # /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/${dst}/sample_table.txt
peakDir=$3 # peakDir="${pre}/output/${dst}/call_peak${dedup}/cfpeak_by_sample/b5_d50_p1"

## half fixed var
tmpDir="tmp/$dst/"
outDir="test/$dst"
logDir="logs/cluster/$dst"
jac_dir="$outDir/jaccard"

# pre="/data2/lulab1/bpf/projects/WCHSU-FTC"
# #pre="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
# dst="lulab_oscc_tissue_diff"
# dedup="_dedup"
# smpTab="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/${dst}/sample_table.txt"
# #header: sample \t group \t ...
# #col1: smp id; col2: grp id; col3: site/organ/cell id

mkdir -p $tmpDir
mkdir -p $outDir/log

colNum=`head -n1 $smpTab | tr '\t' '\n' | wc -l`

if [ $colNum -gt 2 ];then 
echo "use 2nd, 3rd col as groups"
for CELL in `cat ${smpTab} | grep -v 'sample' | cut -f 3 | sort | uniq`
#28
do
	for GRP in `cat ${smpTab} | grep -v 'sample' | cut -f 2 | sort | uniq`
	do

		rm -f $tmpDir/${CELL}_${GRP}.txt
		rm -f $tmpDir/${CELL}_${GRP}.bed
		for SMP in `cat ${smpTab} | awk -v cell=${CELL} -v grp=${GRP} 'BEGIN{OFS="\t";FS="\t"} {if($3==cell && $2==grp) {print $1}}'`
		do
			echo ${SMP}
			bash $codeDir/cfpeakBed2chipr.sh $peakDir/${SMP}.bed
			echo $peakDir/${SMP}.bed.tmp >> $tmpDir/${CELL}_${GRP}.txt
			cat $peakDir/${SMP}.bed.tmp  >> $tmpDir/${CELL}_${GRP}.bed
		done
	done
done
elif [ $colNum == 2 ];then
echo "use 2nd col as groups"
        for GRP in `cat ${smpTab} | grep -v 'sample' | cut -f 2 | sort | uniq`
        do

                rm -f $tmpDir/${GRP}.txt
                rm -f $tmpDir/${GRP}.bed
                for SMP in `cat ${smpTab} | awk -v grp=${GRP} 'BEGIN{OFS="\t";FS="\t"} {if($2==grp) {print $1}}'`
                do
                        echo ${SMP}
                        bash $codeDir/cfpeakBed2chipr.sh $peakDir/${SMP}.bed
                        echo $peakDir/${SMP}.bed.tmp >> $tmpDir/${GRP}.txt
                        cat $peakDir/${SMP}.bed.tmp  >> $tmpDir/${GRP}.bed
                done
        done
fi

mkdir -p $logDir

if [ $colNum -gt 2 ];then
echo "use 2nd, 3rd col as groups"
for CELL in `cat ${smpTab} | grep -v 'sample' | cut -f 3 | sort | uniq`
do
echo $CELL
	for GRP in `cat ${smpTab} | grep -v 'sample' | cut -f 2 | sort | uniq`
	do
		echo $GRP
		
		if [ -s $tmpDir/${CELL}_${GRP}.txt ];
		then
			echo "exist"
		
			#run chipr
			source activate py37 
			chipr --size 10 --alpha 0.05 --rankmethod pvalue --minentries 1 --fragment \
				-i `cat $tmpDir/${CELL}_${GRP}.txt | tr '\n' ' ' ` \
				-o $outDir/${CELL}_${GRP}.bed > $outDir/log/${CELL}_${GRP}.log 2>&1
			bash $codeDir/chipr2cfpeakBed6.sh  $outDir/${CELL}_${GRP}.bed_optimal.bed  $outDir/${CELL}_${GRP}.bed_optimal_cfpeak.bed
			
			sleep 0.5
		fi
	done
done
elif [ $colNum == 2 ];then
echo "use 2nd col as groups"
        for GRP in `cat ${smpTab} | grep -v 'sample' | cut -f 2 | sort | uniq`
        do
                echo $GRP

                if [ -s $tmpDir/${GRP}.txt ];
                then
                        echo "exist"

                        #run chipr
                        source activate py37
                        chipr --size 10 --alpha 0.05 --rankmethod pvalue --minentries 1 --fragment \
                                -i `cat $tmpDir/${GRP}.txt | tr '\n' ' ' ` \
                                -o $outDir/${GRP}.bed > $outDir/log/${GRP}.log 2>&1
                        bash $codeDir/chipr2cfpeakBed6.sh  $outDir/${GRP}.bed_optimal.bed  $outDir/${GRP}.bed_optimal_cfpeak.bed
                        
                        sleep 0.5
                fi
        done
fi
#sbatch --wrap will calculate all `script` witin !!! make sure its execuable first

#calculate overlap
mkdir -p $jac_dir
awk -f $codeDir/combinations.awk <<< `ls $outDir/*_*.bed_optimal_cfpeak.bed | tr "\n" " "` > ${jac_dir}/cmb.txt
bash $codeDir/getBedJaccard.sh ${jac_dir}/cmb.txt  $jac_dir


#merge
mkdir -p $outDir/merge
bash $codeDir/mergeBed.sh "$outDir/*.bed_optimal_cfpeak.bed" $cov_num $chrSize $outDir/merge/cfpeak_smpFreq_${cov_num}.bed 
#op1: ln /data2/lulab1/bpf/projects/WCHSU-FTC/test/$dst/merge/cfpeak_smpFreq_1.bed /data2/lulab1/bpf/projects/WCHSU-FTC/output/GSE62809/call_peak_all/cfpeakCNN/b5_d50_p1.bed
#op2: bash $codeDir/mergeBed.sh "test/GSE254937/merge/cfpeak_smpFreq_1.bed test/GSE62809/merge/cfpeak_smpFreq_1.bed test/lulab_oscc_plasma_diff/merge/cfpeak_smpFreq_1.bed test/lulab_oscc_plasma_total_short_diff/merge/cfpeak_smpFreq_1.bed test/lulab_oscc_tissue_diff/merge/cfpeak_smpFreq_1.bed test/TCGA-HNSC_small_diff/merge/cfpeak_smpFreq_1.bed" 1 $chrSize test/OSCC_mergeDsts/cfpeak_smpFreq_1.bed
# ln test/OSCC_mergeDsts/cfpeak_smpFreq_1.bed /data2/lulab1/bpf/projects/WCHSU-FTC/output/GSE254937/call_peak_all/cfpeakCNN/b5_d50_p1.bed
# ln test/OSCC_mergeDsts/cfpeak_smpFreq_1.bed /data2/lulab1/bpf/projects/WCHSU-FTC/output/GSE62809/call_peak_all/cfpeakCNN/b5_d50_p1.bed
# ln test/OSCC_mergeDsts/cfpeak_smpFreq_1.bed /data2/lulab1/bpf/projects/WCHSU-FTC/output/lulab_oscc_plasma_diff/call_peak_dedup/cfpeakCNN/b5_d50_p1.bed
# ln test/OSCC_mergeDsts/cfpeak_smpFreq_1.bed /data2/lulab1/bpf/projects/WCHSU-FTC/output/lulab_oscc_tissue_diff/call_peak_dedup/cfpeakCNN/b5_d50_p1.bed
# ln test/OSCC_mergeDsts/cfpeak_smpFreq_1.bed /data2/lulab1/bpf/projects/WCHSU-FTC/output/lulab_oscc_plasma_total_short_diff/call_peak_dedup/cfpeakCNN/b5_d50_p1.bed
# ln test/OSCC_mergeDsts/cfpeak_smpFreq_1.bed /data2/lulab1/bpf/projects/WCHSU-FTC/output/TCGA-HNSC_small_diff/call_peak_all/cfpeakCNN/b5_d50_p1.bed
