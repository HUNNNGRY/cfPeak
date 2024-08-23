#!/bin/bash
#SBATCH -J test_job
#SBATCH -p Acluster
#SBATCH -n 1
#SBATCH --cpus-per-task 2
#SBATCH --output=%j.out
#SBATCH --error=%j.err

#eg. clusterA sub cmd:
# sbatch -p Acluster -N 1 -n 1 --job-name=test --output=test.outerr --error=test.outerr get_chipr_consensus_peak_submitA.sh 

# arr=("TCGA-BRCA_small" "TCGA-COAD_small" "TCGA-HNSC_small" "TCGA-LIHC_small" "TCGA-LUAD_small" "TCGA-PRAD_small" "TCGA-STAD_small" "TCGA-THCA_small" "TCGA-LUSC_small" "TCGA-ESCA_small" "TCGA-BLCA_small" "TCGA-CHOL_small" "TCGA-UCEC_small" "TCGA-KIRC_small" "TCGA-KIRP_small" "TCGA-KICH_small")
# for dst in "${arr[@]}"
# do
# 	echo $dst
# 	cat scripts/get_chipr_consensus_peak_submitA.sh | sed "s/TCGA-COAD_small/${dst}/g" > ./get_chipr_consensus_peak_submitA_${dst}.sh
# 	chmod 755 ./get_chipr_consensus_peak_submitA_${dst}.sh
# 	sbatch -p Acluster -N 1 -n 1 --job-name=${dst} --output=${dst}.outerr --error=${dst}.outerr ./get_chipr_consensus_peak_submitA_${dst}.sh
# done
## sleep 30
## rm ./get_chipr_consensus_peak_submitA_*.sh


#eg. local run cmd:
# bash get_chipr_consensus_peak.sh ...

# newest version

## fixed var
export maxJobs=200 # cluster only

# export pre="/data2/lulab1/bpf/projects/WCHSU-FTC"
export pre="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
export dst="TCGA-COAD_small" #"lulab_oscc_tissue_diff"
export dedup="_all" #"_dedup"
export smpTab="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/${dst}/sample_table.txt"
#header: sample \t group \t ...
#col1: smp id; col2: grp id; col3: site/organ/cell id

#smpTab=$1
export chrSize="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/transcriptome_genome_sort_uniq_newTxID"
export cov_cell=1  # num of tissue/cell/group freq.
export cov_threshold=0.1 # num of sample freq. each tissue/cell/group

export peakDir="${pre}/output/${dst}/call_peak${dedup}/cfpeak_by_sample/b5_d50_p1" # TCGA: cfpeak_by_sample; plasma: cfpeakCNN_by_sample
export codeDir="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/scripts"
export tmpDir="tmp/$dst/"
export outDir="test/$dst"
export logDir="logs/cluster/$dst"
export jac_dir="$outDir/jaccard"

mkdir -p $tmpDir
mkdir -p $outDir/log

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
			bash $codeDir/cfpeakBed2chipr.sh $peakDir/${SMP}.bed # not necessary if not run chipr
			echo $peakDir/${SMP}.bed.tmp >> $tmpDir/${CELL}_${GRP}.txt
			cat $peakDir/${SMP}.bed.tmp  >> $tmpDir/${CELL}_${GRP}.bed
		done
	done
done


#Acluster
mkdir -p $logDir
rm $tmpDir/${dst}_job_ids.txt
for CELL in `cat ${smpTab} | grep -v 'sample' | cut -f 3 | sort | uniq`
do
echo $CELL
	export CELL="$CELL"
	for GRP in `cat ${smpTab} | grep -v 'sample' | cut -f 2 | sort | uniq`
	do
		echo $GRP
		export GRP="$GRP"
		if [ -s $tmpDir/${CELL}_${GRP}.txt ];
		then
			# echo "exist"
		
			while [ $(squeue -u baopengfei | wc -l) -gt $maxJobs ]; do sleep 30s; done # prevent too many jobs
			
			declare -i cov_num=$(cat $smpTab| grep -w "$CELL" | grep -w "$GRP" | wc -l | cut -f 1 -d " ")
			# declare -i cov_num=$(echo "scale=0; (${cov_num}*${cov_threshold}+0.5)/1" | bc) # round digit (need install bc)
			declare -i cov_num=$(awk "BEGIN {print int(${cov_num} * ${cov_threshold} + 0.5)}")

			if [ ${cov_num} -lt 1 ]; then
				cov_num=1
			fi
			export cov_num=$cov_num
			echo $cov_num
			
			slurm_log="${logDir}/mergeBed_${CELL}_${GRP}.outerr"
			job_id=$(sbatch --parsable -p Acluster -N 1 -n 1 --exclude node01,node02,node03 \
			--job-name=${CELL}.${GRP} --output=${slurm_log} --error=${slurm_log} \
			--wrap 'cd /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev
			bash $codeDir/mergeBed.sh $tmpDir/${CELL}_${GRP}.txt $cov_num $chrSize $outDir/${CELL}_${GRP}.bed_optimal.bed
			cut -f1-6 $outDir/${CELL}_${GRP}.bed_optimal.bed > $outDir/${CELL}_${GRP}.bed_optimal_cfpeak.bed
			')
			echo $job_id
			echo $job_id >> $tmpDir/${dst}_job_ids.txt
			sleep 2
		fi
	done
done
			#run chipr (tend to underestimate pvalue-power for shorter RNA like miR)
			# source activate py37 
			# chipr --size 10 --alpha 0.05 --rankmethod pvalue --minentries 1 --fragment \
			# 	-i `cat $tmpDir/${CELL}_${GRP}.txt | tr '\n' ' ' ` \
			# 	-o $outDir/${CELL}_${GRP}.bed > $outDir/log/${CELL}_${GRP}.log 2>&1
			# bash $codeDir/chipr2cfpeakBed6.sh  $outDir/${CELL}_${GRP}.bed_optimal.bed  $outDir/${CELL}_${GRP}.bed_optimal_cfpeak.bed

#sbatch --wrap will calculate all `script` witin !!! make sure its execuable first

sleep 120
declare -i jobNum=$(squeue -u baopengfei | grep -f $tmpDir/${dst}_job_ids.txt | wc -l | cut -f 1 -d " ")
if [ $jobNum -gt 0 ]; then
	echo "Waiting for $jobNum jobs to finish..."
	while [ $jobNum -gt 0 ]; do
		sleep 30s
		declare -i jobNum=$(squeue -u baopengfei | grep -f $tmpDir/${dst}_job_ids.txt | wc -l | cut -f 1 -d " ")
	done
fi

#calculate overlap
mkdir -p $jac_dir
awk -f $codeDir/combinations.awk <<< `ls $outDir/*_*.bed_optimal_cfpeak.bed | tr "\n" " "` > ${jac_dir}/cmb.txt
bash $codeDir/getBedJaccard2.sh ${jac_dir}/cmb.txt  $jac_dir


#merge
mkdir -p $outDir/merge
bash $codeDir/mergeBed.sh "$outDir/*.bed_optimal_cfpeak.bed" $cov_cell $chrSize $outDir/merge/cfpeak_smpFreq_${cov_cell}.bed 
#ln /data2/lulab1/bpf/projects/WCHSU-FTC/test/$dst/merge/cfpeak_smpFreq_1.bed /data2/lulab1/bpf/projects/WCHSU-FTC/output/GSE62809/call_peak_all/cfpeakCNN/b5_d50_p1.bed

