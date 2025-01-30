#!/bin/bash
#SBATCH -J MFE
#SBATCH -p Acluster
#SBATCH -n 18
#SBATCH --cpus-per-task 1
#SBATCH --output=%j.out
#SBATCH --error=%j.err

	source activate py37

	
	dst="TCGA_small16" 
	dedup="_all" 
	pre="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
	outpre="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC" 
	cores=16
	method="cfpeakCNN/b5_d50_p1"
	chrSize="${pre}/exSeek-dev/genome/hg38/chrom_sizes/transcriptome_genome_sort_uniq_newTxID"
	tx_fa_ref="$pre/exSeek-dev/genome/hg38/fasta_newTxID/combine19.fa"
	gn_fa_ref="/lulabdata/baopengfei/shared_reference/hg38/genome.fa"
	
cd $pre

for tx in 8DNA_gn 11RNA
do
	i="${outpre}/output/$dst/call_peak${dedup}/${method}_${tx}.bed"
	cat ${i} | awk ' $3-$2 <= 200 {print $0}' | grep -v "," | cut -f1-6 > ${i}.tmp  # #filter too large peak, too slow; bed12 fail to slop block, lead to wrong ext.fa
	bedtools slop -pct -s -l 0.99 -r 0.99  -g ${chrSize} -i ${i}.tmp > ${i}.tmp2
	
	if [[ $tx == "11RNA" ]]
	then bedtools getfasta -split -nameOnly -s -fi $tx_fa_ref -bed ${i}.tmp2 > ${i}.exp.fa
	elif [[ $tx == "8DNA_gn" ]]
	then bedtools getfasta -split -nameOnly -s -fi $gn_fa_ref -bed ${i}.tmp2 > ${i}.exp.fa
	fi
	
	python3 ${pre}/../dsRFinder/scripts/rnafold_dinushuffle_parallel.py  ${i}.exp.fa 50 1234 ${i}.exp.fa.csv ${cores}
	rm ${i}.exp.fa_perm
done
