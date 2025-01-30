## re-run Star snakemake --until mapping_gn_Star sort_gbam with raw fq (smp)
#get Star bam (mk standard true)
## re-run call_peak snakemake --until mapping_tx with raw fq (smp)
#get Bowtie2 tx bam (mk tx true?)

## trim fq
#https://www.biostars.org/p/188878/#240543
#source activate cfpeak
export PATH=/BioII/lulab_b/baopengfei/shared_utils:/BioII/lulab_b/baopengfei/shared_scripts:/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/scripts:/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/bin:$PATH

# for dst in {GSE94533_NCpool,GSE123972_NCpool,GSE110381_NCpool,GSE94582_NCpool,Phospho-RNA-seq_NCpool,AGO2_IP_NCpool};
# do
# echo $dst
# if [ $dst = "AGO2_IP_NCpool" ]; then
# dedup="_dedup"
# else
# dedup="_all"
# fi
# cat config/GSE71008_NCpool2.yaml | sed s/"GSE71008_NCpool2"/"${dst}2"/g > config/${dst}2.yaml
# mkdir -p ../output/${dst}2/trimmed
# mkdir -p logs/cluster/${dst}2
# cp -rp data/$dst  data/${dst}2
# for smp in `cat data/${dst}2/sample_ids.txt`
# do
# echo $smp
# ln -s ../../$dst/trimmed/${smp}.fastq.gz ../output/${dst}2/trimmed/${smp}.fastq.gz
# done
# cp -rp ../output/$dst/tbam ../output/${dst}2/tbam
# mkdir -p ../output/${dst}2/call_peak${dedup}/
# cp -rp ../output/${dst}/call_peak${dedup}/tb* ../output/${dst}2/call_peak${dedup}/
# cp -rp ../output/${dst}/call_peak${dedup}/*_by_sam* ../output/${dst}2/call_peak${dedup}/
# mv ../output/${dst}2/call_peak${dedup}/expeak_by_sample ../output/${dst}2/call_peak${dedup}/cfpeak_by_sample
# mv ../output/${dst}2/call_peak${dedup}/expeakCNN_by_sample ../output/${dst}2/call_peak${dedup}/cfpeakCNN_by_sample
# done

# for dst in TCGA-COAD_small_NC_NCpool CNP0003091_urine_NCpool GSE129255 GSE112343_NCpool GSE56686 snc_pandora_hsa_HeLa_NCpool encode_small_colon_NCpool
# for dst in {GSE71008_NCpool2,GSE94533_NCpool2,GSE123972_NCpool2,GSE110381_NCpool2,GSE94582_NCpool2,Phospho-RNA-seq_NCpool2,AGO2_IP_NCpool2};
for dst in {WSQ_SMARTer_NEB_NCpool2,i-pico_NCpool2}; # ,ipico,GSE278414,FTC_long
do
echo $dst
# if [ $dst = "AGO2_IP_NCpool2" ]; then
if [[ $dst =~ ^(AGO2_IP_NCpool2|FTC_long|ipico|GSE278414)$ ]]; then

dedup="_dedup"
else
dedup="_all"
fi
echo $dedup
outDir=/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/${dst}
CORES=10
cp /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/$dst/sample_ids.txt /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/$dst/sample_ids_raw.txt
for SMP in `cat /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/$dst/sample_ids_raw.txt` #  | grep TCGA| grep 11A
#TCGA-4T-AA8H-01A-11H-A41D-13_mirna_gdc_realn TCGA-A6-2671-01A-01T-1409-13_mirna_gdc_realn
do
    echo $SMP
    FQ=$outDir/trimmed/${SMP}.fastq.gz
    for N in 15
    do
        # N=10
        gzip -d -c $FQ | seqkit subseq -j ${CORES} -r 1:${N} | gzip -c > \
            $outDir/trimmed/${SMP}_${N}.fastq.gz
    done
done
done

#update smp ids (TGIRT-seq not deduped, UMI is tricky in EM)
## re-run StarEM snakemake --until mapping_gn_StarEM call_peaks_blockbuster with trim fq (smp_15)
#get StarEM bam (mk standard false)
## re-run call_peak snakemake --until merge_tbam_EM_19 call_peaks_piranha with trim fq (smp_15)
#get Bowtie2 tx bam (mk tx true?)


## prep bam
source activate blockbuster
mkdir test
CORES=10
# for dst in CNP0003091_urine_NCpool GSE129255 GSE112343_NCpool GSE56686 snc_pandora_hsa_HeLa_NCpool encode_small_colon_NCpool
for dst in {GSE71008_NCpool2,GSE94533_NCpool2,GSE123972_NCpool2,GSE110381_NCpool2,GSE94582_NCpool2,Phospho-RNA-seq_NCpool2}; # AGO2_IP_NCpool2
do
    echo $dst

    if [ $dst = "AGO2_IP_NCpool2" ]; then
        dedup="_dedup"
        bam_dedup_dir="bam-deduped" # output_dir+'/tbam/{wildcards.sample_id}/bam-deduped/{wildcards.rna_type}.bam'  
    else
        dedup="_all"
        bam_dedup_dir="bam"
    fi


    outDir=/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/${dst}
    for smp in `cat /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/${dst}/sample_ids_raw.txt` #  | grep TCGA| grep 11A
    #NCpool
    do
        echo $smp

        ## prep gn standard true (from bam all uniq)
        echo "gn true"
        # inBam="${outDir}/gbamStar_sorted/$smp/bam/genome.bam"
        # outBam="${outDir}/gbamStar_sorted/$smp/bam/genome_uniq.bam"
        inBam="${outDir}/gbamStar_sorted/$smp/genome.bam"
        outBam="${outDir}/gbamStar_sorted/$smp/genome_uniq.bam"
        #keep only those reads which appear exactly once in the SAM file, based on their read names.
        (samtools view -H ${inBam}; samtools view ${inBam} | awk 'BEGIN {FS="\t"; OFS="\t"} 
            /^@/ {print; next} 
            {readname=$1; count[readname]++; if (count[readname] == 1) { reads[readname] = $0 } else { reads[readname] = "" } } 
            END {for (readname in reads) if (reads[readname] != "") print reads[readname]}' ) | samtools view -b -o ${outBam}
        inBam=$outBam
        bamBed="test/${dst}_gbamStarUniq_${smp}.bed"
        bedtools bamtobed -i ${inBam} > ${bamBed}
        cut -f 4 $bamBed | sort | uniq > test/${dst}_sorted_read_names_${smp}.txt # blockBed has pseudo read name due to dedup
        blockBed="test/${dst}_gbamStarUniq_block_${smp}.bed" 
        python scripts/bamBed2blockbusterBed.py -b $bamBed -o $blockBed # standard true: no need to merge. # too many peak cov<3, rm these or also call peak ?
        blockPeak="test/${dst}_gbamStarUniq_block_${smp}.peak"
        blockbuster.x -format 1 -minBlockHeight 3 -minClusterHeight 3 -print 1 -tagFilter 3 <(cat ${blockBed} | LC_COLLATE=C sort -k1,1 -k2,2n) > tmp.${dst}.${smp}.bed
        cat tmp.${dst}.${smp}.bed | \
            grep -E -v '>' | \
            awk '{{ print $2 "\t" $3 "\t" $4 "\t" $1 "\t" $6 "\t" $5}}' | \
            sed s/'>'/''/g | \
            LC_COLLATE=C sort -k1,1 -k2,2n \
            > ${blockPeak}
        bedtools merge -s -c 4,5,6 -o first,mean,distinct -i ${blockPeak} > ${blockPeak}.merged 
        

        ## prep gn standard false (from EM-bam all multi)
        echo "gn false"
        inBam="${outDir}/gbamStarEM/${smp}_15/${bam_dedup_dir}/genome.bam"
        filterBam="${outDir}/gbamStarEM/${smp}_15/${bam_dedup_dir}/genome_filter.bam"
        (samtools view -H ${inBam}; samtools view ${inBam} | grep -F -f test/${dst}_sorted_read_names_${smp}.txt) | samtools view -b -o ${filterBam}
        outBam="${outDir}/gbamStarEM/${smp}_15/bam/genome_filter_multi.bam"
        #keep only those reads which appear > 1 in the SAM file, based on their read names. (not need rm NH:1, since will rm those overlap with uniq )
        (samtools view -H ${filterBam}; samtools view ${filterBam} | awk 'BEGIN {FS="\t"; OFS="\t"} 
            /^@/ {print; next} 
            {readname=$1; count[readname]++; if (count[readname] == 1) { reads[readname] = "" } else { reads[readname] = $0 } } 
            END {for (readname in reads) if (reads[readname] != "") print reads[readname]}' ) | samtools view -b -o ${outBam}
        inBam=$outBam
        bamBed="test/${dst}_gbamStarEMmulti_${smp}_15_filter.bed"
        bedtools bamtobed -i ${inBam} > ${bamBed} # should all belong to sorted_read_names.txt
        blockBed="test/${dst}_gbamStarEMmulti_block_${smp}_15_filter.bed" 
        python scripts/bamBed2blockbusterBed.py -b $bamBed -o $blockBed 
        bedtools intersect -split -s -v -wa -a $blockBed -b test/${dst}_gbamStarUniq_block_${smp}.bed > test/${dst}_gbamStarEMmulti_block_${smp}_15_filter_false.bed # standard false: no need to merge. # too many peak cov<3, rm these ?
        blockPeak="test/${dst}_gbamStarEMmulti_block_${smp}_15_filter.peak"
        blockbuster.x -format 1 -minBlockHeight 3 -minClusterHeight 3 -print 1 -tagFilter 3 <(cat ${blockBed} | LC_COLLATE=C sort -k1,1 -k2,2n) > tmp.${dst}.${smp}.bed
        cat tmp.${dst}.${smp}.bed | \
            grep -E -v '>' | \
            awk '{{ print $2 "\t" $3 "\t" $4 "\t" $1 "\t" $6 "\t" $5}}' | \
            sed s/'>'/''/g | \
            LC_COLLATE=C sort -k1,1 -k2,2n \
            > ${blockPeak}
        bedtools intersect -split -s -v -wa -a $blockPeak -b test/${dst}_gbamStarUniq_block_${smp}.peak > test/${dst}_gbamStarEMmulti_block_${smp}_15_filter_false.peak # standard false: no need to merge. # too many peak cov<3, rm these ?
        bedtools merge -s -c 4,5,6 -o first,mean,distinct -i test/${dst}_gbamStarEMmulti_block_${smp}_15_filter_false.peak > test/${dst}_gbamStarEMmulti_block_${smp}_15_filter_false.peak.merged 

        samtools fastq -@ ${CORES} ${filterBam} | gzip -c > ${outDir}/trimmed/${smp}_15_filter.fastq.gz
        # echo "${smp}_15_filter" >> data/$dst/sample_ids.txt


        ## prep tx standard true (from bam all uniq) To-check: better use bowtie2_def like star_def?
        echo "tx true"
        inBam="${outDir}/tbam/${smp}/${bam_dedup_dir}/merge19_sort.bam"
        outBam="${outDir}/tbam/${smp}/${bam_dedup_dir}/merge19_sort_uniq.bam"
        #keep only those reads which appear exactly once in the SAM file, based on their read names.
        (samtools view -H ${inBam}; samtools view ${inBam} | awk 'BEGIN {FS="\t"; OFS="\t"} 
            /^@/ {print; next} 
            {readname=$1; count[readname]++; if (count[readname] == 1) { reads[readname] = $0 } else { reads[readname] = "" } } 
            END {for (readname in reads) if (reads[readname] != "") print reads[readname]}' ) | samtools view -b -o ${outBam}
        inBam=$outBam
        bamBed="test/${dst}_tbamBowtie2Uniq_${smp}.bed"
        bedtools bamtobed -i ${inBam} > ${bamBed}
        cut -f 4 $bamBed | sort | uniq > test/${dst}_tx_sorted_read_names_${smp}.txt # blockBed has pseudo read name due to dedup
        blockBed="test/${dst}_tbamBowtie2Uniq_block_${smp}.bed" 
        python scripts/bamBed2blockbusterBed.py -b $bamBed -o $blockBed # standard true: no need to merge. # too many peak cov<3, rm these or also call peak ?
        blockPeak="test/${dst}_tbamBowtie2Uniq_block_${smp}.peak"
        blockbuster.x -format 1 -minBlockHeight 3 -minClusterHeight 3 -print 1 -tagFilter 3 <(cat ${blockBed} | LC_COLLATE=C sort -k1,1 -k2,2n) > tmp.${dst}.${smp}.bed
        cat tmp.${dst}.${smp}.bed | \
            grep -E -v '>' | \
            awk '{{ print $2 "\t" $3 "\t" $4 "\t" $1 "\t" $6 "\t" $5}}' | \
            sed s/'>'/''/g | \
            LC_COLLATE=C sort -k1,1 -k2,2n \
            > ${blockPeak}
        bedtools merge -s -c 4,5,6 -o first,mean,distinct -i ${blockPeak} > ${blockPeak}.merged 
        # tbed2gbed /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/bed/long_DNA-long_RNA-tRNA-pri_miRNA-piRNA-rRNA_newTxID.bed \
        # test/${dst}_tbamBowtie2Uniq_block_${smp}.peak.merged test/${dst}_tbamBowtie2Uniq_block_${smp}.peak.merged.gn.bed

        ## prep tx standard false (from EM-bam all multi)
        echo "gn false"
        inBam="${outDir}/tbam/${smp}_15/${bam_dedup_dir}-EM/merge19_sort/merged.sorted.bam"
        filterBam="${outDir}/tbam/${smp}_15/${bam_dedup_dir}-EM/merge19_sort/merged.sorted.filter.bam"
        (samtools view -H ${inBam}; samtools view ${inBam} | grep -F -f test/${dst}_tx_sorted_read_names_${smp}.txt) | samtools view -b -o ${filterBam}
        outBam="${outDir}/tbam/${smp}_15/${bam_dedup_dir}-EM/merge19_sort/genome_filter_multi.bam"
        #keep only those reads which appear > 1 in the SAM file, based on their read names. (not need rm NH:1, since will rm those overlap with uniq )
        (samtools view -H ${filterBam}; samtools view ${filterBam} | awk 'BEGIN {FS="\t"; OFS="\t"} 
            /^@/ {print; next} 
            {readname=$1; count[readname]++; if (count[readname] == 1) { reads[readname] = "" } else { reads[readname] = $0 } } 
            END {for (readname in reads) if (reads[readname] != "") print reads[readname]}' ) | samtools view -b -o ${outBam}
        inBam=$outBam
        bamBed="test/${dst}_tbamBowtie2EMmulti_${smp}_15_filter.bed"
        bedtools bamtobed -i ${inBam} > ${bamBed} # should all belong to sorted_read_names.txt
        blockBed="test/${dst}_tbamBowtie2EMmulti_block_${smp}_15_filter.bed" 
        python scripts/bamBed2blockbusterBed.py -b $bamBed -o $blockBed 
        bedtools intersect -split -s -v -wa -a $blockBed -b test/${dst}_tbamBowtie2Uniq_block_${smp}.bed > test/${dst}_tbamBowtie2EMmulti_block_${smp}_15_filter_false.bed # standard false: no need to merge.
        blockPeak="test/${dst}_tbamBowtie2EMmulti_block_${smp}_15_filter.peak"
        blockbuster.x -format 1 -minBlockHeight 3 -minClusterHeight 3 -print 1 -tagFilter 3 <(cat ${blockBed} | LC_COLLATE=C sort -k1,1 -k2,2n) > tmp.${dst}.${smp}.bed
        cat tmp.${dst}.${smp}.bed | \
            grep -E -v '>' | \
            awk '{{ print $2 "\t" $3 "\t" $4 "\t" $1 "\t" $6 "\t" $5}}' | \
            sed s/'>'/''/g | \
            LC_COLLATE=C sort -k1,1 -k2,2n \
            > ${blockPeak}
        bedtools intersect -split -s -v -wa -a $blockPeak -b test/${dst}_tbamBowtie2Uniq_block_${smp}.peak > test/${dst}_tbamBowtie2EMmulti_block_${smp}_15_filter_false.peak # standard false: no need to merge. # too many peak cov<3, rm these ?
        bedtools merge -s -c 4,5,6 -o first,mean,distinct -i test/${dst}_tbamBowtie2EMmulti_block_${smp}_15_filter_false.peak > test/${dst}_tbamBowtie2EMmulti_block_${smp}_15_filter_false.peak.merged 
        # tbed2gbed /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/bed/long_DNA-long_RNA-tRNA-pri_miRNA-piRNA-rRNA_newTxID.bed \
        # test/${dst}_tbamBowtie2EMmulti_block_${smp}_15_filter_false.peak.merged test/${dst}_tbamBowtie2EMmulti_block_${smp}_15_filter_false.peak.merged.gn.bed
        
        samtools fastq -@ ${CORES} ${filterBam} | gzip -c > ${outDir}/trimmed/${smp}_15_filter_tx.fastq.gz
        # echo "${smp}_15_filter_tx" >> data/$dst/sample_ids.txt
    done
done

#update smp ids
## re-run StarEM snakemake --until call_peaks_blockbuster call_peaks_clam with trim filtered fq (all smp_15_filter) (MAP=tx: smp_15_filter_tx)
#get StarEM + CLAM + blockbuster

## re-run Star snakemake --until call_peaks_blockbuster call_peaks_clipper call_peaks_piranha with trim filtered fq (all smp_15_filter) (MAP=tx: smp_15_filter_tx)
#get Star + CLIPper + Piranha + blockbuster

## re-run call_peak snakemake --until call_peaks_cfpeak CNN call_peaks_piranha call_peaks_clipper call_peaks_clam with trim filtered fq (all smp_15_filter_tx) (MAP=tx: smp_15_filter_tx)
#get Bowtie2 + tx peaks


CORES=6
for MAP in gn tx
do 
echo $MAP

for cfpeakMAP in tx gn
do 
echo $cfpeakMAP

    # #MAP=gn: cfpeak* use tx standard, others use gn standard
    # #MAP=tx: cfpeak* use tx standard, others use tx standard
    # MAP="gn" # tx gn
    # cfpeakMAP="tx" # tx gn

    # for dst in {GSE71008_NCpool2,GSE94533_NCpool2,GSE123972_NCpool2,GSE110381_NCpool2,GSE94582_NCpool2,Phospho-RNA-seq_NCpool2,AGO2_IP_NCpool2}; # TCGA-COAD_small_NC_NCpool
    # for dst in CNP0003091_urine_NCpool GSE129255 GSE112343_NCpool GSE56686 snc_pandora_hsa_HeLa_NCpool encode_small_colon_NCpool
    for dst in {GSE71008_NCpool2,GSE94533_NCpool2,GSE123972_NCpool2,GSE110381_NCpool2,GSE94582_NCpool2,Phospho-RNA-seq_NCpool2,AGO2_IP_NCpool2,CNP0003091_urine_NCpool,GSE129255,GSE112343_NCpool,GSE56686,snc_pandora_hsa_HeLa_NCpool,encode_small_colon_NCpool,TCGA-COAD_small_NC_NCpool}; # TCGA-COAD_small_NC_NCpool
    do
        echo $dst
        outDir=/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/${dst}

        if [ $dst = "AGO2_IP_NCpool2" ]; then
        dedup="_dedup"
        else
        dedup="_all"
        fi
        echo $dedup

        for smp in `cat /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/${dst}/sample_ids_raw.txt`
        # for smp in NCpool TCGA-4T-AA8H-01A-11H-A41D-13_mirna_gdc_realn
        do
            echo $smp
            if [ $MAP = "tx" ];then
                dedup="${dedup}"
                echo $dedup
                parsePeak4benchmark.R -t clam -r T -i ${outDir}/call_peak${dedup}/clam_by_sample/b5_p005/${smp}_15_filter_tx.bed -o test/${dst}_CLAM_${smp}_${MAP}.bed 
                parsePeak4benchmark.R -t piranha -r T -i ${outDir}/call_peak${dedup}/piranha_by_sample/b5_p01/${smp}_15_filter_tx.bed -o test/${dst}_Piranha_${smp}_${MAP}.bed 
                parsePeak4benchmark.R -t clipper -r T -i ${outDir}/call_peak${dedup}/clipper_by_sample/b5_p05/${smp}_15_filter_tx.bed -o test/${dst}_CLIPper_${smp}_${MAP}.bed 

            elif [ $MAP = "gn" ];then
                dedup="_gbamStar"
                echo $dedup
                ## copy STAR-EM bed (win50)
                parsePeak4benchmark.R -t blockbuster -r T -i ${outDir}/call_peak_gbamStarEM/blockbuster_by_sample/min3/${smp}_15_filter_block.bed -o test/${dst}_StarEM_win50_${smp}_${MAP}.bed # EM bed form comparison
                parsePeak4benchmark.R -t clam -r T -i ${outDir}/call_peak_gbamStarEM/clam_by_sample/b5_p005/${smp}_15_filter.bed -o test/${dst}_CLAM_${smp}_${MAP}.bed 

                ## get STAR-EM bed (win50 + minCov100)
                #seem not needed, considered in AUROC of win50 at specific threshold

                ## copy STAR-EM bed (win2000, deprecated)
                # cp ${outDir}/call_peak_gbamStarEM/blockbuster_by_sample/min3/${smp}_15_filter_block.bed test/StarEM_win2000_${smp}.bed

                ## copy STAR bed
                parsePeak4benchmark.R -t piranha -r T -i ${outDir}/call_peak_gbamStar/piranha_by_sample/b5_p01/${smp}_15_filter.bed -o test/${dst}_Piranha_${smp}_${MAP}.bed 
                parsePeak4benchmark.R -t clipper -r T -i ${outDir}/call_peak_gbamStar/clipper_by_sample/b5_p05/${smp}_15_filter.bed -o test/${dst}_CLIPper_${smp}_${MAP}.bed 
                parsePeak4benchmark.R -t blockbuster -r T -i ${outDir}/call_peak_gbamStar/blockbuster_by_sample/min3/${smp}_15_filter_block.bed -o test/${dst}_Star_default_${smp}_${MAP}.bed 

                ## copy cfPeak bed 
                #(not fair ? since gnEM true bed tend to be similar to gn than tx)
                #re-run Bowtie2EM smk with updated sample (filtered bam)  (include call peak)
            fi
        done
    done


    # for dst in {GSE71008_NCpool2,GSE94533_NCpool2,GSE123972_NCpool2,GSE110381_NCpool2,GSE94582_NCpool2,Phospho-RNA-seq_NCpool2}; # AGO2_IP_NCpool2
    # for dst in CNP0003091_urine_NCpool GSE129255 GSE112343_NCpool GSE56686 snc_pandora_hsa_HeLa_NCpool encode_small_colon_NCpool
    # for dst in {GSE71008_NCpool2,GSE94533_NCpool2,GSE123972_NCpool2,GSE110381_NCpool2,GSE94582_NCpool2,Phospho-RNA-seq_NCpool2,AGO2_IP_NCpool2}; # TCGA-COAD_small_NC_NCpool
    for dst in {GSE71008_NCpool2,GSE94533_NCpool2,GSE123972_NCpool2,GSE110381_NCpool2,GSE94582_NCpool2,Phospho-RNA-seq_NCpool2,AGO2_IP_NCpool2,CNP0003091_urine_NCpool,GSE129255,GSE112343_NCpool,GSE56686,snc_pandora_hsa_HeLa_NCpool,encode_small_colon_NCpool,TCGA-COAD_small_NC_NCpool}; # TCGA-COAD_small_NC_NCpool
    do
        echo $dst
        outDir=/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/${dst}

        if [ $dst = "AGO2_IP_NCpool2" ]; then
        dedup="_dedup"
        else
        dedup="_all"
        fi

        for smp in `cat /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/$dst/sample_ids_raw.txt`
        # for smp in NCpool TCGA-4T-AA8H-01A-11H-A41D-13_mirna_gdc_realn
        do
            echo $smp
            if [ $cfpeakMAP = "tx" ];then
                echo "tx standard"
                dedup="${dedup}"
                ## opt2: get cfPeak bed ( use tx standard)
                parsePeak4benchmark.R -t cfpeak -r T -i ${outDir}/call_peak${dedup}/cfpeak_by_sample/b5_d50_p1/${smp}_15_filter_tx.bed -o test/${dst}_bowtie2_cfpeak_${smp}_${cfpeakMAP}.bed # same
                # parsePeak4benchmark.R -t cfpeakCNN -r T -i ${outDir}/call_peak${dedup}/cfpeakCNN_by_sample/b5_d50_p1/${smp}_15_filter_tx.bed -o test/${dst}_bowtie2_cfpeakCNN_${smp}.bed
            elif [ $cfpeakMAP = "gn" ];then
                echo "gn standard"
                ## opt1: get cfPeak bed ( all use gn standard, not fair since gnEM true bed tend to be similar to gn than tx, but reasonable, though undermines cfpeak)
                #cfPeak
                tbed2gbed /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/bed/long_DNA-long_RNA-tRNA-pri_miRNA-piRNA-rRNA_newTxID.bed \
                ${outDir}/call_peak${dedup}/cfpeak_by_sample/b5_d50_p1/${smp}_15_filter.bed ${outDir}/call_peak${dedup}/cfpeak_by_sample/b5_d50_p1/${smp}_15_filter_${cfpeakMAP}.bed
                parsePeak4benchmark.R -t cfpeak -r T -i ${outDir}/call_peak${dedup}/cfpeak_by_sample/b5_d50_p1/${smp}_15_filter_${cfpeakMAP}.bed -o test/${dst}_bowtie2_cfpeak_${smp}_${cfpeakMAP}.bed
                # #cfPeakCNN
                # tbed2gbed /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/bed/long_DNA-long_RNA-tRNA-pri_miRNA-piRNA-rRNA_newTxID.bed \
                # ${outDir}/call_peak${dedup}/cfpeakCNN_by_sample/b5_d50_p1/${smp}_15_filter.bed ${outDir}/call_peak${dedup}/cfpeakCNN_by_sample/b5_d50_p1/${smp}_15_filter_gn.bed
                # parsePeak4benchmark.R -t cfpeakCNN -r T -i ${outDir}/call_peak${dedup}/cfpeakCNN_by_sample/b5_d50_p1/${smp}_15_filter_gn.bed -o test/${dst}_bowtie2_cfpeakCNN_${smp}.bed
            fi
        done
    done



    ## test AUROC
    # for dst in {GSE71008_NCpool2,GSE94533_NCpool2,GSE123972_NCpool2,GSE110381_NCpool2,GSE94582_NCpool2,Phospho-RNA-seq_NCpool2}; # AGO2_IP_NCpool2
    # for dst in CNP0003091_urine_NCpool GSE129255 GSE112343_NCpool GSE56686 snc_pandora_hsa_HeLa_NCpool encode_small_colon_NCpool
    # for dst in {GSE71008_NCpool2,GSE94533_NCpool2,GSE123972_NCpool2,GSE110381_NCpool2,GSE94582_NCpool2,Phospho-RNA-seq_NCpool2,AGO2_IP_NCpool2}; # TCGA-COAD_small_NC_NCpool
    for dst in {GSE71008_NCpool2,GSE94533_NCpool2,GSE123972_NCpool2,GSE110381_NCpool2,GSE94582_NCpool2,Phospho-RNA-seq_NCpool2,AGO2_IP_NCpool2,CNP0003091_urine_NCpool,GSE129255,GSE112343_NCpool,GSE56686,snc_pandora_hsa_HeLa_NCpool,encode_small_colon_NCpool,TCGA-COAD_small_NC_NCpool}; # TCGA-COAD_small_NC_NCpool
    do
        echo $dst
        outDir=/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/${dst}

        if [ $dst = "AGO2_IP_NCpool2" ]; then
        dedup="_dedup"
        else
        dedup="_all"
        fi

        for smp in `cat /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/$dst/sample_ids_raw.txt`
        # for smp in NCpool TCGA-4T-AA8H-01A-11H-A41D-13_mirna_gdc_realn
        do
            echo $smp
            if [ $MAP = "tx" ];then
                arr=("bowtie2_cfpeak" "Piranha" "CLIPper" "CLAM") # bowtie2_cfpeakCNN
            elif [ $MAP = "gn" ];then
                arr=("StarEM_win50" "Star_default" "bowtie2_cfpeak" "Piranha" "CLIPper" "CLAM") # StarEM_win2000,Star_default,Star_k100, bowtie2_cfpeakCNN
            fi

            for BED in "${arr[@]}" 
            do
                echo "start $BED `date`"
                if [ $BED = "bowtie2_cfpeak" -o $BED = "bowtie2_cfpeakCNN" ];then
                    echo "cfpeak mode"
                    MAP2=$cfpeakMAP
                    cat test/${dst}_${BED}_${smp}_${cfpeakMAP}.bed | grep -v "," | cut -f1-6 > test/${dst}_${BED}_${smp}_${MAP2}.bed.tmp
                    if [ $cfpeakMAP = "tx" ];then
                        echo "cfpeak* use tx standard & tx coord"
                        # standard_bed="test/${dst}_gbamStarUniq_block_${smp}.peak.merged" # Pos: merge not change result
                        # standard_False_bed="test/${dst}_gbamStarEMmulti_block_${smp}_15_filter_false.peak.merged" # Neg: merge not change result
                        standard_bed="test/${dst}_tbamBowtie2Uniq_block_${smp}.peak.merged" # Pos: merge not change result
                        standard_False_bed="test/${dst}_tbamBowtie2EMmulti_block_${smp}_15_filter_false.peak.merged" # Neg: merge not change result
                    elif [ $cfpeakMAP = "gn" ];then
                        echo "cfpeak* use gn standard & gn coord"
                        standard_bed="test/${dst}_gbamStarUniq_block_${smp}.peak.merged"
                        standard_False_bed="test/${dst}_gbamStarEMmulti_block_${smp}_15_filter_false.peak.merged"
                    fi
                else
                    echo "other mode"
                    MAP2=$MAP
                    cat test/${dst}_${BED}_${smp}_${MAP}.bed | grep -v "," | cut -f1-6 > test/${dst}_${BED}_${smp}_${MAP2}.bed.tmp
                    if [ $MAP = "tx" ];then
                        echo "others use tx standard & tx coord"
                        standard_bed="test/${dst}_tbamBowtie2Uniq_block_${smp}.peak.merged" # Pos
                        standard_False_bed="test/${dst}_tbamBowtie2EMmulti_block_${smp}_15_filter_false.peak.merged" # Neg
                    elif [ $MAP = "gn" ];then
                        echo "others use gn standard & gn coord"
                        standard_bed="test/${dst}_gbamStarUniq_block_${smp}.peak.merged" # Pos
                        standard_False_bed="test/${dst}_gbamStarEMmulti_block_${smp}_15_filter_false.peak.merged" # Neg
                    fi
                fi
                
                cat test/${dst}_${BED}_${smp}_${MAP2}.bed.tmp | LC_COLLATE=C sort -k1,1 -k2,2n | bedtools merge -s -c 4,5,6 -o first,max,distinct -i stdin > test/${dst}_${BED}_${smp}.bed.tmp2 # max --> mean

                # ### opt1: rm not called standard for each method (normal but AUC~0.5, normal AUPR)
                # bedtools intersect -split -s -u -wa -a ${standard_bed} -b test/${dst}_${BED}_${smp}.bed.tmp2 > test/${dst}_gbamStarUniq_block_${smp}_${BED}_called.bed        
                # bedtools intersect -split -s -u -wa -a ${standard_False_bed} -b test/${dst}_${BED}_${smp}.bed.tmp2 > test/${dst}_gbamStarEMmulti_block_${smp}_15_filter_false_${BED}_called.bed
                # cat test/${dst}_gbamStarEMmulti_block_${smp}_15_filter_false_${BED}_called.bed test/${dst}_gbamStarUniq_block_${smp}_${BED}_called.bed > tmp.${dst}.${smp}.bed

                # bedtools intersect -split -s -u -wa -a test/${dst}_${BED}_${smp}.bed.tmp2 -b tmp.${dst}.${smp}.bed > test/${dst}_${BED}_${smp}_called.bed

                # TYPE="called"
                # echo "start $BED `date`"
                # standard_bed2="test/${dst}_gbamStarUniq_block_${smp}_${BED}_called.bed" # Pos
                # standard_False_bed2="test/${dst}_gbamStarEMmulti_block_${smp}_15_filter_false_${BED}_called.bed" # Neg
                # output=test/${dst}_${BED}_${smp}_called.bed # rank by 4th column
                # output_tbl="test/${dst}_${BED}_${smp}_${MAP2}_${TYPE}.txt" #
                # cal_confusion_mat2_EM.sh  -i $output  -s ${standard_bed2} -m ${BED} -r replicate -t ${CORES} -f ${standard_False_bed2} >  $output_tbl  

                # ### opt1.2: use all
                # TYPE="all3"
                # echo "start $BED `date`"
                # standard_bed2="$standard_bed" # Pos
                # standard_False_bed2="$standard_False_bed" # Neg
                # output=test/${dst}_${BED}_${smp}.bed.tmp2 # rank by 4th column
                # output_tbl="test/${dst}_${BED}_${smp}_${MAP2}_${TYPE}.txt" #
                # cal_confusion_mat2_EM.sh  -i $output  -s ${standard_bed2} -m ${BED} -r replicate -t ${CORES} -f ${standard_False_bed2} >  $output_tbl  

                ### opt1.3: use all
                TYPE="all4"
                echo "start $BED `date`"
                standard_bed2="$standard_bed" # Pos
                standard_False_bed2="$standard_False_bed" # Neg
                output=test/${dst}_${BED}_${smp}.bed.tmp2 # rank by 4th column
                output_tbl="test/${dst}_${BED}_${smp}_${MAP2}_${TYPE}.txt" #
                cal_confusion_mat2_EM2.sh  -i $output  -s ${standard_bed2} -m ${BED} -r replicate -t ${CORES} -f ${standard_False_bed2} >  $output_tbl  

                # # opt2: use all standard ( low AUC & AUPR, only compare F1 and precision due to redundant peaks)
                # TYPE="all"
                # output="test/${dst}_${BED}_${smp}.bed.tmp2" # rank by 4th column
                # output_tbl="test/${dst}_${BED}_${smp}_${MAP2}_${TYPE}.txt" # why not much changed
                # cal_confusion_mat2.sh  -i $output  -s ${standard_bed} -m ${BED} -r replicate -t ${CORES} >  $output_tbl  

                # # opt2 v3: use all standard ( low AUC & AUPR, only compare F1 and precision due to redundant peaks)
                # TYPE="all2"
                # output="test/${dst}_${BED}_${smp}.bed.tmp2" # rank by 4th column
                # output_tbl="test/${dst}_${BED}_${smp}_${MAP2}_${TYPE}.txt" # why not much changed
                # cal_confusion_mat3.sh  -i $output  -s ${standard_bed} -m ${BED} -r replicate -t ${CORES} >  $output_tbl  
            done
        done
    done
done
done
