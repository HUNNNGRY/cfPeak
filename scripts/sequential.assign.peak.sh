#!/bin/bash

# last 2406 by b.p.f@qq.com
# Usage: bash sequential.assign.peak.sh <reads.bed> <outputdir> <beddir> <regions>
# support bed6 (not bed12 ?)

bed=$1
outputdir=$2
beddir=$3
#regions=`cat ${beddir}/priority.txt`
regions=`cat $4`
chrSize=/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/transcriptome_genome_sort_uniq 
fracA=0.5 #Default is 1E-9 (i.e., 1bp).

mkdir -p $outputdir
#watch out for bedtools coverage -s/-S -split
#v2.30.0
bedtools="/BioII/lulab_b/baopengfei/biosoft/bedtools"

#/BioII/lulab_b/jinyunfan/software/bedops/bin/sort-bed ${bed} > ${outputdir}/sorted.reads.bed
cat ${bed} | bedtools sort | LC_COLLATE=C sort -k1,1 -k2,2n > $outputdir/sorted.reads.bed

regionalCounts=""
for region in ${regions}
do
echo ${region}
regionalCounts="${regionalCounts} ${outputdir}/${region}.count.txt"
if [ ! -s ${outputdir}/${region}.count.txt ];then
echo "Counting ${region} ..." #>> ${outputdir}/log.txt
## PE bed (2 rows each read pair) pico reverse version
#cat ${outputdir}/sorted.reads.bed \
#| bedtools coverage -counts -s -g /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/transcriptome_genome_sort -a - -b ${beddir}/${region}.bed \
#| sort -k 4 \
#|awk 'BEGIN{{IFS="\t";OFS="\t";}}{{if(NR%2==1){{former=$7;}}else{{if(former>0||$7>0){{print 1;former=0;}}else{{print 0;former=0;}}}}}}'  > ${outputdir}/${region}.count.txt

#all_transcript_id
#transcriptome_genome_sort
cat ${outputdir}/sorted.reads.bed \
| bedtools coverage -split -f $fracA -counts -s -g $chrSize -a - -b ${beddir}/${region}.bed | awk 'BEGIN{{IFS="\t";OFS="\t";}} {{ if($7>0) {{print 1}} else{{print 0}} }}'  > ${outputdir}/${region}.count.txt
echo "Done!" >> ${outputdir}/log.txt
else 
echo "${region} already counted!" >> ${outputdir}/log.txt
fi
done
paste ${regionalCounts} > ${outputdir}/result.txt

# only assign each peak to one region: counts.txt
cat ${outputdir}/result.txt  | awk -v regions="${regions}" 'BEGIN{OFS="\t";split(regions,regionList," ");nRegions=length(regionList);for(i=1;i<=nRegions;i++){regionCounts[i]=0;}}{split($0,data," ");for(i=1;i<=nRegions;i++){if(data[i]>0){regionCounts[i]+=1;break;}}}END{for(i=1;i<=nRegions;i++){print regionList[i],regionCounts[i];}}' > ${outputdir}/counts.txt
