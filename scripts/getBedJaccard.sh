#calculate overlap
cmb_file=$1 # cmb.txt
jac_dir=$2 #"test/$dst/jaccard"
mkdir -p $jac_dir

#awk -f scripts/combinations.awk <<< `ls test/${dst}/*_*.bed_optimal_cfpeak.bed | tr "\n" " "` > ${jac_dir}/cmb.txt

rm -f ${jac_dir}/jac.txt
echo -e "from\tto\tfrom_n\tto_n\tlen_intersection\tlen_union\tjaccard\tn_intersections" >> ${jac_dir}/jac.txt
for tmp in `cat $cmb_file | awk '{print $1 "," $2}'`
do
	i=`echo $tmp | cut -d ',' -f 1`
	j=`echo $tmp | cut -d ',' -f 2`
	echo -e "${i}:${j}"
	i_n=`wc -l $i | cut -d ' ' -f 1`
	j_n=`wc -l $j | cut -d ' ' -f 1`
	
	bedtools jaccard -s -f 0.1 -F 0.1 -a  ${i} -b ${j} | grep -v 'jaccard' | awk -v i=$i -v j=$j -v i_n=$i_n -v j_n=$j_n 'BEGIN {FS=OFS="\t"} {print i,j,i_n,j_n,$0}' >> ${jac_dir}/jac.txt
done

