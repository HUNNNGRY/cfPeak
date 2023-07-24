#!/usr/bin/bash

$i=$1
pre="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008/call_domain_withRepeats_all/domains_localmax_significant_EM/b5_d05_p01/"
echo $i
mkdir -p ${pre}/intersect
cat ${pre}/${i}.bed | grep -v "_pos" | grep -v "_neg" | grep -v "NR_" > ${pre}/intersect/${i}.bed
/usr/bin/Rscript bin/domain_intersect_G4iM_RBP.R -i ${pre}/intersect/${i}.bed -o ${pre}/intersect/${i}.intersect.bed
