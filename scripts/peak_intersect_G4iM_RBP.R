#! /usr/bin/env Rscript

#domain intersect with RBPs, G4/i-motif
# last 2207 by bpf 
# b.p.f@qq.com

suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser(description='intersect domain bed with G4-iM, RBPs, AGO2, RBPhotspot')
parser$add_argument('-i', '--inputFile', type='character', default='-',
                    help='input domain BED file') #, default: -
parser$add_argument('-o', '--outputFile', type='character', default='-',
                    help='output BED file')
parser$add_argument('--bedtoolsPath', type='character', default="/BioII/lulab_b/baopengfei/biosoft",
                    help='bedtools path, default: /BioII/lulab_b/baopengfei/biosoft')
parser$add_argument('--chrSize', type='character', default="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/transcriptome_genome_sort_uniq_newTxID", 
                    help='chrSize file path, used for shuffle random regions, default: /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/10RNA_newTxID')
# /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/10RNA_newTxID
parser$add_argument('--seed', type="integer", default=1234,
                    help='start seed for shuffle random regions with the same length, default: 1234')
parser$add_argument('--RNAtype', type='character', default="RNA",
                    help='RNA shuffle type range, default: RNA, means shuf across all 11 RNA types')
# parser$add_argument('--cores', type="integer", #default=4,
#                     help='multi-cores processing, default: half of auto detect')
args <- parser$parse_args()

for(i in 1:length(args)){
  v.name <- names(args)[i]
  assign(v.name,args[[i]])
  message(paste0(v.name,": ",get(v.name)))
}

# eg. running code -----
## get second structure p value
#GSE71008
#GSE133684
# bed="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008/call_domain_withRepeats_all/domains_localmax_EM/b5_d05_p01.bed"
# bed="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008/call_domain_withRepeats_all/domains_localmax_EM/b5_d05_p01.bed.shuffle"
# bed="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE133684/call_domain_withRepeats_dedupByPos/domains_localmax_EM/b5_d05_p01.bed"
# bed="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE133684/call_domain_withRepeats_dedupByPos/domains_localmax_EM/b5_d05_p01.bed.shuffle"
# bedtools slop -s -l 20 -r 20 \
#   -g /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/transcriptome_genome_sort_uniq_newTxID \
#   -i $bed > \
#   ${bed}.tmp
# bedtools getfasta -name \
#   -fi /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/fasta_newTxID/combine18.fa \
#   -bed ${bed}.tmp > \
#   ${bed}.fa
# /BioII/lulab_b/baopengfei/anaconda3/envs/exvariance/bin/python3 /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/bin/rnafold_dinushuffle_parallel.py \
#   ${bed}.fa 20 1234 ${bed}.fa.csv
# rm ${bed}.tmp ${bed}.fa_perm



# # run
#ssh cnode
# dst="GSE71008" # GSE71008,GSE94533,GSE110381
# pre="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/${dst}/call_domain_withRepeats_all/domains_localmax_significant/b5_d05_p10"
# mkdir -p ${pre}/intersect
# for i in `cat data/$dst/sample_ids.txt`
# do
# echo $i
# cat ${pre}/${i}.bed | grep -v "_pos" | grep -v "_neg" > ${pre}/intersect/${i}.bed # | grep -v "NR_"
# /usr/bin/Rscript bin/domain_intersect_G4iM_RBP.R -i ${pre}/intersect/${i}.bed -o ${pre}/intersect/${i}.intersect.bed
# done &
# 
# pre="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/${dst}/call_domain_withRepeats_all/domains_localmax_by_sample_EM2/b5_d05_p10"
# mkdir -p ${pre}/intersect
# for i in `cat data/$dst/sample_ids.txt`
# do
# echo $i
# cat ${pre}/${i}.bed | grep -v "_pos" | grep -v "_neg" > ${pre}/intersect/${i}.bed # | grep -v "NR_"
# /usr/bin/Rscript bin/domain_intersect_G4iM_RBP.R -i ${pre}/intersect/${i}.bed -o ${pre}/intersect/${i}.intersect.bed
# done &
# 
# pre="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/${dst}/call_domain_withRepeats_all/domains_by_sample/b5_p10"
# mkdir -p ${pre}/intersect
# for i in `cat data/$dst/sample_ids.txt`
# do
# echo $i
# cat ${pre}/${i}.bed | grep -v "_pos" | grep -v "_neg" > ${pre}/intersect/${i}.bed # | grep -v "NR_"
# /usr/bin/Rscript bin/domain_intersect_G4iM_RBP.R -i ${pre}/intersect/${i}.bed -o ${pre}/intersect/${i}.intersect.bed
# done &
# 
# pre="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/${dst}/call_domain_withRepeats_all/domains_clipper_by_sample/b5_p10"
# mkdir -p ${pre}/intersect
# for i in `cat data/$dst/sample_ids.txt`
# do
# echo $i
# cat ${pre}/${i}.bed | grep -v "_pos" | grep -v "_neg" > ${pre}/intersect/${i}.bed # | grep -v "NR_"
# /usr/bin/Rscript bin/domain_intersect_G4iM_RBP.R -i ${pre}/intersect/${i}.bed -o ${pre}/intersect/${i}.intersect.bed
# done &
# 
# pre="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/${dst}/call_domain_withRepeats_all/domains_clam_by_sample/b5_p10"
# mkdir -p ${pre}/intersect
# for i in `cat data/$dst/sample_ids.txt`
# do
# echo $i
# cat ${pre}/${i}.bed | grep -v "_pos" | grep -v "_neg" > ${pre}/intersect/${i}.bed # | grep -v "NR_"
# /usr/bin/Rscript bin/domain_intersect_G4iM_RBP.R -i ${pre}/intersect/${i}.bed -o ${pre}/intersect/${i}.intersect.bed
# done &





# # test
# setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/")
# bedtoolsPath <- "/BioII/lulab_b/baopengfei/biosoft"
# chrSize <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/transcriptome_genome_sort_uniq_newTxID"
# inputFile <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008/call_domain_withRepeats_all/domains_by_sample/b5_p10//intersect/SRR2105127.bed" # b5_d05_p01.bed
# outputFile <- "./test.intersect.bed"
# seed <- 1234
# RNAtype <- "RNA"

# formal running code ------------------------------
## library
options(stringsAsFactors = F)
options(bedtools.path = bedtoolsPath) # /BioII/lulab_b/baopengfei/anaconda3/envs/exvariance/bin
suppressPackageStartupMessages(library(bedtoolsr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(dplyr))


## get AGO2 tx bed (no DNA)
ago2 <- data.table::fread("/BioII/lulab_b/baopengfei/shared_reference/RBP/splitByRBP/AGO2_tx.bed",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)
ago2 <- ago2[,1:6]
colnames(ago2) <- c("chr","start","end","name","score","strand")


## get plasma otherRBPs tx bed (no DNA)
otherRBPs <- data.table::fread("/BioII/lulab_b/baopengfei/shared_reference/RBP/splitByRBP/otherRBPs_tx.bed",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)
otherRBPs <- otherRBPs[,1:6]
colnames(otherRBPs) <- c("chr","start","end","name","score","strand")


## get RBPhotspot tx bed (no DNA)
hotspot <- data.table::fread("/BioII/lulab_b/baopengfei/shared_reference/RBP/POSTAR3_hotspot/all_merge_newTxID.bed",sep = '\t',header = F,stringsAsFactors = F)
hotspot <- hotspot[,1:6]
colnames(hotspot) <- c("chr","start","end","name","score","strand")


## G4,i-motif (no DNA)
G4 <- data.table::fread("/BioII/lulab_b/baopengfei/shared_reference/structure/G4iMGrinder/bedg4im_merge_tx.bed6",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)
G4 <- G4[,1:6]
colnames(G4) <- c("chr","start","end","name","score","strand")

k <- list(ago2,otherRBPs,hotspot,G4)
l <- c("AGO2","otherRBPs","RBPhotspot","G4iM")



## get tx type
ref <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)
#ref <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length.txt",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)
# table(ref$transcript_type)
# ref <- ref[is.na(ref$tx.length),]
if(RNAtype=="RNA"){
  ref <- ref[grepl("RNA",ref$transcript_type),]
}
#todo: tx_len=NA?

## read domain bed
domain <- as.data.frame(data.table::fread(inputFile,data.table = F,sep = '\t',check.names = F,stringsAsFactors = F))
domain <- domain[,1:6]
colnames(domain) <- c("chr","start","end","name","score","strand")
domain$name <- paste0(domain$name,"--",domain$chr,"--",domain$start,"--",domain$end)
if(RNAtype=="RNA"){
  domain <- domain[!grepl("__chr",domain$chr),] # remove DNA
}


chrSize <-  data.table::fread(chrSize)
chrSize <- chrSize[chrSize$V1 %in% ref$transcript_id,]
## shuf domain
#domain2 <- bedtoolsr::bt.shuffle(seed = seed, chrom = T, i = domain, g = chrSize) # noOverlapping=F, excl=NULL, 
domain2 <- bedtoolsr::bt.shuffle(seed = seed, noOverlapping=T, i = domain, g = chrSize) # chrom = F, noOverlapping=F, excl=NULL, 
write.table(domain2,paste0(inputFile,".shuffle"),sep = '\t',row.names = F,col.names = F,quote = F)

## flank domain
domain3 <- bedtoolsr::bt.shift(p=-1, m=0, pct=T, i = domain, g = chrSize)
domain3$V3[domain3$V2<0] <- 1
domain3$V2[domain3$V2<0] <- 0
write.table(domain3,paste0(inputFile,".flank"),sep = '\t',row.names = F,col.names = F,quote = F)

domains <- list(domain,domain2,domain3)
outfiles <- c(outputFile,paste0(outputFile,".shuffle"),paste0(outputFile,".flank"))

## run bedtools intersect
for (m in 1:length(domains)){
  #print(m)
  domain.tmp <- domains[[m]]
  outfile <- outfiles[[m]]
  
  res <- list()
  for (i in 1:length(l)){
    # x <- "AGO2"
      x <- l[i]
    #  i <- 1
    print(x)
    
    domain.ago2 <- bedtoolsr::bt.intersect(a = domain.tmp,wao = T,s = T,b = k[i])
    domain.ago2 <- as_tibble(domain.ago2) %>% 
      dplyr::mutate(overlap_ratio_product = V13/(V3-V2) * (V13/(V9-V8)), site = x) %>% 
      dplyr::arrange(desc(overlap_ratio_product)) %>% 
      dplyr::distinct(V4, .keep_all = TRUE)   # only keep max overlap_ratio_product each peak
      # dplyr::rename_all(funs(paste0(., "e", "E")))
      #dplyr::rename()
    colnames(domain.ago2) <- c( c("chr","start","end","name","score","strand"),paste0("site","_",c("chr","start","end","name","score","strand","overlap_base","overlap_ratio_product") ) ,"site" ) 
    res[[x]] <- domain.ago2
  }
  res.df <- do.call("rbind",res)
  #res.df <- na.omit(res.df)
  #dim( na.omit(domain.recenter) )
  #length(unique(res.df$chr))
  
  res.df$RNA <- ref$transcript_type[match(res.df$chr,ref$transcript_id)]
  # res.df$RNA <- gsub("_for|_rev|\\.for|\\.rev","",res.df$RNA,perl = T)
  # rna <- c("pri_miRNA","lncRNA", "mRNA", "piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA')
  # dna <- c("intron","promoter", "enhancer","repeats") # 
  # res.df$RNA <- factor(res.df$RNA,levels = c(rna,dna))
  # table(res.df$RNA)
  
  #get tx.len
  res.df$tx_len <- ref$tx.length[match(res.df$chr,ref$transcript_id)]
  
  # res.df <- res.df[order(res.df$site_overlap_ratio_product,decreasing = T),]
  # res.df <- na.omit(res.df)
  # res.df[1:3,]
  
  ## write table
  write.table(res.df,outfile,sep = '\t',row.names = F,col.names = T,quote = F)
}



