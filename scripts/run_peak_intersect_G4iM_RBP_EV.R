#! /usr/bin/env Rscript

#domain intersect with RBPs, G4/i-motif, Alu et al.
# last 2406 by bpf 
# b.p.f@qq.com

#todo: add 8DNA tx reference bed, and change 8DNA intersect 


suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser(description='intersect domain bed with G4-iM, RBPs, AGO2, RBPhotspot')
parser$add_argument('-i', '--inputFile', type='character', default='-',
                    help='input domain BED file') #, default: -
parser$add_argument('-o', '--outputFile', type='character', default='-',
                    help='output BED file')
parser$add_argument('--bedtoolsPath', type='character', default="/BioII/lulab_b/baopengfei/biosoft",
                    help='bedtools path, default: /BioII/lulab_b/baopengfei/biosoft')
parser$add_argument('--chrSize', type='character', default="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/transcriptome_genome_sort_uniq_newTxID", 
                    help='chrSize file path, used for shuffle random regions, contains both 11+8 tx and 23 chr. default: /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/transcriptome_genome_sort_uniq_newTxID')
# /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/10RNA_newTxID
parser$add_argument('--seed', type="integer", default=1234,
                    help='start seed for shuffle random regions with the same length, default: 1234')
parser$add_argument('--RNAtype', type='character', default="RNA",
                    help='RNA shuffle type range, default: RNA, means shuf across all 11 RNA types; option: DNA, means shuf across all 8 DNA types')
parser$add_argument('--coord', type='character', default="tx",
                    help='peak coord type, default: tx, ENST*** 1 20; op: gn, chr1 10001 10020')
parser$add_argument('--shufSameRNATxInPriority', type='logical', default=FALSE,
                    help='random shuffle region in the same RNA species (may get err if one species has little tx), default: FALSE; op: TRUE')
parser$add_argument('--OverlapTypes', type='character', default='AGO2,otherRBPs,RBPhotspot,G4,iM,EV,Alu', help='OverlapTypes to run intersect, default: AGO2,otherRBPs,RBPhotspot,G4,iM,EV,Alu')


# parser$add_argument('--cores', type="integer", #default=4,
#                     help='multi-cores processing, default: half of auto detect')
args <- parser$parse_args()

#todo: add param for shuf in all tx of same RNA type

for(i in 1:length(args)){
  v.name <- names(args)[i]
  assign(v.name,args[[i]])
  message(paste0(v.name,": ",get(v.name)))
}

# eg. running code -----
# # run
#see onenote lulab exRNA fragment/peak



# # test
# setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/")
# bedtoolsPath <- "/BioII/lulab_b/baopengfei/biosoft"
# chrSize <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/transcriptome_genome_sort_uniq_newTxID"
# inputFile <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008/call_peak_all/piranha_by_sample/b5_p01/intersect/SAMN03863400.bed6" # b5_d05_p01.bed
# outputFile <- "./test.intersect.bed"
# seed <- 1234
# RNAtype <- "RNA"
# coord <- "tx"
# shufSameRNATxInPriority <- TRUE

# setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/")
# bedtoolsPath <- "/BioII/lulab_b/baopengfei/biosoft"
# chrSize <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/transcriptome_genome_sort_uniq_newTxID"
# inputFile <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008/call_domain_withRepeats_all/domains_by_sample/b5_p05/intersect_8DNA/SRR2105340_gn.bed" # b5_d05_p01.bed
# outputFile <- "./test.intersect.bed"
# seed <- 1234
# coord <- "gn"



# formal running code ------------------------------
## library
options(stringsAsFactors = F)
options(bedtools.path = bedtoolsPath) # /BioII/lulab_b/baopengfei/anaconda3/envs/exvariance/bin
suppressPackageStartupMessages(library(bedtoolsr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(dplyr))


## get AGO2 tx bed (no DNA)
if(coord=="tx"){
  AGO2 <- data.table::fread("/BioII/lulab_b/baopengfei/shared_reference/RBP/splitByRBP/AGO2_tx.bed",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)
}else if(coord=="gn"){
  AGO2 <- data.table::fread("/BioII/lulab_b/baopengfei/shared_reference/RBP/splitByRBP/AGO2.bed",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)
}
AGO2 <- AGO2[,1:6]
colnames(AGO2) <- c("chr","start","end","name","score","strand")


## get plasma otherRBPs tx bed (no DNA)
if(coord=="tx"){
  otherRBPs <- data.table::fread("/BioII/lulab_b/baopengfei/shared_reference/RBP/splitByRBP/otherRBPs_tx.bed",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)
}else if(coord=="gn"){
  otherRBPs <- data.table::fread("/BioII/lulab_b/baopengfei/shared_reference/RBP/splitByRBP/otherRBPs.bed",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)
}
otherRBPs <- otherRBPs[,1:6]
colnames(otherRBPs) <- c("chr","start","end","name","score","strand")


## get RBPhotspot tx bed (no DNA)
if(coord=="tx"){
  RBPhotspot <- data.table::fread("/BioII/lulab_b/baopengfei/shared_reference/RBP/POSTAR3_hotspot/all_merge_newTxID.bed",sep = '\t',header = F,stringsAsFactors = F)
}else if(coord=="gn"){
  RBPhotspot <- data.table::fread("/BioII/lulab_b/baopengfei/shared_reference/RBP/POSTAR3_hotspot/all_merge_newTxID_gn.bed",sep = '\t',header = F,stringsAsFactors = F)
}
RBPhotspot <- RBPhotspot[,1:6]
colnames(RBPhotspot) <- c("chr","start","end","name","score","strand")


## G4 (no DNA)
if(coord=="tx"){
  # G4 <- data.table::fread("/BioII/lulab_b/baopengfei/shared_reference/structure/G4iMGrinder/bedg4im_merge_tx.bed6",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)
  G4 <- data.table::fread("/BioII/lulab_b/baopengfei/shared_reference/structure/quadratlas_g4grinder_tx.bed",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)
}else if(coord=="gn"){
  G4 <- data.table::fread("/BioII/lulab_b/baopengfei/shared_reference/structure/quadratlas_g4grinder.bed",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)
}
G4 <- G4[,1:6]
colnames(G4) <- c("chr","start","end","name","score","strand")


## i-motif (no DNA)
if(coord=="tx"){
  iM <- data.table::fread("/BioII/lulab_b/baopengfei/shared_reference/structure/G4iMGrinder/bedim_tx.bed",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)
}else if(coord=="gn"){
  iM <- data.table::fread("/BioII/lulab_b/baopengfei/shared_reference/structure/G4iMGrinder/bedim.bed",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)
}
iM <- iM[,1:6]
colnames(iM) <- c("chr","start","end","name","score","strand")


## EV
if(coord=="tx"){
  EV <- data.table::fread("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/lulab/11RNAsmallDomain_EVenrich.bed",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)
}else if(coord=="gn"){
  EV <- data.table::fread("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/lulab/8DNAsmallDomain_EVenrich_gn.bed",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)
} #8DNA with gn
EV <- EV[,1:6]
colnames(EV) <- c("chr","start","end","name","score","strand")


## TO TEST:
## Alu
if(coord=="tx"){
  Alu <- data.table::fread("/BioII/lulab_b/baopengfei/projects/RP-RNA/bed/Alu_tx.bed",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)
}else if(coord=="gn"){
  Alu <- data.table::fread("/BioII/lulab_b/baopengfei/projects/RP-RNA/bed/Alu.bed",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)
} #8DNA with gn
Alu <- Alu[,1:6]
colnames(Alu) <- c("chr","start","end","name","score","strand")


l <- c("AGO2","otherRBPs","RBPhotspot","G4","iM","EV","Alu")
k <- list(AGO2,otherRBPs,RBPhotspot,G4,iM,EV,Alu)
names(k) <- l
OverlapTypes <- strsplit(OverlapTypes,",")[[1]]
for(i in l){
  if(!(i %in% OverlapTypes)){
    message("rm ",i)
    k[[i]] <- NULL
  }
}
l <- names(k)


## get tx type (no chr, only ENST)
ref <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)


if(coord=="tx"){
  if(RNAtype=="RNA"){
    ref <- ref[grepl("RNA",ref$transcript_type),]
  }else if(RNAtype=="DNA"){
    ref <- ref[!grepl("RNA",ref$transcript_type),]
  }
}else if(coord=="gn"){
  # ref <- ref[!grepl("RNA",ref$transcript_type),]
  print("gn coord, not use ref")
}



## read domain bed
domain <- as.data.frame(data.table::fread(inputFile,data.table = F,sep = '\t',check.names = F,stringsAsFactors = F))
domain <- domain[,1:6] # this version bed12 not considered
colnames(domain) <- c("chr","start","end","name","score","strand")
domain$name <- paste0(domain$name,"--",domain$chr,"--",domain$start,"--",domain$end)
if(coord=="tx"){
  if(RNAtype=="RNA"){
    domain <- domain[!grepl("__chr",domain$chr),] # remove 8DNA
  } else if (RNAtype=="DNA"){
    domain <- domain[grepl("__chr",domain$chr),] # remove 11RNA
  }
}

## shuf only in domain related ? 


## filter chrsize
chrSize <-  data.table::fread(chrSize) # contains both chr, ENST
if(coord=="tx"){
  chrSize <- chrSize[chrSize$V1 %in% ref$transcript_id,] # 11 RNA or 8DNA
} else if(coord=="gn"){
  chrSize <- chrSize[chrSize$V1 %in% unique(domain$chr),] # 23 chr (only domain related, too long, not affected)
}

## shuf domain
if (!shufSameRNATxInPriority){ # mainly for 11RNA
#domain2 <- bedtoolsr::bt.shuffle(seed = seed, chrom = T, i = domain, g = chrSize) # noOverlapping=F, excl=NULL, 
domain2 <- bedtoolsr::bt.shuffle(seed = seed, noOverlapping=T, i = domain, g = chrSize) # chrom = F, noOverlapping=F, excl=NULL, 
write.table(domain2,paste0(inputFile,".shuffle"),sep = '\t',row.names = F,col.names = F,quote = F)
} else {
  #table(ref$transcript_type)
  message("shufSameRNATxInPriority !")
  rnaType <- ref$transcript_type[match(domain$chr,ref$transcript_id)]
  chrType <- ref$transcript_type[match(chrSize$V1,ref$transcript_id)]
  tmpList <- list()
  for (rna in unique(rnaType)){
    # message(rna)
    #rna <- "piRNA"
    domain.tmp <- domain[rnaType==rna,]
    chrSize.tmp <- chrSize[chrType==rna,]
    domain2.tmp <- bedtoolsr::bt.shuffle(seed = seed, chromFirst=T, noOverlapping=F, i = domain.tmp, g = chrSize.tmp) # noOverlapping=T: error for longest full-tx length peak (also need chromFirst ?)
    tmpList[[rna]] <- domain2.tmp
  }
  domain2 <- as.data.frame(do.call(rbind,tmpList))
  domain2 <- domain2[match(domain$name,domain2$V4),] #arrange order
  write.table(domain2,paste0(inputFile,".shuffle2"),sep = '\t',row.names = F,col.names = F,quote = F)
}

## flank domain
### 1st: use left flank 
domain3 <- bedtoolsr::bt.shift(p=-1, m=0, pct=T, i = domain, g = chrSize)
idx <- which((domain3$V3-domain3$V2)<(domain$end-domain$start))

### 2nd: use right flank if left flank shorter than peak
domain3.right <- bedtoolsr::bt.shift(p=+1, m=0, pct=T, i = domain[idx,], g = chrSize)
idx2 <- which((domain3$V3[idx]-domain3$V2[idx])<(domain3.right$V3-domain3.right$V2))
domain3[idx,][idx2,] <- domain3.right[idx2,]

### 3rd: use 11RNA random shuffle if both flank shorter than peak
idx3 <- which((domain3$V3-domain3$V2)<(domain$end-domain$start))
idx4 <- which((domain3$V3[idx3]-domain3$V2[idx3])<(domain2$V3[idx3]-domain2$V2[idx3])) # actually not needed
domain3[idx3,][idx4,] <- domain2[idx3,][idx4,]

# table((domain3.right$V3-domain3.right$V2)<=1) 
table((domain3$V3-domain3$V2)<(domain$end-domain$start)) # should be all FALSE !!!!!!!!!!!!!!!!!!!!
# table((domain3.right$V3-domain3.right$V2)<(domain$end[idx]-domain$start[idx]))
# domain3$V3[domain3$V2<0] <- 1
# domain3$V2[domain3$V2<0] <- 0
#table(domain3$V2<0)
write.table(domain3,paste0(inputFile,".flank"),sep = '\t',row.names = F,col.names = F,quote = F)

domains <- list(domain,domain2,domain3)
if (!shufSameRNATxInPriority){
  outfiles <- c(outputFile,paste0(outputFile,".shuffle"),paste0(outputFile,".flank"))
}else{
  outfiles <- c(outputFile,paste0(outputFile,".shuffle2"),paste0(outputFile,".flank"))
}






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
    
    domain.ago2 <- bedtoolsr::bt.intersect(a = domain.tmp,wao = T,s = T,b = k[i]) # this version bed12 not considered, need add -split for bed12 input
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
  
  if(coord=="tx"){
    res.df$RNA <- ref$transcript_type[match(res.df$chr,ref$transcript_id)]
    # res.df$RNA <- gsub("_for|_rev|\\.for|\\.rev","",res.df$RNA,perl = T)
    # rna <- c("pri_miRNA","lncRNA", "mRNA", "piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA')
    # dna <- c("intron","promoter", "enhancer","repeats") # 
    # res.df$RNA <- factor(res.df$RNA,levels = c(rna,dna))
    # table(res.df$RNA)
    
    #get tx.len
    res.df$tx_len <- ref$tx.length[match(res.df$chr,ref$transcript_id)]
  }else if(coord=="gn"){
    res.df$RNA <- "8DNA"
    res.df$tx_len <- -1
  }

  
  # res.df <- res.df[order(res.df$site_overlap_ratio_product,decreasing = T),]
  # res.df <- na.omit(res.df)
  # res.df[1:3,]
  
  ## write table
  write.table(res.df,outfile,sep = '\t',row.names = F,col.names = T,quote = F)
}



