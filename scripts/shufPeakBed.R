#! /usr/bin/env Rscript

# shuffle peak bed6 
# last 2207 by bpf 
# b.p.f@qq.com

suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser(description='shuf domain bed randomly')
parser$add_argument('-i', '--inputFile', type='character', default='-',
                    help='input domain BED file') #, default: -
# parser$add_argument('-o', '--outputFile', type='character', default='-',
#                     help='output BED file')
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
# parser$add_argument('--cores', type="integer", #default=4,
#                     help='multi-cores processing, default: half of auto detect')
parser$add_argument('--region', type='character', default="flank,shuffle",
                    help='output region type, default: flank,shuffle (means output both flank & shuffle), opt: flank or shuffle')
parser$add_argument('--shuffleSurfix', type='character', default="shuffle",
                    help='output shuffle surfix of appended to input path')
parser$add_argument('--flankSurfix', type='character', default="flank",
                    help='output flank surfix of appended to input path')
args <- parser$parse_args()

for(i in 1:length(args)){
  v.name <- names(args)[i]
  assign(v.name,args[[i]])
  message(paste0(v.name,": ",get(v.name)))
}
regions <- strsplit(region,",")[[1]]

# eg. running code -----
# # run
#see onenote lulab exRNA fragment/peak



# # test
# setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/")
# bedtoolsPath <- "/BioII/lulab_b/baopengfei/biosoft"
# chrSize <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/transcriptome_genome_sort_uniq_newTxID"
# inputFile <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008/call_domain_withRepeats_all/domains_by_sample/b5_p05/intersect/SRR2105340.bed" # b5_d05_p01.bed
# outputFile <- "./test.intersect.bed"
# seed <- 1234
# RNAtype <- "RNA"
# coord <- "tx"

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


## get tx type (no chr, only ENST)
ref <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)
rna <- c("rRNA","pri_miRNA","lncRNA", "mRNA", "piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA')
dna <- c("intron","promoter", "enhancer", "repeats") 


if(coord=="tx"){
  if(RNAtype=="RNA"){
    ref <- ref[grepl(paste(rna,collapse = "|"),ref$transcript_type,perl = T),]
  }else if(RNAtype=="DNA"){
    # ref <- ref[!grepl("RNA",ref$transcript_type),]
    ref <- ref[grepl(paste(dna,collapse = "|"),ref$transcript_type,perl = T),]
  }
}else if(coord=="gn"){
  # ref <- ref[!grepl("RNA",ref$transcript_type),]
  print("gn coord, not use ref")
}


## read domain bed
domain <- as.data.frame(data.table::fread(inputFile,data.table = F,sep = '\t',check.names = F,stringsAsFactors = F))
domain <- domain[,1:6]
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


if("shuffle" %in% regions){
  ## shuf/rand domain/peak
  if (!shufSameRNATxInPriority){ # mainly for 11RNA
    #domain2 <- bedtoolsr::bt.shuffle(seed = seed, chrom = T, i = domain, g = chrSize) # noOverlapping=F, excl=NULL, 
    domain2 <- bedtoolsr::bt.shuffle(seed = seed, noOverlapping=T, i = domain, g = chrSize) # chrom = F, noOverlapping=F, excl=NULL, 
    write.table(domain2,paste0(inputFile,".",shuffleSurfix),sep = '\t',row.names = F,col.names = F,quote = F)
  }else{
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
    write.table(domain2,paste0(inputFile,".",shuffleSurfix,".pri"),sep = '\t',row.names = F,col.names = F,quote = F)
  }
}

if("flank" %in% regions){
  ## flank domain
  ### 1st: use left flank 
  domain3 <- bedtoolsr::bt.shift(p=-1, m=0, pct=T, i = domain, g = chrSize)
  idx <- which((domain3$V3-domain3$V2)<(domain$end-domain$start))
  
  ### 2nd: use right flank if left flank shorter than peak
  domain3.right <- bedtoolsr::bt.shift(p=+1, m=0, pct=T, i = domain[idx,], g = chrSize)
  idx2 <- which((domain3$V3[idx]-domain3$V2[idx])<(domain3.right$V3-domain3.right$V2))
  domain3[idx,][idx2,] <- domain3.right[idx2,]
  
  ### 3rd: use 11RNA random shuffle if both flank shorter than peak
  if("shuffle" %in% regions){
  idx3 <- which((domain3$V3-domain3$V2)<(domain$end-domain$start))
  idx4 <- which((domain3$V3[idx3]-domain3$V2[idx3])<(domain2$V3[idx3]-domain2$V2[idx3])) # actually not needed
  domain3[idx3,][idx4,] <- domain2[idx3,][idx4,]
  }
  
  # table((domain3.right$V3-domain3.right$V2)<=1) 
  table((domain3$V3-domain3$V2)<(domain$end-domain$start)) # should be all FALSE !!!!!!!!!!!!!!!!!!!!
  # table((domain3.right$V3-domain3.right$V2)<(domain$end[idx]-domain$start[idx]))
  # domain3$V3[domain3$V2<0] <- 1
  # domain3$V2[domain3$V2<0] <- 0
  #table(domain3$V2<0)
  write.table(domain3,paste0(inputFile,".",flankSurfix),sep = '\t',row.names = F,col.names = F,quote = F)
}



