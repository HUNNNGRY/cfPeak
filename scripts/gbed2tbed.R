#! /usr/bin/env Rscript

# convert/map gbed to tbed
# last 2204 by bpf 
# b.p.f@qq.com

# no change checked at 220607


suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser(description='convert/map gbed to tbed')
parser$add_argument('-i', '--inputFile', type='character', default='-',
                    help='input BED file with covearge in the 5th column, default: -')
# parser$add_argument('-b','--rawReadBed', type='character', 
#                     help='rawread bed file for calculating count, canbe output from bedtools bamtobed, required')
parser$add_argument('-o', '--outputFile', type='character', default='-',
                    help='output BED file with adjusted p-value in the last column, default: -')
# parser$add_argument('-z', '--zeroTruncDistribution', type='character', default='negbin',
#                     help='model distribution, dist = c("poisson", "negbin", "geometric"), default: negbin')
# parser$add_argument('-p', '--pvalue', type='double', default=0.1,
#                     help='adjusted p-value threshold for defining significant domains, default: 0.1', metavar='NUMBER')
# parser$add_argument('-b', '--background', type='double', default=0.99,
#                     help='set coverage quantile below this quantile as background regions')
parser$add_argument('--bedtoolsPath', type='character', default="/BioII/lulab_b/baopengfei/biosoft",
                    help='bedtools path, default: /BioII/lulab_b/baopengfei/biosoft')
parser$add_argument('--mode', type='character', default="1toN",
                    help='default: 1toN (convert one gn peak to multi RNA type multi tx name); optional: 1to1 (convert one gn peak to one RNA type one tx name, multi-block still can exist)')
# parser$add_argument('--chrSize', type='character', default="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/all_transcript_id",
#                     help='chrSize file path, used for shuffle random regions, default: /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/all_transcript_id')
# parser$add_argument('--peakLengthBin', type="integer", default=10,
#                     help='bin for peak length treated as having the same distribution, default: 10')
# parser$add_argument('--minFitRegionNum', type="integer", default=100,
#                     help='number of regions with count>0 used for fit model, default: 100')
# parser$add_argument('--startSeed', type="integer", default=1234,
#                     help='start seed for shuffle random regions with the same length, default: 1234')
parser$add_argument('--cores', type="integer", default=4,
                    help='multi-cores processing, default: half of auto detect')
args <- parser$parse_args()

for(i in 1:length(args)){
  v.name <- names(args)[i]
  assign(v.name,args[[i]])
  message(paste0(v.name,": ",get(v.name)))
}

## test 
# inputFile <- "/BioII/lulab_b/baopengfei/shared_reference/RBP/tmp/JYF_top500.bed" #"/BioII/lulab_b/baopengfei/shared_reference/RBP/splitByRBP/AGO2_test.bed"
# outputFile <- "./test.bed"
# bedtoolsPath <- "/BioII/lulab_b/baopengfei/biosoft"



## library
options(stringsAsFactors = F)
#suppressPackageStartupMessages(library(countreg))
options(bedtools.path = bedtoolsPath) # /BioII/lulab_b/baopengfei/anaconda3/envs/exvariance/bin
suppressPackageStartupMessages(library(bedtoolsr))
suppressPackageStartupMessages(library(IRanges))
suppressPackageStartupMessages(library(BiocGenerics))


## def func. 
GRanges2bed <- function(gr=ds){
  #1 based to 0 based
  out <- data.frame(chr=gr@seqnames, start=gr@ranges@start-1, end=gr@ranges@start+gr@ranges@width-1, name=gr$name, score=gr$score, strand=gr@strand )
  out$chr <- as.character(out$chr)
  out$name <- as.character(out$name)
  out$strand <- as.character(out$strand)
  out <- as.data.frame(out)
  return(out)
}



# convert gbed to tbed ----------------------------------------------------
## read peak bed
#peak <- rtracklayer::import.bed(inputFile)
#peak.df <- GRanges2bed(peak)
peak.df <- data.table::fread(inputFile,data.table = F,header = F,sep = "\t",stringsAsFactors = F)
peak.df <- as.data.frame(peak.df)
peak.df <- peak.df[,1:6]
colnames(peak.df) <- c("chr","start","end","name","score","strand")
peak.df$chr <- as.character(peak.df$chr)
peak.df$name <- as.character(peak.df$name)
peak.df$strand <- as.character(peak.df$strand)

peak.df <- peak.df[peak.df$strand == "+" | peak.df$strand == "-" ,] # need contain strand info
#peak.df[1:3,]
#peak.df$chr <- paste0("chr",peak.df$chr)
#peak.df$chr <- gsub("chrMT","chrM",peak.df$chr)
#peak.df <- peak.df[,1:6]
#table(peak.df$chr)
#head(peak.df)
#str(peak.df)

# peak.df <- data.table::fread("/BioII/lulab_b/baopengfei/shared_reference/RBP/splitByRBP/AGO2.bed",data.table = F,header = F,sep = "\t",stringsAsFactors = F)
# peak.df[1:3,]
# str(peak.df)


## read tx table
chrSize.df <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/tbed/11RNA_map_newTxID.txt",data.table = F,header = F,sep = "\t",stringsAsFactors = F)
#chrSize.df <- chrSize.df[chrSize.df$V1 %in% unique(peak.cov$V1),]
#head(chrSize.df)
#table(duplicated(chrSize.df$transcript_id))
#chrSize.df <- chrSize.df[!duplicated(chrSize.df$transcript_id),]
colnames(chrSize.df) <- c("transcript_id","s","e","name","score","strand","transcript_type")
#table(chrSize.df$transcript_type)
priority <- c("rRNA","pri_miRNA","lncRNA","mRNA","piRNA","snoRNA","snRNA","srpRNA","tRNA","tucpRNA","Y_RNA") # ,"intron.for","intron.rev","promoter.for","promoter.rev","enhancer.for","enhancer.rev"
chrSize.df <- chrSize.df[chrSize.df$transcript_type %in% priority,]
#head(chrSize.df[grep("ENST",chrSize.df$transcript_id),],4)
#table(chrSize.df$chrom==chrSize.df$transcript_id)
# chrSize.df$start <- 0
# chrSize.df <- chrSize.df[,c(1,3,2)]
# chrSize.df$name <- "."
# chrSize.df$score <- "."
# chrSize.df$strand <- "+"

# long_RNA_table <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/transcript_table/long_RNA.txt",data.table = F,header = T,sep = "\t",stringsAsFactors = F)
# table(duplicated(long_RNA_table$transcript_id))
# head(long_RNA_table,3)
# table(long_RNA_table$chrom)



## read tx block bed
tx.ref <- rtracklayer::import.bed(file('/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/bed/11RNA_map_bed12_newTxID.bed'))
#table(tx.ref$name %in% chrSize.df$transcript_id)
tx.ref <- tx.ref[tx.ref$name %in% chrSize.df$transcript_id,]
#tx.ref <- tx.ref[sample(1:length(tx.ref),0.1*length(tx.ref),replace = F)] # test small sample
tx.ref.block0 <- rtracklayer::blocks(tx.ref)

#pre-computed block txt (10 tx bed)
tx.ref.block <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/bed/11RNA_map_block_newTxID.txt",sep = "\t",data.table = F,header = T,)
#tx.ref.block[1:3,]
tx.ref.block <- tx.ref.block[!grepl("alt|random|chrUn",tx.ref.block$seqnames,perl=T),] # filter non-canonical chr


#bedtools intersect -a A.bed -b B.bed -wao
ts <- bedtoolsr::bt.intersect(wao=T, s = T, a = peak.df, b = tx.ref.block) # -f, sorted = T, g = chr.size seems not working
ts <- ts[ts$V15 > 0 & !is.na(ts$V15),]
ts$a.ratio <- ts$V15/(ts$V3-ts$V2)
# ts <- ts[ts$a.ratio==1,] # only convert those peak with full-range within exon ?
#ts[1:3,]
ts$transcript_type <- chrSize.df$transcript_type[match(ts$V10,chrSize.df$transcript_id)]
#table(ts$V10 %in% chrSize.df$transcript_id) # should be all true
ts$transcript_type <- factor(ts$transcript_type,levels = priority)
#table(ts$transcript_type)
if(mode == '1to1'){
ts <- ts[order(paste0(ts$V4,ts$V1,ts$V2,ts$V3,ts$V6),as.numeric(ts$transcript_type),-as.numeric(ts$V13),ts$V10,ts$V14,decreasing = F),] # select longest tx of the same tx type, if multi ENST overlap. (v10: name; v13: block number; v14: block idx)
#table(duplicated(paste0(ts$V4,ts$V1,ts$V2,ts$V3,ts$V6))) # ,ts$transcript_type,ts$V13,ts$V10
ts <- ts[!duplicated(paste0(ts$V4,ts$V1,ts$V2,ts$V3,ts$V6,ts$V14)),] # only keep one gn peak with highest priority and overlap base num, if multi exist 
}
ts$tx.strand <- "+"


ts.list <- list()
getTxPos <- function(i){
#i <- 4
if(i%%1000 == 0){
  message(i)
}
tmp <- tx.ref.block0[[ts$V10[i]]] # v10: tx name
#which(ts$V10=="ENST00000609909.1")
#tmpn <- as.data.frame(tx.ref.block0[["ENST00000609909.1"]])

idx <- as.numeric(ts$V14[i]) # v14: block idx
if(ts$V12[i]=="+"){
  block.start <- max(0,(ts$V2[i] - ts$V8[i]))
  if(idx-1 >= 1){
    ts.list[[i]] <- c(tx.start=sum(width(tmp)[1:(idx-1)],na.rm = T) + block.start,
                      tx.end=sum(width(tmp)[1:(idx-1)],na.rm = T) + block.start + ts$V15[i])
  }else{
    ts.list[[i]] <- c(tx.start=block.start,
                      tx.end=block.start + ts$V15[i])
  }
  #ts$tx.start[i] <- sum(width(tmp)[1:(idx-1)]) + (ts$V2[i] - ts$V8[i])
  #ts$tx.end[i] <- ts$tx.start[i] + ts$V15[i]
} else if (ts$V12[i]=="-"){
  block.end <- max(0,(ts$V9[i] - ts$V3[i]))
  if(idx+1 <= length(tmp)){
    ts.list[[i]] <- c(tx.start=sum(width(tmp)[(idx+1):length(tmp)],na.rm = T) + block.end,
                    tx.end=sum(width(tmp)[(idx+1):length(tmp)],na.rm = T) + block.end + ts$V15[i])
  }else{
    ts.list[[i]] <- c(tx.start=block.end,
                      tx.end=block.end + ts$V15[i])   
  }
  #ts$tx.start[i] <- sum(width(tmp)[(idx+1):length(width(tmp))]) + (ts$V9[i] - ts$V3[i])
  #ts$tx.end[i] <- ts$tx.start[i] + ts$V15[i]
}
}

library(parallel)
# cores <- 0.5*detectCores() 
ts.1 <- mclapply(1:nrow(ts),getTxPos,mc.cores = args$cores)
#ts.1 <- lapply(1:nrow(ts),getTxPos)
#head(ts.1,4)

ts.1 <- do.call(rbind,ts.1)
ts.1 <- as.data.frame(ts.1)
#table((ts.1$tx.end-ts.1$tx.start)==(ts$V3-ts$V2))
#tx.ref.block0[["ENST00000640238.1"]]
ts <- cbind(ts,ts.1)

#head(ts,3)
#colnames(ts)
ts$orignal.id <- paste(ts$V4,ts$V1,ts$V2,ts$V3,ts$V6,sep = "|")
table((ts$tx.end-ts$tx.start)==(ts$V3-ts$V2))             # should be all TRUE if filter with (ts <- ts[ts$a.ratio==1,])
ts <- ts[,c("V10","tx.start","tx.end","orignal.id","V5","tx.strand")]


## write outfile
message('successfully mapped intervals: ', nrow(ts))
if(outputFile == '-'){
  data.table::fwrite(ts, row.names=FALSE, col.names=FALSE, sep='\t', quote = F)
} else {
  data.table::fwrite(ts, outputFile, row.names=FALSE, col.names=FALSE, sep='\t', quote = F)
}

#peak.df[peak.df$name=="AGO2-chr1:100007036-100007072",]

