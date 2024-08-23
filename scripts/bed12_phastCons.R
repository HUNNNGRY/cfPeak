#! /usr/bin/env Rscript

# last 2023.12 by bpf 
# b.p.f@qq.com


suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser(description='calculate conservation phastCons score for bed12 input')
parser$add_argument('--inBed12', type='character', required=TRUE,
                    help='input bed12 file, 4th name column need to be unique each row')
parser$add_argument('--outFile', type='character', required=TRUE,
                    help='output outFile is bed6: score column means average phastCons score for each bed12 record')
parser$add_argument('--cores', type="integer", default=-1,
                    help='multi-cores processing, default: half of auto detect')
args <- parser$parse_args()

for(i in 1:length(args)){
  v.name <- names(args)[i]
  assign(v.name,args[[i]])
  message(paste0(v.name,": ",get(v.name)))
}



# # test
# pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
# outpre <- "/data2/lulab1/bpf/projects/WCHSU-FTC"
# dst <- "SLE"
# dedup <- "dedup"
# 
# inBed12 <- paste0(outpre,"/output/",dst,"/call_peak_",dedup,"/cfpeakCNN/b5_d50_p1_11RNA_gn.bed")
# 
# outFile <- paste0(outpre,"/output/",dst,"/call_peak_",dedup,"/cfpeakCNN/b5_d50_p1.phastCons")
# cores <- 5




options(stringsAsFactors = F)
suppressPackageStartupMessages(library(phastCons100way.UCSC.hg38))
#suppressPackageStartupMessages(library(IRanges))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(IRanges))
suppressPackageStartupMessages(library(GenomicScores))
suppressPackageStartupMessages(library(parallel))



####################
# phastCons (no strand info)
####################
#gr <- GenomicRanges::makeGRangesFromDataFrame(df = bed2[674:676,7:18], seqnames.field = "gn", start.field = "gn_start", end.field = "gn_end", keep.extra.columns = TRUE)
rtrack_obj <- rtracklayer::import.bed(inBed12)
#rtrack_obj2 <- rtracklayer::import.bed(inBed.8DNA.gn)
#rtrack_obj <- c(rtrack_obj1,rtrack_obj2)
#tmp <- rtracklayer::asBED(gr)
tmp <-  as.data.frame(rtracklayer::blocks(rtrack_obj))
#table(duplicated(tmp$group_name))
#ids <- tmp$group_name[duplicated(tmp$group_name)]
#tmp <- tmp[tmp$group_name %in% ids,]
tmp$pos <- paste0(tmp$seqnames,":",tmp$start,"-",tmp$end)
#tmp$phastCons100way <- sapply(X = tmp$pos, FUN = phastCons100way)

if(cores==-1){
  cores <- as.integer(0.5*detectCores())
  message("use half of auto-detect cores:",cores)
}else{
  message("use predefined cores:",cores)
}
cl <- makeCluster(cores) 
clusterExport(cl=cl, varlist=c("tmp"))
matrix_of_sums <- parLapply(cl, 1:nrow(tmp), function(i) 
  #i <- 1
  #By default, parLapply does not guarantee any specific order of the results.
  tryCatch(
    {
      #mean: pos <- "chr1:1000000-1000010"
      #single base: pos <- c("chr1:1000000-1000001", "chr1:1000001-1000002")
      GenomicScores::gscores(phastCons100way.UCSC.hg38::phastCons100way.UCSC.hg38, GenomicRanges::GRanges(tmp$pos[i]))$default # summaryFunFunction: mean
    },
    error = function(cond) {
      message(paste("non-auto chr positions not suport in phastCons100way.UCSC.hg38: ", tmp$pos[i]))
      #non-auto chr not suport!
      # Choose a return value in case of error
      -1
    }
  )
)


#i <- 1
stopCluster(cl)
tmp$phastCons100way <- (Reduce("cbind", matrix_of_sums)[1,])
tmp$phastCons100way[tmp$phastCons100way==-1] <- NA
#380 peaks: parLapply: 1.2 min
#380 peaks: sapply: 4.6 min

tmp2 <- tmp %>% 
  group_by(group_name) %>%  # seqnames,start,end,group_name,strand
  summarise(phastCons100way=mean(phastCons100way,na.rm=T))
colnames(tmp2)[1] <- "name"


# keep originial order
#old.order <- paste(rtrack_obj$name, rtrack_obj@seqnames, rtrack_obj@ranges@start, rtrack_obj@ranges@start+rtrack_obj@ranges@width, rtrack_obj@strand)
tmp2 <- tmp2[match(rtrack_obj$name,tmp2$name),]

# write table
write.table(tmp2,outFile,quote = F,sep = "\t",row.names = F,col.names = T)