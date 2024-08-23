#! /usr/bin/env Rscript

# parse different forms of peak caller results into bed6, 5th score is >0 interger and represents quality
# last 2408 by bpf 
# b.p.f@qq.com


suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser(description='intersect domain bed with G4-iM, RBPs, AGO2, RBPhotspot')
parser$add_argument('-i', '--inputFile', type='character', # default='-',
                    help='input domain BED file') #, default: -
parser$add_argument('-o', '--outputFile', type='character',#  default='-',
                    help='output BED file')
parser$add_argument('-t', '--tool', type='character',  required=T,
                    help='opts: cfpeak, cfpeakcnn, clipper, clam, piranha')
parser$add_argument('-r', '--rename', type="logical", default=F,
                    help='wheather to rename left records by order, default: F')
args <- parser$parse_args()

#todo: add param for shuf in all tx of same RNA type

for(i in 1:length(args)){
  v.name <- names(args)[i]
  assign(v.name,args[[i]])
  message(paste0(v.name,": ",get(v.name)))
}



## functions

extractClamScore <- function(textVec) {
  # Extract the numeric values using a regular expression
  # textVec <- domain$V4[1:3]
  numeric_values <- as.numeric(sub(".*:(.*)", "\\1", textVec)) # only last values if multi exist
  # min_value <- min(numeric_values, na.rm = TRUE) #
  return(numeric_values)
}
#


# ### test blockbuster
# inputFile="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/TCGA-COAD_small_NC_NCpool/call_peak_gbamStarEM/blockbuster_by_sample/min3/NCpool_15_filter_block.bed"
# outputFile="test.bed"
# tool="blockbuster"
# rename=T
# ### test piranha
# inputFile="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/TCGA-COAD_small_NC_NCpool/call_peak_gbamStar/piranha_by_sample/b5_p01/NCpool_15_filter.bed"
# outputFile="test.bed"
# tool="piranha"
# rename=T
# ### test clipper
# inputFile="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/TCGA-COAD_small_NC_NCpool/call_peak_gbamStar/clipper_by_sample/b5_p05/NCpool_15_filter.bed"
# outputFile="test.bed"
# tool="clipper"
# rename=T
# ### test clam
# inputFile="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/TCGA-COAD_small_NC_NCpool/call_peak_gbamStarEM/clam_by_sample/b5_p005/NCpool_15_filter.bed"
# outputFile="test.bed"
# tool="clam"
# rename=T
# ### test cfpeak
# inputFile="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/TCGA-COAD_small_NC_NCpool/call_peak_all/cfpeak_by_sample/b5_d50_p1/NCpool_15_filter.bed"
# outputFile="test.bed"
# tool="cfpeak"
# rename=T
# ### test cfpeakCNN
# inputFile="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/TCGA-COAD_small_NC_NCpool/call_peak_all/cfpeakCNN_by_sample/b5_d50_p1/NCpool_15_filter.bed"
# outputFile="test.bed"
# tool="cfpeakCNN"
# rename=T


#domain <- as.data.frame(data.table::fread(inputFile,data.table = F,sep = '\t',check.names = F,stringsAsFactors = F))
domain <- (read.table(inputFile, sep = '\t',check.names = F,stringsAsFactors = F))

if(tolower(tool)=="blockbuster"){
  message("format: ",tool)
  colnames(domain) <- c("chr","start",'end','name','score','strand')
}else if(tolower(tool)=="piranha"){
  message("format: ",tool)
  colnames(domain) <- c("chr","start",'end','name','score','strand')
}else if(tolower(tool)=="clipper"){
  message("format: ",tool)
  domain <- domain[,1:6]
  colnames(domain) <- c("chr","start",'end','name','score','strand')
  minScore <- domain$score
  minScore <- minScore[minScore>0]
  minScore <- min(minScore)
  domain$score[domain$score==0] <- minScore
  domain$score <- -log10(domain$score) # [0,]
  domain$score <- round(domain$score,digits = 0)+1 # [1,]
  # summary(domain$score)
}else if(tolower(tool)=="clam"){
  message("format: ",tool)
  colnames(domain) <- c("chr","start",'end','name','score','strand')
  domain$score <- extractClamScore(domain$name)
  domain$score <- -log10(domain$score) # [0,]
  domain$score <- round(domain$score,digits = 0)+1 # [1,]
}else if(tolower(tool)=="cfpeak"){
  message("format: ",tool)
  domain <- domain[,1:6]
  colnames(domain) <- c("chr","start",'end','name','score','strand')
}else if(tolower(tool)=="cfpeakcnn"){
  message("format: ",tool)
  colnames(domain) <- c("chr","start",'end','name','score','strand')
  # summary(domain$score)
  minScore <- domain$score
  minScore <- minScore[minScore>0]
  minScore <- min(minScore)
  maxScore <- domain$score
  maxScore <- maxScore[maxScore>0]
  maxScore <- max(maxScore)
  # table(domain$score==0)
  domain$score[domain$score==0] <- minScore
  domain$score[domain$score<0] <- maxScore
  domain$score <- -log10(domain$score) # [0,]
  domain$score <- round(domain$score,digits = 0)+1 # [1,]
}else{
  stop("not valid input")
}



if(rename){
  message("renaming.")
  domain$name <- paste0(tolower(tool),"_",1:nrow(domain))
}

#data.table::fwrite(domain,outputFile,quote = F,sep = '\t',row.names = F,col.names = F)
write.table(domain,outputFile,quote = F,sep = '\t',row.names = F,col.names = F)

