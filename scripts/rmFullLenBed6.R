#! /usr/bin/env Rscript

# rm (near) full-length bed6
# last 2407 by bpf 
# b.p.f@qq.com


suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser(description='intersect domain bed with G4-iM, RBPs, AGO2, RBPhotspot')
parser$add_argument('-i', '--inputFile', type='character', default='-',
                    help='input domain BED file') #, default: -
parser$add_argument('-o', '--outputFile', type='character', default='-',
                    help='output BED file')
parser$add_argument('--chrSize', type='character', default="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/transcriptome_genome_sort_uniq_newTxID", 
                    help='chrSize file path, used for shuffle random regions, contains both 11+8 tx and 23 chr. default: /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/transcriptome_genome_sort_uniq_newTxID')
parser$add_argument('--deltaLen', type="integer", default=5,
                    help='rm those length less than a given integer relative full-length tx, default: 5')
parser$add_argument('--rename', type="logical", default=F,
                    help='wheather to rename left records by order, default: F')
args <- parser$parse_args()

#todo: add param for shuf in all tx of same RNA type

for(i in 1:length(args)){
  v.name <- names(args)[i]
  assign(v.name,args[[i]])
  message(paste0(v.name,": ",get(v.name)))
}


chrSize <-  data.table::fread(chrSize) # contains both chr, ENST
domain <- as.data.frame(data.table::fread(inputFile,data.table = F,sep = '\t',check.names = F,stringsAsFactors = F))
message(paste0(dim(domain),collapse = "\t"))
tx <- chrSize$V2[match(domain$V1,chrSize$V1)]
len <- as.integer(domain$V3)-as.integer(domain$V2)
domain <- domain[len<=tx-deltaLen,]
message(paste0(dim(domain),collapse = "\t"))
if(rename){
  message("renaming.")
  domain$V4 <- paste0("peak_",1:nrow(domain))
}
data.table::fwrite(domain,outputFile,quote = F,sep = '\t',row.names = F,col.names = F)

