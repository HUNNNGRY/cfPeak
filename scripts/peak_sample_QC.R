# quality control for cfRNA peak analysis
# last 2208 by bpf
# b.p.f@qq.com
# todo: add duplication ratio as metric, for call peak will be affected by whether or not to dedup

suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser(description='sample map stat QC for peak')
parser$add_argument('-i', '--inputLogDir', type='character', # default='-',
                    help='input bowtie2 tmap dir: /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/FTC_small/tbam')
parser$add_argument('-o', '--outputFile', type='character', #default='-',
                    help='output map stats txt file with QC status: /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/FTC_small/qc.txt')
args <- parser$parse_args()

for(i in 1:length(args)){
  v.name <- names(args)[i]
  assign(v.name,args[[i]])
  message(paste0(v.name,": ",get(v.name)))
}
inputLogDir <- args$inputLogDir
outputFile <- args$outputFile

## test
#inputLogDir <- "../output/GSE71008_diff/tbam"
#outputFile <- "data/GSE71008_diff/qc.txt" 


# small/SE map ratio (main QC module) ------------------------------------------------------
read.SE.map <- function(x){
  #x <- "output/FTC_small/log/tbam/csFTA-10_1/enhancer.for.log"
  #x <- "output/FTC_long/log/clFTA-10/tbam/enhancer.for.log"
  #x <- "output/WSQ_SMARTer_NEB/tbam/NC_ChQ-21_smart_1/log/promoter_for.log"
  # print(x)
  sample <- unlist(lapply(strsplit(x,"/",fixed = T), function(y) tail(y,3)[1]))
  rna <- gsub(".log","",basename(x)) 
  map0 <- readLines(x)
  in.num <- as.numeric(gsub(" reads; of these:","",map0[grepl("reads; of these:",map0)]))
  not.map.num <- as.numeric(unlist(sapply(strsplit(map0[grepl("aligned 0 times",map0)]," (", fixed = T), "[" ,1)))
  uniq.map.num <- as.numeric(unlist(sapply(strsplit(map0[grepl("aligned exactly 1 time",map0)]," (", fixed = T), "[" ,1)))
  multi.map.num <- as.numeric(unlist(sapply(strsplit(map0[grepl("aligned >1 time",map0)]," (", fixed = T), "[" ,1)))
  #map.num <- uniq.map.num + multi.map.num
  #in.num==not.map.num+uniq.map.num+multi.map.num # TURE
  map <- as.numeric(gsub("% overall alignment rate","",map0[6]))
  map <- data.frame(sample=sample,type=rna,overall.align.ratio=map,in.num=in.num, uniq.map.num=uniq.map.num, multi.map.num=multi.map.num,not.map.num=not.map.num)
  return(map)
}

map.df <- lapply(Sys.glob( paste0(inputLogDir,"/*/log/*.log") ),read.SE.map)
map.df <- do.call(rbind,map.df)
map.df$overall.align.ratio <- map.df$overall.align.ratio/100
#table(map.df$type)
library(tidyr)
#map.df.a <- map.df %>% tidyr::pivot_wider(names_from = type, values_from = value)
#map.df.a2 <- map.df.a
#map.df.a2[,2:ncol(map.df.a2)] <- 1-map.df.a2[,2:ncol(map.df.a2)]
#map.df.a$unmap <- 1-rowSums(map.df.a[,2:ncol(map.df.a)])
#rowSums(map.df.a[,2:ncol(map.df.a)])
#map.df2 <- map.df.a %>% tidyr::pivot_longer(cols = 2:ncol(map.df.a), names_to = "type", values_to = "value")

dat <- as_tibble(map.df)
dat <- dat %>% dplyr::group_by(sample) %>% 
  dplyr::mutate(map.num=uniq.map.num+multi.map.num, 
                total.uniq.map.num=sum(uniq.map.num),
                all.num=max(in.num),
                unmap.num=min(not.map.num),
                map.ratio=map.num/all.num,
                unmap.ratio=unmap.num/all.num) %>% 
  dplyr::select(sample,type,all.num,map.ratio,total.uniq.map.num,unmap.ratio)
dat.wide <- dat %>% tidyr::pivot_wider(names_from = "type", values_from = "map.ratio")
colnames(dat.wide)[1:4] <- c("sample","read_num","uniq_num","unmap")
colnames(dat.wide) <- gsub(".","_",colnames(dat.wide),fixed = T)
colnames(dat.wide) <- gsub("miRNA","pri_miRNA",colnames(dat.wide),fixed = T)
colnames(dat.wide) <- gsub("pri_pri_miRNA","pri_miRNA",colnames(dat.wide),fixed = T)
dat.wide <-  dat.wide[,c("sample","read_num","uniq_num","unmap","rRNA","spikein","univec","pri_miRNA","lncRNA", "mRNA", "piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA', "intron_for", "intron_rev", "promoter_for", "promoter_rev", "enhancer_for", "enhancer_rev", "repeats_for", "repeats_rev")]
rowSums(dat.wide[,4:ncol(dat.wide)]) # check: should be all ones !!!

fqNum_cutoff <- 5000000
spikeinRatio_cutoff <- 0.3
univecRatio_cutoff <- 0.3
txUnmapRatio_cutoff <- 0.9
txUniqMapNum_cutoff <- fqNum_cutoff*(1-txUnmapRatio_cutoff)*0.1 # aver 10% mapped are unique mapped
dat.wide$passQC <- FALSE
dat.wide$passQC[dat.wide$read_num>=txUnmapRatio_cutoff & 
                  dat.wide$uniq_num>=txUniqMapNum_cutoff &
                  dat.wide$unmap<=txUnmapRatio_cutoff & 
                  dat.wide$spikein<=spikeinRatio_cutoff &
                  dat.wide$univec<=univecRatio_cutoff] <- TRUE
#dat.wide$QCrank <- order(dat.wide$unmap,dat.wide$univec,dat.wide$spikein,-dat.wide$uniq_num,decreasing = F) # smaller is better
#order(c(1,3,4,2,3),c(1,3,4,2,4),decreasing = F)





# small/SE (miR) dedup ratio (supp QC module, need pre-run bam-deduped-samtools)  ------------------------------------------------------

# # need run below in bash
# dst="GSE110381_diff" # GSE110381_diff, GSE71008_diff
# pre="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
# cores=8
# for i in `cat $pre/exSeek-dev/data/$dst/sample_ids.txt` #  | tail -n119
# do
# echo $i;
# mkdir -p $pre/output/$dst/tbam/$i/bam-deduped-samtools/log/;
# for rna in rRNA miRNA piRNA mRNA repeats_for
# do
# echo $rna;
# #samtools rmdup -s
# /BioII/lulab_b/baopengfei/anaconda3/envs/exvariance/bin/samtools markdup -s -r -m t -@ ${cores} \
#   $pre/output/$dst/tbam/$i/bam/$rna.bam \
#   $pre/output/$dst/tbam/$i/bam-deduped-samtools/$rna.bam \
#   >  $pre/output/$dst/tbam/$i/bam-deduped-samtools/log/$rna.log 2>&1
# done
# done
# ##try picard?

read.SE.dedup <- function(x){
  #x <- "../output/GSE71008_diff/tbam/SAMN03863387/bam-deduped-samtools/log/miRNA.log"
  # print(x)
  sample <- unlist(lapply(strsplit(x,"/",fixed = T), function(y) tail(y,4)[1]))
  rna <- gsub(".log","",basename(x))
  dedup.stat <- readLines(x) # samtools rmdup
  before.dup <- dedup.stat[grepl("^EXAMINED:",dedup.stat,perl = T)]
  before.dup <- as.numeric(gsub("EXAMINED:","",before.dup))
  after.dup <- dedup.stat[grepl("^DUPLICATE TOTAL:",dedup.stat,perl = T)]
  after.dup <- as.numeric(gsub("DUPLICATE TOTAL:","",after.dup))
  dup.ratio <- after.dup/before.dup
  dedup.mir.num <- before.dup-after.dup # miR lib size ?
  map <- data.frame(sample=sample,mir_dup_ratio=dup.ratio,mir_rmdup_num=dedup.mir.num) #type=rna,
  return(map)
}
#
dedup.df <- lapply( Sys.glob( paste0(inputLogDir,"/*/bam-deduped-samtools/log/miRNA.log") ) ,read.SE.dedup)
dedup.df <- do.call(rbind,dedup.df)




# sum mat ------------------------------------------------------
dat.wide <- merge(dat.wide,dedup.df,by="sample")
dat.wide$QCrankSum <- order(dat.wide$unmap,decreasing = F) + order(dat.wide$univec,-dat.wide$uniq_num,decreasing = F) + order(-dat.wide$uniq_num,decreasing = F) + order(dat.wide$mir_dup_ratio,decreasing = F) + order(-dat.wide$mir_rmdup_num,decreasing = F) #before.dup-after.dup

message(paste0("total: ",nrow(dat.wide)))
message(paste0("passQC: ",sum(dat.wide$passQC)))

write.table(x = dat.wide, file = args$outputFile, quote = F, sep = "\t", row.names = F, col.names = T)

# cfRNA pico:
# All sample	Condition	247
# Raw Filter	10000000.00 	246
# Clean Filter	5000000.00 	246
# Spike in 	0.50 	246
# Genome Aligned Reads	500000.00 	241
# Long RNA ratio	0.20 	239
# Unclassified	0.30 	228
# Intron Spanning	100000.00 	228

# cf small:
# (1) datasets are required to have at least 100,000 reads that overlap with any annotated RNA transcript in the host genome
# (2) over 50% of the reads that map to the host genome also align to any RNA annotation.



# # sum qc
# dst <- "GSE110381_diff"
# dst.NC <- "GSE110381"
# qc <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/",dst,"/qc.txt")
# qc2 <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/",dst,"/sample_table.txt")
# qc3 <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/GSE110381/SraRunTable.txt"
# passQC <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/",dst.NC,"/sample_ids_NCpool_test15.txt")
# #diff.passQC <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/",dst.NC,"/sample_ids_diff_passQC.txt")
# diff.passQC2 <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/",dst.NC,"/sample_ids_diff_passQC2.txt")
# diff.passQC.batch <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/",dst.NC,"/sample_ids_diff_passQC2_rmBatch.txt")

#
# df <- read.table(qc,header = T,check.names = F)
# df2 <- read.table(qc2,header = T,check.names = F)
# df <- df[df$sample %in% df2$sample[df2$group=="NC"],]
#
# table(df$uniq_num>=500000)
# table(df$unmap<=.3)
# table(df$mir_rmdup_num>=2500)
#
# df <- df[df$uniq_num>=800000 & df$unmap<=.3 & df$mir_rmdup_num>=2800,]
# write.table(file = passQC, x = df$sample, quote = F, sep = "\t", row.names = F, col.names = F)

# 
# 
# 
# 
# df <- read.table(qc,header = T,check.names = F)
# df2 <- read.table(qc2,header = T,check.names = F)
# #df <- df[df$sample %in% df2$sample[df2$group=="NC"],]
# table(df$uniq_num>=500000)
# table(df$unmap<=.4)
# table(df$mir_rmdup_num>=2500)
# df <- df[df$uniq_num>=500000 & df$unmap<=.4 & df$mir_rmdup_num>=2500,]
# write.table(file = diff.passQC2, x = df$sample, quote = F, sep = "\t", row.names = F, col.names = F)
# table(df2$group[match(df$sample,df2$sample)])
#
#
#
# df <- read.table(qc,header = T,check.names = F)
# df2 <- read.table(qc2,header = T,check.names = F)
# sra <- read.table(qc3,check.names = F,sep = "\t",header = T)
# table(df2$sample %in% sra$Run)
# table(sra$Run %in% df2$sample)
# sra2 <- sra[sra$Run %in% df2$sample,]
# table(sra$Sequencing_replicate) # seems Hongke/me? only downloaded samples with seq_rep==A
# #hist(sra2$Bases,breaks = 100, xlim = c(0,3e+09))
# hist(sra2$Age) # 1.filter AGE: !!!!!!!!!!!!!
# table(sra2$Age>=45 & sra2$Age<=70)
# #FALSE  TRUE
# #46   186
# sra2 <- sra2[sra2$Age>=45 & sra2$Age<=70,]
# table(sra2$`self-identified_ancestry`)
# sra2 <- sra2[sra2$`self-identified_ancestry` %in% c("AA","C"),] # 2.filter ancestry !!!!!!!!
# table(sra2$finding)
# 
# df2 <- df2[df2$sample %in% sra2$Run,]
# 
# passQCfile <- read.table(diff.passQC2,header = F)$V1
# table(df2$sample %in% passQCfile)
# df2 <- df2[df2$sample %in% passQCfile,]
# table(df2$group)
# write.table(file = diff.passQC.batch, x = df2$sample, quote = F, sep = "\t", row.names = F, col.names = F)





# # sum qc
# dst <- "GSE71008_diff" # seem no need to qc 
# dst.NC <- "GSE71008"
# qc <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/",dst,"/qc.txt")
# qc2 <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/",dst,"/sample_table.txt")
# passQC <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/",dst.NC,"/sample_ids_NCpool_test15.txt2")
# df <- read.table(qc,header = T,check.names = F)
# df2 <- read.table(qc2,header = T,check.names = F)
# df <- df[df$sample %in% rownames(df2)[df2$disease_type=="Control"],]
# 
# table(df$uniq_num>=500000)
# table(df$unmap<=.3)
# table(df$mir_rmdup_num>=2000)
# 
# df <- df[df$uniq_num>=500000 & df$unmap<=.3 & df$mir_rmdup_num>=2800,]
# write.table(file = passQC, x = df$sample, quote = F, sep = "\t", row.names = F, col.names = F)
