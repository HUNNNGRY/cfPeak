# transmit form local to bioii at 220604


# add miR ENST ------------------------------------------------------------
#setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/")
setwd("/Users/baopengfei/Desktop/lulab/tmp/projects/WCHSU-FTC/exSeek-dev/")
options(stringsAsFactors = F)

lncrna <- read.table("./genome/hg38/transcript_table/lncRNA.txt", header = T)
lncrna[1:3,]
colnames(lncrna)
summary(lncrna$score)

tmp <- read.table("genome/hg38/bed_by_biotype/miRNA.bed", header = F)
colnames(tmp) <- c("chrom","start","end","name","score","strand","thickStart","thichEnd","itemRgb","blockCount","blockSizes","blcokStarts")
tmp <- tmp[,c("chrom","start","end","name","score","strand")]
rownames(tmp) <- tmp$name
head(tmp$name)

## read miR gtf
gtf <- read.table("genome/hg38/gtf_by_biotype/miRNA.gtf",header = F,sep = "\t")
gtf <- gtf[gtf$V3=="transcript",]
gtf[1:3,]
l <- c("gene_id","transcript_id","gene_type","gene_name")
for (i in 1:length(l)){
  #i <- 5
  gtf[[l[i]]] <- unlist(sapply(strsplit(gtf$V9,";",fixed = T),"[",i))
  #gtf[[l[i]]] <- unlist(sapply(strsplit(gtf[[l[i]]]," "),"[",2))
  gtf[[l[i]]] <- gsub("gene_id |transcript_id |gene_type |gene_name |transcript_name ","",gtf[[l[i]]],perl = T)
}
gtf$transcript_name <- unlist(sapply(strsplit(gtf$V9,";",fixed = T),"[",6))
gtf$transcript_name <- gsub("transcript_name ","",gtf$transcript_name)
gtf$gene_type <- "pri_miRNA"
gtf$transcript_type <- "pri_miRNA"
gtf$source <- gtf$V2
gtf <- gtf[,colnames(gtf) %in% colnames(lncrna)]
#head((gtf$transcript_name)) # transcript_id, gene_name, transcript_name
gtf$transcript_id <- gsub(" ","",gtf$transcript_id)
gtf$gene_name <- gsub(" ","",gtf$gene_name)
gtf$transcript_name <- gsub(" ","",gtf$transcript_name)
#rownames(gtf) <- gtf$gene_id
rownames(gtf) <- gtf$transcript_id
gtf <- gtf[rownames(tmp),]
#table(rownames(gtf) %in% rownames(tmp)) 
#table(rownames(tmp) %in% rownames(gtf)) 

#head(tmp,3)
all(gtf$transcript_id == tmp$name) # should be all TRUE !
tmp <- cbind(tmp,gtf)
tmp <- tmp[,colnames(lncrna)]
write.table(tmp, "genome/hg38/transcript_table/pri_miRNA.txt",sep="\t",quote=F,row.names=F,col.names=T)
tmp <- tmp[,c(4,2,3,5,6)]
tmp$len <- tmp$end-tmp$start
tmp$start <- 0
tmp$end <- NULL
tmp$ID <- "."
tmp$strand <- "+"
write.table(tmp[,c(1,2,5,6,3,4)], "genome/hg38/tbed/pri_miRNA.bed",sep="\t",quote=F,row.names=F,col.names=F)
#head(tmp)
#??read.gtf






# prepare exseek tx.table: add DNA elements -----------------------------------------------------------------
#setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/")
setwd("/Users/baopengfei/Desktop/lulab/tmp/projects/WCHSU-FTC/exSeek-dev/")
options(stringsAsFactors = F)

lncrna <- read.table("./genome/hg38/transcript_table/lncRNA.txt", header = T)
lncrna[1:3,]

dat <- list()
for (i in c("promoter.rev","promoter.for","enhancer.rev","enhancer.for","repeats.rev","repeats.for","intron.rev","intron.for")){
  #i <- "intron.for"
  print(i)
  x <- paste0("./genome/hg38/bed/",i,".bed")
  tmp <- read.table(x, header = F)
  l <- gsub(".bed","",basename(x))
  tmp[1:3,]
  
  colnames(tmp) <- c("chrom", "start","end","gene_type","score","strand")
  tmp$name <- paste0(tmp$gene_type,"::",tmp$chrom,":",tmp$start,"-",tmp$end,"(",tmp$strand,")")
  tmp$gene_id <- tmp$name
  tmp$transcript_id <-  tmp$name
  tmp$gene_name <- tmp$name
  tmp$transcript_name <- tmp$name
  tmp$transcript_type <- l
  tmp$gene_type <- l
  tmp$source <- "encodev27"
  tmp <- tmp[,colnames(lncrna)]
  dat[[i]] <- tmp
  #head(tmp)
  #print( hist((as.numeric(tmp$end)-as.numeric(tmp$start)),breaks = 1000000,xlim = c(0,10000)) )
  #return(tmp)
}


head(dat[[7]])
dat.df <- as.data.frame(do.call("rbind",dat))

all <- read.table("genome/hg38/transcript_table/all_old2.txt", header = T)
all[1:3,]
all <- all[all$chrom!="chrom",]
all <- rbind(all,dat.df)
tail(all)
#all[all$transcript_id=="(AATGG)n::chr16:13155048-13155369(+)",]
write.table(all, "genome/hg38/transcript_table/all.txt",sep="\t",quote=F,row.names=F,col.names=T)






# gbed2tbed ---------------------------------------------------------------
# get tbed (10 rna types
#miRNA,mRNA,lncRNA,piRNA,snoRNA,snRNA,srpRNA,tRNA,tucpRNA,Y_RNA  # srpRNA has one record duplicated with mRNA
#transcript_table/all.txt中多数rna类型都是redundant，需要过滤只选择longest的ID，也就是map的ID （tbed/*RNA.bed）

fs <- Sys.glob("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/tbed/*RNA.bed") # should be 10 ?
chrSize.df <- list()
read.RNA <- function(x){
  #x <- fs[1]
  print(x)
  tmp.bed <- data.table::fread(x,data.table = F,header = F,sep = "\t",stringsAsFactors = F)
  type <- basename(x)
  tmp.bed$type <- gsub(".bed","",type)
  chrSize.df[[type]] <- tmp.bed
  #tmp.bed[1:3,]
}
chrSize.df <- do.call(rbind,lapply(fs,read.RNA))
#table(duplicated(chrSize.df$V1))
chrSize.df <- chrSize.df[!duplicated(chrSize.df$V1),]
data.table::fwrite(chrSize.df[,1:6],"/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/tbed/11RNA_map.bed",sep = "\t",quote = F,col.names = F,row.names = F)
data.table::fwrite(chrSize.df,"/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/tbed/11RNA_map.txt",sep = "\t",quote = F,col.names = F,row.names = F)

chrSize.df2 <- chrSize.df
chrSize.df2$V1 <- gsub("(-)","_neg",chrSize.df2$V1,fixed = T)
chrSize.df2$V1 <- gsub("(+)","_pos",chrSize.df2$V1,fixed = T)
chrSize.df2$V1 <- gsub("::","__",chrSize.df2$V1,fixed = T)
chrSize.df2$V1 <- gsub(":","___",chrSize.df2$V1,fixed = T)
chrSize.df2$V1 <- gsub("-","____",chrSize.df2$V1,fixed = T)
chrSize.df2$V1 <- gsub(".","_____",chrSize.df2$V1,fixed = T)
chrSize.df2[1:3,]
data.table::fwrite(chrSize.df2,"/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/tbed/11RNA_map_newTxID.txt",sep = "\t",quote = F,col.names = F,row.names = F)


long_RNA <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/bed/long_RNA.bed",data.table = F,header = F,sep = "\t",stringsAsFactors = F)
long_RNA[1:2,]
long_RNA2 <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/bed/tRNA.bed",data.table = F,header = F,sep = "\t",stringsAsFactors = F)
long_RNA2[1:2,]
long_RNA3 <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/bed/pri_miRNA.bed",data.table = F,header = F,sep = "\t",stringsAsFactors = F)
long_RNA3[1:2,]
long_RNA4 <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/bed/piRNA.bed",data.table = F,header = F,sep = "\t",stringsAsFactors = F)
long_RNA4[1:2,]
long_RNA5 <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/bed/rRNA.bed",data.table = F,header = F,sep = "\t",stringsAsFactors = F)
long_RNA5[1:2,]
#exSeek-dev/genome/hg38/bed/{long_RNA,tRNA,pri_miRNA,piRNA}_newTxID.bed

long_RNA <- rbind(long_RNA,long_RNA2,long_RNA3,long_RNA4,long_RNA5)
#table(duplicated(long_RNA$V4)) # 108 dup
#head(long_RNA[duplicated(long_RNA$V4),])
#long_RNA[long_RNA$V4=="T018925",]
long_RNA <- long_RNA[!duplicated(long_RNA$V4),]
#table(long_RNA$V4 %in% chrSize.df$V1)
#table(chrSize.df$V1 %in% long_RNA$V4)
RNAs <- c("rRNA","pri_miRNA","lncRNA","mRNA","piRNA","snoRNA","snRNA","srpRNA","tRNA","tucpRNA","Y_RNA")
long_RNA <- long_RNA[long_RNA$V4 %in% chrSize.df$V1,]
data.table::fwrite(long_RNA,"/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/bed/11RNA_map_bed12.bed",sep = "\t",quote = F,col.names = F,row.names = F)

long_RNA$V4 <- gsub("(-)","_neg",long_RNA$V4,fixed = T)
long_RNA$V4 <- gsub("(+)","_pos",long_RNA$V4,fixed = T)
long_RNA$V4 <- gsub("::","__",long_RNA$V4,fixed = T)
long_RNA$V4 <- gsub(":","___",long_RNA$V4,fixed = T)
long_RNA$V4 <- gsub("-","____",long_RNA$V4,fixed = T)
long_RNA$V4 <- gsub(".","_____",long_RNA$V4,fixed = T)
long_RNA[1:2,]
data.table::fwrite(long_RNA,"/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/bed/11RNA_map_bed12_newTxID.bed",sep = "\t",quote = F,col.names = F,row.names = F)


# get tbed block df
#only consider longest tx records
# tx.ref <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/bed/long_DNA.bed",data.table = F,header = F,sep = "\t",stringsAsFactors = F)
# tx.ref <- tx.ref[!grepl("alt|random|chrUn",tx.ref$V1,perl=T),] # filter non-canonical chr
# tx.ref[1:3,]
# str(tx.ref)

# ## create bed12 for rRNA (failed)
# id.conv <- read.table("genome/hg38/source/refSeq_rRNA.gene_names.txt")
# colnames(id.conv) <- c("ref_seq_id","gene_name")
# 
# rrna <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/bed/11RNA_map_bed12.bed")
# rrna <- rrna[rrna$V4 %in% id.conv$ref_seq_id,]
# write.table(rrna,file = "rRNA_map_bed12.bed",quote = F,sep = "\t",row.names = F,col.names = F)

tx.ref <- rtracklayer::import.bed(file('/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/bed/11RNA_map_bed12.bed'))
#table(tx.ref$name %in% chrSize.df$transcript_id)
#tx.ref <- tx.ref[tx.ref$name %in% chrSize.df$transcript_id,]
#tx.ref <- tx.ref[sample(1:length(tx.ref),0.01*length(tx.ref),replace = F)] # test small sample
head(tx.ref)
tx.ref.block0 <- rtracklayer::blocks(tx.ref)
#head(tx.ref.block)

tx.ref.block.df <- list()
#for (i in 1:length(tx.ref.block0)){
getBlock <- function(i){
  if(i%%1000 == 0){
    print(i)
  }
  enst <- names(tx.ref.block0)[i]
  tmp <- as.data.frame(tx.ref.block0[[i]])
  tmp$enst <- enst
  #tmp$length <- tmp$end-tmp$start+1 # still 1 based GRanges
  tmp$total <- length(tx.ref.block0[[i]])
  tmp$nth <- 1:length(tx.ref.block0[[i]])
  tx.ref.block.df[[i]] <- tmp
}
library(parallel)
cores <- 6
tx.ref.block.df <- mclapply(1:length(tx.ref.block0),getBlock,mc.cores = cores)
#tx.ref.block.df <- lapply(1:length(tx.ref.block0),getBlock)

tx.ref.block <- do.call(rbind,tx.ref.block.df)
tx.ref.block.df <- NULL
tx.ref.block <- as.data.frame(tx.ref.block,stringsAsFactors=F) #rtracklayer::blocks
tx.ref.block[1:3,]
#head(tx.ref.block,3)
tx.ref.block$start <- tx.ref.block$start-1
tx.ref.block$end <- tx.ref.block$start + tx.ref.block$width
tx.ref.block$score <- "."
tx.ref.block <- tx.ref.block[,c("seqnames","start","end","enst","score","strand","total","nth")] # ,"group"
tx.ref.block$strand <- as.character(tx.ref.block$strand)
tx.ref.block$seqnames <- as.character(tx.ref.block$seqnames)
data.table::fwrite(tx.ref.block,"/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/bed/10RNA_map_block.txt",sep = "\t",quote = F,col.names = T,row.names = F)
tx.ref.block[1:3,]

tx.ref.block$enst <- gsub("(-)","_neg",tx.ref.block$enst,fixed = T)
tx.ref.block$enst <- gsub("(+)","_pos",tx.ref.block$enst,fixed = T)
tx.ref.block$enst <- gsub("::","__",tx.ref.block$enst,fixed = T)
tx.ref.block$enst <- gsub(":","___",tx.ref.block$enst,fixed = T)
tx.ref.block$enst <- gsub("-","____",tx.ref.block$enst,fixed = T)
tx.ref.block$enst <- gsub(".","_____",tx.ref.block$enst,fixed = T)
data.table::fwrite(tx.ref.block,"/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/bed/10RNA_map_block_newTxID.txt",sep = "\t",quote = F,col.names = T,row.names = F)




# add hg38_tx reference to CLIPper ----------------------------------------
#setwd("/Users/baopengfei/Desktop/lulab/tmp/projects/WCHSU-FTC/exSeek-dev/")
setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/")
options(stringsAsFactors = F)

library(tidyverse)
library(dplyr)

# read ref 
tx <- data.table::fread("genome/hg38/chrom_sizes/tx_gn_length.txt",stringsAsFactors = F,sep = "\t",header = T)
table(tx$transcript_type)
table(duplicated(tx$transcript_id))
tx <- tx[!duplicated(tx$transcript_id),]
tx <- tx %>%
  # dplyr::filter(!(transcript_type %in% c("rRNA"))) %>% 
  dplyr::select(transcript_id,transcript_type)
table(tx$transcript_type) # no repeats
length(unique(tx$transcript_type)) # 18,19 (rRNA)
# enhancer.for enhancer.rev   intron.for   intron.rev       lncRNA         mRNA        piRNA    pri_miRNA promoter.for promoter.rev  repeats.for  repeats.rev 
# 397307       397307         4211         5707        27908       146890        23431         3762        89275        89275      3158443      2361574 
# snoRNA        snRNA       srpRNA         tRNA      tucpRNA        Y_RNA 
# 955         1900          679          649        11550          756 

# read all bed: 11 RNAs + 8 DNAs 
#ENST00000607096
bed <- data.table::fread("genome/hg38/tbed/all.bed",stringsAsFactors = F,sep = "\t",header = F)
pi.bed <- data.table::fread("genome/hg38/tbed/piRNA.bed",stringsAsFactors = F,sep = "\t",header = F)
pi.bed[1:3,]
mi.bed <- data.table::fread("genome/hg38/tbed/pri_miRNA.bed",stringsAsFactors = F,sep = "\t",header = F)
mi.bed[1:3,]
r.bed <- data.table::fread("genome/hg38/tbed/rRNA.bed",stringsAsFactors = F,sep = "\t",header = F)
r.bed[1:3,]
bed <- rbind(as.data.frame(bed),as.data.frame(pi.bed),as.data.frame(mi.bed),as.data.frame(r.bed))
bed[1:3,]
#grep("ENST00000607096",bed$V1)
#nrow(bed)==nrow(tx) # should be TRUE ?
table(bed$V1 %in% tx$transcript_id)

#table(tx$transcript_type[tx$transcript_id %in% bed$V1])
length(unique(tx$transcript_type[tx$transcript_id %in% bed$V1])) # 19
colnames(bed)[1] <- "transcript_id"
tx[1:3,]
dim(bed)
#6572662
#grep("ENST00000607096",tx$transcript_id)
#bed <- as_tibble(bed)
#tx <- as_tibble(tx[tx$transcript_id,])
bed <- bed %>% 
  left_join(x = bed, y = tx, by = "transcript_id")
dim(bed)
#6572662
table(bed$transcript_type)
length(unique((bed$transcript_type)))
bed$V4 <- bed$transcript_id
bed$V5 <- 0
#6574543-6572662=1881


bed2 <- bed %>% 
  # filter(transcript_type != "repeats.for", transcript_type != "repeats.rev") %>% # hg38txNoRepeats
  dplyr::filter(transcript_type %in% c("rRNA","pri_miRNA","lncRNA", "mRNA", "piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA') ) %>% # hg38txNoDNA
  dplyr::select(-transcript_type)
# data.table::fwrite(bed2,"/BioII/lulab_b/baopengfei/gitsoft/clipper/clipper/data/regions/hg38tx_exons.bed",quote = F,sep = "\t",row.names = F,col.names = F)
# data.table::fwrite(bed2,"/BioII/lulab_b/baopengfei/gitsoft/clipper/clipper/data/regions/hg38tx_genes.bed",quote = F,sep = "\t",row.names = F,col.names = F)
# data.table::fwrite(bed2,"/BioII/lulab_b/baopengfei/gitsoft/clipper/clipper/data/regions/hg38txNoRepeats_exons.bed",quote = F,sep = "\t",row.names = F,col.names = F)
# data.table::fwrite(bed2,"/BioII/lulab_b/baopengfei/gitsoft/clipper/clipper/data/regions/hg38txNoRepeats_genes.bed",quote = F,sep = "\t",row.names = F,col.names = F)
data.table::fwrite(bed2,"/BioII/lulab_b/baopengfei/gitsoft/clipper/clipper/data/regions/hg38txNoDNA_exons.bed",quote = F,sep = "\t",row.names = F,col.names = F)
data.table::fwrite(bed2,"/BioII/lulab_b/baopengfei/gitsoft/clipper/clipper/data/regions/hg38txNoDNA_genes.bed",quote = F,sep = "\t",row.names = F,col.names = F)

gtf <- bed2
gtf$V2 <- gtf$V2+1
gtf[1:2, ]
#chr1    .       ncRNA   120725  120932  .       -       .       gene_id "ENSG00000238009.6";gene_name "ENSG00000238009";biotype "lncRNA";
gtf <- gtf %>% 
  dplyr::mutate(source="AS_STRUCTURE",type="gene",phase=".",score=".",meta=paste0("gene_id=",transcript_id,";mrna_length=",V3-V2+1,";premrna_length=",V3-V2+1),x="") %>% 
  dplyr::select(transcript_id,source,type,V2,V3,score,V6,phase,meta,x)
# data.table::fwrite(gtf,"/BioII/lulab_b/baopengfei/gitsoft/clipper/clipper/data/hg38tx.AS.STRUCTURE.COMPILED.gff",quote = F,sep = "\t",row.names = F,col.names = F)
#data.table::fwrite(gtf,"/BioII/lulab_b/baopengfei/gitsoft/clipper/clipper/data/hg38txNoRepeats.AS.STRUCTURE.COMPILED.gff",quote = F,sep = "\t",row.names = F,col.names = F)
data.table::fwrite(gtf,"/BioII/lulab_b/baopengfei/gitsoft/clipper/clipper/data/hg38txNoDNA.AS.STRUCTURE.COMPILED.gff",quote = F,sep = "\t",row.names = F,col.names = F)


## newTxID
#ENST00000607781.1
#(-), "", $3); gsub(/\(\+\)/, "_pos", $3); gsub(/::/, "__", $3); gsub(/:/, "___", $3); gsub(/\-/, "____", $3); gsub(/\./, "_____", $3); gsub(/\(-\)/, "_neg", $7); gsub(/\(\+\)/, "_pos", $7); gsub(/::/, "__", $7); gsub(/:/, "___", $7); gsub(/\-/, "____", $7); gsub(/\./, "_____"
#bed2[bed2$transcript_id=="ENST00000607781.1",]
bed2$transcript_id <- gsub("(-)","_neg",bed2$transcript_id,fixed = T)
bed2$transcript_id <- gsub("(+)","_pos",bed2$transcript_id,fixed = T)
bed2$transcript_id <- gsub("::","__",bed2$transcript_id,fixed = T)
bed2$transcript_id <- gsub(":","___",bed2$transcript_id,fixed = T)
bed2$transcript_id <- gsub("-","____",bed2$transcript_id,fixed = T)
bed2$transcript_id <- gsub(".","_____",bed2$transcript_id,fixed = T)
bed2$V4 <- bed2$transcript_id
head(bed2)
# data.table::fwrite(bed2,"/BioII/lulab_b/baopengfei/gitsoft/clipper/clipper/data/regions/hg38txnewTxID_exons.bed",quote = F,sep = "\t",row.names = F,col.names = F)
# data.table::fwrite(bed2,"/BioII/lulab_b/baopengfei/gitsoft/clipper/clipper/data/regions/hg38txnewTxID_genes.bed",quote = F,sep = "\t",row.names = F,col.names = F)
data.table::fwrite(bed2,"/BioII/lulab_b/baopengfei/gitsoft/clipper/clipper/data/regions/hg38txNoDNAnewTxID_exons.bed",quote = F,sep = "\t",row.names = F,col.names = F)
data.table::fwrite(bed2,"/BioII/lulab_b/baopengfei/gitsoft/clipper/clipper/data/regions/hg38txNoDNAnewTxID_genes.bed",quote = F,sep = "\t",row.names = F,col.names = F)
# data.table::fwrite(bed2,"/BioII/lulab_b/baopengfei/gitsoft/clipper/clipper/data/regions/hg38txNoRepeatsnewTxID_exons.bed",quote = F,sep = "\t",row.names = F,col.names = F)
# data.table::fwrite(bed2,"/BioII/lulab_b/baopengfei/gitsoft/clipper/clipper/data/regions/hg38txNoRepeatsnewTxID_genes.bed",quote = F,sep = "\t",row.names = F,col.names = F)
 
gtf <- bed2
gtf$V2 <- gtf$V2+1
gtf <- gtf %>% 
  dplyr::mutate(source="AS_STRUCTURE",type="gene",phase=".",score=".",meta=paste0("gene_id=",transcript_id,";mrna_length=",V3-V2+1,";premrna_length=",V3-V2+1),x="") %>% 
  dplyr::select(transcript_id,source,type,V2,V3,score,V6,phase,meta,x)
# data.table::fwrite(gtf,"/BioII/lulab_b/baopengfei/gitsoft/clipper/clipper/data/hg38txnewTxID.AS.STRUCTURE.COMPILED.gff",quote = F,sep = "\t",row.names = F,col.names = F)
data.table::fwrite(gtf,"/BioII/lulab_b/baopengfei/gitsoft/clipper/clipper/data/hg38txNoDNAnewTxID.AS.STRUCTURE.COMPILED.gff",quote = F,sep = "\t",row.names = F,col.names = F)
# data.table::fwrite(gtf,"/BioII/lulab_b/baopengfei/gitsoft/clipper/clipper/data/hg38txNoRepeatsnewTxID.AS.STRUCTURE.COMPILED.gff",quote = F,sep = "\t",row.names = F,col.names = F)


#SPECIES.AS.STRUCTURE.COMPILED.gff
#chrX AS_STRUCTURE gene 41333284 41364472 . + . gene_id=ENSG00000215301.10;mrna_length=15892;premrna_length=3118

# #get gtf for CLAM.permutation_peakcallerpy (filter $3=="gene" records with ' gene_id "XX" ')
# gff <- data.table::fread("/BioII/lulab_b/baopengfei/gitsoft/clipper/clipper/data/hg38txnewTxID.AS.STRUCTURE.COMPILED.gff",stringsAsFactors = F,sep = "\t",header = F)
# gff[1:3,]
# gff <- gff[,1:9]
# gff$V9 <- unlist(sapply(strsplit(gff$V9,";"),"[",1))
# gff$V9 <- gsub("gene_id=","gene_id ",gff$V9)
# gff$V9 <- gsub("gene_id ","gene_id \"",gff$V9)
# gff$V9 <- paste0(gff$V9,"\"")
# data.table::fwrite(gff,"/BioII/lulab_b/baopengfei/gitsoft/clipper/clipper/data/hg38txnewTxID.AS.STRUCTURE.COMPILED.gtf",quote = F,sep = "\t",row.names = F,col.names = F)

gff <- data.table::fread("/BioII/lulab_b/baopengfei/gitsoft/clipper/clipper/data/hg38txNoDNAnewTxID.AS.STRUCTURE.COMPILED.gff",stringsAsFactors = F,sep = "\t",header = F)
gff[1:3,]
gff <- gff[,1:9]
gff$V9 <- unlist(sapply(strsplit(gff$V9,";"),"[",1))
gff$V9 <- gsub("gene_id=","gene_id ",gff$V9)
gff$V9 <- gsub("gene_id ","gene_id \"",gff$V9)
gff$V9 <- paste0(gff$V9,"\"")
data.table::fwrite(gff,"/BioII/lulab_b/baopengfei/gitsoft/clipper/clipper/data/hg38txNoDNAnewTxID.AS.STRUCTURE.COMPILED.gtf",quote = F,sep = "\t",row.names = F,col.names = F)





# archive -----
# # prepare gtf for PEKA motif analysis tool (failed)
# gtf <- bed2
# gtf$V2 <- gtf$V2+1
# gtf[1:2, ]
# #chr1    .       ncRNA   120725  120932  .       -       .       gene_id "ENSG00000238009.6";gene_name "ENSG00000238009";biotype "lncRNA";
# gtf <- gtf %>% 
#   mutate(source=".",type="ncRNA",phase=".",score=".",meta=paste0("gene_id \"",transcript_id,"\";gene_name \"",transcript_id,"\";biotype \"txref\"")) %>% 
#   select(transcript_id,source,type,V2,V3,score,V6,phase,meta)
# data.table::fwrite(gtf,"/BioII/lulab_b/baopengfei/gitsoft/peka/TestData/inputs/hg38_10RNA_sorted.regions.gtf",quote = F,sep = "\t",row.names = F,col.names = F)
# #data.table::fwrite(gtf,"/BioII/lulab_b/baopengfei/gitsoft/clipper/clipper/data/hg38txNoRepeats.AS.STRUCTURE.COMPILED.gff",quote = F,sep = "\t",row.names = F,col.names = F)
# # data.table::fwrite(gtf,"/BioII/lulab_b/baopengfei/gitsoft/clipper/clipper/data/hg38txNoDNA.AS.STRUCTURE.COMPILED.gff",quote = F,sep = "\t",row.names = F,col.names = F)
# 
# 
# 
# tmp <- read.table("/BioII/lulab_b/baopengfei/gitsoft/peka/TestData/inputs/peak_gn.bed")
# tmp <- tmp[,1:6]
# tmp[1:3,]
# pos <- as.integer((tmp[,2]+tmp[,3]) * 0.5)
# pos_left <- pos - 10
# # pos_right <- pos + 10
# 
# res <- list()
# for (i in 1:nrow(tmp)){
#   rec <- tmp[i,]
#   s <- pos_left[i]
#   e <- s+1 #pos_right[i]
#   df <- data.frame("chr"=rec$V1, "s"=s:(s+19), "e"=e:(e+19), "name"=rec$V4, "score"=rec$V5, "strand"=rec$V6)
#   res[[i]] <- df
# }
# res[[2]]
# res.df <- do.call(rbind,res)
# # tmp[,2] <- pos
# # tmp[,3] <- pos+1
# write.table(x=res.df,file = "/BioII/lulab_b/baopengfei/gitsoft/peka/TestData/inputs/peak_gn_center.bed",quote = F, row.names = F, col.names = F, sep = '\t')
# 




### archive
## add genomic stranded element and miR (ENST ID) to chromsome size
# cp genome/hg38/chrom_sizes/transcriptome_genome_old genome/hg38/chrom_sizes/transcriptome_genome
# for i in intron.rev promoter.rev enhancer.rev repeats.rev ;
# do
# cat genome/hg38/bed/${i}.bed | awk '{print($4"::"$1":"$2"-"$3"(-)" "\t" $3-$2)}' > genome/hg38/chrom_sizes/${i}
# cat genome/hg38/chrom_sizes/${i} >> genome/hg38/chrom_sizes/transcriptome_genome
# done
# for i in intron.for promoter.for enhancer.for repeats.for ;
# do
# cat genome/hg38/bed/${i}.bed | awk '{print($4"::"$1":"$2"-"$3"(+)" "\t" $3-$2)}' > genome/hg38/chrom_sizes/${i}
# cat genome/hg38/chrom_sizes/${i} >> genome/hg38/chrom_sizes/transcriptome_genome
# done
# i="miRNA" 
# cat genome/hg38/bed_by_biotype/${i}.bed | awk '{print($4 "\t" $3-$2)}' > genome/hg38/chrom_sizes/pri_${i}
# cat genome/hg38/chrom_sizes/pri_${i} >> genome/hg38/chrom_sizes/transcriptome_genome

## also add genomic stranded element and miR (ENST ID) to transcript table
#/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/scripts/addDNAelementToTxMapping.R

## need consider strand when extract fasta 
# for i in intron.for intron.rev promoter.for promoter.rev enhancer.for enhancer.rev repeats.for repeats.rev ;
# do 
# bedtools getfasta -s -name -fi genome/hg38/fasta/genome.fa -bed genome/hg38/bed/${i}.bed > genome/hg38/fasta/${i}.fa; 
# done

# for i in intron.for intron.rev promoter.for promoter.rev enhancer.for enhancer.rev repeats.for repeats.rev ;
# do
# samtools faidx genome/hg38/fasta/${i}.fa
# done

# for i in intron.for intron.rev promoter.for promoter.rev enhancer.for enhancer.rev repeats.for repeats.rev ;
# do
# bowtie2-build --threads 8  genome/hg38/fasta/${i}.fa genome/hg38/genome_index/bowtie2/${i};
# done

# for i in intron.for intron.rev promoter.for promoter.rev enhancer.for enhancer.rev repeats.for repeats.rev ;
# do
# for j in 1 2 3 4 rev.1 rev.2 ;
# do
# ln -s ../../genome_index/bowtie2/${i}.${j}.bt2 genome/hg38/rsem_index/bowtie2/${i}.${j}.bt2 ; 
# done
# done

# create bed12 --------------------------------------------------------
## create bed12 for genomics intervals
cat exSeek-dev/genome/hg38/bed/{promoter,enhancer,intron,repeats}.{for,rev}.bed | awk 'BEGIN{{OFS=FS="\t"}}{print $1,$2,$3,$4"::"$1":"$2"-"$3"("$6")",100,$6,$2,$3,0,1,$3-$2",",0}' > exSeek-dev/genome/hg38/bed/long_DNA.bed
cat exSeek-dev/genome/hg38/bed/{promoter,enhancer,intron}.{for,rev}.bed | awk 'BEGIN{{OFS=FS="\t"}}{print $1,$2,$3,$4"::"$1":"$2"-"$3"("$6")",100,$6,$2,$3,0,1,$3-$2",",0}' > exSeek-dev/genome/hg38/bed/long_DNA_noRepeats.bed
#-> % head genome/hg38/bed/promoter.rev.bed
#chr1	36274	36474	promoter	0	-
#chr1	37674	41274	promoter	0	-


## create bed12 for rRNA (failed)
id.conv <- read.table("genome/hg38/source/refSeq_rRNA.gene_names.txt")
colnames(id.conv) <- c("ref_seq_id","gene_name")
# 


#bash
gffread  genome/hg38/gff3/refseq_sorted.gff3 -F -T -o genome/hg38/gff3/refseq_sorted.gtf

#R
gtf <- rtracklayer::import.gff("genome/hg38/gff3/refseq_sorted.gtf")
gtf.df <- as.data.frame(gtf) # 1632
gtf.df$start <- gtf.df$start - 1 # convert to bed coord
head(gtf.df)
gtf.df <- gtf.df[gtf.df$type == "exon",] # 544 = 1632 / 3
gtf.df <- gtf.df[!duplicated(gtf.df[,c("seqnames","start","end","gene_name","transcript_id")]),]
# # length(unique(gtf.df$gene_name))
table(gtf.df$gene_name %in% (id.conv$gene_name))
table((id.conv$gene_name) %in% gtf.df$gene_name)
gtf.df <- gtf.df[(gtf.df$gene_name %in% (id.conv$gene_name)),]
length(unique(gtf.df$gene_name))
length(unique(gtf.df$transcript_id))
table(gtf.df$transcript_id %in% id.conv$ref_seq_id)
bed <- gtf.df[,c("seqnames","start","end","transcript_id","width","strand")]
bed$width <- 100
#table(gtf.df$strand)
write.table(bed,file = "tmp.bed",quote = F,sep = "\t",row.names = F,col.names = F)

#bash
cat tmp.bed | awk 'BEGIN{{OFS=FS="\t"}}{print $1,$2,$3,$4,100,$6,$2,$3,0,1,$3-$2",",0}' > genome/hg38/bed_by_biotype/rRNA.bed
ln -s genome/hg38/bed_by_biotype/rRNA.bed genome/hg38/bed/rRNA.bed
bed="genome/hg38/bed/rRNA.bed"
outbed=$(echo $bed | sed s/"\.bed"/""/g)"_newTxID.bed"
(cat ${bed} |  awk 'BEGIN{FS=OFS="\t"} {gsub(/\(-\)/, "_neg", $4); gsub(/\(\+\)/, "_pos", $4); gsub(/::/, "__", $4); gsub(/:/, "___", $4); gsub(/\-/, "____", $4); gsub(/\./, "_____", $4)} 1 '  )  >  ${outbed}




# create tbed --------------------------------------------------------
## add rRNA tbed (pri_miR, piRNA the same?)
cd /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev
cat genome/hg38/chrom_sizes/rRNA | awk 'BEGIN{{OFS=FS="\t"}} {print $1,0,$2,"X",".","+"}' > genome/hg38/tbed/rRNA.bed





