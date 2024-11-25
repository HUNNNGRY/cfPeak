#! /usr/bin/env Rscript

# append annotation to bed6 sncRNA peak

suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser(description='filter intra-stable records for RNAfold-permutation result')
parser$add_argument('--inputDir', type='character', required=TRUE,
                    help='inputDir of cfpeak(CNN) contains multple pre-calculated files')
parser$add_argument('--outFile', type='character', required=TRUE,
                    help='outfile path to save peak summary')
parser$add_argument('--online', type='character', default="FALSE",
                    help='online conversion of biomart ?')
parser$add_argument('--fasta', type='character', default="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/fasta_newTxID/combine19.fa",
                    help='gn/tx fasta file path')
parser$add_argument('--chrSize', type='character', default="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt",
                    help='gn/tx chrSize file path')
# parser$add_argument('--dedup', type='character', required=TRUE,
#                     help='dedup')

args <- parser$parse_args()

for(i in 1:length(args)){
  v.name <- names(args)[i]
  assign(v.name,args[[i]])
  message(paste0(v.name,": ",get(v.name)))
}
online <- ifelse(toupper(online) %in% c("TRUE","T"), TRUE, FALSE)
message(online)

# load all func.
source("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/scripts/util.R")
#


# # test OSCC
# pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
# outpre <- "/data2/lulab1/bpf/projects/WCHSU-FTC"
# dst <- "lulab_oscc_tissue_diff" # lulab_oscc_plasma_diff
# dedup <- "_dedup"
# #paste0(outpre,"/output/",dst,"/call_peak",dedup,"/cfpeakCNN/b5_d50_p1_8DNA_gn.bed")
# online <- F

# # # # test SLE
# pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
# outpre <- "/data2/lulab1/bpf/projects/WCHSU-FTC"
# dst <- "SLE"
# dedup <- "_dedup"
# online <- F

# # # test TCGA (2406)
# pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
# outpre <- pre
# dst <- "TCGA_small16"
# dedup <- "_all"
# online <- F


txFasta="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/fasta_newTxID/combine19.fa"
txChrSize <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt"

inputDir <- paste0(outpre,"/output/",dst,"/call_peak",dedup,"/cfpeakCNN")

bedtoolsPath <- "/BioII/lulab_b/baopengfei/biosoft"#inBed <- paste0(inputDir,"/b5_d50_p1.bed")
inBed.11RNA.tx <- paste0(inputDir,"/b5_d50_p1_11RNA.bed")
inBed.8DNA.tx <- paste0(inputDir,"/b5_d50_p1_8DNA.bed")
inBed.11RNA.gn <- paste0(inputDir,"/b5_d50_p1_11RNA_gn.bed")
inBed.8DNA.gn <- paste0(inputDir,"/b5_d50_p1_8DNA_gn.bed")
inIntersect.11RNA.tx <- paste0(inputDir,"/b5_d50_p1_11RNA.intersect.bed")
inIntersect.8DNA.gn <- paste0(inputDir,"/b5_d50_p1_8DNA_gn.intersect.bed")
inIntersect.AluIR.RNA.tx <- paste0(inputDir,"/b5_d50_p1_11RNA.intersect-Alu_IR.bed") # op: Alu
inIntersect.AluIR.DNA.gn <- paste0(inputDir,"/b5_d50_p1_8DNA_gn.intersect-Alu_IR.bed") # op: Alu
inIntersect.AluNonIR.RNA.tx <- paste0(inputDir,"/b5_d50_p1_11RNA.intersect-Alu_nonIR.bed") # op: Alu
inIntersect.AluNonIR.DNA.gn <- paste0(inputDir,"/b5_d50_p1_8DNA_gn.intersect-Alu_nonIR.bed") # op: Alu
nearest.11RNA.gn <- paste0(inputDir,"/b5_d50_p1_11RNA_gn.bed.nearest")
nearest.8DNA.gn <- paste0(inputDir,"/b5_d50_p1_8DNA_gn.bed.nearest")
inFasta <- paste0(inputDir,"/b5_d50_p1.bed.fa")
inFasta.CPC <- paste0(inputDir,"/b5_d50_p1.bed.fa.CPC.txt")
inFdr.11RNA.tx <- paste0(inputDir,"/b5_d50_p1_11RNA.bed.exp.fa.csv")
inFdr.8DNA.gn <- paste0(inputDir,"/b5_d50_p1_8DNA_gn.bed.exp.fa.csv")
inPhastCon.11RNA.gn <- paste0(inputDir,"/b5_d50_p1_11RNA_gn.bed.phastCons")
inPhastCon.11RNA.gn.flank <- paste0(inputDir,"/b5_d50_p1_11RNA_gn.bed.phastCons.flank")
inPhastCon.8DNA.gn <- paste0(inputDir,"/b5_d50_p1_8DNA_gn.bed.phastCons")
inPhastCon.8DNA.gn.flank <- paste0(inputDir,"/b5_d50_p1_8DNA_gn.bed.phastCons.flank")
outFile <- paste0(inputDir,"/b5_d50_p1.summary")


bedtoolsPath <- "/BioII/lulab_b/baopengfei/biosoft"
options(bedtools.path = bedtoolsPath) # /BioII/lulab_b/baopengfei/anaconda3/envs/exvariance/bin
options(stringsAsFactors = F)
suppressPackageStartupMessages(library(bedtoolsr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(IRanges))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(GenomicScores))
suppressPackageStartupMessages(library(biomaRt))



# read ref
ref <- data.table::fread(txChrSize,data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)


# read bed
bed.11RNA.tx <- readBedDf(inBed.11RNA.tx)
bed.8DNA.tx <- readBedDf(inBed.8DNA.tx)
bed.11RNA.gn <- readBedDf(inBed.11RNA.gn)
bed.8DNA.gn <- readBedDf(inBed.8DNA.gn)
bed.11RNA <- cbind(bed.11RNA.tx,bed.11RNA.gn)
bed.8DNA <- cbind(bed.8DNA.tx,bed.8DNA.gn)
bed2 <- as.data.frame(rbind(bed.11RNA,bed.8DNA))




####################
# gn re-assign
####################
# bed2 <- data.table::fread(outFile,sep='\t',header = T,data.table = F,check.names = F,stringsAsFactors = F)
# dim(bed2)
sncRNA.table <- data.table::fread( paste0(outpre,"/output/",dst,"/call_peak",dedup,"/cfpeakCNN/gnAssign/snc_assign.txt"), header = T,check.names = F,stringsAsFactors = F,sep = "\t")
sncRNA.table <- sncRNA.table[order(as.integer(gsub("peak_","",sncRNA.table$name)),decreasing = F),]
dim(sncRNA.table)
#sncRNA.table[1:13,]  # filterd [10,200]
s <- sncRNA.table$snc_type[match(bed2$name,sncRNA.table$name)]
table(s)
bed2$gn_type <- s
# colnames(bed2)
# bed2 <- bed2[,c(1:18,ncol(bed2),19:(ncol(bed2)-1))]
# write.table(bed2, outFile,quote=F,row.names=F,col.names=T,sep="\t")




####################
# precursor: pri_miR...intron
####################
#class
#symbol, ENSG, function
#R --inBed
#inBed="/data2/lulab1/bpf/projects/WCHSU-FTC/output/SLE/call_peak_dedup/cfpeakCNN/b5_d50_p1.bed"
#ref[1:3,]
bed2$tx_type <- ref$transcript_type[match(bed2$tx,ref$transcript_id)]
bed2$tx_len <- ref$tx.length[match(bed2$tx,ref$transcript_id)]
bed2$tx_type <- gsub("_rev|_for|\\.for|\\.rev","",bed2$tx_type)
bed2$tx_ENST <- new2oldTxID(bed2$tx)
bed2$tx_ENST <- unlist(sapply(strsplit(bed2$tx_ENST,".",fixed=T),"[",1))

q <- unique(bed2$tx_ENST[grepl("^ENST",bed2$tx_ENST)])
#length(q)
#rm(list = c("mart"))
ann <- enst2ensgEntrezSymbolDescritptionBiotypeDF(query = q, online = online)
colnames(ann) <- paste0("tx_",colnames(ann) )
#table(duplicated(ann$tx_ENST))
ann <- ann[!duplicated(ann$tx_ENST),]
nrow(ann)<=length(q)
bed2 <- left_join(bed2,ann, by="tx_ENST")
dim(bed2)




####################
# nearest_gene: pri_miR...ENSG***
####################
nearestTx <- read.table(nearest.11RNA.gn,sep = "\t",header = F)
nearestTx2 <- read.table(nearest.8DNA.gn,sep = "\t",header = F)
nearestTx <- rbind(nearestTx,nearestTx2)
nearestTx <- nearestTx[,c("V4","V25","V16")]
#dim(nearestTx)
colnames(nearestTx) <- c("name","nearestCodingTxDist","nearestCodingTx_ENST")
nearestTx <- nearestTx[!duplicated(nearestTx$name),]
nearestTx$nearestCodingTx_ENST <- unlist(sapply(strsplit(nearestTx$nearestCodingTx_ENST,".",fixed=T),"[",1))
q <- unique(nearestTx$nearestCodingTx_ENST[grepl("^ENST",nearestTx$nearestCodingTx_ENST)])
#length(q)
ann <- enst2ensgEntrezSymbolDescritptionBiotypeDF(query = q, online = online)
colnames(ann) <- paste0( "nearestCodingTx_",colnames(ann) )
#table(duplicated(ann$nearestCodingTx_ENST))
ann <- ann[!duplicated(ann$nearestCodingTx_ENST),]
nrow(ann)<=length(q)
nearestTx <- left_join(nearestTx,ann, by="nearestCodingTx_ENST")
#nearestTx.bak <- nearestTx

# append to bed
bed2 <- left_join(bed2,nearestTx, by="name")
#bed2.bak <- bed2



####################
# peak_target_gene: top5 mRNA (holding)
####################
#symbol, ENSG, function
#ORA mode only:
#miR # use R-multiMiR pre-calculated reference file: hsa-miR-1-5p --> ENSG**
#IntaRNA: predict all types, but only use those not miR. holding, too much calculation, and directly use mRNA UTR3 is not enough. there are many potential regulation approach for sncRNA other than AGO-dependent or -independent base-pairing with target or by 3-D aptamer structure: like compete with AGO-miRNA to indirectly regulation mRNA expr (Fish et al. 2018, Nat Med)

####################
# region_annotation
####################
#RBP
#EV
#RGS
#(Alu)
tmp.list <- list()
for (inIntersect in c(inIntersect.11RNA.tx, inIntersect.8DNA.gn) ){
  message(inIntersect)
  # inIntersect <- inIntersect.8DNA.gn
  res.df <- data.table::fread(inIntersect, data.table = F,sep = "\t",header = T,fill = T,check.names = F,stringsAsFactors = F)
  res.df <- res.df[!is.na(res.df$site) & res.df$site!="",]
  # res.df[268609:268612,]
  # res.df$RNA <- ref$transcript_type[match(res.df$chr,ref$transcript_id)]
  #table(is.na(res.df$site_overlap_ratio_product))
  res.df$site_overlap_ratio_product[is.na(res.df$site_overlap_ratio_product)] <- 0
  res.df$site <- factor(res.df$site,levels = c("AGO2","otherRBPs","RBPhotspot","EV","G4","iM","Alu")) 
  res.df <- dplyr::as_tibble(res.df) %>% 
    dplyr::arrange(desc(site),desc(site_overlap_ratio_product),desc(site_overlap_base)) %>% 
    dplyr::distinct(name, site, .keep_all = TRUE) # each peak, each site only keep one most-overlpped record
  # table(res.df$site) # i-motif not enough ?
  res.df$name <- unlist(sapply(strsplit(res.df$name,"--"),"[",1))
  #colnames(res.df)
  res.df <- res.df[,c("name","site_name","site_overlap_base","site_overlap_ratio_product","site")] # "site_start","site_end",
  # res.df2 <- reshape2::acast(formula = , value.var = , id.var = c(''))
  #dplyr::across()
  res.df2 <- tidyr::pivot_wider(data = res.df, id_cols = "name", names_from = "site", names_prefix = "", names_sep = "_",
                                values_from = c("site_name","site_overlap_base","site_overlap_ratio_product") # "site_start","site_end",
                                )
  #table(duplicated(res.df2))
  #rownames(res.df2) <- res.df2$name
  tmp.list[[inIntersect]] <- res.df2
}
tmp.df <- as.data.frame(do.call(rbind,tmp.list))
rownames(tmp.df) <- tmp.df$name


{
#Alu_nonIR (op)
tmp.list <- list()
for (inIntersect in c(inIntersect.AluNonIR.RNA.tx, inIntersect.AluNonIR.DNA.gn) ){
  message(inIntersect)
  # inIntersect <- inIntersect.8DNA.gn
  res.df <- data.table::fread(inIntersect, data.table = F,sep = "\t",header = T,fill = T,check.names = F,stringsAsFactors = F)
  res.df <- res.df[!is.na(res.df$site) & res.df$site!="",]
  # res.df[268609:268612,]
  # res.df$RNA <- ref$transcript_type[match(res.df$chr,ref$transcript_id)]
  #table(is.na(res.df$site_overlap_ratio_product))
  res.df$site_overlap_ratio_product[is.na(res.df$site_overlap_ratio_product)] <- 0
  # res.df$site <- factor(res.df$site,levels = c("AGO2","otherRBPs","RBPhotspot","EV","G4","iM")) 
  res.df <- dplyr::as_tibble(res.df) %>% 
    dplyr::arrange(desc(site),desc(site_overlap_ratio_product),desc(site_overlap_base)) %>% 
    dplyr::distinct(name, site, .keep_all = TRUE) # each peak, each site only keep one most-overlpped record
  # table(res.df$site) # i-motif not enough ?
  res.df$name <- unlist(sapply(strsplit(res.df$name,"--"),"[",1))
  #colnames(res.df)
  res.df <- res.df[,c("name","site_name","site_overlap_base","site_overlap_ratio_product","site")] # "site_start","site_end",
  # res.df2 <- reshape2::acast(formula = , value.var = , id.var = c(''))
  #dplyr::across()
  res.df2 <- tidyr::pivot_wider(data = res.df, id_cols = "name", names_from = "site", names_prefix = "", names_sep = "_",
                                values_from = c("site_name","site_overlap_base","site_overlap_ratio_product") # "site_start","site_end",
  )
  #table(duplicated(res.df2))
  #rownames(res.df2) <- res.df2$name
  tmp.list[[inIntersect]] <- res.df2
}
tmp.df2 <- as.data.frame(do.call(rbind,tmp.list))
rownames(tmp.df2) <- tmp.df2$name

#Alu_IR (op)
tmp.list <- list()
for (inIntersect in c(inIntersect.AluIR.RNA.tx, inIntersect.AluIR.DNA.gn) ){
  message(inIntersect)
  # inIntersect <- inIntersect.8DNA.gn
  res.df <- data.table::fread(inIntersect, data.table = F,sep = "\t",header = T,fill = T,check.names = F,stringsAsFactors = F)
  res.df <- res.df[!is.na(res.df$site) & res.df$site!="",]
  # res.df[268609:268612,]
  # res.df$RNA <- ref$transcript_type[match(res.df$chr,ref$transcript_id)]
  #table(is.na(res.df$site_overlap_ratio_product))
  res.df$site_overlap_ratio_product[is.na(res.df$site_overlap_ratio_product)] <- 0
  # res.df$site <- factor(res.df$site,levels = c("AGO2","otherRBPs","RBPhotspot","EV","G4","iM")) 
  res.df <- dplyr::as_tibble(res.df) %>% 
    dplyr::arrange(desc(site),desc(site_overlap_ratio_product),desc(site_overlap_base)) %>% 
    dplyr::distinct(name, site, .keep_all = TRUE) # each peak, each site only keep one most-overlpped record
  # table(res.df$site) # i-motif not enough ?
  res.df$name <- unlist(sapply(strsplit(res.df$name,"--"),"[",1))
  #colnames(res.df)
  res.df <- res.df[,c("name","site_name","site_overlap_base","site_overlap_ratio_product","site")] # "site_start","site_end",
  # res.df2 <- reshape2::acast(formula = , value.var = , id.var = c(''))
  #dplyr::across()
  res.df2 <- tidyr::pivot_wider(data = res.df, id_cols = "name", names_from = "site", names_prefix = "", names_sep = "_",
                                values_from = c("site_name","site_overlap_base","site_overlap_ratio_product") # "site_start","site_end",
  )
  #table(duplicated(res.df2))
  #rownames(res.df2) <- res.df2$name
  tmp.list[[inIntersect]] <- res.df2
}
tmp.df3 <- as.data.frame(do.call(rbind,tmp.list))
rownames(tmp.df3) <- tmp.df3$name

tmp.df <- left_join(x = tmp.df, y = tmp.df2, by="name")
tmp.df <- left_join(x = tmp.df, y = tmp.df3, by="name")
#tmp.df <- tmp.df[,order(colnames(tmp.df))]
}


# append to bed
#bed2 <- merge(bed2,res.df2, by.x="name", by.y="name") # ?
bed2 <- left_join(x = bed2, y = tmp.df, by="name")





####################
# 2ndStructure: extend peak to 20bp (or one-fold length?) each side till tx boundary
####################
#MFE
#MFE.zscore
#MFE.FDR
# add MFE
tmp.list <- list()
for (inFdr in c(inFdr.11RNA.tx, inFdr.8DNA.gn) ){
  message(inFdr)
#inFdr <- "/data2/lulab1/bpf/projects/WCHSU-FTC/output/SLE/call_peak_dedup/cfpeakCNN/archive/b5_d50_p1_11RNA.bed.exp.fa.csv"
fdr <- read.csv(inFdr,sep = ",",header = F)
fdr$V1 <- gsub(">","",fdr$V1,fixed = T)
fdr$V1 <- gsub("(+)","",fdr$V1,fixed = T)
fdr$V1 <- gsub("(-)","",fdr$V1,fixed = T)
fdr$V1 <- unlist(sapply(strsplit(fdr$V1,"::",fixed = T),"[",1)) # if not bedtools getfasta -nameOnly
#table(fdr$V2<=intraFdr) # ~50%
#below seem to handle vec input
query.1 <- lengths(regmatches(fdr$V4, gregexpr("(", fdr$V4, fixed = T))) 
query.2 <- lengths(regmatches(fdr$V4, gregexpr(")", fdr$V4, fixed = T))) 
query.total <- nchar(fdr$V4)
#all(query.1==query.2)
query.pair.frac <- (query.1+query.2)/query.total
fdr$intraPairLen <- query.1+query.2
fdr$intraPairFrac <- query.pair.frac
colnames(fdr) <- c("name","MFE_FDR","MFE","MFE_structure","MFE_fasta","MFE_intraPairLen","MFE_intraPairFrac")

tmp.list[[inFdr]] <- fdr
}
tmp.df <- as.data.frame(do.call(rbind,tmp.list))
rownames(tmp.df) <- tmp.df$name

# append to bed
#bed2 <- merge(bed2,res.df2, by.x="name", by.y="name") # ?
bed2 <- left_join(x = bed2, y = tmp.df, by="name")





####################
# CodingPotential (CPC)
####################
# add peak fa
fasta_sequences <- Biostrings::readDNAStringSet(inFasta)
fasta_dataframe <- data.table::data.table(id = names(fasta_sequences), sequence = as.character(fasta_sequences))
fasta_dataframe$id <- unlist(sapply(strsplit(fasta_dataframe$id,"::"),"[",1))
colnames(fasta_dataframe) <- c("name","peak_fasta")
fasta_dataframe$name <- gsub("(+)","",fasta_dataframe$name,fixed = T)
fasta_dataframe$name <- gsub("(-)","",fasta_dataframe$name,fixed = T)
bed2 <- left_join(x = bed2, y = fasta_dataframe, by="name")


# add CPC
cpc <- data.table::fread(inFasta.CPC, data.table = F,sep = "\t",header = T,fill = T,check.names = F,stringsAsFactors = F)
cpc$`#ID` <- gsub("(+)","",cpc$`#ID`,fixed = T)
cpc$`#ID` <- gsub("(-)","",cpc$`#ID`,fixed = T)
cpc <- cpc[,c("#ID","peptide_length","coding_probability","label")]
colnames(cpc)[1] <- "name"
colnames(cpc)[4] <- "coding_label"
bed2 <- left_join(x = bed2, y = cpc, by="name")




####################
# phastCons (no strand info)
####################
# saveRDS(object = bed2, file = "test.rds")
tmp1 <- read.table(inPhastCon.11RNA.gn,header = T,sep = "\t")
tmp2 <- read.table(inPhastCon.8DNA.gn,header = T,sep = "\t")
tmp1.flank <- read.table(inPhastCon.11RNA.gn.flank,header = T,sep = "\t")
tmp2.flank <- read.table(inPhastCon.8DNA.gn.flank,header = T,sep = "\t")

tmp12 <- rbind(tmp1,tmp2)
#tmp12.flank[1:3,]
tmp12.flank <- rbind(tmp1.flank,tmp2.flank)
tmp12.flank$name <- unlist(sapply(strsplit(tmp12.flank$name,"--"),"[",1))
colnames(tmp12.flank)[2] <- "phastCons100way_flank"
#table(is.na(tmp12.flank$name ))
#tmp <- tmp12.flank[is.na(tmp12.flank$name ),]
all(tmp12.flank$name==tmp12$name) # all T

#add UntraConsIsland + UntraConsIsland.phastCons + UntraConsIsland.flank.phastCons columns
add UntraConsIsland (flank) phastCons
...

bed2 <- left_join(x = bed2, y = tmp12, by="name")
bed2 <- left_join(x = bed2, y = tmp12.flank, by="name")

write.table(bed2, outFile,quote=F,row.names=F,col.names=T,sep="\t")
# dim(bed2) # 62,71 columns
#
# bed2 <- data.table::fread(outFile,sep='\t',header = T,data.table = F,check.names = F,stringsAsFactors = F)
#colnames(bed2)
#dim(bed2)

# add orphan type label

# add GWAS leading SNP,gene





# bed2 <- data.table::fread(outFile,data.table = F,header = T,sep = '\t',check.names = F,stringsAsFactors = F)
# table( bed2$tx_type %in% dna)
# bed2.rna <- bed2[!bed2$tx_type %in% dna,]
# bed2.dna <- bed2[bed2$tx_type %in% dna,]
# dna.fa <- data.table::fread(inFdr.8DNA.gn,data.table = F,header = F,sep = ',',check.names = F,stringsAsFactors = F)
# fa <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/TCGA_small16/call_peak_all/cfpeakCNN/b5_d50_p1_8DNA_gn.bed.exp.fa",data.table = F,header = F,sep = ',',check.names = F,stringsAsFactors = F)
# rownames(dna.fa) <- dna.fa$V1
# rownames(dna.fa) <- gsub("(+)","",rownames(dna.fa),fixed = T)
# rownames(dna.fa) <- gsub("(-)","",rownames(dna.fa) ,fixed = T)
# rownames(dna.fa) <- gsub(">","",rownames(dna.fa) ,fixed = T)
# dim(dna.fa)
# dim(bed2.dna)
# bed2.dna$MFE_fasta <- dna.fa$V5[match(bed2.dna$name,rownames(dna.fa))]
# bed2.dna$MFE_fasta[1:3]
# bed2 <- as.data.frame(rbind(bed2.rna,bed2.dna))
# #