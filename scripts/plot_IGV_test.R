# plot peak cov,bed IGV
# last 221003 by bpf

# test single sample bw + peak.bed (diff exPeak params) ------------------------------------------------------------
meta.list4 <- list()
# meta.list4[["WSQ_SMARTer_NEB"]] <- c("NC_PKU-2392860_smart_PNK_1")
meta.list4[["GSE71008"]] <- c("SRR2105340") # ,"SRR2105335"

pre="call_domain_withRepeats_all" # call_domain_withRepeats_all

genes.list4 <- list()
genes.list4[[1]] <- data.frame(chr="NR_146144_____1",start=0,width=13315)
genes.list4[[2]] <- data.frame(chr="ENST00000602385_____1",start=350,width=150)
genes.list4[[3]] <- data.frame(chr="ENST00000362104_____3",start=0,width=92)
genes.list4[[4]] <- data.frame(chr="ENST00000602385_____1",start=0,width=541)
genes.list4[[5]] <- data.frame(chr="NR_023363_____1",start=0,width=150)

# chr <- "ENST00000362104_____3"
# start <- 1
# width <- 92
# chr <- "ENST00000602385_____1"
# start <- 1
# width <- 541
# chr <- "NR_023363_____1"
# start <- 1
# width <- 150
# chr <- "ENST00000602385_____1"
# start <- 350
# width <- 150

# chr <- "NR_146144_____1"
# start <- 1
# width <- 13315


for (dst in names(meta.list4)){
  print(dst)
  #dst <- "GSE71008" #GSE71008,WSQ_SMARTer_NEB
  smps <- meta.list4[[dst]]
  for (smp in smps){
    #smp <- "SRR2105340" #SRR2105340,NC_PKU-2392860_smart_PNK_1
    print(smp)
    for (i in 1:length(genes.list4)){
      #i <- 3
      print(i)
      chr <- genes.list4[[i]]$chr
      start <- genes.list4[[i]]$start
      width <- genes.list4[[i]]$width
      print(paste0(chr,".",start,".",width))
      res.list2 <- list()
      
      read_bw <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/",pre,"/tbigwig_long_RNA_EM/",smp,".transcriptome.bigWig")
      # x <- rtracklayer::import.bw(read_bw, which=GRanges(c(chr), IRanges(start = start, width = width))) # 1-based ?
      #x <- as.data.frame(x)
      bw.res.df <- readBW(bw=read_bw,chr=chr,start=start,width=width)
      res.list2[[paste0(dst,"_",smp,"_",chr,".",start,".",width)]] <- plotBW(bw.res.df,chr=chr,start=start,width=width,ylab="Depth",title=paste0(dst,"_",smp,"\n",chr,"_",start,"_",width),plotXbreak = T,plotYbreak = T,plotYtext = F)
      
      # peak
      expeak_paths <- Sys.glob(paste0("*multi3*_*_decay*_pval*_t*p*_",smp,".bed"))
      for (expeak_path in expeak_paths){
        # expeak_path <- "multi3global_background_decay0.1_pval0.05_t8p1_NC_PKU-2392860_smart_PNK_1.bed"
        # expeak_path <- "multi3global_background_decay0.1_pval0.05_t8p1_SRR2105340.bed"
        print(expeak_path)
        #expeak_path <- "multi3local_localmaxdecay_decay0.9_pval0.05_t8p1_NC_PKU-2392860_smart_PNK_1.bed"
        recur <- ifelse(grepl("notrecursive_",expeak_path),"r0","r1")
        mode <- ifelse(grepl("global",expeak_path),"global","local")
        boundary <- ifelse(grepl("localmaxdecay",expeak_path),"decay","bg")
        if(grepl("decay0.9",expeak_path)){
          d <- "d90"
        } else if (grepl("decay0.5",expeak_path)){
          d <- "d50"
        } else {
          d <- "d10"
        }
        if(grepl("pval0.05",expeak_path)){
          p <- "p05"
        } else {
          p <- "p1"
        }
        ## localma peak bed
        localmaxEM2 <- data.table::fread(expeak_path,sep = "\t",header = F,stringsAsFactors = F,check.names = F) # 1-based ?
        #localmaxEM2[1:3,]
        localmaxEM2 <- localmaxEM2[,1:6]
        colnames(localmaxEM2) <- c("seqnames","start","end","name","type","strand")
        localmaxEM2$score <- 1
        localmaxEM2 <- localmaxEM2[,c("seqnames","start","end","name","score","strand")]
        #head(localmaxEM2,3)
        localmaxEM2.res.df <- readPeak(bed = localmaxEM2,chr=chr,start=start,width=width)
        bed.localmaxEM2 <- plotPeak(bed = localmaxEM2.res.df,chr=chr,start=start,width=width,
                                    ylab = paste0(recur,"_",mode,"_",boundary,"_",d,"_",p), 
                                    color = "grey30",fill = "grey30")
        # bed.localmaxEM2
        res.list2[[expeak_path]] <- bed.localmaxEM2
      }
  
      ## collect
      tmp2 <- ggpubr::ggarrange(plotlist = res.list2, ncol = 1,heights = c(0.5*length(expeak_paths),rep(1,length(expeak_paths))), nrow = length(expeak_paths)+1, align = "v") # 
      ggsave(plot = tmp2, filename = paste0("expeaks_params_",dst,"_",smp,"-",chr,".",start,".",width,".pdf"),width = 10, height = 16)
    }
  }
}
# x <- 0.00017652042929768404
# x <- 0.9
# 3+2/(0.05+x)
# 3+40*(1-log2(x+1))



# test NC_pool bw + peak.bed ------------------------------------------------------------
setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev")

meta.list2 <- list()
meta.list2[["GSE71008_NCpool"]] <- c("NCpool") # ,"SRR2105335"


genes.list2 <- list()
genes.list2[[1]] <- data.frame(chr="ENST00000635274_____1",start=50,width=230)


for (i in 1:length(genes.list2)){
  #i <- 1
  chr <- genes.list2[[i]]$chr
  start <- genes.list2[[i]]$start
  width <- genes.list2[[i]]$width
  print(paste0(chr,".",start,".",width))
  
  for (dst in names(meta.list2)){
    # dst <- "sim"
    print(dst)
    if (grepl("AGO2_IP|FTC_small",dst,perl=T)){
      pre="call_peak_dedup" 
    } else {
      pre="call_peak_all"
    }
    #dst <- "GSE50676"
    #meta.list2[["GSE50676"]] <- c("GM12878")
    smps <- meta.list2[[dst]]
    for (smp in smps){
      #smp <- "tissue05_blood95"
      res.list2 <- list()
      print(smp)
   
      read_bw <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/",pre,"/tbigwig_11RNA_primary/",smp,".transcriptome.bigWig")
      read_bw2 <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/",pre,"/tbigwig_RNA_EM/",smp,".transcriptome.bigWig")
      
      # x <- rtracklayer::import.bw(read_bw) # 1-based ?
      #x <- as.data.frame(x)
      bw.res.df <- readBW(read_bw,chr = chr,start = start,width = width) 
      res.list2[[paste0(dst,"_",smp,"bw1")]] <- plotBW(bw.res.df,chr = chr,start = start,width = width,ylab = "",xlab = "",title = paste0(dst,"_",smp,"\n",chr,"_",start,"_",width))  
      bw.res.df2 <- readBW(read_bw2,chr = chr,start = start,width = width) 
      res.list2[[paste0(dst,"_",smp,"bw2")]] <- plotBW(bw.res.df2,chr = chr,start = start,width = width,ylab = "",xlab = "",title = paste0(dst,"_",smp,"\n",chr,"_",start,"_",width))  
      
      
      # peak
      bed_clipper <-  paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/",pre,"/clipper_by_sample/b5_p05/",smp,".bed")
      bed_localmaxEM2 <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/",pre,"/expeak_by_sample/b5_d05_p05/",smp,".bed")

      ## CLIPper peak bed
      CLIPper <- data.table::fread(bed_clipper,sep = "\t",header = F,stringsAsFactors = F,check.names = F) # 1-based ?
      if(nrow(CLIPper)==0){
        CLIPper <- data.frame("V1"="","V2"="","V3"="","V4"="","V5"="","V6"="","V7"="","V8"="")
      }
      colnames(CLIPper) <- c("seqnames","start","end","name","score","strand","V7","V8") # score: pval
      CLIPper <- CLIPper[,c("seqnames","start","end","name","score","strand")]
      CLIPper.res.df <- readPeak(CLIPper,chr = chr,start = start,width = width)
      bed.clipper <- plotPeak(bed = CLIPper.res.df,chr = chr,start = start,width = width,ylab = "CLIPper")
      # bed.clipper
      res.list2[[paste0(dst,"_",smp,"CLIPper")]] <- bed.clipper
      
      ## localmax peak bed
      localmaxEM2 <- data.table::fread(bed_localmaxEM2,sep = "\t",header = F,stringsAsFactors = F,check.names = F) # 1-based ?
      if(nrow(localmaxEM2)==0){
        localmaxEM2 <- data.frame("V1"="","V2"="","V3"="","V4"="","V5"="","V6"="")
      } 
      localmaxEM2 <- localmaxEM2[,1:6]
      colnames(localmaxEM2) <- c("seqnames","start","end","name","type","strand")
      localmaxEM2$score <- 1
      localmaxEM2 <- localmaxEM2[,c("seqnames","start","end","name","score","strand")]
      localmaxEM2.res.df <- readPeak(localmaxEM2,chr = chr,start = start,width = width)
      bed.localmaxEM2 <- plotPeak(bed = localmaxEM2.res.df,chr = chr,start = start,width = width,ylab = "exPeak")
      # bed.localmaxEM2
      res.list2[[paste0(dst,"_",smp,"exPeak")]] <- bed.localmaxEM2
      
      ## collect
      tmp2 <- ggpubr::ggarrange(plotlist = res.list2, ncol = 1,heights = c(4,4,1,1),align = "v") 
      ggsave(plot = tmp2, filename = paste0("peaks_2bw_bed_",dst,"_",smp,"-",chr,".",start,".",width,".pdf"),width = 25, height = 8)
    }
  }
}







# test plot multi-tx exons ------
plotMultiTxTrack <- function(txID_new,gtf,fill="grey50"){
#txID <- "NR_023363_____1"
#gtf <- gtf0
# fill <- "black"
txID <- txID_new
txID.old <- new2oldTxID(txID)
is.rRNA <- grepl("NR_",txID)
if(is.rRNA){
  rRNA <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/source/refSeq_rRNA.gene_names.txt")
  txID2 <- rRNA$V2[rRNA$V1==txID.old]
  txID.old2 <- txID2
}
table(gtf$gene_name==gID)
head(gtf,4)
gtf <- as.data.frame(gtf)
if(is.rRNA){
  gID <- txID.old2
  gtf <- gtf[gtf$gene_name==gID & !is.na(gtf$gene_name),]
}else{
  gID <- unique(na.omit(gtf$gene_id[gtf$transcript_id==txID.old]))
  gtf <- gtf[gtf$gene_id==gID,]
}
# gtf[1:3,]
gtf <- gtf[gtf$type=="exon",]
gtf$ID <- paste0("bed_",1:nrow(gtf))
gtf <- gtf[,c("seqnames","start","end","ID","score","strand","type","gene_id","gene_name","gene_type","transcript_id","transcript_type","exon_number","exon_id")]

# peak.df <- data.table::fread(inputFile,data.table = F,header = F,sep = "\t",stringsAsFactors = F)
peak.df <- as.data.frame(gtf)
peak.df <- peak.df[,c("seqnames","start","end","ID","score","strand")]
peak.df$start <- peak.df$start-1
colnames(peak.df) <- c("chr","start","end","name","score","strand")
peak.df$chr <- as.character(peak.df$chr)
peak.df$name <- as.character(peak.df$name)
peak.df$strand <- as.character(peak.df$strand)
peak.df <- peak.df[peak.df$strand == "+" | peak.df$strand == "-" ,] # need contain strand info

## read tx table
chrSize.df <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/tbed/11RNA_map.txt",data.table = F,header = F,sep = "\t",stringsAsFactors = F)
colnames(chrSize.df) <- c("transcript_id","s","e","name","score","strand","transcript_type")
chrSize.df <- chrSize.df[chrSize.df$transcript_id == txID.old,] # chrSize.df[chrSize.df$transcript_type %in% priority,]

## read tx block bed
tx.ref <- rtracklayer::import.bed(file('/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/bed/11RNA_map_bed12.bed'))
tx.ref <- tx.ref[tx.ref$name %in% chrSize.df$transcript_id,]
tx.ref.block0 <- rtracklayer::blocks(tx.ref)

#pre-computed block txt (10 tx bed)
tx.ref.block <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/bed/11RNA_map_block.txt",sep = "\t",data.table = F,header = T,)
#tx.ref.block[1:3,]
tx.ref.block <- tx.ref.block[!grepl("alt|random|chrUn",tx.ref.block$seqnames,perl=T),] # filter non-canonical chr
tx.ref.block <- tx.ref.block[tx.ref.block$enst == txID.old,]

#bedtools intersect -a A.bed -b B.bed -wao
ts <- bedtoolsr::bt.intersect(wao=T, s = T, a = peak.df, b = tx.ref.block) # -f, sorted = T, g = chr.size seems not working
ts <- ts[ts$V15 > 0 & !is.na(ts$V15),]
ts$a.ratio <- ts$V15/(ts$V3-ts$V2)
ts$transcript_type <- chrSize.df$transcript_type[match(ts$V10,chrSize.df$transcript_id)]
ts$tx.strand <- "+"

ts.list <- list()
getTxPos <- function(i){
  #i <- 1
  # if(i%%1000 == 0){
  #   message(i)
  # }
  tmp <- tx.ref.block0[[as.character(ts$V10)[i]]] # v10: tx name
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
ts.1 <- mclapply(1:nrow(ts),getTxPos,mc.cores = 4) # args$cores
ts.1 <- do.call(rbind,ts.1)
ts.1 <- as.data.frame(ts.1)
ts <- cbind(ts,ts.1)

ts <- ts[,c("V10","tx.start","tx.end","V4","V5","tx.strand")]
colnames(ts)[4] <- "ID"
ts <- as_tibble(ts)
gtf <- as_tibble(gtf)

res.tbl <- left_join(ts,gtf)
colnames(res.tbl)[7:11] <- c("gn.seqnames","gn.start","gn.end","gn.score","gn.strand")
colnames(res.tbl)[1:6] <- c("seqnames","start","end","ID","score","strand")


## plot
library(magrittr)
library(dplyr)
library(ggplot2)
library(ggtranscript)
res.tbl <- as_tibble(res.tbl)
# res.tbl %>% head()

# extract exons
res.tbl.exons <- res.tbl %>% dplyr::filter(type == "exon")
display_tx_num <- 5
if(!is.rRNA){
  if (length(table(res.tbl.exons$transcript_id)) > display_tx_num){
    res.tbl.exons <- res.tbl.exons[res.tbl.exons$transcript_id %in% c(txID.old,sample(unique(res.tbl.exons$transcript_id[res.tbl.exons$transcript_id!=txID.old]),display_tx_num-1,replace = F)),] # only select  tx (include longest)
  }
  res.tbl.exons$transcript_id <- factor(res.tbl.exons$transcript_id,levels = c(txID.old,unique(res.tbl.exons$transcript_id[res.tbl.exons$transcript_id!=txID.old])))
}else{
  if (length(table(res.tbl.exons$transcript_id)) > display_tx_num){
    res.tbl.exons <- res.tbl.exons[res.tbl.exons$transcript_id %in% c(sample(unique(res.tbl.exons$transcript_id[res.tbl.exons$transcript_id!=txID.old]),display_tx_num,replace = F)),] # only select  tx (include longest)
  }
  res.tbl.exons$transcript_id <- factor(res.tbl.exons$transcript_id,levels = c(unique(res.tbl.exons$transcript_id[res.tbl.exons$transcript_id!=txID.old])))
}

# str(res.tbl.exons)
# table(res.tbl.exons$transcript_id)
tmp <- res.tbl.exons %>%
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_id
  )) +
  geom_range(
    # aes(fill = transcript_type)
    fill=fill
  ) +
  geom_intron(
    data = to_intron(res.tbl.exons, "transcript_id"),
    aes(strand = strand)
  ) +
  bed_theme
return(tmp)
}


inputFile <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/gtf_by_biotype/11RNA.gtf"
gtf0 <- rtracklayer::import(inputFile)

# txID <- "ENST00000377315_____4"
# txID <- "NR_023363_____1"
# txID <- "26266"
txID <- "ENST00000564646_____1"
test <- plotMultiTxTrack(txID_new = txID, gtf = gtf0, fill="salmon")
test


# test ggtranscript --------
library(magrittr)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(ggplot2)
library(ggtranscript)

sod1_annotation %>% head()

# extract exons
sod1_exons <- sod1_annotation %>% dplyr::filter(type == "exon")

sod1_exons %>%
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_name
  )) +
  geom_range(
    aes(fill = transcript_biotype)
  ) +
  geom_intron(
    data = to_intron(sod1_exons, "transcript_name"),
    aes(strand = strand)
  )


# # reduce gap
# sod1_rescaled <- shorten_gaps(
#   sod1_exons, 
#   to_intron(sod1_exons, "transcript_name"), 
#   group_var = "transcript_name"
# )
# 
# sod1_rescaled %>%
#   dplyr::filter(type == "exon") %>%
#   ggplot(aes(
#     xstart = start,
#     xend = end,
#     y = transcript_name
#   )) +
#   geom_range(
#     aes(fill = transcript_biotype)
#   ) +
#   geom_intron(
#     data = sod1_rescaled %>% dplyr::filter(type == "intron"), 
#     arrow.min.intron.length = 200
#   )

# # distinguish CDS from UTR
# # filter for only exons from protein coding transcripts
# sod1_exons_prot_cod <- sod1_exons %>%
#   dplyr::filter(transcript_biotype == "protein_coding")
# 
# # obtain cds
# sod1_cds <- sod1_annotation %>% dplyr::filter(type == "CDS")
# 
# tmp <- sod1_exons_prot_cod %>%
#   ggplot(aes(
#     xstart = start,
#     xend = end,
#     y = transcript_name
#   )) +
#   geom_range(
#     fill = "white",
#     height = 0.25
#   ) +
#   geom_range(
#     data = sod1_cds
#   ) +
#   geom_intron(
#     data = to_intron(sod1_exons_prot_cod, "transcript_name"),
#     aes(strand = strand),
#     arrow.min.intron.length = 500,
#   )
# 
# ggpubr::ggarrange(plotlist = res.list2, ncol = 1,heights = c(7,1,1,1,1,1), nrow = 6,  align = "v")














# test sim data bw ------------------------------
setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev")
# read_bw <- "/BioII/lulab_b/baopengfei/biosoft/art_bin_MountRainier/ENST00000602385.sort.+.bw"
source <- "plasma.t50b50" # tissue,blood,plasma.t50b50

# read_bw <- "test.bw"
# chr <- "ENST00000602385_____1"
# start <- 0
# width <- 541
pre <- "ENST00000385247_t50_b50"
chr <- "ENST00000385247_____1"
start <- 0
width <- 97
# pre <- "ENST00000364849"
# chr <- "ENST00000364849_____1"
# start <- 0
# width <- 89

read_bw <- paste0("/BioII/lulab_b/baopengfei/biosoft/art_bin_MountRainier/archive/",pre,"/gold.",source,".all.read.sort.bw")
bw.res.df <- readBW(bw=read_bw,chr=chr,start=start,width=width)
plotXbreak <- FALSE
title <- ifelse(paste0(chr,":",start,"-",start+width))
p1 <- plotBW(single_base_bw=bw.res.df,chr=chr,start=start,width=width,ylab="",title=chr,plotXbreak = plotXbreak,plotYbreak = T,plotYtext = F)
peak_bed <- paste0("/BioII/lulab_b/baopengfei/biosoft/art_bin_MountRainier/archive/",pre,"/gold.","plasma",".bed")
# peak_bed <- "/BioII/lulab_b/baopengfei/biosoft/art_bin_MountRainier/gold.bed"
bed <- read.table(peak_bed)
bed <- readPeak(bed,chr=chr,start=start,width=width)
p2 <- plotPeak(bed = bed, chr=chr,start=start,width=width,ylab="")
ggpubr::ggarrange(plotlist = list(p1,p2),ncol = 1,heights = c(4,1),  align = "hv")
#ggbio::autoplot(object = "/BioII/lulab_b/baopengfei/biosoft/art_bin_MountRainier/ENST00000602385.sort.+.bed")
#Sushi::plotBed(chrom = chr, chromstart = start, chromend = end, beddata = "/BioII/lulab_b/baopengfei/biosoft/art_bin_MountRainier/ENST00000602385.sort.+.bed")
# tmp <- rtracklayer::import("gold.tissue.all.read.sort.bw")


# test sim data barplot (cfPeak Fig3) ------------------------------
df <- read.table("/BioII/lulab_b/baopengfei/biosoft/art_bin_MountRainier/res.txt",header = T)
df$method <- factor(df$method,levels = c("Piranha","CLIPper","CLAM","exPeak"))
df$ratio <- factor(df$ratio,levels = c(0.5,0.05,0.005))
ggplot(df,aes(x=method,y=value,fill=ratio))+ # method
  geom_bar(stat = "identity", width = 0.6, position = position_dodge(width = 0.8),color="black")+
  # ggsci::scale_color_d3()+
  # ggsci::scale_fill_jama()+
  scale_fill_manual(values = c(alpha("#7C0D0E",0.8),alpha("firebrick",0.7),"salmon")) +
  # scale_fill_brewer(palette = "OrRd",direction = -1) +
  theme_minimal()+
  theme(aspect.ratio = 0.6,
        plot.title = element_text(size = 24,color="black",hjust = 0.5,face="bold"),
        axis.title = element_text(size = 24,color ="black"), 
        axis.title.y = element_text(size = 24,,color ="black"),  # ,angle = 90
        axis.text = element_text(size= 24,color = "black"),
        #panel.grid=element_blank(),
        # panel.grid.major.x=element_blank(),
        # panel.grid.minor.x = element_blank(),
        # panel.grid.major.y = element_line(color = "black"), #size= 1,
        #panel.grid.minor.y = element_blank(),
        # panel.border = element_blank(),
        axis.text.x = element_text(size= 24, angle = 45, vjust = 1.3, hjust=1), # angle = 45, hjust = 1
        # axis.text.y = element_blank(),
        legend.position = "right" , # c(0.2,0.8),#"right",
        legend.text = element_text(size= 24),
        legend.title= element_text(size= 24),
        strip.text = element_text(size= 24) ) + 
  facet_grid(~type)
ggsave("bar.pdf",width = 20,height = 10)



# tes select best tx ------------------------
# setwd("")
#read ref
ref <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt",header = T,sep = "\t", stringsAsFactors = F)
ref <- ref[!is.na(ref$tx.length),]
#ref[1:3,]



meta.list2 <- list()
meta.list2[["GSE71008"]] <- c("SAMN03863524","SAMN03863398","SAMN03863475","SAMN03863560")  # "SAMN03863475","SAMN03863524"; "SAMN03863398","SAMN03863560";
meta.list2[["GSE94582"]] <- c("NEBNext_Lab1","NEBNext_Lab2","N4_B_Lab5","N4_A_Lab5") # "NEBNext_Lab1","NEBNext_Lab2"; "N4_B_Lab5";
meta.list2[["Phospho-RNA-seq"]] <- c("ULMC157_none","ULMC123_none","ULMC157_T4PNK","ULMC123_T4PNK") #"ULMC157_none"; "ULMC157_T4PNK", "ULMC148_T4PNK", "ULMC123_T4PNK"; 

#df[1:3,]
l <- list()
for(dst in names(meta.list2)){
  print(dst)
  for (smp in meta.list2[[dst]]){
    print(smp)
    df <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_all/expeak_by_sample/b5_d50_p1/intersect_8DNA/",smp,".bed6"))
    df$smp <- smp
    l[[smp]] <- df
  }
}

df <- do.call(rbind,l)

df$RNA <- ref$transcript_type[match(df$V1,ref$transcript_id)]
table(df$RNA)
table(df$smp)
#df <- df[df$V5>=20,] # localmax

table(duplicated(df$V1))
df <- df[df$V1 %in% df$V1[duplicated(df$V1)],]
df$txLength <- ref$tx.length[match(df$V1,ref$transcript_id)]

df2 <- as.data.frame(table(df$V1,df$smp))
df2 <- df2[df2$Freq>=2,]
df2$txLength <- ref$tx.length[match(df2$Var1,ref$transcript_id)]
df3 <- as.data.frame(table(df2$Var1))
df3 <- df3[df3$Freq>=2,]
df3$txLength <- ref$tx.length[match(df3$Var1,ref$transcript_id)]
#G033992__chr17___22520709____22523110_pos

enhancer__chr7___69059541____69062941_pos,2500,1000
enhancer__chr9___5092200____5092800_pos,0,600
promoter__chr7___148941042____148941642_pos,0,600
promoter__chr7___148941042____148941642_pos,400,200
SSU____rRNA_Hsa__chr22___11629545____11631288_pos,0,2000
SSU____rRNA_Hsa__chr5___175114739____175115174_neg,0,435
G033992__chr17___22520709____22523110_pos,0,2400
G087360__chrM___14742____15537_neg,0,795

df2 <- df[df$RNA=="mRNA",]
tmp <- as.data.frame(table(df2$V1))
#df <- df[order(df$V5,decreasing = T),]

#novel miR
G004258__chr1___145982791____146044553_pos,52896,52933
G004258__chr1___145982791____146044553_pos,52940,52973

G033992__chr17___8100551____8128648_pos,25191,25226
G033992__chr17___8100551____8128648_pos,25245,25267




# test peak specific
# dst <- "GSE71008"
# smp <- "SAMN03863396"
# df <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_all/expeak_by_sample/b5_d50_p1/intersect_8DNA/",smp,".bed6"))
dst <- "GSE71008_NCpool"
dedup <- "all"
method <- "clam" #"expeak","clipper","clam" {piranha_by_sample/b5_p01,clipper_by_sample/b5_p05}
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
bed <- Sys.glob(paste0(pre,"/output/",dst,"/call_peak_",dedup,"/",method,"_by_sample/b5*/intersect/",smp,"_",method,"Only.bed6"))
df <- read.table(bed)

df$RNA <- ref$transcript_type[match(df$V1,ref$transcript_id)]
table(df$RNA)
df$txLen <- ref$tx.length[match(df$V1,ref$transcript_id)]

df2 <- as.data.frame(table(df$V1))
colnames(df2) <- c("V1","V2")
df2$RNA <- ref$transcript_type[match(df2$V1,ref$transcript_id)]
table(df2$RNA)
df2$txLen <- ref$tx.length[match(df2$V1,ref$transcript_id)]


# #expeak
# NR_146154_____1,0,5055
# NR_146152_____1,0,500
# NR_146148_____1,0,5054
# ENST00000636484_____1,0,328
# ENST00000600213_____3,0,1049

# #clipper
# ENST00000571127_____1,2,mRNA,791
# ENST00000545075_____2,3,mRNA,1530
# 
# ENST00000319763_____1,1,lncRNA,2081
# ENST00000391625_____2,1,lncRNA,3894
# ENST00000414890_____1,1,lncRNA,737

#clam
ENST00000521141_____1,1,lncRNA,1630
T004607,1,tucpRNA,1634
T234700,2,tucpRNA,10167
30032,1,tRNA,72
ENST00000310125_____4,1,mRNA,1680

# #niche of clipper ?
# df <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_all/expeak_by_sample/b5_d50_p1/intersect/",smp,".bed6"))
# df$RNA <- ref$transcript_type[match(df$V1,ref$transcript_id)]
# table(df$RNA)
# df$txLen <- ref$tx.length[match(df$V1,ref$transcript_id)]
# table(df$V3==df$txLen)
