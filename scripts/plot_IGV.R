# plot peak cov,bed IGV
# last 221003 by bpf

library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggbreak)
# install.packages("ggbreak")
library(gg.gap)
library(IRanges)
library(GenomicRanges)
# install.packages ('gg.gap') 
#ggbreak::scale_x_break(c(7, 17)) 
options(stringsAsFactors = F)
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("which", "base")
setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev")
# setwd("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/")

# def func. ------------------------------------------------------
old2newTxID <- function(x){
  x <- gsub("(-)","_neg",x,fixed = T)
  x <- gsub("(+)","_pos",x,fixed = T)
  x <- gsub("::","__",x,fixed = T)
  x <- gsub(":","___",x,fixed = T)
  x <- gsub("-","____",x,fixed = T)
  x <- gsub(".","_____",x,fixed = T)
  return(x)
}
new2oldTxID <- function(x){
  # x <- txID
  x <- gsub("_neg","(-)",x,fixed = T)
  x <- gsub("_pos","(+)",x,fixed = T)
  x <- gsub("_____",".",x,fixed = T)
  x <- gsub("____","-",x,fixed = T)
  x <- gsub("___",":",x,fixed = T)
  x <- gsub("__","::",x,fixed = T)
  return(x)
}

GRange2bed <- function(gr){
  # gr <- hg38[[1]] # 1-base to 0-base
  bed <- data.frame("chr"=gr@seqnames,"start"=gr@ranges@start-1,"end"=gr@ranges@start,"name"=gr$name,"score"=gr$score,"strand"=gr@strand)
  return(bed)
}

bw_theme <- theme_void() + 
  theme(plot.title = element_text(size = 20,color="black",hjust = 0.5,face="bold"),
        axis.title = element_text(size = 20,color ="black"), 
        axis.title.y = element_text(size = 20,color ="black"),  #, vjust = 1,angle = 90
        axis.line.y.left = element_blank(),, # colour, size, linetype, lineend, color
        axis.line.x.bottom = element_blank(), 
        axis.ticks = element_line(color = "black",linewidth = 1, size = 1),
        #panel.grid=element_blank(),
        # panel.grid.major.x=element_blank(),
        # panel.grid.minor.x = element_blank(),
        # panel.grid.major.y = element_line(color = "black"), #size= 1,
        #panel.grid.minor.y = element_blank(),
        # panel.border = element_blank(),
        axis.text.x = element_text(size= 24,color = "black"), # angle = 45, hjust = 1 
        axis.text.y = element_text(size= 24,color = "black"),#, angle=90
        legend.position = "right",
        legend.text = element_text(size= 16),
        legend.title= element_text(size= 16),
        strip.text = element_text(size= 16) ) #+
# facet_grid(sample~.,scales = "free_y")
bed_theme <-  theme_void() + 
  theme(plot.title = element_text(size = 20,color="black",hjust = 0.5,face="bold"),
        axis.title = element_text(size = 20,color ="black"), 
        axis.title.y = element_text(size = 20,color ="black"),  # ,angle = 90
        #panel.grid=element_blank(),
        # panel.grid.major.x=element_blank(),
        # panel.grid.minor.x = element_blank(),
        # panel.grid.major.y = element_line(color = "black"), #size= 1,
        #panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_blank(), # angle = 45, hjust = 1
        axis.text.y = element_blank(),
        legend.position = "right",
        legend.text = element_text(size= 16),
        legend.title= element_text(size= 16),
        strip.text = element_text(size= 16) )#+
# facet_grid(sample~.,scales = "free_y")

readBW <- function(bw,chr,start,width){
  #bw <- read_bw
  bg <- rtracklayer::import.bw(bw, which=GRanges(c(chr), IRanges(start = start+1, width = width))) # warning if no record found, IRanges is 1-based
  bg <- as.data.frame(bg)[,1:6]
  bg$start <- bg$start-1 # convert to 0-based bed: no err if no rec; start can be 0 !
  left <- bg$start
  wid <- bg$width
  score <- bg$score
  
  bg <- bg[,c("seqnames","start","end","width","score","strand")]
  # bg$width <- "X"
  colnames(bg)[4] <- "name"
  # head(bg,3)
  # bg <- bg %>% 
  #   dplyr::filter(seqnames==chr) # %>% 
  # bg <- bg[bg$seqnames==chr & bg$start<=start+width & bg$end>=start,] # rtracklayer::import.bw: GRanges
   
  if (nrow(bg)==0){
    bw.res.df <- data.frame("seqnames"=chr, "start"=start:(start+width-1), "end"=(start+1):(start+width), "name"="X", "score"=0, "strand"=".")
    return(bw.res.df)
  }
  
  # bw.res <- list()
  df <- data.frame("seqnames"=chr, "start"=start:(start+width-1), "end"=(start+1):(start+width), "name"="X", "score"=0, "strand"=".")
  for (i in 1:nrow(bg)){
    #rec <- bg[i,]
    # i <- 1
    # print(i)
    s <- left[i] # rec$start
    w <- wid[i]
    df[which(df$start==s):which(df$start==s+w-1),"score"] <- score[i] # 1-based
    # df <- data.frame("seqnames"=rec$seqnames, "start"=s:(s+w-1), "end"=(s+1):(s+w), "name"="X", "score"=rec$score, "strand"=rec$strand)
    # bw.res[[i]] <- df
  }
  # bw.res.df <- do.call(rbind,bw.res)
  return(df)
}

plotBW <- function(single_base_bw,chr,start,width,ylab,xlab="",title="",color="steelblue4",fill="steelblue4",plotXbreak=TRUE,plotYbreak=TRUE,plotYtext=FALSE,annotate.size=8,maxY=-1,libsize=0){
  #single_base_bw <- bw.res.df
  #bw.res.df$score <- 1
  # ylab="test"
  # xlab=""
  # title=""
  # color="steelblue4"
  # fill="steelblue4"
  # plotXbreak=TRUE
  # plotYbreak=FALSE
  # plotYtext=TRUE
  # start <- min(x$start)
  # width <- width-1
  # single_base_bw <- single_base_bw[single_base_bw$seqnames==chr & single_base_bw$start<=start+width & single_base_bw$end>=start,]
  single_base_bw <- as.data.frame(single_base_bw[,1:6])
  colnames(single_base_bw) <- c("seqnames","start","end","width","score","strand")
  if (nrow(single_base_bw)==0){
    single_base_bw <- data.frame("seqnames"=chr, "start"=start:(start+width-1), "end"=(start+1):(start+width), "name"="X", "score"=0, "strand"=".") # 0-based
    # color <- "white"
    fill <- "white"
  }
  if(libsize>0){
    # single_base_bw$score <- as.double(single_base_bw$score/libsize*1000000) # normalized by external lib.size
    single_base_bw$score <- (single_base_bw$score/libsize*1000000) # round(s,digits = 1)
  }

  if(plotXbreak){
    x_scale <- scale_x_continuous(breaks = c(start,start+width), limits = c(start,start+width) ) # 1-based
  }else{
    x_scale <- scale_x_continuous(breaks = NULL, limits = c(start,start+width) )
  }
  

  if(maxY==-1){ # max value
    if(max(single_base_bw$score)==0){
      maxY <- 1 # enough even for cpm
    }else{
      maxY <- max(single_base_bw$score) # round(,digits = 1)
    }
  }
  
    #use given max yaxis
  if(plotYbreak){
    y_scale <- scale_y_continuous(breaks = c(0,maxY), limits = c(0,maxY), labels = c("0",format(maxY,digits = 1,nsmall=0)) ) # 1-based
  }else{
    y_scale <- scale_y_continuous(breaks = NULL, limits = c(0,maxY), labels = c("0",format(maxY,digits = 1,nsmall=0)) )
  }
  

  if(plotYtext){
    anno <- annotate(geom="text",x=(start+1),y=(0.8*maxY), label=paste0("[0-",format(maxY,digits = 1,nsmall=0),"]"), color="grey30",size = annotate.size,hjust = 0) #  hjust = 0
      #annotate(geom="text",xmin=(start+1),xmax=(start+1+width*0.2), ymin=(0.8*max(x$score-1,1)),ymax=max(x$score-1,1), label=paste0("[0-",max(x$score,1),"]"), color="black")
  }else{
    anno <- annotate(geom="text",x=(start+1),y=(0.8*maxY), label="", color="grey30",size = annotate.size,hjust = 0) #
  }

  tmp <- ggplot(single_base_bw, aes(x=start,y=(score))) +  # , group=sample
    # geom_density(stat = "identity", color="steelblue4", fill="steelblue4")+ # salmon
    geom_bar(stat="identity",color=color,fill=fill) + # position = 'dodge', "steelblue4"
    # geom_rect(aes(xmin=start,xmax=start+1,ymin=0,ymax=score))+
    # ggsci::scale_color_d3()+
    # ggsci::scale_fill_d3()+
    ylab(ylab)+ # (log2)
    xlab(xlab)+
    labs(title=title) +
    x_scale +
    y_scale + # trans = ""
    anno +
    geom_hline(yintercept = c(0),color=color)+
    # ylim(c(0,5))+
    # xlim(obj = )+
    bw_theme
  return(tmp)
}


readPeak <- function(bed,chr,start,width){
  #bed <- localmaxEM2
  #tbl only keep exist position, not like bw plot (every pos) !!!
  
  bed <- as.data.frame(bed[,1:6])
  colnames(bed) <- c("seqnames","start","end","width","score","strand")
  bed <- bed[bed$seqnames==chr & bed$start<=start+width & bed$end>=start,]
  
  if (nrow(bed)==0){
    bed.res.df <-  data.frame("seqnames"=".", "start"=0, "end"=0, "name"="X", "score"=0, "strand"=".") # 0-based
    return(bed.res.df)
  }
  
  left <- bed$start
  wid <- bed$end-bed$start
  score <- bed$score 
  strand <- bed$strand
  
  # pivot_longer(cols = start:end, names_to = "type", values_to = "pos") %>% 
  # filter(!duplicated(pos)) #%>% 
  # mutate(sample=id,group=paste0(g1,g2))
  # s <- bed$start
  # w <- bed$width
  bed.res <- list()
  for (i in 1:nrow(bed)){
    # rec <- bed[i,]
    s <- left[i] # rec$start
    w <- wid[i] # rec$width
    df <- data.frame("seqnames"=chr, "start"=s:(s+w-1), "end"=(s+1):(s+w), "name"="X", "score"=1, "strand"=".") # get single base bed:  ignore socre col for peak bed
    bed.res[[i]] <- df
  }
  #res[[2]]
  bed.res.df <- do.call(rbind,bed.res)
  return(bed.res.df)
}
plotPeak <- function(bed,chr,start,width,ylab,color="grey30",fill="grey30"){
  #x <- expeak_path
  # bed <- RBPs.res.df
  bed <- as.data.frame(bed[,1:6])
  colnames(bed) <- c("seqnames","start","end","width","score","strand")
  bed <- bed[bed$seqnames==chr & bed$start<=start+width & bed$end>=start,]
  
  if (nrow(bed)==0){
    bed <- data.frame("seqnames"=".", "start"=0, "end"=0, "name"="X", "score"=0, "strand"=".")
    color <- "white"
    fill <- "white"
  }
  tmp <- ggplot(bed) +  # , aes(x=start+1), group=sample
    # geom_line() +
    # geom_segment(bed2, aes(x = start, y = score, xend = end, yend = score), fill = "blue",colour = "blue",color = "blue") + #aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment")
    # geom_bar() +
    # geom_density(color="white", fill="steelblue4")+ # salmon
    geom_rect(aes(xmin=start,xmax=start+1,ymin=0,ymax=score), color=color, fill=fill)+
    ylab(ylab) +
    xlab("") +
    geom_hline(yintercept = c(-0.1,1.1),color="grey90") +
    scale_y_continuous(breaks=c(-0.1,1.1), limits=c(-0.1,1.1)) + # trans = ""
    scale_x_continuous(breaks=c(start, start+width), limits = c(start,start+width)) +
    # ylim(c(0,100))+
    # xlim(obj = c(start,start+width))+
    bed_theme
  return(tmp)
}

readAnno <- function(bed,chr,start,width) {
  # x <- bed_bed
  bed <- data.table::fread(bed,data.table = F,sep = '\t',check.names = F,stringsAsFactors = F) # diff from readPeak: read direct from filename
  bed.res.df <- readPeak(bed = bed, chr = chr, start = start, width = width)
  return(bed.res.df)
}
plotAnno <- plotPeak
# plotAnno <- function(x,lab){
#   tmp <- ggplot(x) +  # , aes(x=start+1), group=sample
#     # geom_line() +
#     # geom_segment(bed2, aes(x = start, y = score, xend = end, yend = score), fill = "blue",colour = "blue",color = "blue") + #aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment")
#     # geom_bar() +
#     # geom_density(color="white", fill="steelblue4")+ # salmon
#     geom_rect(aes(xmin=start,xmax=start+1,ymin=0,ymax=score),color="salmon", fill="salmon")+
#     ylab(lab)+
#     xlab("")+
#     geom_hline(yintercept = c(-0.5,1.5),color="grey90")+
#     # ylim(c(0,100))+
#     xlim(obj = c(start,start+width))+
#     bed_theme
#   return(tmp)
# }


plotMultiTxTrack <- function(txID_new,gtf,fill="grey50"){
  #txID <- "NR_146154_____1"
  # txID_new = chr, gtf = gtf0
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
  # table(gtf$gene_name==gID)
  # head(gtf,4)
  gtf <- as.data.frame(gtf)
  if(is.rRNA){
    gID <- txID.old2
    #table(gtf$gene_name==gID)
    #table(gtf$transcript_type)
    #gtf2 <- gtf[grepl("RNA5S",gtf$gene_name),]
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
    bed_theme +
    ylab("")
  
  return(tmp)
}


bam2bw <- function(bamFile,outBwPath){
  # a little slow in R (mem constriant ?)
  library(GenomicAlignments)
  library(rtracklayer)
  # bamFile <- system.file("extdata", "ex1.bam", package="Rsamtools", mustWork=TRUE)
  alignment <- readGAlignments(bamFile)
  reads_coverage <- coverage(alignment)
  # return(reads_coverage)
  export.bw(reads_coverage, con = outBwPath)
}



## get rev complement seq
seq_rev <- function(char) {
  alphabets <- strsplit(char, split = "")[[1]]
  return(rev(alphabets))
}

seq_compl <- function(seq) {
  # Check if there's "T" in the sequence
  RNA <- Reduce(`|`, seq == "U")
  cmplvec <- sapply(seq, function(base) {
    # This makes DNA the default
    # As long as there's no U, the sequence is treated as DNA
    if (RNA) {
      switch(base, "A" = "U", "C" = "G", "G" = "C", "U" = "A")
    } else {
      switch(base, "A" = "T", "C" = "G", "G" = "C", "T" = "A")
    }
  })
  return(paste(cmplvec, collapse = ""))
}



# 1.plot multi sample bigwig only ------------------------------------------
setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev")
meta.list <- list()
#meta.list <- meta.list2

#meta.list[["i-pico"]] <- c("L38-4-BC5","L39-5-BC2","L39-5-BC3","L41-7-BC5","L40-6-BC2","L38-4-BC3") # i-pico: DT1
meta.list[["GSE50676"]] <- c("GM12878","GM12891") # CLIP
# meta.list[["GSE52600"]] <- c("SRR1035214","SRR1035216" # RIP-seq
#                              # "SRR1035213","SRR1035215")
meta.list[["TCGA_small"]] <- c("TCGA-A7-A0CE-11A-21R-A090-13_mirna_gdc_realn", # TCGA_small_BRCA: BRCA,adj_NC
                               "TCGA-A7-A0D9-11A-53R-A090-13_mirna_gdc_realn"
                               #"TCGA-A7-A0CE-01A-11R-A010-13_mirna_gdc_realn","TCGA-A7-A0D9-01A-31R-A057-13_mirna_gdc_realn",
                               ) 

meta.list[["GSE148861_GSE148862"]] <- c("SRR11563522","SRR11563463") # WBC,sncRNA

# meta.list[["DT2_short"]] <- c("DT2-CRC-2_2","DT2-CRC-3_2") # detector2
#meta.list[["TGIRT-seq_2020_eLife"]] <- c("SRR12047433","SRR12047420") # SRR12047433,SRR12047420 # PNK,SRR12047409 # smart olig dT
meta.list[["TGIRT-seq_2020_eLife_short"]] <- c("SRR12047433","SRR12047432","SRR12047420","SRR12047419") # SRR12047433,SRR12047432,SRR12047419(PNK),SRR12047420(PNK),SRR12047409(smart olig dT, seem little signal)
meta.list[["GSE110381"]] <- c("SRR6700252","SRR6700253") #
meta.list[["GSE123972"]] <- c("SRR8507737","SRR8507736") # TC
# meta.list[["AGO2_IP"]] <- c("AGO2_IP_SUP_s1","AGO2_IP_SUP_s2") # exRNA-AGO2-IP
# meta.list[["WSQ_SMARTer_NEB"]] <- c("NC_PKU-2392860_smart_1","NC_PKU-2397512_smart_1","NC_PKU-2392860_smart_PNK_1","NC_PKU-2397512_smart_PNK_1","NC_PKU-2392860_1","NC_PKU-2397512_1")
meta.list[["Phospho-RNA-seq"]] <- c("ULMC123_none","ULMC134_none","ULMC123_T4PNK","ULMC134_T4PNK") #ULMC123_T4PNK (GSE126051)
meta.list[["GSE71008"]] <- c("SAMN03863523","SAMN03863521")   #c("SRR2105340","SRR2105335")
# GSE110381 seem has no signal in miR and mRNA eg tx
meta.list[["GSE94582"]] <- c("NEBNext_Lab1","NEBNext_Lab2","TruSeq_Lab1","TruSeq_Lab2","N4_A_Lab5","N4_B_Lab5","CleanTag_Lab5") # 

meta.col <- list() 
meta.col[["i-pico"]] <- c("#2C68A9") 
meta.col[["GSE50676"]] <- c("#FD6905") # ,"","steelblue","","steelblue3","steelblue4"
meta.col[["TCGA_small"]] <- c("#843692")
meta.col[["GSE148861_GSE148862"]] <- c("#843692")
meta.col[["GSE110381"]] <- c("#2C68A9") 
meta.col[["GSE123972"]] <- c("#2C68A9") 
#meta.col[["TGIRT-seq_2020_eLife"]] <- c("#2C68A9")
meta.col[["TGIRT-seq_2020_eLife_short"]] <- c("#2C68A9")
# meta.col[["AGO2_IP"]] <- c("steelblue4")
meta.col[["Phospho-RNA-seq"]] <- c("#2C68A9")
meta.col[["GSE71008"]] <- c("#2C68A9") 
meta.col[["GSE94582"]] <- c("#2C68A9")
# meta.col[["DT2_short"]] <- c("#2C68A9")
# meta.col[["GSE50676"]] <- c("steelblue4")
#RColorBrewer::display.brewer.all()


genes.list <- list()
##Fig1
#genes.list <- genes.list2
# genes.list[[1]] <- data.frame(chr="ENST00000385025_____1",start=0,width=91) # 4 contribution overlap MIR
# genes.list <- genes.list2
# genes.list[[1]] <- data.frame(chr="ENST00000385247_____1",start=0,width=97)
#genes.list[[1]] <- data.frame(chr="ENST00000362309_____3",start=0,width=100) # let7i
genes.list[[1]] <- data.frame(chr="ENST00000385227_____1",start=0,width=89) # miR/MIR30C1, RBP + EV
genes.list[[2]] <- data.frame(chr="NR_146151_____1",start=8800,width=500) # rRNA, EV + G4，13309
#genes.list[[4]] <- data.frame(chr="ENST00000624233_____1",start=800,width=400) # mRNA, EV,1426
genes.list[[4]] <- data.frame(chr="ENST00000622269_____1",start=200,width=200) # mRNA, EV
genes.list[[3]] <- data.frame(chr="ENST00000315707_____3",start=400,width=300) # lncRNA, RBP,2000
# genes.list[[4]] <- data.frame(chr="12969",start=0,width=73) # tRNA, EV
# genes.list[[5]] <- data.frame(chr="ENST00000217398_____3",start=0,width=979) # mRNA, G4 (high heter)
# genes.list[[1]] <- data.frame(chr="T227251",start=2000,width=500) # mRNA, G4:  ENST00000257359_____6   4010


# #select EV eg in Fig2
# #bedtools intersect -u -s -a /BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/lulab/11RNAsmallDomain_EVenrich.bed -b ../output/GSE123972/call_peak_all/expeakCNN_by_sample/b5_d50_p1/SRR8507737.bed > tmp/EVoverlap_SRR8507737.bed # overlap with exRNA
# #bedtools intersect -v -s -a /BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/lulab/11RNAsmallDomain_EVenrich.bed -b ../output/GSE50676/call_peak_all/expeakCNN_by_sample/b5_d50_p1/GM12891.bed > tmp/EVoverlap_GM12891.bed # not overlap with CLIP
# #bedtools intersect -v -s -a /BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/lulab/11RNAsmallDomain_EVenrich.bed -b /BioII/lulab_b/baopengfei/shared_reference/RBP/splitByRBP/RBPs_tx.bed > tmp/EVoverlap_RBP.bed # not overlap with RBP
# #bedtools intersect -v -s -a /BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/lulab/11RNAsmallDomain_EVenrich.bed -b /BioII/lulab_b/baopengfei/shared_reference/structure/quadratlas_g4grinder_tx.bed > tmp/EVoverlap_G4.bed # not overlap with G4
# inter.list <- list()
# #"CleanTag_Lab5","N4_A_Lab5","N4_B_Lab5","NEBNext_Lab1","NEBNext_Lab2","TruSeq_Lab1","TruSeq_Lab2","SAMN03863521","SAMN03863523","SRR8507737","SRR8507736","ULMC123_none","ULMC134_none","ULMC123_T4PNK","ULMC134_T4PNK","GM12878","GM12891","G4","RBP"
# for (i in c("CleanTag_Lab5","N4_A_Lab5","NEBNext_Lab1","TruSeq_Lab1","SAMN03863521","SRR8507737","ULMC123_none","ULMC123_T4PNK","GM12878","GM12891","G4","RBP")){
#   #i <- "CleanTag_Lab5"
#   print(i)
#   tmp <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/tmp/EVoverlap_",i,".bed"),header = F)
#   inter.list[[i]] <- tmp$V4
# }
# l <- Reduce(intersect,inter.list)
# tmp <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/tmp/EVoverlap_","CleanTag_Lab5",".bed"),header = F)
# tmp <- tmp[tmp$V4 %in% l,]
# 
# #ref <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt") #,header=T
# tmp$RNA <- ref$transcript_type[match(tmp$V1,ref$transcript_id)]
# tmp$txLen <- ref$tx.length[match(tmp$V1,ref$transcript_id)]
# #choose EV peak in mRNA in priority




# #EV
# # NR_146151_____1,4465,4572
# # NR_146151_____1,9293,9318
# # NR_146151_____1,11889,12006
# #G4
# # NR_146151_____1,9042,9101
# # NR_146151_____1,11731,11862
# 
# genes.list[[1]] <- data.frame(chr="NR_146151_____1",start=4000,width=1000) 
# genes.list[[2]] <- data.frame(chr="NR_146151_____1",start=8500,width=1000) 
# genes.list[[3]] <- data.frame(chr="NR_146151_____1",start=11500,width=1000) 



# genes.list[[2]] <- data.frame(chr="NR_023363_____1",start=0,width=150)
# genes.list[[2]] <- data.frame(chr="26266",start=0,width=75)
# genes.list[[4]] <- data.frame(chr="ENST00000564646_____1",start=0,width=2873)
# genes.list[[1]] <- data.frame(chr="ENST00000613283_____2",start=0,width=1365) # MYC (poor cov)
# genes.list[[1]] <- data.frame(chr="ENST00000623130_____1",start=0,width=84332) # many peak on lnc
# genes.list[[2]] <- data.frame(chr="ENST00000626826_____1",start=0,width=205012) # many peak on lnc


## fig4 (all 4 types eg)
meta.list <- list()
# meta.list[["DT2"]] <- c("DT2-CRC-2_2","DT2-CRC-3_2") # detector2
meta.list[["TGIRT-seq_2020_eLife_short"]] <- c("SRR12047433","SRR12047432","SRR12047420","SRR12047419") # SRR12047433,SRR12047432,SRR12047419(PNK),SRR12047420(PNK),SRR12047409(smart olig dT, seem little signal)
meta.list[["GSE71008"]] <- c("SAMN03863524","SAMN03863398","SAMN03863475","SAMN03863560")  # "SAMN03863475","SAMN03863524"; "SAMN03863398","SAMN03863560";
meta.list[["GSE94582"]] <- c("NEBNext_Lab1","NEBNext_Lab2","N4_B_Lab5","N4_A_Lab5") # "NEBNext_Lab1","NEBNext_Lab2"; "N4_B_Lab5";
meta.list[["Phospho-RNA-seq"]] <- c("ULMC157_none","ULMC123_none","ULMC157_T4PNK","ULMC123_T4PNK") #"ULMC157_none"; "ULMC157_T4PNK", "ULMC148_T4PNK", "ULMC123_T4PNK";
meta.col <- list()
meta.col[["TGIRT-seq_2020_eLife_short"]] <- c("#2C68A9")
meta.col[["GSE71008"]] <- c("#2C68A9")
meta.col[["GSE94582"]] <- c("#2C68A9")
meta.col[["Phospho-RNA-seq"]] <- c("#2C68A9")
# meta.col[["DT2"]] <- c("#2C68A9")

genes.list <- list()
#genes.list[[1]] <- data.frame(chr="SSU____rRNA_Hsa__chr22___11629545____11631288_pos",start=500,width=200)
genes.list[[1]] <- data.frame(chr="G033992__chr17___333488____338834_pos",start=450,width=200) # Fig4: miR
#genes.list[[2]] <- data.frame(chr="SSU____rRNA_Hsa__chr5___175114739____175115174_neg",start=0,width=435)
genes.list[[2]] <- data.frame(chr="enhancer__chr7___69059541____69062941_pos",start=2800,width=600) # deprecated ?
#genes.list[[4]] <- data.frame(chr="enhancer__chr9___5092200____5092800_pos",start=0,width=600)
genes.list[[1]] <- data.frame(chr="promoter__chr7___148941042____148941642_pos",start=300,width=300) # SuppFig7: RSS
#genes.list[[6]] <- data.frame(chr="G087360__chrM___14742____15537_neg",start=0,width=795)
#genes.list[[4]] <- data.frame(chr="G033992__chr17___22520709____22523110_pos",start=1200,width=1200)
#genes.list[[5]] <- data.frame(chr="AluJb__chr15___68202323____68202601_pos",start=0,width=278)
#genes.list[[6]] <- data.frame(chr="G____rich__chr6___7146362____7146415_pos",start=0,width=53)
genes.list[[4]] <- data.frame(chr="G____rich__chr6___7146362____7146415_pos",start=0,width=103) # SuppFig7: exp G4
genes.list[[5]] <- data.frame(chr="enhancer__chr11___62853176____62855776_neg",start=2321,width=192) # Fig4: exp C/D snoRNA (3fold)
genes.list[[1]] <- data.frame(chr="enhancer__chr2___202619187____202620387_neg",start=346,width=198) # SuppFig7: exp tRNA pseu
genes.list[[1]] <- data.frame(chr="enhancer__chr6___27752463____27755463_neg",start=1918,width=219) # Fig4: exp tRNA standard

#genes.list[[1]] <- data.frame(chr="G033992__chr17___22520709____22523110_pos",start=2139,width=51) # k-turn (3fold)
#genes.list[[1]] <- data.frame(chr="G033992__chr17___22520709____22523110_pos",start=2022,width=285) # k-turn (5fold)
# genes.list[[1]] <- data.frame(chr="enhancer__chr11___62853176____62855776_neg",start=2365,width=110) # k-turn (5fold)
genes.list[[2]] <- data.frame(chr="G078490__chr7___52897807____52917741_neg",start=2720,width=114) # k-turn


# ## fig4 (miRDeep2)
# meta.list <- list()
# meta.list[["GSE71008"]] <- c("SAMN03863524","SAMN03863398","SAMN03863475","SAMN03863560")  # "SAMN03863475","SAMN03863524"; "SAMN03863398","SAMN03863560";
# meta.list[["GSE94582"]] <- c("NEBNext_Lab1","NEBNext_Lab2","N4_B_Lab5","N4_A_Lab5") # "NEBNext_Lab1","NEBNext_Lab2"; "N4_B_Lab5";
# meta.list[["Phospho-RNA-seq"]] <- c("ULMC157_none","ULMC123_none","ULMC157_T4PNK","ULMC123_T4PNK") #"ULMC157_none"; "ULMC157_T4PNK", "ULMC148_T4PNK", "ULMC123_T4PNK";
# meta.col <- list()
# meta.col[["GSE71008"]] <- c("#2C68A9")
# meta.col[["GSE94582"]] <- c("#2C68A9")
# meta.col[["Phospho-RNA-seq"]] <- c("#2C68A9")
# genes.list <- list()
# #genes.list[[1]] <- data.frame(chr="G078490__chr7___128168336____128234552_neg",start=26400,width=400)
# genes.list[[1]] <- data.frame(chr="G033992__chr17___27274477____27280481_pos",start=2950,width=200)
# genes.list[[2]] <- data.frame(chr="G033992__chr17___333488____338834_pos",start=450,width=200)





## fig7 MvsL diff eg
#diff3 <- diff[diff$log2FoldChange>0,]
meta.list <- list()
meta.list[["lulab_oscc_plasma_diff"]] <- c("W2","W4","W5","Y5","Y6","Y7")
meta.col <- list()
meta.col[["lulab_oscc_plasma_diff"]] <- c("#2C68A9")
genes.list <- list()
# T031004_14945_14959_+|tucpRNA|T031004|peak_30622|T031004|14945|14959
#T031004,14945,14, 
genes.list[[1]] <- data.frame(chr="T031004",start=14900,width=150) # canditate 1
# ENST00000606662_____1_1514_1531_+|lncRNA|ENST00000606662_____1|peak_17884|ENST00000606662_____1|1514|1531
genes.list[[1]] <- data.frame(chr="ENST00000606662_____1",start=1500,width=150) # canditate 2
#We constructed a signature composed of nine ferroptosis-related lncRNAs AC034236.2 (ENSG00000271918.1)
#ferroptosis is involved in tumor-host interactions, modulates the tumor microenvironment, and serves as an antimetastatic mechanism
#overlapped with TARDBP RBP: major component of ubiquitinated inclusions found in neurons and glial cells in ALS patients

# ENST00000526482_____1_1520_1534_+|mRNA|ENST00000526482_____1|peak_10921|ENST00000526482_____1|1520|1534
genes.list[[1]] <- data.frame(chr="ENST00000526482_____1",start=1500,width=150) # canditate 3 !!!!!
#overlapped with G4,  ENSG00000140090 (symbol: NCKX4,SHEP6,SLC24A2) associate with fasting plasma glucose

# G078490__chr7___86440584____86773221_neg_129137_129155_+|intron_rev|G078490__chr7___86440584____86773221_neg|peak_27474|G078490__chr7___86440584____86773221_neg|129137|129155
genes.list[[1]] <- data.frame(chr="G078490__chr7___86440584____86773221_neg",start=129000,width=200) # canditate 4 !!!
#gn: chr7:86644021-86644221    overlped with RNA polII, nearest gene GRM3 dysregulate cAMP signaling, cAMP signaling has been implicated in melanoma progression and drug resistance

# promoter__chr12___59464285____59466885_pos_403_417_+|promoter_for|promoter__chr12___59464285____59466885_pos|peak_36636|promoter__chr12___59464285____59466885_pos|403|417
# genes.list[[1]] <- data.frame(chr="promoter__chr12___59464285____59466885_pos",start=300,width=120) # 
#gn: chr12:59464585-59464705   no much gn signal

# enhancer__chr13___33743263____33752263_pos_413_431_+|enhancer_for|enhancer__chr13___33743263____33752263_pos|peak_32641|enhancer__chr13___33743263____33752263_pos|413|431
# genes.list[[1]] <- data.frame(chr="enhancer__chr13___33743263____33752263_pos",start=300,width=150) # 
#gn:  chr13:33743563-33752713, no sig in ucsc

# MIRb__chr17___34133669____34133797_neg_58_74_+|repeats_rev|MIRb__chr17___34133669____34133797_neg|peak_29152|MIRb__chr17___34133669____34133797_neg|58|74
genes.list[[1]] <- data.frame(chr="MIRb__chr17___34133669____34133797_neg",start=0,width=150) # candidate 5 !!!
#gn chr17:34133647-34133797  overlap with gene RP11-17M24.1 (ENSG00000265356, minus strand) m6A-induced lncRNA RP11 triggers the dissemination of colorectal cancer ...

#limit to up-regulated
# ENST00000420315_____1_338_349_+|lncRNA|ENST00000420315_____1|peak_4127|ENST00000420315_____1|338|349
#genes.list[[1]] <- data.frame(chr="ENST00000420315_____1",start=300,width=100) #
#Gene: RP11-255A11.2 (ENSG00000228072)

# ENST00000635227_____1_316_328_+|lncRNA|ENST00000635227_____1|peak_23159|ENST00000635227_____1|316|328
genes.list[[1]] <- data.frame(chr="ENST00000635227_____1",start=300,width=100) # candidate 
# Gene: LINC01108 (ENSG00000226673)

#enhancer__chr14___92232950____92236950_neg_1303_1317_+|enhancer_rev|enhancer__chr14___92232950____92236950_neg|peak_32833|enhancer__chr14___92232950____92236950_neg|1303|1317
#genes.list[[1]] <- data.frame(chr="enhancer__chr14___92232950____92236950_neg",start=1300,width=100)
#gn: chr14:92235550-92235650 no much ucsc signal

# promoter__chr5___2684486____2691286_pos_4370_4384_+|promoter_for|promoter__chr5___2684486____2691286_pos|peak_37997|promoter__chr5___2684486____2691286_pos|4370|4384
# genes.list[[1]] <- data.frame(chr="promoter__chr5___2684486____2691286_pos",start=4300,width=100)  # candidate
#gn:  chr5:2688786-2688886  overlap with gencode sncRNA, no nearest gene ?
#

#MLT1K__chr20___61804800____61805198_neg_185_199_+|repeats_rev|MLT1K__chr20___61804800____61805198_neg|peak_29266|MLT1K__chr20___61804800____61805198_neg|185|199
genes.list[[1]] <- data.frame(chr="MLT1K__chr20___61804800____61805198_neg", start=150, width=100)  # candidate !!!!
#gn:  chr20:61804948-61805048   no much ucsc signal, overlapping with CDH4: inhibits ferroptosis in oral squamous cell carcinoma cells
libsizes <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/lulab_oscc_plasma_diff/call_peak_dedup/count_matrix/EM.libSize",header = F)$V1
smpids <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/lulab_oscc_plasma_diff/call_peak_dedup/count_matrix/cfpeakCNN_b5_d50_p1.txt",check.names = F,header = T)
smpids <- smpids[,2:ncol(smpids)]
smpids <- colnames(smpids)
names(libsizes) <- smpids
#CPM coverage
for (i in 1:length(genes.list)){
  # i <- 2
  chr <- genes.list[[i]]$chr
  start <- genes.list[[i]]$start
  width <- as.numeric(genes.list[[i]]$width)
  print(paste0(chr,".",start,".",width))

  res.list <- list()
  for (dst in names(meta.list)){
    # dst <- "GSE50676"
    print(dst)
    smps <- meta.list[[dst]]
    for (smp in smps){
      # smp <- "GM12878"
      print(smp)
      if (grepl("AGO2_IP|FTC_small|TGIRT|DT2|i-pico|oscc",dst,perl=T)){
        pre="call_peak_dedup" 
      } else {
        pre="call_peak_all"
      }
      read_bw <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/",pre,"/tbigwig_RNA_EM/",smp,".transcriptome.bigWig") # tbigwig_RNA_EM, tbigwig_11RNA_primary
      # x <- rtracklayer::import.bw(read_bw) # 1-based ?
      # bw.res.df <- readBW(x)
      bw.res.df <- readBW(bw=read_bw,chr=chr,start=start,width=width)
      bw.res.df[1:3,]
      plotXbreak <- ifelse(dst==names(meta.list)[length(names(meta.list))] & smp==smps[length(smps)], TRUE, FALSE)
      # title <- ifelse(dst==names(meta.list)[1] & smp==smps[1], 
      #                 paste0(chr,":",start,"-",start+width),
      #                 "")
      title <- ""
      res.list[[paste0(dst,"_",smp)]] <- plotBW(single_base_bw=bw.res.df,chr=chr,start=start,width=width,ylab="",title=title,plotXbreak = plotXbreak,plotYbreak = T,plotYtext = F,annotate.size = 6, color = meta.col[[dst]], fill = meta.col[[dst]], maxY = 5, libsize = as.numeric(libsizes[smp])) # color = col, 
    }
  }

  ## collect
  tmp <- ggpubr::ggarrange(plotlist = res.list,ncol = 1,  align = "hv") #heights = c(7,1,1,1,1,1,1,1,1,1,1), nrow = length(res.list),
  ggsave(plot = tmp,filename = paste0("bw","-",chr,".",start,".",width,".pdf"),width = 15, height = 1.3*length(res.list)) # 15
}


#


# 1.2 plot multi sample bigwig (fig6,CRC,NC) ------------------
## fig6.CRC+NC,plasma+TCGA (fig6)
#meta.list <- list()
meta.col <- list()
genes.list0 <- list()
annotation_row_merge <- read.table(file = "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/tmp/PRNNA-TCGA3_annotation_row_merge.txt",header = T,sep = "\t",check.names = F,stringsAsFactors = F)

for(i in 1:nrow(annotation_row_merge)){
  genes.list0[[i]] <- data.frame(chr=annotation_row_merge[i,"chr"],start=annotation_row_merge[i,"start"],width=annotation_row_merge[i,"width"],peak=annotation_row_merge[i,"peak"]) #mir
}

dsts <- c("PRJNA540919 plasma","TCGA tissue") # GSE110381

# #GSE110381
# {
# genes.list <- genes.list0[c(8,69)]
# #ENST00000512398_____1.183.15,ENSG00000248790 (novelLnc)
# #genes.list <- genes.list[8]
# genes.list[[1]]$start <- 153
# genes.list[[1]]$width <- 75 # 5-fold
# 
# #enhancer__chr6___148273878____148276878_neg.1355.15
# #exp gn coord: chr6:148275478-148275553 # 5-fold
# #no signal in UCSC
# #genes.list <- genes.list[69]
# genes.list[[2]]$start <- 1325
# genes.list[[2]]$width <- 75 # 5-fold

# bed
#ENST00000512398_____1,153,228,peak_14332,1,+
#enhancer__chr6___148273878____148276878_neg,1325,1400,peak_72888,1,+
#

# genes.list <- genes.list[8]
# genes.list[[1]]$start <- 69
# genes.list[[1]]$width <- 45
# maxY.list <- list()
# maxY.list[["GSE110381_diff"]] <- 6
# maxY.list[["TCGA_small_diff"]] <- 1

# genes.list <- genes.list[c(3)]
# genes.list[[1]]$start <- 0
# genes.list[[1]]$width <- 83 # full
# #MIR762: Circulating microRNA-762 upregulates colorectal cancer may be accompanied by Wnt-1/β-catenin signaling
# #ENST00000390236_____1    47    18 peak_6344
}




#PRJNA
{
#miR eg
genes.list <- genes.list0[c(6)]
genes.list[[1]]$start <- 0
genes.list[[1]]$width <- 65 # full
#ENST00000615997_____1,MIR6803: exosomal miR-6803-5p was associated with poor prognosis in CRC independent of other confounding factors. 

# genes.list <- genes.list0[c(9)]
# genes.list[[1]]$start <- 1900
# genes.list[[1]]$width <- 100 # full
# #MIR194-2HG: 
# #3647

#lnc eg1
genes.list <- genes.list0[c(89)]
genes.list[[1]]$start <- 2100
genes.list[[1]]$width <- 200 # full
#G033992__chr17___18931287____18939853_pos,2201,17
#chr17:18933488-18933505

#lnc eg2
genes.list <- genes.list0[c(92)]
genes.list[[1]]$start <- 6300
genes.list[[1]]$width <- 200
#G090327__chrX___144234365____144242179_neg,6357,6373,16
#gn: chrX:144235806-144235822
}


#gn:  chr20:61804948-61805048   no much ucsc signal, overlapping with CDH4: inhibits ferroptosis in oral squamous cell carcinoma cells

dsts
libsizes <- list() 
libsizes[[dsts[[1]]]] <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/PRJNA540919_diff/call_peak_all/count_matrix/EM.libSize",header = F)$V1
libsizes[[dsts[[2]]]] <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/TCGA_small_diff3/call_peak_all/count_matrix/EM.libSize",header = F)$V1

#smpids <- NULL
smpids <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/PRJNA540919_diff/call_peak_all/count_matrix/cfpeakCNN_b5_d50_p1_GSE110381.txt",check.names = F,header = T)
smpids <- smpids[,2:ncol(smpids)]
smpids <- colnames(smpids)
names(libsizes[[dsts[[1]]]]) <- smpids
smpids <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/TCGA_small_diff3/call_peak_all/count_matrix/cfpeakCNN_b5_d50_p1_GSE110381.txt",check.names = F,header = T)
smpids <- smpids[,2:ncol(smpids)]
smpids <- colnames(smpids)
names(libsizes[[dsts[[2]]]]) <- smpids



for (i in 1:length(genes.list)){
  #i <- 1 "ENST00000385214_____1_30_49_+|pri_miRNA|ENST00000385214_____1|peak_6181|ENST00000385214_____1|30|49"
  chr <- genes.list[[i]]$chr
  start <- as.numeric(genes.list[[i]]$start)
  width <- as.numeric(genes.list[[i]]$width)
  peakID <- genes.list[[i]]$peak
  print(paste0(chr,".",start,".",width))
  
  res.list <- list()
  for (dst in dsts){ # names(meta.list)
    #dst <- "TCGA tissue" #"in-house plasma"
    print(dst)
    if (grepl("AGO2_IP|FTC_small|TGIRT",dst,perl=T)){
      pre="call_peak_dedup" 
    } else {
      pre="call_peak_all"
    }
    #smps <- meta.list[[dst]]
    smps <- as.character(annotation_row_merge[annotation_row_merge$peak==peakID,grepl(dst,colnames(annotation_row_merge))]) #[,paste0(dst,"_",group)]
    
    grps <- rep(c(rep("CRC",3),rep("NC",3)),times=1) # CRC,CRC,CRC,NC,NC,NC,  CRC,CRC,CRC,NC,NC,NC,
    
    for (j in 1:length(smps)){
      #j <- 6
      smp <- smps[j]
      grp <- grps[j]
      # print(smp)
      # print(grp)
      if (grepl("CRC",grp,perl=T)){
        bw.color <- "firebrick"
      } else if (grepl("NC",grp,perl=T)){
        bw.color <- "grey30" #"seagreen"
      }
      plotXbreak <- ifelse(dst==dsts[length(dsts)] & smp==smps[length(smps)], TRUE, FALSE)
      
      if ( dst == "in-house plasma"){
        dst.lab <- "WSQ_SMARTer_NEB_diff"
      } else if (dst == "GSE110381 plasma"){
        dst.lab <- "GSE110381_diff"
      } else if (dst == "PRJNA540919 plasma"){
        dst.lab <- "PRJNA540919_diff"
      } else if (dst == "TCGA tissue"){
        dst.lab <- "TCGA_small_diff3"
      }
      read_bw <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst.lab,"/",pre,"/tbigwig_RNA_EM/",smp,".transcriptome.bigWig")
      # x <- rtracklayer::import.bw(read_bw) # 1-based ?
      # bw.res.df <- readBW(x)
      bw.res.df <- readBW(bw=read_bw,chr=chr,start=start,width=width)
      # title <- ifelse(dst==names(meta.list)[1] & smp==smps[1], 
      #                 paste0(chr,":",start,"-",start+width),
      #                 "")
      title <- "" #paste0(chr,":",start,"-",start+width)#""
      res.list[[paste0(dst.lab,"_",smp)]] <- plotBW(single_base_bw=bw.res.df,chr=chr,start=start,width=width,ylab="",title=title,plotXbreak = plotXbreak, plotYbreak = T, plotYtext = F, color = bw.color, fill = bw.color,annotate.size = 8,libsize = libsizes[[dst]][[smp]], maxY = 1.5) # , maxY = maxY.list[[dst.lab]] 
    }
  }
  ## collect
  tmp <- ggpubr::ggarrange(plotlist = res.list, ncol = 1,  align = "hv") #heights = c(7,1,1,1,1,1,1,1,1,1,1), nrow = length(res.list),
  ggsave(plot = tmp,filename = paste0("bw","-",chr,".",start,".",width,".pdf"),width = 15, height = 1.3*length(res.list)) # 15
}



#plot other igv eg. (not dcb)
#ENST00000362125_____3: hsa-miR-340-5p	logfc<0
genes.list <- list()
genes.list[[1]] <- data.frame(chr="ENST00000362125_____3",start=0,width=95,peak="peak_4049") # miR
# ENST00000362125_____3   15      36      peak_4049       238     +
# ENST00000362125_____3   57      79      peak_4050       219.273 +
dsts <- "PRJNA540919_diff"
smp.tbl <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/",dsts[1],"/sample_table.txt"),header = T)
smp.tbl$group <- factor(smp.tbl$group,levels = c("CRC","NC"))
smp.tbl <- smp.tbl[order(smp.tbl$group),]

meta.list <- list()
meta.list[["PRJNA540919_diff"]] <- smps
meta.col <- list()
meta.col[["PRJNA540919_diff"]] <- "steelblue"

for (i in 1:length(genes.list)){
  #i <- 1 "ENST00000385214_____1_30_49_+|pri_miRNA|ENST00000385214_____1|peak_6181|ENST00000385214_____1|30|49"
  chr <- genes.list[[i]]$chr
  start <- as.numeric(genes.list[[i]]$start)
  width <- as.numeric(genes.list[[i]]$width)
  peakID <- genes.list[[i]]$peak
  print(paste0(chr,".",start,".",width))
  
  res.list <- list()
  for (dst in dsts){ # names(meta.list)
    #dst <- "TCGA tissue" #"in-house plasma"
    print(dst)
    if (grepl("AGO2_IP|FTC_small|TGIRT",dst,perl=T)){
      pre="call_peak_dedup" 
    } else {
      pre="call_peak_all"
    }
    #smps <- meta.list[[dst]]
    # smps <- as.character(annotation_row_merge[annotation_row_merge$peak==peakID,grepl(dst,colnames(annotation_row_merge))]) #[,paste0(dst,"_",group)]
    smps <- smp.tbl$sample
    grps <- smp.tbl$group # rep(c(rep("CRC",3),rep("NC",3)),times=1) # CRC,CRC,CRC,NC,NC,NC,  CRC,CRC,CRC,NC,NC,NC,
    
    for (j in 1:length(smps)){
      #j <- 6
      smp <- smps[j]
      grp <- grps[j]
      # print(smp)
      # print(grp)
      if (grepl("CRC",grp,perl=T)){
        bw.color <- "firebrick"
      } else if (grepl("NC",grp,perl=T)){
        bw.color <- "grey30" #"seagreen"
      }
      plotXbreak <- ifelse(dst==dsts[length(dsts)] & smp==smps[length(smps)], TRUE, FALSE)
      
      if ( dst == "in-house plasma"){
        dst.lab <- "WSQ_SMARTer_NEB_diff"
      } else if (dst == "GSE110381 plasma"){
        dst.lab <- "GSE110381_diff"
      } else if (dst == "PRJNA540919 plasma"){
        dst.lab <- "PRJNA540919_diff"
      } else if (dst == "TCGA tissue"){
        dst.lab <- "TCGA_small_diff3"
      }
      read_bw <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst.lab,"/",pre,"/tbigwig_RNA_EM/",smp,".transcriptome.bigWig")
      # x <- rtracklayer::import.bw(read_bw) # 1-based ?
      # bw.res.df <- readBW(x)
      bw.res.df <- readBW(bw=read_bw,chr=chr,start=start,width=width)
      # title <- ifelse(dst==names(meta.list)[1] & smp==smps[1], 
      #                 paste0(chr,":",start,"-",start+width),
      #                 "")
      title <- "" #paste0(chr,":",start,"-",start+width)#""
      res.list[[paste0(dst.lab,"_",smp)]] <- plotBW(single_base_bw=bw.res.df,chr=chr,start=start,width=width,ylab="",title=title,plotXbreak = plotXbreak, plotYbreak = T, plotYtext = F, color = bw.color, fill = bw.color, maxY=5200) # , maxY = maxY.list[[dst.lab]] 
    }
  }
  ## collect
  tmp <- ggpubr::ggarrange(plotlist = res.list, ncol = 1,  align = "hv") #heights = c(7,1,1,1,1,1,1,1,1,1,1), nrow = length(res.list),
  ggsave(plot = tmp,filename = paste0("bw","-",chr,".",start,".",width,".pdf"),width = 15, height = 1.3*length(res.list)) # 15
}

#fig7


# 2.plot RBP,G4,EV ref ------------------------------------------------------------
setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev")
# bed_ago2 <- "/BioII/lulab_b/baopengfei/shared_reference/RBP/splitByRBP/AGO2_tx.bed"
# bed_otherRBPs <- "/BioII/lulab_b/baopengfei/shared_reference/RBP/splitByRBP/otherRBPs_tx.bed"
bed_RBPs <- "/BioII/lulab_b/baopengfei/shared_reference/RBP/splitByRBP/RBPs_tx.bed"
# bed_hotspot <- "/BioII/lulab_b/baopengfei/shared_reference/RBP/POSTAR3_hotspot/all_merge_newTxID.bed"
bed_EV <- "/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/lulab/11RNAsmallDomain_EVenrich.bed"
# bed_G4 <- "/BioII/lulab_b/baopengfei/shared_reference/structure/G4iMGrinder/bedg4im_merge_tx.bed6"
bed_G4 <- "/BioII/lulab_b/baopengfei/shared_reference/structure/quadratlas_g4grinder_tx.bed"
#bed_gold <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/tmp/4dst_filter2_11RNA.bed"
bed_gold <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/tmp/6dst_filter3_11RNA.bed"

# mirdeep_mature <- Sys.glob("/BioII/lulab_b/baopengfei/gitsoft/mirdeep2/test/output/mirna_results*/novel_mature_*_*_2023_t_*_*_*_score-50_to_na.bed6")
# mirdeep_pre <- Sys.glob("/BioII/lulab_b/baopengfei/gitsoft/mirdeep2/test/output/mirna_results*/novel_pres_*_*_2023_t_*_*_*_score-50_to_na.bed6")


# #Fig6: CRCvsNC
# genes.list <- list()
# genes.list[[1]] <- data.frame(chr="ENST00000365690_____4",start=0,width=78) #mir
# genes.list[[2]] <- data.frame(chr="G078490__chr7___76260333____76318258_neg",start=31850,width=100) #intron

    


#genes.list  # same as bw
for (i in 1:length(genes.list)){
  # i <- 1
  chr <- genes.list[[i]]$chr
  start <- genes.list[[i]]$start
  width <- as.numeric(genes.list[[i]]$width)
  print(paste0(chr,".",start,".",width))
  
  # res.list <- list()
  
  # # RBP, G4 annotation
  # ago2.res.df <- readAnno(bed_ago2)
  # bed.ago2 <- plotAnno(ago2.res.df,"AGO2")
  # # bed.ago2
  # 
  # otherRBPs.res.df <- readAnno(bed_otherRBPs)
  # bed.otherRBPs <-  plotAnno(otherRBPs.res.df,"otherRBPs")
  # # bed.otherRBPs

  RBPs.res.df <- readAnno(bed = bed_RBPs, chr = chr, start = start, width = width)
  bed.RBPs <-  plotAnno(bed = RBPs.res.df, chr = chr, start = start, width = width, ylab = "") # RBPs

  # hotspot.res.df <- readAnno(bed_hotspot)
  # bed.hotspot <- plotAnno(hotspot.res.df,"hotspot")
  # # bed.hotspot

  EV.res.df <- readAnno(bed = bed_EV, chr = chr, start = start, width = width)
  bed.EV <- plotAnno(bed = EV.res.df, chr = chr, start = start, width = width, ylab = "") # EV
  # bed.EV


  G4.res.df <- readAnno(bed = bed_G4, chr = chr, start = start, width = width)
  bed.G4 <- plotAnno(bed = G4.res.df, chr = chr, start = start, width = width, ylab = "") # G4iM
  # bed.G4

  gold.res.df <- readAnno(bed = bed_gold, chr = chr, start = start, width = width)
  bed.gold <- plotAnno(bed = gold.res.df, chr = chr, start = start, width = width, ylab = "") #

  ## get list
  res.list3 <- list()
  # res.list3[["bed.ago2"]] <- bed.ago2
  # res.list3[["bed.otherRBPs"]] <- bed.otherRBPs
  res.list3[["bed.RBPs"]] <- bed.RBPs
  # res.list3[["bed.hotspot"]] <- bed.hotspot
  res.list3[["bed.EV"]] <- bed.EV
  res.list3[["bed.G4"]] <- bed.G4
  res.list3[["bed.gold"]] <- bed.gold
  
  ## miRDeep2
  # bed.mirdeep.mature.res.df <- readAnno(bed = mirdeep_mature, chr = chr, start = start, width = width)
  # bed.mirdeep.mature <- plotAnno(bed = bed.mirdeep.mature.res.df, chr = chr, start = start, width = width,color = "salmon",fill = "salmon", ylab = "") # 
  # res.list3[["bed_mirdeep_mature"]] <- bed.mirdeep.mature
  # bed.mirdeep.pre.res.df <- readAnno(bed = mirdeep_pre, chr = chr, start = start, width = width)
  # bed.mirdeep.pre <- plotAnno(bed = bed.mirdeep.pre.res.df, chr = chr, start = start, width = width, ylab = "") 
  # res.list3[["bed_mirdeep_pre"]] <- bed.mirdeep.pre
  # 
  
  ## collect
  tmp3 <- ggpubr::ggarrange(plotlist = res.list3, ncol = 1,heights = rep(1,length(res.list3)), nrow = length(res.list3),  align = "hv") # 
  ggsave(plot = tmp3, filename = paste0("ref_bed-",chr,".",start,".",width,".pdf"),width = 10, height = 0.5*length(res.list3))
}

#






# 3.plot multi-tx exons ------------------------------
setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev")

inputFile <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/gtf_by_biotype/11RNA.gtf"
gtf0 <- rtracklayer::import(inputFile)
#lnc seem has little multi-exon tx
#45S rRNA seem has no gtf record


#genes.list <- genes.list2
#genes.list  # same as bw
for (i in 1:length(genes.list)){
  # i <- 1
  chr <- genes.list[[i]]$chr
  start <- genes.list[[i]]$start
  width <- genes.list[[i]]$width
  #chr <- "NR_146151_____1"
  #start <- 0
  #width <- 3000
  print(paste0(chr,".",start,".",width))
  print(chr)
  # txID <- "ENST00000377315_____4"
  # txID <- "NR_023363_____1"
  # txID <- "26266"
  # txID <- "ENST00000564646_____1"
  tmp3 <- plotMultiTxTrack(txID_new = chr, gtf = gtf0, fill="grey50") # 
  ggsave(plot = tmp3, filename = paste0("exons_bed-",chr,".",start,".",width,".pdf"),width = 10, height = 0.8*length(genes.list))
}







# 4.plot single sample bw + 4 peak + gold peak (fig3) ------------------------------------------------------------
setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev")
conflict_prefer("which", "base")

#bed_gold <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/tmp/4dst_filter2_11RNA.bed"
bed_gold <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/tmp/6dst_filter3_11RNA.bed"

meta.list2 <- list()
# meta.list[["WSQ_SMARTer_NEB"]] <- c("NC_PKU-2392860_smart_1","NC_PKU-2397512_smart_1","NC_PKU-2392860_smart_PNK_1","NC_PKU-2397512_smart_PNK_1","NC_PKU-2392860_1","NC_PKU-2397512_1")
meta.list2[["GSE71008"]] <- "SAMN03863396"
# meta.list2[["GSE71008_NCpool"]] <- "NCpool"  # SAMN03863523,SAMN03863523 c("SRR2105340","SRR2105341") # ,"SRR2105335"
# meta.list[["GSE94582"]] <- c("NEBNext_Lab1","TruSeq_Lab1","N4_A_Lab5","N4_B_Lab2","CleanTag_Lab5")
# meta.list2[["GSE50676"]] <- c("GM12878") # CLIP
# meta.list2[["GSE52600"]] <- c("SRR1035214", # RIP-seq
#                              "SRR1035213")
# meta.list2[["sim"]] <- c("tissue50_blood50","tissue005_blood995") # sim plasma, tissue05_blood95, tissue50_blood50, tissue005_blood995
# pre="call_domain_withRepeats_all" # call_domain_withRepeats_all
# meta.list2 <- meta.list2[[1]]

genes.list2 <- list()
# #fig3
genes.list2[[1]] <- data.frame(chr="ENST00000385227_____1",start=0,width=89) #miR
# genes.list2[[2]] <- data.frame(chr="ENST00000564646_____1",start=0,width=2873) #
genes.list2[[3]] <- data.frame(chr="ENST00000538098_____2",start=0,width=548) # mRNA
genes.list2[[4]] <- data.frame(chr="ENST00000554988_____1",start=0,width=450) # lncRNA,638
genes.list2[[2]] <- data.frame(chr="NR_023363_____1",start=0,width=121)
# genes.list2[[1]] <- data.frame(chr="ENST00000315707_____3",start=0,width=2000)
# genes.list2[[5]] <- data.frame(chr="ENST00000536684_____2",start=0,width=1303) # mRNA


# genes.list2[[1]] <- data.frame(chr="ENST00000602385_____1",start=0,width=541) # sim: terc gene
# genes.list2[[2]] <- data.frame(chr="ENST00000385247_____1",start=0,width=97) # sim: miR gene
# # genes.list2[[3]] <- data.frame(chr="ENST00000364849_____1",start=0,width=89) # sim: snoRNA gene
# genes.list2[[3]] <- data.frame(chr="ENST00000521127_____1",start=0,width=560)	# sim: snoRNA gene map to lnc
# 


# #peak specific
##expeak
# genes.list2[[1]] <- data.frame(chr="NR_146154_____1",start=0,width=5055)
genes.list2[[1]] <- data.frame(chr="NR_146152_____1",start=0,width=600) # 250
# genes.list2[[3]] <- data.frame(chr="NR_146148_____1",start=0,width=5054)
# genes.list2[[4]] <- data.frame(chr="ENST00000636484_____1",start=0,width=328)
##clipper
# genes.list2[[1]] <- data.frame(chr="ENST00000571127_____1",start=0,width=791)
# genes.list2[[2]] <- data.frame(chr="ENST00000545075_____2",start=0,width=1530)
# genes.list2[[3]] <- data.frame(chr="ENST00000319763_____1",start=0,width=2081)
# genes.list2[[4]] <- data.frame(chr="ENST00000391625_____2",start=0,width=3894)
# genes.list2[[5]] <- data.frame(chr="ENST00000414890_____1",start=0,width=737)
##clam
# genes.list2[[1]] <- data.frame(chr="ENST00000521141_____1",start=0,width=1630)
# genes.list2[[2]] <- data.frame(chr="T004607",start=0,width=1634)
# genes.list2[[3]] <- data.frame(chr="T234700",start=0,width=10167)
# genes.list2[[4]] <- data.frame(chr="30032",start=0,width=72)
# genes.list2[[5]] <- data.frame(chr="ENST00000310125_____4",start=0,width=1680)



#meta.list2 <- meta.list
#genes.list2 <- genes.list
for (i in 1:length(genes.list2)){
  #i <- 1
  chr <- genes.list2[[i]]$chr
  start <- genes.list2[[i]]$start
  width <- genes.list2[[i]]$width
  print(paste0(chr,".",start,".",width))
  
  for (dst in names(meta.list2)){
    # dst <- "GSE71008_NCpool"
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
      #smp <- "NCpool"
      res.list2 <- list()
      print(smp)

      read_bw <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/",pre,"/tbigwig_RNA_EM/",smp,".transcriptome.bigWig")
      # x <- rtracklayer::import.bw(read_bw) # 1-based ?
      #x <- as.data.frame(x)
      bw.res.df <- readBW(read_bw,chr = chr,start = start,width = width) 
      res.list2[[paste0(dst,"_",smp,"_EM")]] <- plotBW(bw.res.df,chr = chr,start = start,width = width,ylab = "",xlab = "",title = "") # paste0(dst,"_",smp,"\n",chr,"_",start,"_",width)  
      # res.list2[[paste0(dst,"_",smp)]] <- gg.gap(plot=res.list2[[1]],rel_heights = c(0.2,1), segments=c(20,120),ylim = c(0,140)) # gg.gap truncate y axis for specific tx, # + ggbreak::scale_y_break(breaks = c(50, 300)), ggbreak seem not fit ggarrange ?
      # res.list2[[paste0(dst,"_",smp)]] <- gg.gap(plot=res.list2[[1]], rel_heights = c(0.2,1), segments=c(40,220),ylim = c(0,260)) # gg.gap truncate y axis for specific tx, # + ggbreak::scale_y_break(breaks = c(50, 300)), ggbreak seem not fit ggarrange ?
      # ggsave(plot = res.list2[[paste0(dst,"_",smp)]], filename = paste0("peaks_bw_bed_",dst,"_",smp,"-",chr,".",start,".",width,".2.pdf"),width = 13, height = 3)
      
      # read_bw <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/",pre,"/tbigwig_11RNA_primary/",smp,".transcriptome.bigWig")
      # # x <- rtracklayer::import.bw(read_bw) # 1-based ?
      # #x <- as.data.frame(x)
      # bw.res.df <- readBW(read_bw,chr = chr,start = start,width = width) 
      # res.list2[[paste0(dst,"_",smp,"_pri")]] <- plotBW(bw.res.df,chr = chr,start = start,width = width,ylab = "",xlab = "",title = "") # paste0(dst,"_",smp,"\n",chr,"_",start,"_",width)  
    
        
      # peak
      bed_piranha <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/",pre,"/piranha_by_sample/b5_p01/",smp,".bed")
      bed_clipper <-  paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/",pre,"/clipper_by_sample/b5_p05/",smp,".bed")
      bed_clam <-  paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/",pre,"/clam_by_sample/b5_p005/",smp,".bed")
      # bed_localmax <-  paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/",pre,"/domains_localmax_significant/b5_d05_p05/",smp,".bed")
      bed_localmaxEM2 <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/",pre,"/expeakCNN_by_sample/b5_d50_p1/",smp,".bed")
      
      ## piranha peak bed
      piranha <- data.table::fread(bed_piranha,sep = "\t",header = F,stringsAsFactors = F,check.names = F) # 1-based ?
      if(nrow(piranha)==0){
        piranha <- data.frame("V1"="","V2"="","V3"="","V4"="","V5"="","V6"="","V7"="")
      }
      colnames(piranha) <- c("seqnames","start","end","name","score","strand","pval")
      piranha <- piranha[,c("seqnames","start","end","name","score","strand")]
      piranha.res.df <- readPeak(piranha,chr = chr,start = start,width = width)
      bed.piranha <- plotPeak(bed = piranha.res.df,chr = chr,start = start,width = width,ylab = "")
      # bed.piranha
      res.list2[[paste0(dst,"_",smp,"Piranha")]] <- bed.piranha
      
      ## CLIPper peak bed
      CLIPper <- data.table::fread(bed_clipper,sep = "\t",header = F,stringsAsFactors = F,check.names = F) # 1-based ?
      if(nrow(CLIPper)==0){
        CLIPper <- data.frame("V1"="","V2"="","V3"="","V4"="","V5"="","V6"="","V7"="","V8"="")
      }
      colnames(CLIPper) <- c("seqnames","start","end","name","score","strand","V7","V8") # score: pval
      CLIPper <- CLIPper[,c("seqnames","start","end","name","score","strand")]
      CLIPper.res.df <- readPeak(CLIPper,chr = chr,start = start,width = width)
      bed.clipper <- plotPeak(bed = CLIPper.res.df,chr = chr,start = start,width = width,ylab = "")
      # bed.clipper
      res.list2[[paste0(dst,"_",smp,"CLIPper")]] <- bed.clipper
      
      ## CLAM peak bed
      CLAM <- data.table::fread(bed_clam,sep = "\t",header = F,stringsAsFactors = F,check.names = F) # 1-based ?
      if(nrow(CLAM)==0){
        CLAM <- data.frame("V1"="","V2"="","V3"="","V4"="","V5"="","V6"="")
      }      
      colnames(CLAM) <- c("seqnames","start","end","name","type","strand")
      CLAM$score <- 1
      CLAM <- CLAM[,c("seqnames","start","end","name","score","strand")]
      clam.res.df <- readPeak(CLAM,chr = chr,start = start,width = width)
      bed.clam <- plotPeak(bed = clam.res.df,chr = chr,start = start,width = width,ylab = "")
      # bed.clam
      res.list2[[paste0(dst,"_",smp,"CLAM")]] <- bed.clam 

      ## expeak peak bed
      localmaxEM2 <- data.table::fread(bed_localmaxEM2,sep = "\t",header = F,stringsAsFactors = F,check.names = F) # 1-based ?
      if(nrow(localmaxEM2)==0){
        localmaxEM2 <- data.frame("V1"="","V2"="","V3"="","V4"="","V5"="","V6"="")
      } 
      localmaxEM2 <- localmaxEM2[,1:6]
      colnames(localmaxEM2) <- c("seqnames","start","end","name","type","strand")
      localmaxEM2$score <- 1
      localmaxEM2 <- localmaxEM2[,c("seqnames","start","end","name","score","strand")]
      localmaxEM2.res.df <- readPeak(localmaxEM2,chr = chr,start = start,width = width)
      bed.localmaxEM2 <- plotPeak(bed = localmaxEM2.res.df,chr = chr,start = start,width = width,ylab = "") # , color = "firebrick", fill = "firebrick"
      # bed.localmaxEM2
      res.list2[[paste0(dst,"_",smp,"exPeak")]] <- bed.localmaxEM2
      
      
      # gold peak bed
      gold <- data.table::fread(bed_gold,sep = "\t",header = F,stringsAsFactors = F,check.names = F) # 1-based ?
      if(nrow(gold)==0){
        gold <- data.frame("V1"="","V2"="","V3"="","V4"="","V5"="","V6"="")
      } 
      gold <- gold[,1:6]
      colnames(gold) <- c("seqnames","start","end","name","type","strand")
      gold$score <- 1
      gold <- gold[,c("seqnames","start","end","name","score","strand")]
      gold.res.df <- readPeak(gold,chr = chr,start = start,width = width)
      bed.gold <- plotPeak(bed = gold.res.df,chr = chr,color = "salmon",fill = "salmon",start = start,width = width, ylab = "")
      # bed.gold
      res.list2[[paste0(dst,"_",smp,"gold")]] <- bed.gold
      
      
      # # method specific peak bed
      # bed_piranha <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/",pre,"/piranha_by_sample/b5_p01/intersect/",smp,"_piranhaOnly.bed6")
      # bed_clipper <-  paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/",pre,"/clipper_by_sample/b5_p05/intersect/",smp,"_clipperOnly.bed6")
      # bed_clam <-  paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/",pre,"/clam_by_sample/b5_p005/intersect/",smp,"_clamOnly.bed6")
      # bed_localmaxEM2 <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/",pre,"/expeak_by_sample/b5_d50_p1/intersect/",smp,"_expeakOnly.bed6")
      # 
      # piranha <- data.table::fread(bed_piranha,sep = "\t",header = F,stringsAsFactors = F,check.names = F) # 1-based ?
      # if(nrow(piranha)==0){
      #   piranha <- data.frame("V1"="","V2"="","V3"="","V4"="","V5"="","V6"="","V7"="")
      # }
      # colnames(piranha) <- c("seqnames","start","end","name","score","strand")
      # # piranha <- piranha[,c("seqnames","start","end","name","score","strand")]
      # piranha.res.df <- readPeak(piranha,chr = chr,start = start,width = width)
      # bed.piranha <- plotPeak(bed = piranha.res.df,chr = chr,start = start,width = width,ylab = "")
      # # bed.piranha
      # res.list2[[paste0(dst,"_",smp,"PiranhaOnly")]] <- bed.piranha
      # 
      # CLIPper <- data.table::fread(bed_clipper,sep = "\t",header = F,stringsAsFactors = F,check.names = F) # 1-based ?
      # if(nrow(CLIPper)==0){
      #   CLIPper <- data.frame("V1"="","V2"="","V3"="","V4"="","V5"="","V6"="","V7"="","V8"="")
      # }
      # colnames(CLIPper) <- c("seqnames","start","end","name","score","strand") # score: pval
      # #CLIPper <- CLIPper[,c("seqnames","start","end","name","score","strand")]
      # CLIPper.res.df <- readPeak(CLIPper,chr = chr,start = start,width = width)
      # bed.clipper <- plotPeak(bed = CLIPper.res.df,chr = chr,start = start,width = width,ylab = "")
      # # bed.clipper
      # res.list2[[paste0(dst,"_",smp,"CLIPperOnly")]] <- bed.clipper
      # 
      # CLAM <- data.table::fread(bed_clam,sep = "\t",header = F,stringsAsFactors = F,check.names = F) # 1-based ?
      # if(nrow(CLAM)==0){
      #   CLAM <- data.frame("V1"="","V2"="","V3"="","V4"="","V5"="","V6"="")
      # }      
      # colnames(CLAM) <- c("seqnames","start","end","name","type","strand")
      # CLAM$score <- 1
      # CLAM <- CLAM[,c("seqnames","start","end","name","score","strand")]
      # clam.res.df <- readPeak(CLAM,chr = chr,start = start,width = width)
      # bed.clam <- plotPeak(bed = clam.res.df,chr = chr,start = start,width = width,ylab = "")
      # # bed.clam
      # res.list2[[paste0(dst,"_",smp,"CLAMOnly")]] <- bed.clam 
      # 
      # localmaxEM2 <- data.table::fread(bed_localmaxEM2,sep = "\t",header = F,stringsAsFactors = F,check.names = F) # 1-based ?
      # if(nrow(localmaxEM2)==0){
      #   localmaxEM2 <- data.frame("V1"="","V2"="","V3"="","V4"="","V5"="","V6"="")
      # } 
      # localmaxEM2 <- localmaxEM2[,1:6]
      # colnames(localmaxEM2) <- c("seqnames","start","end","name","type","strand")
      # localmaxEM2$score <- 1
      # localmaxEM2 <- localmaxEM2[,c("seqnames","start","end","name","score","strand")]
      # localmaxEM2.res.df <- readPeak(localmaxEM2,chr = chr,start = start,width = width)
      # bed.localmaxEM2 <- plotPeak(bed = localmaxEM2.res.df,chr = chr,start = start,width = width,ylab = "") # , color = "firebrick", fill = "firebrick"
      # # bed.localmaxEM2
      # res.list2[[paste0(dst,"_",smp,"exPeakOnly")]] <- bed.localmaxEM2
      
      
      ## collect
      # tmp2 <- ggpubr::ggarrange(plotlist = res.list2, ncol = 1,heights = c(5,5, 1,1,1,1, 1 ),align = "v")  # 1,1,1,1
      ggsave(plot = tmp2, filename = paste0("peaks_bw_bed_",dst,"_",smp,"-",chr,".",start,".",width,".pdf"),width = 12, height = 10)
      tmp2 <- ggpubr::ggarrange(plotlist = res.list2, ncol = 1,heights = c(5, 1,1,1,1, 1 ),align = "v")  # 1,1,1,1
      ggsave(plot = tmp2, filename = paste0("peaks_bw_bed_",dst,"_",smp,"-",chr,".",start,".",width,".pdf"),width = 12, height = 4) 
    }
  }
}

#



# 4.2 plot single sample bw + 4 peak + POSTAR3 peak (supp fig3: CLIP-seq) ------------------------------------------------------------
setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev")
#bed_gold <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/tmp/4dst_filter2_11RNA.bed"
bed_gold <- "/lulabdata/baopengfei/shared_reference/RBP/splitByRBP/AGO2_tx_merge.bed"
bed_ribosnitch <- "/lulabdata/baopengfei/shared_reference/2014_nature_riboSNitches/RiboSNitches_hg38_tx.bed"

meta.list2 <- list()
meta.list2[["GSE50676"]] <- c("GM12878","GM12891","GM12892") # CLIP
# meta.list2[["GSE50676_NCpool"]] <- c("NCpool") # CLIP
# meta.list2[["GSE52600"]] <- c("SRR1035214", # RIP-seq
#                              "SRR1035213")


genes.list2 <- list()

# # #fig3 supp
genes.list2[[1]] <- data.frame(chr="ENST00000522153_____1",start=0,width=556)
# genes.list2[[2]] <- data.frame(chr="ENST00000369278_____4",start=700,width=700)

#genes.list2[[1]] <- data.frame(chr="ENST00000416016_____2",start=900,width=400)
#genes.list2[[1]] <- data.frame(chr="ENST00000341446_____8",start=800,width=200)
genes.list2[[1]] <- data.frame(chr="ENST00000456707_____1",start=0,width=609) # PARS_26|chr1|201876303|201876304|+
#genes.list2[[1]] <- data.frame(chr="ENST00000435041_____2",start=800,width=200)
#genes.list2[[1]] <- data.frame(chr="ENST00000485763_____1",start=400,width=200)

# ENST00000485763_____1,511,558 (616)
# ENST00000416016_____2,1097,1178 (1336)
# ENST00000341446_____8,931,979   (6648)
# ENST00000456707_____1,523,591   (609 )
# ENST00000435041_____2,912,999   (4290)

#meta.list2 <- meta.list
#genes.list2 <- genes.list
for (i in 1:length(genes.list2)){
  #i <- 1
  chr <- genes.list2[[i]]$chr
  start <- genes.list2[[i]]$start
  width <- genes.list2[[i]]$width
  print(paste0(chr,".",start,".",width))
  
  for (dst in names(meta.list2)){
    # dst <- "GSE71008_NCpool"
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
      #smp <- "NCpool"
      res.list2 <- list()
      print(smp)
      
      read_bw <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/",pre,"/tbigwig_RNA_EM/",smp,".transcriptome.bigWig")
      # x <- rtracklayer::import.bw(read_bw) # 1-based ?
      #x <- as.data.frame(x)
      bw.res.df <- readBW(read_bw,chr = chr,start = start,width = width) 
      res.list2[[paste0(dst,"_",smp,"_EM")]] <- plotBW(bw.res.df,chr = chr,start = start,width = width,ylab = "",xlab = "",title = "") # paste0(dst,"_",smp,"\n",chr,"_",start,"_",width)  
      # res.list2[[paste0(dst,"_",smp)]] <- gg.gap(plot=res.list2[[1]],rel_heights = c(0.2,1), segments=c(20,120),ylim = c(0,140)) # gg.gap truncate y axis for specific tx, # + ggbreak::scale_y_break(breaks = c(50, 300)), ggbreak seem not fit ggarrange ?
      # res.list2[[paste0(dst,"_",smp)]] <- gg.gap(plot=res.list2[[1]], rel_heights = c(0.2,1), segments=c(40,220),ylim = c(0,260)) # gg.gap truncate y axis for specific tx, # + ggbreak::scale_y_break(breaks = c(50, 300)), ggbreak seem not fit ggarrange ?
      # ggsave(plot = res.list2[[paste0(dst,"_",smp)]], filename = paste0("peaks_bw_bed_",dst,"_",smp,"-",chr,".",start,".",width,".2.pdf"),width = 13, height = 3)
      
      # read_bw <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/",pre,"/tbigwig_11RNA_primary/",smp,".transcriptome.bigWig")
      # # x <- rtracklayer::import.bw(read_bw) # 1-based ?
      # #x <- as.data.frame(x)
      # bw.res.df <- readBW(read_bw,chr = chr,start = start,width = width) 
      # res.list2[[paste0(dst,"_",smp,"_pri")]] <- plotBW(bw.res.df,chr = chr,start = start,width = width,ylab = "",xlab = "",title = "") # paste0(dst,"_",smp,"\n",chr,"_",start,"_",width)  
      # 
      
      # peak
      bed_piranha <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/",pre,"/piranha_by_sample/b5_p01/",smp,".bed")
      bed_clipper <-  paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/",pre,"/clipper_by_sample/b5_p05/",smp,".bed")
      bed_clam <-  paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/",pre,"/clam_by_sample/b5_p005/",smp,".bed")
      # bed_localmax <-  paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/",pre,"/domains_localmax_significant/b5_d05_p05/",smp,".bed")
      bed_localmaxEM2 <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/",pre,"/expeakCNN_by_sample/b5_d50_p1/",smp,".bed")
      
      ## piranha peak bed
      piranha <- data.table::fread(bed_piranha,sep = "\t",header = F,stringsAsFactors = F,check.names = F) # 1-based ?
      if(nrow(piranha)==0){
        piranha <- data.frame("V1"="","V2"="","V3"="","V4"="","V5"="","V6"="","V7"="")
      }
      colnames(piranha) <- c("seqnames","start","end","name","score","strand","pval")
      piranha <- piranha[,c("seqnames","start","end","name","score","strand")]
      piranha.res.df <- readPeak(piranha,chr = chr,start = start,width = width)
      bed.piranha <- plotPeak(bed = piranha.res.df,chr = chr,start = start,width = width,ylab = "")
      # bed.piranha
      res.list2[[paste0(dst,"_",smp,"Piranha")]] <- bed.piranha
      
      ## CLIPper peak bed
      CLIPper <- data.table::fread(bed_clipper,sep = "\t",header = F,stringsAsFactors = F,check.names = F) # 1-based ?
      if(nrow(CLIPper)==0){
        CLIPper <- data.frame("V1"="","V2"="","V3"="","V4"="","V5"="","V6"="","V7"="","V8"="")
      }
      colnames(CLIPper) <- c("seqnames","start","end","name","score","strand","V7","V8") # score: pval
      CLIPper <- CLIPper[,c("seqnames","start","end","name","score","strand")]
      CLIPper.res.df <- readPeak(CLIPper,chr = chr,start = start,width = width)
      bed.clipper <- plotPeak(bed = CLIPper.res.df,chr = chr,start = start,width = width,ylab = "")
      # bed.clipper
      res.list2[[paste0(dst,"_",smp,"CLIPper")]] <- bed.clipper
      
      ## CLAM peak bed
      CLAM <- data.table::fread(bed_clam,sep = "\t",header = F,stringsAsFactors = F,check.names = F) # 1-based ?
      if(nrow(CLAM)==0){
        CLAM <- data.frame("V1"="","V2"="","V3"="","V4"="","V5"="","V6"="")
      }      
      colnames(CLAM) <- c("seqnames","start","end","name","type","strand")
      CLAM$score <- 1
      CLAM <- CLAM[,c("seqnames","start","end","name","score","strand")]
      clam.res.df <- readPeak(CLAM,chr = chr,start = start,width = width)
      bed.clam <- plotPeak(bed = clam.res.df,chr = chr,start = start,width = width,ylab = "")
      # bed.clam
      res.list2[[paste0(dst,"_",smp,"CLAM")]] <- bed.clam 
      
      ## expeak peak bed
      localmaxEM2 <- data.table::fread(bed_localmaxEM2,sep = "\t",header = F,stringsAsFactors = F,check.names = F) # 1-based ?
      if(nrow(localmaxEM2)==0){
        localmaxEM2 <- data.frame("V1"="","V2"="","V3"="","V4"="","V5"="","V6"="")
      } 
      localmaxEM2 <- localmaxEM2[,1:6]
      colnames(localmaxEM2) <- c("seqnames","start","end","name","type","strand")
      localmaxEM2$score <- 1
      localmaxEM2 <- localmaxEM2[,c("seqnames","start","end","name","score","strand")]
      localmaxEM2.res.df <- readPeak(localmaxEM2,chr = chr,start = start,width = width)
      bed.localmaxEM2 <- plotPeak(bed = localmaxEM2.res.df,chr = chr,start = start,width = width,ylab = "") # , color = "firebrick", fill = "firebrick"
      # bed.localmaxEM2
      res.list2[[paste0(dst,"_",smp,"exPeak")]] <- bed.localmaxEM2
      
      
      # gold peak bed
      gold <- data.table::fread(bed_gold,sep = "\t",header = F,stringsAsFactors = F,check.names = F) # 1-based ?
      if(nrow(gold)==0){
        gold <- data.frame("V1"="","V2"="","V3"="","V4"="","V5"="","V6"="")
      } 
      gold <- gold[,1:6]
      colnames(gold) <- c("seqnames","start","end","name","type","strand")
      gold$score <- 1
      gold <- gold[,c("seqnames","start","end","name","score","strand")]
      gold.res.df <- readPeak(gold,chr = chr,start = start,width = width)
      bed.gold <- plotPeak(bed = gold.res.df,chr = chr,color = "salmon",fill = "salmon",start = start,width = width, ylab = "")
      # bed.gold
      res.list2[[paste0(dst,"_",smp,"gold")]] <- bed.gold
      
      
      # RiboSNitch peak bed
      ribosnitch <- data.table::fread(bed_ribosnitch,sep = "\t",header = F,stringsAsFactors = F,check.names = F) # 1-based ?
      if(nrow(ribosnitch)==0){
        ribosnitch <- data.frame("V1"="","V2"="","V3"="","V4"="","V5"="","V6"="")
      } 
      ribosnitch <- ribosnitch[,1:6]
      colnames(ribosnitch) <- c("seqnames","start","end","name","type","strand")
      ribosnitch$score <- 1
      ribosnitch <- ribosnitch[,c("seqnames","start","end","name","score","strand")]
      ribosnitch.res.df <- readPeak(ribosnitch,chr = chr,start = start,width = width)
      bed.ribosnitch <- plotPeak(bed = ribosnitch.res.df,chr = chr,color = "salmon",fill = "salmon",start = start,width = width, ylab = "")
      # bed.ribosnitch
      res.list2[[paste0(dst,"_",smp,"ribosnitch")]] <- bed.ribosnitch
      
      
      
      # # method specific peak bed
      # bed_piranha <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/",pre,"/piranha_by_sample/b5_p01/intersect/",smp,"_piranhaOnly.bed6")
      # bed_clipper <-  paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/",pre,"/clipper_by_sample/b5_p05/intersect/",smp,"_clipperOnly.bed6")
      # bed_clam <-  paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/",pre,"/clam_by_sample/b5_p005/intersect/",smp,"_clamOnly.bed6")
      # bed_localmaxEM2 <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/",pre,"/expeakCNN_by_sample/b5_d50_p1/intersect/",smp,"_expeakCNNOnly.bed6")
      # 
      # piranha <- data.table::fread(bed_piranha,sep = "\t",header = F,stringsAsFactors = F,check.names = F) # 1-based ?
      # if(nrow(piranha)==0){
      #   piranha <- data.frame("V1"="","V2"="","V3"="","V4"="","V5"="","V6"="","V7"="")
      # }
      # colnames(piranha) <- c("seqnames","start","end","name","score","strand")
      # # piranha <- piranha[,c("seqnames","start","end","name","score","strand")]
      # piranha.res.df <- readPeak(piranha,chr = chr,start = start,width = width)
      # bed.piranha <- plotPeak(bed = piranha.res.df,chr = chr,start = start,width = width,ylab = "")
      # # bed.piranha
      # res.list2[[paste0(dst,"_",smp,"PiranhaOnly")]] <- bed.piranha
      # 
      # CLIPper <- data.table::fread(bed_clipper,sep = "\t",header = F,stringsAsFactors = F,check.names = F) # 1-based ?
      # if(nrow(CLIPper)==0){
      #   CLIPper <- data.frame("V1"="","V2"="","V3"="","V4"="","V5"="","V6"="","V7"="","V8"="")
      # }
      # colnames(CLIPper) <- c("seqnames","start","end","name","score","strand") # score: pval
      # #CLIPper <- CLIPper[,c("seqnames","start","end","name","score","strand")]
      # CLIPper.res.df <- readPeak(CLIPper,chr = chr,start = start,width = width)
      # bed.clipper <- plotPeak(bed = CLIPper.res.df,chr = chr,start = start,width = width,ylab = "")
      # # bed.clipper
      # res.list2[[paste0(dst,"_",smp,"CLIPperOnly")]] <- bed.clipper
      # 
      # CLAM <- data.table::fread(bed_clam,sep = "\t",header = F,stringsAsFactors = F,check.names = F) # 1-based ?
      # if(nrow(CLAM)==0){
      #   CLAM <- data.frame("V1"="","V2"="","V3"="","V4"="","V5"="","V6"="")
      # }      
      # colnames(CLAM) <- c("seqnames","start","end","name","type","strand")
      # CLAM$score <- 1
      # CLAM <- CLAM[,c("seqnames","start","end","name","score","strand")]
      # clam.res.df <- readPeak(CLAM,chr = chr,start = start,width = width)
      # bed.clam <- plotPeak(bed = clam.res.df,chr = chr,start = start,width = width,ylab = "")
      # # bed.clam
      # res.list2[[paste0(dst,"_",smp,"CLAMOnly")]] <- bed.clam 
      # 
      # localmaxEM2 <- data.table::fread(bed_localmaxEM2,sep = "\t",header = F,stringsAsFactors = F,check.names = F) # 1-based ?
      # if(nrow(localmaxEM2)==0){
      #   localmaxEM2 <- data.frame("V1"="","V2"="","V3"="","V4"="","V5"="","V6"="")
      # } 
      # localmaxEM2 <- localmaxEM2[,1:6]
      # colnames(localmaxEM2) <- c("seqnames","start","end","name","type","strand")
      # localmaxEM2$score <- 1
      # localmaxEM2 <- localmaxEM2[,c("seqnames","start","end","name","score","strand")]
      # localmaxEM2.res.df <- readPeak(localmaxEM2,chr = chr,start = start,width = width)
      # bed.localmaxEM2 <- plotPeak(bed = localmaxEM2.res.df,chr = chr,start = start,width = width,ylab = "") # , color = "firebrick", fill = "firebrick"
      # # bed.localmaxEM2
      # res.list2[[paste0(dst,"_",smp,"exPeakOnly")]] <- bed.localmaxEM2
      
      
      ## collect
      tmp2 <- ggpubr::ggarrange(plotlist = res.list2, ncol = 1,heights = c(5, 1,1,1,1, 1,1 ),align = "v")  # 1,1,1,1
      ggsave(plot = tmp2, filename = paste0("CLIPpeaks_bw_bed_",dst,"_",smp,"-",chr,".",start,".",width,".pdf"),width = 12, height = 10)  # height=4
    }
  }
}


for (i in 1:length(genes.list2)){
  #i <- 1
  chr <- genes.list2[[i]]$chr
  start <- genes.list2[[i]]$start
  width <- genes.list2[[i]]$width
  print(paste0(chr,".",start,".",width))
  
  for (dst in names(meta.list2)){
    # dst <- "GSE71008_NCpool"
    print(dst)
    if (grepl("AGO2_IP|FTC_small",dst,perl=T)){
      pre="call_peak_dedup" 
    } else {
      pre="call_peak_all"
    }
    #dst <- "GSE50676"
    #meta.list2[["GSE50676"]] <- c("GM12878")
    smps <- meta.list2[[dst]]
    res.list2 <- list()
    for (smp in smps){
      #smp <- "NCpool"
      print(smp)
      
      read_bw <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/",pre,"/tbigwig_RNA_EM/",smp,".transcriptome.bigWig")
      # x <- rtracklayer::import.bw(read_bw) # 1-based ?
      #x <- as.data.frame(x)
      bw.res.df <- readBW(read_bw,chr = chr,start = start,width = width) 
      res.list2[[paste0(dst,"_",smp,"_EM")]] <- plotBW(bw.res.df,chr = chr,start = start,width = width,ylab = "",xlab = "",title = "") # paste0(dst,"_",smp,"\n",chr,"_",start,"_",width)
    }
    
    # gold peak bed
    gold <- data.table::fread(bed_gold,sep = "\t",header = F,stringsAsFactors = F,check.names = F) # 1-based ?
    if(nrow(gold)==0){
      gold <- data.frame("V1"="","V2"="","V3"="","V4"="","V5"="","V6"="")
    } 
    gold <- gold[,1:6]
    colnames(gold) <- c("seqnames","start","end","name","type","strand")
    gold$score <- 1
    gold <- gold[,c("seqnames","start","end","name","score","strand")]
    gold.res.df <- readPeak(gold,chr = chr,start = start,width = width)
    bed.gold <- plotPeak(bed = gold.res.df,chr = chr,color = "salmon",fill = "salmon",start = start,width = width, ylab = "")
    # bed.gold
    res.list2[[paste0(dst,"_","gold")]] <- bed.gold
    
    ## collect
    tmp2 <- ggpubr::ggarrange(plotlist = res.list2, ncol = 1,heights = c(5,5,5, 1),align = "v") 
    ggsave(plot = tmp2, filename = paste0("CLIPpeaks_bw_",dst,"_","-",chr,".",start,".",width,".pdf"),width = 12, height = 5)  # height=4
  }
}

## select best eg.
ref <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)
table(ref$transcript_type)
# bt="/BioII/lulab_b/baopengfei/biosoft/bedtools"
# $bt intersect -u -s \
#   -a /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE50676/call_peak_all/expeak_by_sample/b5_d50_p1/intersect/GM12891_expeakOnly.bed6 \
#   -b /lulabdata/baopengfei/shared_reference/RBP/splitByRBP/AGO2_tx.bed > tmp/GM12891_expeakOnly_intersectRBP.bed6
tmp <- read.table("tmp/GM12891_expeakOnly_intersectRBP.bed6")
tmp$RNA <- ref$transcript_type[match(tmp$V1,ref$transcript_id)]
tmp$RNA.len <- ref$tx.length[match(tmp$V1,ref$transcript_id)]
#tmp <- tmp[tmp$V5>=5,]

#limit to mRNA peak, overlap with ribosnitch
# $bt intersect -u -s \
#   -a /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE50676/call_peak_all/expeak_by_sample/b5_d50_p1/intersect/GM12891_expeakOnly.bed6 \
#   -b /lulabdata/baopengfei/shared_reference/2014_nature_riboSNitches/RiboSNitches_hg38_tx.bed > tmp/GM12891_expeakOnly_intersectRiboSNitch.bed6
tmp2 <- read.table("tmp/GM12891_expeakOnly_intersectRiboSNitch.bed6")$V4
tmp <- tmp[tmp$V4 %in% tmp2,] # 5 mRNA left

# ENST00000416016_____2,1097,1178 (1336)
# ENST00000341446_____8,931,979   (6648)
# ENST00000456707_____1,523,591   (609 )
# ENST00000435041_____2,912,999   (4290)



# ribosnitch <- read.table("/lulabdata/baopengfei/shared_reference/2014_nature_riboSNitches/PARS_RiboSNitches.txt", header = T,sep = "\t", check.names = F, stringsAsFactors = F)
# ribosnitch$start <- ribosnitch$position-1
# ribosnitch$score <- 1
# ribosnitch$strand <- "+"
# ribosnitch.bed1 <- ribosnitch[,c("chr","start","position","SNPID","score","strand")]
# ribosnitch$strand <- "-"
# ribosnitch.bed2 <- ribosnitch[,c("chr","start","position","SNPID","score","strand")]
# ribosnitch.bed <- rbind(ribosnitch.bed1,ribosnitch.bed2)
# table(ribosnitch.bed$strand)
# 
# #utils::download.file('ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz',"/lulabdata/baopengfei/shared_reference/ucsc/liftover/hg19ToHg38.over.chain.gz")
# chain <- rtracklayer::import.chain("/lulabdata/baopengfei/shared_reference/ucsc/liftover/hg19ToHg38.over.chain")
# ribosnitch.gr <- GRanges(seqnames = ribosnitch.bed$chr, name = ribosnitch.bed$SNPID, score=ribosnitch.bed$score, strand = ribosnitch.bed$strand, 
#                          ranges = IRanges(start = ribosnitch.bed$start+1,width = 1)) # 0-base to 1-base
# ribosnitch.gr[1]@ranges
# hg38 <- rtracklayer::liftOver(x = ribosnitch.gr, chain = chain)
# length(ribosnitch.gr) # 3814
# length(hg38) # 3814
# 
# tmp <- lapply(hg38,GRange2bed)
# length(tmp) # 3814
# hg38.bed <- as.data.frame(do.call("rbind",tmp))
# nrow(hg38.bed) # 3808
# write.table(hg38.bed,"/lulabdata/baopengfei/shared_reference/2014_nature_riboSNitches/PARS_RiboSNitches_hg38_gn.bed",quote = F,sep = "\t",row.names = F, col.names = F)
# 
# # fail <- ribosnitch.bed$SNPID[!(ribosnitch.bed$SNPID %in% hg38.bed$name)]
# # ribosnitch.bed.fail <- ribosnitch.bed[ribosnitch.bed$SNPID %in% fail,]
# # ribosnitch.gr.fail <- GRanges(seqnames = ribosnitch.bed.fail$chr, name = ribosnitch.bed.fail$SNPID, score=ribosnitch.bed.fail$score, strand = ribosnitch.bed.fail$strand, 
# #                          ranges = IRanges(start = ribosnitch.bed.fail$start+1,width = 1)) 
# # rtracklayer::liftOver(x = ribosnitch.gr.fail, chain = chain)
# #seem some sequence fail to liftover
# #
# 
# 
# 
# ribosnitch <- read.table("/lulabdata/baopengfei/shared_reference/2014_nature_riboSNitches/eQTL_RiboSNitches.txt", header = T,sep = "\t", check.names = F, stringsAsFactors = F)
# colnames(ribosnitch) <- tolower(colnames(ribosnitch))
# table(duplicated(ribosnitch$snpid))
# table(ribosnitch$snpid==ribosnitch$`eqtl snpid`)
# table(duplicated(ribosnitch$`eqtl snpid`))
# ribosnitch$start <- ribosnitch$position-1
# ribosnitch$score <- 1
# ribosnitch$strand <- "+"
# ribosnitch.bed1 <- ribosnitch[,c("chr","start","position","snpid","score","strand")]
# ribosnitch$strand <- "-"
# ribosnitch.bed2 <- ribosnitch[,c("chr","start","position","snpid","score","strand")]
# ribosnitch.bed <- rbind(ribosnitch.bed1,ribosnitch.bed2)
# table(ribosnitch.bed$strand)
# 
# #utils::download.file('ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz',"/lulabdata/baopengfei/shared_reference/ucsc/liftover/hg19ToHg38.over.chain.gz")
# chain <- rtracklayer::import.chain("/lulabdata/baopengfei/shared_reference/ucsc/liftover/hg19ToHg38.over.chain")
# ribosnitch.gr <- GRanges(seqnames = ribosnitch.bed$chr, name = ribosnitch.bed$snpid, score=ribosnitch.bed$score, strand = ribosnitch.bed$strand, 
#                          ranges = IRanges(start = ribosnitch.bed$start+1,width = 1)) # 0-base to 1-base
# 
# hg38 <- rtracklayer::liftOver(x = ribosnitch.gr, chain = chain)
# length(ribosnitch.gr) # 420
# length(hg38) # 420
# tmp <- lapply(hg38,GRange2bed)
# length(tmp) # 420
# hg38.bed <- as.data.frame(do.call("rbind",tmp))
# nrow(hg38.bed) # 420
# write.table(hg38.bed,"/lulabdata/baopengfei/shared_reference/2014_nature_riboSNitches/eQTL_RiboSNitches_hg38_gn.bed",quote = F,sep = "\t",row.names = F, col.names = F)
# 
# #




# 4.3 plot single sample bw + 4 peak + gold peak (i-pico, RBP e.g.) ------------------------------------------------------------
setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev")
conflict_prefer("which", "base")

# bed_gold <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/tmp/6dst_filter3_11RNA.bed"
RBPs <- c("ABCF1","SFPQ","SUB1", "METAP2","RPS11","FKBP4","DDX3X",
          "BCLAF1","HNRNPU","APOBEC3C","SRSF9","NPM1","YWHAG","SRSF3","GRWD1",
          "all")
bed_RBP <- list()
for(i in 1:length(RBPs)){
  bed_RBP[[RBPs[i]]] <- paste0("/BioII/lulab_b/baopengfei/shared_reference/RBP/splitByRBP/",RBPs[i],"_tx.bed")
}


meta.list2 <- list()
meta.list2[["i-pico"]] <- c("L38-4-BC5","L39-5-BC2","L39-5-BC3","L41-7-BC5","L40-6-BC2","L38-4-BC3") # i-pico: DT1


tx.tbl <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/transcript_table/all_newTxID.txt",sep = "\t")
tx.tbl[1:3,]
table(tx.tbl$transcript_type)
tx.tbl <- tx.tbl[tx.tbl$transcript_type %in% c("mRNA","lncRNA"),]
tx.tbl$gene_id2 <- unlist(sapply(strsplit(tx.tbl$gene_id,".",fixed=T),"[",1))

tx.len <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt",sep = "\t")
tx.len[1:2,]
table(tx.len$transcript_type)
tx.len <- tx.len[tx.len$transcript_type %in% c("mRNA","lncRNA"),]
tx.len <- tx.len[!is.na(tx.len$tx.length),]

tx.tbl <- tx.tbl[tx.tbl$transcript_id %in% tx.len$transcript_id,]

res$raw.enst <- tx.tbl$transcript_id[match(res$ensg,tx.tbl$gene_id2)]
res$raw.enst2 <- gsub(".","_____",res$raw.enst,fixed = T)
table(is.na(res$raw.enst ))
res <- res[!is.na(res$raw.enst),]

res$tx.length <- tx.len$tx.length[match(res$raw.enst,tx.len$transcript_id)]

# #table(is.na(tx.len$tx.length))
# res[1:3,]
# table(res$raw.enst2 %in% unique(tx.tbl$transcript_id))
# table(res$raw.enst2 %in% unique(tx.len$transcript_id))
# 
# res <- res[res$raw.enst2 %in% unique(tx.tbl$transcript_id),]
# table(res$ensg %in% tx.tbl$gene_id2)
# res <- res[res$ensg %in% tx.tbl$gene_id2,]
# res$raw.enst.expeak <- tx.tbl$transcript_id[match(res$ensg,tx.tbl$gene_id2)]
# 
# res$raw.enst.expeak.len <- tx.len$tx.length[match(res$raw.enst.expeak,tx.len$transcript_id)]
# res[1:3,]
# table(is.na(res$raw.enst.expeak.len))
# table(res$raw.enst2==res$raw.enst.expeak)
# #seems most version of tx are updated in zhanqing
# #FALSE  TRUE 
# #7474   365 
#res <- res[!is.na(res$raw.enst.expeak.len),]
top <-10
genes.list2 <- list()
for (i in 1:top){
  # i <- 1
  genes.list2[[i]] <- data.frame(chr=res[res$logFC<0 & res$FDR<=0.1,]$raw.enst[i],start=0,width=res[res$logFC<0 & res$FDR<=0.1,]$tx.length[i])
  genes.list2[[i+top]] <- data.frame(chr=res[res$logFC>0 & res$FDR<=0.1,]$raw.enst[i+top],start=0,width=res[res$logFC>0 & res$FDR<=0.1,]$tx.length[i+top])
}


# #genes.list[[1]] <- data.frame(chr="ENST00000315707.3",start=400,width=300) # reveiw of ipico: lnc
# genes.list[[1]] <- data.frame(chr="ENST00000602361.1",start=0,width=268) # reveiw of ipico: top EV
# genes.list[[1]] <- data.frame(chr="ENST00000530124.5",start=0,width=1946) # reveiw of ipico: top EV

# #kept:
#   ENST00000440597_____1.0.600
#   ENST00000549592_____1.0.3000
#   ENST00000319041_____6.0.876
#   ENST00000484216_____1.0.827
#   ENST00000602361_____1.0.268
genes.list2 <- list()
genes.list2[[1]] <- data.frame(chr="ENST00000440597_____1",start=0,width=600)
genes.list2[[2]] <- data.frame(chr="ENST00000549592_____1",start=0,width=3000)
genes.list2[[3]] <- data.frame(chr="ENST00000319041_____6",start=0,width=876)
genes.list2[[4]] <- data.frame(chr="ENST00000484216_____1",start=0,width=827)
genes.list2[[5]] <- data.frame(chr="ENST00000602361_____1",start=0,width=268)

for (i in 1:length(genes.list2)){
  #i <- 1
  chr <- genes.list2[[i]]$chr
  start <- as.numeric(genes.list2[[i]]$start)
  width <- as.numeric(genes.list2[[i]]$width)
  print(paste0(chr,".",start,".",width))
  
  for (dst in names(meta.list2)){
    # dst <- "GSE71008_NCpool"
    print(dst)
    if (grepl("AGO2_IP|FTC_small|pico",dst,perl=T)){
      pre="call_peak_dedup" 
    } else {
      pre="call_peak_all"
    }
    #dst <- "GSE50676"
    #meta.list2[["GSE50676"]] <- c("GM12878")
    smps <- meta.list2[[dst]]
    
    # plot multi smps bw
    res.list2 <- list()
    for (smp in smps){
      #smp <- "NCpool"
      print(smp)
      read_bw <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/",pre,"/tbigwig_11RNA_primary/",smp,".transcriptome.bigWig") # tbigwig_RNA_EM, tbigwig_11RNA_primary
      bw.res.df <- readBW(read_bw,chr = chr,start = start,width = width) 
      res.list2[[paste0(dst,"_",smp)]] <- plotBW(bw.res.df,chr = chr,start = start,width = width,ylab = smp,xlab = "",title = "") # paste0(dst,"_",smp,"\n",chr,"_",start,"_",width)  
    }
    #ggpubr::ggarrange(plotlist = res.list2, ncol = 1,heights = c(1,1,1,1,1,1 ),align = "v")  # 1,1,1,1
  
    # plot multi RBPs bed
    res.list2.2 <- list()
    for(i in 1:length(RBPs)){
      # gold peak bed
      gold <- data.table::fread(bed_RBP[[RBPs[i]]],sep = "\t",header = F,stringsAsFactors = F,check.names = F) # 1-based ?
      if(nrow(gold)==0){
        gold <- data.frame("V1"="","V2"="","V3"="","V4"="","V5"="","V6"="")
      } 
      gold <- gold[,1:6]
      colnames(gold) <- c("seqnames","start","end","name","type","strand")
      gold$score <- 1
      gold <- gold[,c("seqnames","start","end","name","score","strand")]
      gold.res.df <- readPeak(gold,chr = chr,start = start,width = width)
      bed.gold <- plotPeak(bed = gold.res.df,chr = chr,color = "grey30",fill = "grey30",start = start,width = width, ylab = RBPs[i])
      # bed.gold
      res.list2.2[[paste0(dst,"_",RBPs[i])]] <- bed.gold
    }
    # ggpubr::ggarrange(plotlist = res.list2.2, ncol = 1,heights = c(1,1,1,1,1,1,1 ),align = "v")  # 1,1,1,1
    
    ## collect
    #ggsave(plot = tmp2, filename = paste0("peaks_bw_bed_",dst,"_",smp,"-",chr,".",start,".",width,".pdf"),width = 12, height = 10)
    tmp2 <- ggpubr::ggarrange(plotlist = c(res.list2,res.list2.2), ncol = 1,heights = c(rep(3,length(res.list2)),rep(1,length(res.list2.2)) ),align = "v")  # 1,1,1,1
    ggsave(plot = tmp2, filename = paste0("ipico_peaks_bw_bed_",dst,"_",chr,".",start,".",width,".pdf"),width = 12, height = 16) 
  }
}

#



# 5.EM/raw double bw + 4 peak (Fig5, suppl Fig2?)  ------------------------------------
setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev")



#EM effect not significant, and only expeak has peak record, thus we use part 1-3 for fig5 instead
#TCGA_small also seem no signal (even sample-wise also no signal)


meta.list2 <- list()
meta.list2[["GSE71008"]] <- "SAMN03863396" #c("SAMN03863524","SAMN03863398","SAMN03863475","SAMN03863560")  # "SAMN03863475","SAMN03863524"; "SAMN03863398","SAMN03863560";
meta.list2[["GSE94582"]] <- c("NEBNext_Lab1","NEBNext_Lab2","N4_B_Lab5","N4_A_Lab5") # "NEBNext_Lab1","NEBNext_Lab2"; "N4_B_Lab5";
meta.list2[["Phospho-RNA-seq"]] <- c("ULMC157_none","ULMC123_none","ULMC157_T4PNK","ULMC123_T4PNK") #"ULMC157_none"; "ULMC157_T4PNK", "ULMC148_T4PNK", "ULMC123_T4PNK"; 
#meta.list2[["TCGA_small"]] <- c("TCGA-A7-A0CE-11A-21R-A090-13_mirna_gdc_realn") 
# pre="call_domain_withRepeats_all" # call_domain_withRepeats_all

genes.list2 <- list()
# genes.list2[[1]] <- data.frame(chr="SSU____rRNA_Hsa__chr22___11629545____11631288_pos",start=400,width=300)
# genes.list2[[2]] <- data.frame(chr="SSU____rRNA_Hsa__chr5___175114739____175115174_neg",start=0,width=435)
# 
# genes.list2[[1]] <- data.frame(chr="NR_023363_____1",start=0,width=121) # rRNA
# # genes.list2[[1]] <- data.frame(chr="ENST00000384875_____3",start=0,width=90) # miRNA
# 
# genes.list2[[1]] <- data.frame(chr="NR_046235_____3",start=6000,width=500) # rRNA
genes.list2[[1]] <- data.frame(chr="ENST00000635274_____1",start=130,width=100) # 

# #genes.list2[[1]] <- data.frame(chr="NR_146148_____1",start=2400,width=200) # rRNA
# 
# 
# genes.list2[[3]] <- data.frame(chr="enhancer__chr7___69059541____69062941_pos",start=2500,width=1000)
# genes.list2[[4]] <- data.frame(chr="enhancer__chr9___5092200____5092800_pos",start=0,width=600)
# genes.list2[[5]] <- data.frame(chr="promoter__chr7___148941042____148941642_pos",start=0,width=600)
# genes.list2[[9]] <- data.frame(chr="G087360__chrM___14742____15537_neg",start=0,width=795)
# genes.list2[[10]] <- data.frame(chr="G033992__chr17___22520709____22523110_pos",start=0,width=2400)
# genes.list2[[1]] <- data.frame(chr="NR_146117_____1",start=6500,width=200) # rRNA
# genes.list2[[2]] <- data.frame(chr="NR_146117_____1",start=7830,width=200) # rRNA
# genes.list2[[3]] <- data.frame(chr="NR_146117_____1",start=8432,width=200) # rRNA
# genes.list2[[1]] <- data.frame(chr="NR_023371_____1",start=0,width=121) # rRNA
# genes.list2[[4]] <- data.frame(chr="G078490__chr7___63480989____63561170_neg",start=2052,width=200) # intron
# genes.list2[[2]] <- data.frame(chr="G078490__chr7___41705834____41907280_neg",start=3826,width=200) # intron

# NR_146117_____1,7830,200
# NR_146117_____1,8432,200
# NR_023371_____1,0,121
# G078490__chr7___63480989____63561170_neg,2052,200
# G078490__chr7___41705834____41907280_neg,3826,200




# #run in bash
# {
# preFix="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
# pre="call_peak_all"
# dst="GSE71008" # "Phospho-RNA-seq"
# 
# for smp in SAMN03863396  #{SAMN03863524,SAMN03863398,SAMN03863475,SAMN03863560} # {ULMC157_none,ULMC123_none,ULMC157_T4PNK,ULMC123_T4PNK}
# do
# echo $smp
#   #smp="TCGA-A7-A0CE-11A-21R-A090-13_mirna_gdc_realn"
#   dedup="" # -deduped
#   threads=10
# 
#   # 19 primary
#   tempDir="${preFix}/output/${dst}/${pre}/tbigwig_19RNA_primary"
#   mkdir -p $tempDir
#   mergeBamFile="${preFix}/output/$dst/tbam/${smp}/bam${dedup}/merge19_sort.bam" #"${tempDir}/${smp}.bam"
#   #samtools merge -f ${mergeBamFile} ${preFix}/output/$dst/tbam/${smp}/bam${dedup}/{intron,promoter,enhancer,repeats}_{for,rev}.bam
#   #bam2bed
#   bedFile="${tempDir}/${smp}_19RNA_primary.bed"
#   samtools view -F 256 -bh ${mergeBamFile} | bedtools bamtobed -i stdin | pigz -c -p ${threads} > ${bedFile}
#   #bed2bg
#   bgFile="${tempDir}/${dst}_${smp}_19RNA_primary.bg"
#   chrom_sizes="${preFix}/exSeek-dev/genome/hg38/chrom_sizes/transcriptome_genome_sort_uniq_newTxID"
#   bedtools genomecov -i ${bedFile} -g ${chrom_sizes} -bg -split | LC_COLLATE=C sort -T ${tempDir} -k1,1 -k2,2n > ${bgFile}
#   #bg2bw
#   bwFile="${tempDir}/${smp}.transcriptome.bigWig"
#   tempFile="${tempDir}/${smp}_tmp_bg"
#   cat ${bgFile} | bedtools sort > ${tempFile}
#   /BioII/lulab_b/baopengfei/biosoft/bedGraphToBigWig ${tempFile} ${chrom_sizes} ${bwFile}
#   #rm  ${bedFile} ${bgFile} ${tempFile}
# 
#   # 19 uniq
#   tempDir="${preFix}/output/${dst}/${pre}/tbigwig_19RNA_uniq"
#   mkdir -p $tempDir
#   mergeBamFile="${preFix}/output/$dst/tbam/${smp}/bam${dedup}-EM/merge19_sort/unique.sorted.revertFullLengthReads.sorted.bam" #"${tempDir}/${smp}.bam"
#   #samtools merge -f ${mergeBamFile} ${preFix}/output/$dst/tbam/${smp}/bam${dedup}/{intron,promoter,enhancer,repeats}_{for,rev}.bam
#   #bam2bed
#   bedFile="${tempDir}/${smp}_19NA_uniq.bed"
#   samtools view -bh ${mergeBamFile} | bedtools bamtobed -i stdin | pigz -c -p ${threads} > ${bedFile}
#   #bed2bg
#   bgFile="${tempDir}/${dst}_${smp}_19RNA_uniq.bg"
#   chrom_sizes="${preFix}/exSeek-dev/genome/hg38/chrom_sizes/transcriptome_genome_sort_uniq_newTxID"
#   bedtools genomecov -i ${bedFile} -g ${chrom_sizes} -bg -split | LC_COLLATE=C sort -T ${tempDir} -k1,1 -k2,2n > ${bgFile}
#   #bg2bw
#   bwFile="${tempDir}/${smp}.transcriptome.bigWig"
#   tempFile="${tempDir}/${smp}_tmp_bg"
#   cat ${bgFile} | bedtools sort > ${tempFile}
#   /BioII/lulab_b/baopengfei/biosoft/bedGraphToBigWig ${tempFile} ${chrom_sizes} ${bwFile}
#   #rm ${bedFile} ${bgFile} ${tempFile}
# done
# }



#EM + primary + uniq bw 
for (i in 1:length(genes.list2)){
  #i <- 1
  chr <- genes.list2[[i]]$chr
  start <- genes.list2[[i]]$start
  width <- genes.list2[[i]]$width
  print(paste0(chr,".",start,".",width))
  
  for (dst in names(meta.list2)){
    # dst <- "GSE71008_NCpool"
    print(dst)
    if (grepl("AGO2_IP|FTC_small|TGIRT",dst,perl=T)){
      pre="call_peak_dedup" 
    } else {
      pre="call_peak_all"
    }
    #dst <- "GSE50676"
    #meta.list2[["GSE50676"]] <- c("GM12878")
    smps <- meta.list2[[dst]]
    for (smp in smps){
      #smp <- "NCpool"
      res.list2 <- list()
      print(smp)
      
      read_bw1 <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/",pre,"/tbigwig_RNA_EM/",smp,".transcriptome.bigWig")
      # x <- rtracklayer::import.bw(read_bw1) # 1-based ?
      #x <- as.data.frame(x)
      bw.res.df1 <- readBW(read_bw1,chr = chr,start = start,width = width) 
      res.list2[[paste0(dst,"_",smp,"_","EM")]] <- plotBW(bw.res.df1,chr = chr,start = start,width = width,ylab = "",xlab = "",title = "") # paste0(dst,"_",smp,"\n",chr,"_",start,"_",width)  
      # res.list2[[paste0(dst,"_",smp)]] <- gg.gap(plot=res.list2[[1]],rel_heights = c(0.05,1), segments=c(20,120),ylim = c(0,140)) # gg.gap truncate y axis for specific tx, # + ggbreak::scale_y_break(breaks = c(50, 300)), ggbreak seem not fit ggarrange ?
      # ggsave(plot = res.list2[[paste0(dst,"_",smp)]], filename = paste0("peaks_bw_bed_",dst,"_",smp,"-",chr,".",start,".",width,".2.pdf"),width = 13, height = 3)
      
      read_bw2 <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/",pre,"/tbigwig_19RNA_primary/",smp,".transcriptome.bigWig")
      # x <- rtracklayer::import.bw(read_bw2) # 1-based ?
      #x <- as.data.frame(x)
      bw.res.df2 <- readBW(read_bw2,chr = chr,start = start,width = width) 
      res.list2[[paste0(dst,"_",smp,"_","primary")]] <- plotBW(bw.res.df2,chr = chr,start = start,width = width,ylab = "",xlab = "",title = "") # paste0(dst,"_",smp,"\n",chr,"_",start,"_",width)  
      # res.list2[[paste0(dst,"_",smp)]] <- gg.gap(plot=res.list2[[1]],rel_heights = c(0.05,1), segments=c(20,120),ylim = c(0,140)) # gg.gap truncate y axis for specific tx, # + ggbreak::scale_y_break(breaks = c(50, 300)), ggbreak seem not fit ggarrange ?
      # ggsave(plot = res.list2[[paste0(dst,"_",smp)]], filename = paste0("peaks_bw_bed_",dst,"_",smp,"-",chr,".",start,".",width,".2.pdf"),width = 13, height = 3)
 
      read_bw3 <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/",pre,"/tbigwig_19RNA_uniq/",smp,".transcriptome.bigWig")
      # x <- rtracklayer::import.bw(read_bw2) # 1-based ?
      #x <- as.data.frame(x)
      bw.res.df3 <- readBW(read_bw3,chr = chr,start = start,width = width) 
      res.list2[[paste0(dst,"_",smp,"_","uniq")]] <- plotBW(bw.res.df3,chr = chr,start = start,width = width,ylab = "",xlab = "",title = "") # paste0(dst,"_",smp,"\n",chr,"_",start,"_",width)  
      # res.list2[[paste0(dst,"_",smp)]] <- gg.gap(plot=res.list2[[1]],rel_heights = c(0.05,1), segments=c(20,120),ylim = c(0,140)) # gg.gap truncate y axis for specific tx, # + ggbreak::scale_y_break(breaks = c(50, 300)), ggbreak seem not fit ggarrange ?
      # ggsave(plot = res.list2[[paste0(dst,"_",smp)]], filename = paste0("peaks_bw_bed_",dst,"_",smp,"-",chr,".",start,".",width,".2.pdf"),width = 13, height = 3)
      
      
      # peak
      bed_piranha <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/",pre,"/piranha_by_sample/b5_p01/",smp,".bed")
      bed_clipper <-  paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/",pre,"/clipper_by_sample/b5_p05/",smp,".bed")
      bed_clam <-  paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/",pre,"/clam_by_sample/b5_p005/",smp,".bed")
      # bed_localmax <-  paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/",pre,"/domains_localmax_significant/b5_d05_p05/",smp,".bed")
      bed_localmaxEM2 <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/",pre,"/expeak_by_sample/b5_d50_p1/",smp,".bed")
      ## piranha peak bed
      piranha <- data.table::fread(bed_piranha,sep = "\t",header = F,stringsAsFactors = F,check.names = F) # 1-based ?
      if(nrow(piranha)==0){
        piranha <- data.frame("V1"="","V2"="","V3"="","V4"="","V5"="","V6"="","V7"="")
      }
      colnames(piranha) <- c("seqnames","start","end","name","score","strand","pval")
      piranha <- piranha[,c("seqnames","start","end","name","score","strand")]
      piranha.res.df <- readPeak(piranha,chr = chr,start = start,width = width)
      bed.piranha <- plotPeak(bed = piranha.res.df,chr = chr,start = start,width = width,ylab = "")
      # bed.piranha
      res.list2[[paste0(dst,"_",smp,"Piranha")]] <- bed.piranha
      
      
      ## CLIPper peak bed
      CLIPper <- data.table::fread(bed_clipper,sep = "\t",header = F,stringsAsFactors = F,check.names = F) # 1-based ?
      if(nrow(CLIPper)==0){
        CLIPper <- data.frame("V1"="","V2"="","V3"="","V4"="","V5"="","V6"="","V7"="","V8"="")
      }
      colnames(CLIPper) <- c("seqnames","start","end","name","score","strand","V7","V8") # score: pval
      CLIPper <- CLIPper[,c("seqnames","start","end","name","score","strand")]
      CLIPper.res.df <- readPeak(CLIPper,chr = chr,start = start,width = width)
      bed.clipper <- plotPeak(bed = CLIPper.res.df,chr = chr,start = start,width = width,ylab = "")
      # bed.clipper
      res.list2[[paste0(dst,"_",smp,"CLIPper")]] <- bed.clipper
      
      
      ## CLAM peak bed
      CLAM <- data.table::fread(bed_clam,sep = "\t",header = F,stringsAsFactors = F,check.names = F) # 1-based ?
      if(nrow(CLAM)==0){
        CLAM <- data.frame("V1"="","V2"="","V3"="","V4"="","V5"="","V6"="")
      }      
      colnames(CLAM) <- c("seqnames","start","end","name","type","strand")
      CLAM$score <- 1
      CLAM <- CLAM[,c("seqnames","start","end","name","score","strand")]
      clam.res.df <- readPeak(CLAM,chr = chr,start = start,width = width)
      bed.clam <- plotPeak(bed = clam.res.df,chr = chr,start = start,width = width,ylab = "")
      # bed.clam
      res.list2[[paste0(dst,"_",smp,"CLAM")]] <- bed.clam 
      
      # ## localma peak bed
      # localmax <- data.table::fread(bed_localmax,sep = "\t",header = F,stringsAsFactors = F,check.names = F) # 1-based ?
      # #localmax[1:3,]
      # localmax <- localmax[,1:6]
      # colnames(localmax) <- c("seqnames","start","end","name","type","strand")
      # localmax$score <- 1
      # localmax <- localmax[,c("seqnames","start","end","name","score","strand")]
      # localmax.res.df <- readPeak(localmax,chr = chr,start = start,width = width)
      # bed.localmax <- plotPeak(bed = localmax.res.df,chr = chr,start = start,width = width,ylab = "LocalMax")
      # # bed.localmax
      # res.list2[[paste0(dst,"_",smp,"LocalMax")]] <- bed.localmax
      
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
      bed.localmaxEM2 <- plotPeak(bed = localmaxEM2.res.df,chr = chr,start = start,width = width, ylab = "")
      # bed.localmaxEM2
      res.list2[[paste0(dst,"_",smp,"exPeak")]] <- bed.localmaxEM2
      
      
      # anno. rmsk (11RNA only)
      rmsk <- "/lulabdata/baopengfei/shared_reference/hg38/backup/rmsk_tx.bed6"
      rmsk <- data.table::fread(rmsk,sep = "\t",header = F,stringsAsFactors = F,check.names = F) # 1-based ?
      if(nrow(rmsk)==0){
        rmsk <- data.frame("V1"="","V2"="","V3"="","V4"="","V5"="","V6"="")
      } 
      rmsk <- rmsk[,1:6]
      colnames(rmsk) <- c("seqnames","start","end","name","type","strand")
      rmsk$score <- 1
      rmsk <- rmsk[,c("seqnames","start","end","name","score","strand")]
      rmsk.res.df <- readPeak(rmsk,chr = chr,start = start,width = width)
      bed.rmsk <- plotPeak(bed = rmsk.res.df,chr = chr,start = start,width = width, fill = "salmon", color = "salmon", ylab = "")
      # bed.rmsk
      res.list2[[paste0(dst,"_",smp,"rmsk")]] <- bed.rmsk
      
      # anno. rmsk.alu (11RNA only)
      rmsk <- "/lulabdata/baopengfei/shared_reference/hg38/backup/rmsk_Alu_tx.bed6"
      rmsk <- data.table::fread(rmsk,sep = "\t",header = F,stringsAsFactors = F,check.names = F) # 1-based ?
      if(nrow(rmsk)==0){
        rmsk <- data.frame("V1"="","V2"="","V3"="","V4"="","V5"="","V6"="")
      } 
      rmsk <- rmsk[,1:6]
      colnames(rmsk) <- c("seqnames","start","end","name","type","strand")
      rmsk$score <- 1
      rmsk <- rmsk[,c("seqnames","start","end","name","score","strand")]
      rmsk.res.df <- readPeak(rmsk,chr = chr,start = start,width = width)
      bed.rmsk <- plotPeak(bed = rmsk.res.df,chr = chr,start = start,width = width, fill = "salmon", color = "salmon", ylab = "")
      # bed.rmsk
      res.list2[[paste0(dst,"_",smp,"rmsk_Alu")]] <- bed.rmsk
      
      
      ## collect
      tmp2 <- ggpubr::ggarrange(plotlist = res.list2, ncol = 1,heights = c(5,5,5,1,1,1,1,1,1),align = "v") 
      ggsave(plot = tmp2, filename = paste0("peaks_bw_bed_",dst,"_",smp,"-",chr,".",start,".",width,".pdf"),width = 12, height = 8)
    }
  }
}



#macs2 bdgcmp -p 1 -m logFE -t /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008/call_peak_all/tbedgraph_RNA_EM/SAMN03863396.transcriptome.bedGraph -c /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008/call_peak_all/tbigwig_19RNA_primary/GSE71008_SAMN03863396_19RNA_primary.bg -o macs2_bdgcmp.bed
# macs2 bdgdiff -l 10 -g 3 --cutoff 3 \
#   --t1 /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008/call_peak_all/tbedgraph_RNA_EM/SAMN03863396.transcriptome.bedGraph \
#   --c1 /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008/call_peak_all/tbigwig_19RNA_uniq/GSE71008_SAMN03863396_19RNA_uniq.bg \
#   --t2 /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008/call_peak_all/tbigwig_19RNA_primary/GSE71008_SAMN03863396_19RNA_primary.bg \
#   --c2 /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008/call_peak_all/tbigwig_19RNA_uniq/GSE71008_SAMN03863396_19RNA_uniq.bg \
#   --o-prefix macs2_bdgdiff

#diff <- read.table("macs2_bdgcmp.bed",header = F)
#diff <- diff[diff$V1=="NR_146148_____1" & diff$V2>=2300 & diff$V2<=2600,]
diff <- read.table("macs2_bdgdiff_c3.0_cond1.bed",sep = "\t", skip = 1,header = F)
#diff <- diff[(diff$V3-diff$V2)<=100 & (diff$V3-diff$V2)>=10,]
diff$strand <- "+"

overlapAlu <- bedtoolsr::bt.intersect(u = T,wa = T,
                                      a="../output/GSE71008/call_peak_all/expeak_by_sample/b5_d50_p1/intersect_8DNA/SAMN03863396_gn.bed", 
                                      b="/lulabdata/baopengfei/shared_reference/hg38/backup/rmsk_Alu.bed6")
expeak <- read.table("../output/GSE71008/call_peak_all/expeak_by_sample/b5_d50_p1/intersect_8DNA/SAMN03863396.bed")
expeak <- expeak[expeak$V4 %in% overlapAlu$V4,]

expeak <- bedtoolsr::bt.intersect(u = T,wa = T,a=expeak[,1:6], b=diff[,1:6])



# # #peaks
# # AluJb__chr15___68202323____68202601_pos,135,155 # 278 (loose: >=1)
# # G041621__chr19___53051683____53057401_neg,550,564 # 5718
# # #tRNA____Cys____TGY__chr6___149027805____149027840_neg,0,23 # 35
# # G____rich__chr6___7146362____7146415_pos,29,49 # 53
# # #promoter__chr12___121806128____121809528_neg,914,934 # 3400
# # #promoter__chr16___30486277____30489677_pos,2926,2944 # 3400
# # #enhancer__chr9___2824800____2829600_neg,1639,1652 # 4800
# # enhancer__chr9___128906479____128909279_pos,568,586 # 2800

# enhancer__chr7___69059541____69062941_pos,,1000
# enhancer__chr9___5092200____5092800_pos,0,600
# promoter__chr7___148941042____148941642_pos,0,600
# SSU____rRNA_Hsa__chr22___11629545____11631288_pos,0,2000
# SSU____rRNA_Hsa__chr5___175114739____175115174_neg,0,435
# G033992__chr17___22520709____22523110_pos,0,2400
# G087360__chrM___14742____15537_neg,0,795
# # genes.list2[[1]] <- data.frame(chr="AluJb__chr15___68202323____68202601_pos",start=0,width=278)
# # genes.list2[[2]] <- data.frame(chr="G041621__chr19___53051683____53057401_neg",start=0,width=1000)
# # # genes.list2[[3]] <- data.frame(chr="tRNA____Cys____TGY__chr6___149027805____149027840_neg",start=0,width=35)
# # genes.list2[[4]] <- data.frame(chr="G____rich__chr6___7146362____7146415_pos",start=0,width=53)
# # #genes.list2[[5]] <- data.frame(chr="promoter__chr12___121806128____121809528_neg",start=500,width=1000)
# # #genes.list2[[6]] <- data.frame(chr="promoter__chr16___30486277____30489677_pos",start=2500,width=1000)
# # #genes.list2[[7]] <- data.frame(chr="enhancer__chr9___2824800____2829600_neg",start=1000,width=1000)
# # genes.list2[[8]] <- data.frame(chr="enhancer__chr9___128906479____128909279_pos",start=0,width=1000)
# 









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
