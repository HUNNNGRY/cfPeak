# func. for cfPeak
# last 2305 by bpf 
# b.p.f@qq.com

options(stringsAsFactors = F) # ,scipen = 4,digits = 5

# load lib ahead
suppressPackageStartupMessages(library(ggplot2))

# R common -----
detach_package <- function(pkg, character.only = FALSE)
{
  if(!character.only)
  {
    pkg <- deparse(substitute(pkg))
  }
  search_item <- paste("package", pkg, sep = ":")
  while(search_item %in% search())
  {
    detach(search_item, unload = TRUE, character.only = TRUE)
  }
}

tocapital <- function(my_string){
  capital_string <- paste(toupper(substring(my_string, 1, 1)), substring(my_string, 2), sep = "")
  return(capital_string)
}

splitVec <- function(vec,delim,fixed=F,perl=F,nth){
  #return(unlist(sapply(strsplit(vec,delim,fixed = TRUE,per=perl),function(x) x[nth])))
  return(as.character(lapply(strsplit(vec,delim,fixed = TRUE,per=perl),function(x) x[nth])))
}


#' Wrapper around mclapply to track progress (works on mac and linux)
#' 
#' Based on http://stackoverflow.com/questions/10984556
#' 
#' @param X         a vector (atomic or list) or an expressions vector. Other
#'                  objects (including classed objects) will be coerced by
#'                  ‘as.list’
#' @param FUN       the function to be applied to
#' @param ...       optional arguments to ‘FUN’
#' @param mc.preschedule see mclapply
#' @param mc.set.seed see mclapply
#' @param mc.silent see mclapply
#' @param mc.cores see mclapply
#' @param mc.cleanup see mclapply
#' @param mc.allow.recursive see mclapply
#' @param mc.progress track progress?
#' @param mc.style    style of progress bar (see txtProgressBar)
#'
#' @examples
#' x <- mclapply2(1:1000, function(i, y) Sys.sleep(0.01))
#' x <- mclapply2(1:3, function(i, y) Sys.sleep(1), mc.cores=1)
#' 
#' dat <- lapply(1:10, function(x) rnorm(100)) 
#' func <- function(x, arg1) mean(x)/arg1 
#' mclapply2(dat, func, arg1=10, mc.cores=2)

#library(parallel)
mclapply2 <- function(X, FUN, ..., 
                      mc.preschedule = TRUE, mc.set.seed = TRUE,
                      mc.silent = FALSE, mc.cores = getOption("mc.cores", 2L),
                      mc.cleanup = TRUE, mc.allow.recursive = TRUE,
                      mc.progress=TRUE, mc.style=3) 
{
  if (!is.vector(X) || is.object(X)) X <- as.list(X)
  
  if (mc.progress) {
    f <- fifo(tempfile(), open="w+b", blocking=T)
    p <- parallel:::mcfork()
    pb <- txtProgressBar(0, length(X), style=mc.style)
    setTxtProgressBar(pb, 0) 
    progress <- 0
    if (inherits(p, "masterProcess")) {
      while (progress < length(X)) {
        readBin(f, "double")
        progress <- progress + 1
        setTxtProgressBar(pb, progress) 
      }
      cat("\n")
      parallel:::mcexit()
    }
  }
  tryCatch({
    result <- mclapply(X, ..., function(...) {
      res <- FUN(...)
      if (mc.progress) writeBin(1, f)
      res
    }, 
    mc.preschedule = mc.preschedule, mc.set.seed = mc.set.seed,
    mc.silent = mc.silent, mc.cores = mc.cores,
    mc.cleanup = mc.cleanup, mc.allow.recursive = mc.allow.recursive
    )
    
  }, finally = {
    if (mc.progress) close(f)
  })
  result
}

#

# Id conversion func. ---------------------------------------------------------------
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

splitPeakFeatureID <- function(peakFeatureID){
  tmp.list1 <- strsplit(peakFeatureID,"|",fixed=T)
  tmp.df1 <- data.frame( feature=peakFeatureID,
                         RNA=unlist(sapply(tmp.list1,"[",2)),
                         TxNewID=unlist(sapply(tmp.list1,"[",3)),
                         # TxID=new2oldTxID(TxNewID),
                         PeakID=unlist(sapply(tmp.list1,"[",4)),
                         PeakStart=unlist(sapply(tmp.list1,"[",6)),
                         PeakEnd=unlist(sapply(tmp.list1,"[",7)),
                         RNA2=unlist(sapply(tmp.list1,"[",8))
                         
                         # PeakWidth=PeakEnd-PeakStart
  )
  tmp.df1$TxID <- new2oldTxID(tmp.df1$TxNewID)
  tmp.df1$TxID2 <- unlist(sapply(strsplit(tmp.df1$TxID,".",fixed=T),"[",1))
  tmp.df1$PeakWidth=as.numeric(tmp.df1$PeakEnd)-as.numeric(tmp.df1$PeakStart)
  rownames(tmp.df1) <- peakFeatureID
  return(tmp.df1)
}

ensg2symbolEntrezDF <- function(enst,returnall=T){
  library(mygene)
  #returnall=TRUE, you will get both duplicate or missing query terms
  tmp <- as.data.frame(mygene::queryMany(enst, scopes="ensembl.gene",fields=c("symbol","entrezgene"),species="human", returnall=returnall)) #$entrezgene
  rownames(tmp) <- tmp$query
  # print(duplicated())
  return(tmp)
}
enst2ensgEntrezDF <- function(enst,returnall=T){
  library(mygene)
  #returnall=TRUE, you will get both duplicate or missing query terms
  tmp <- as.data.frame(mygene::queryMany(enst, scopes="ensembl.transcript",fields=c("ensembl.gene","symbol","entrezgene"),species="human", returnall=returnall)) #$entrezgene
  rownames(tmp) <- tmp$query
  # print(duplicated())
  return(tmp)
}
enst2ensgDF <- function(enst,returnall=T){
  library(mygene)
  #returnall=TRUE, you will get both duplicate or missing query terms
  tmp <- as.data.frame(mygene::queryMany(enst, scopes="ensembl.transcript",fields=c("ensembl.gene"),species="human", returnall=returnall)) #$entrezgene
  rownames(tmp) <- tmp$query
  # print(duplicated())
  return(tmp)
}
enst2entrezDF <- function(enst,returnall=T){
  library(mygene)
  #returnall=TRUE, you will get both duplicate or missing query terms
  tmp <- as.data.frame(mygene::queryMany(enst, scopes="ensembl.transcript",fields=c("entrezgene"),species="human", returnall=returnall)) #$entrezgene
  rownames(tmp) <- tmp$query
  # print(duplicated())
  return(tmp)
}
enst2symbolDF <- function(enst,returnall=T){
  library(mygene)
  #returnall=TRUE, you will get both duplicate or missing query terms
  tmp <- as.data.frame(mygene::queryMany(enst, scopes="ensembl.transcript",fields=c("symbol"),species="human", returnall=returnall)) #$entrezgene
  rownames(tmp) <- tmp$query
  # print(duplicated())
  return(tmp)
}
symbol2ensgEntrezDF <- function(hgnc_symbol,returnall=T){
  library(mygene)
  #returnall=TRUE, you will get both duplicate or missing query terms
  tmp <- as.data.frame(mygene::queryMany(hgnc_symbol, scopes="symbol",fields=c("ensembl.gene","entrezgene"),species="human", returnall=returnall)) #$entrezgene
  rownames(tmp) <- tmp$query
  return(tmp)
}
# get1stEnsgMygene <- function(x){
#   # tmp[3,"ensembl"]
#   return(as.character(unlist(x))[1])
# }
# tmp$ensembl2 <- lapply(tmp$ensembl,get1stEnsgMygene)


# better use offline-mode BioMart
enst2ensgEntrezSymbolDescritptionBiotypeDF <- function(query, online=F){
  if(online){
    message("try online mart conversion")
    this.date <- unlist(sapply(strsplit(x = as.character(Sys.time()), split = " ",fixed = T),"[",1))
    tryCatch(
      {
        message("try main site: www.ensembl.org")
        mart <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl") # need online network
        saveRDS(object = mart, file = paste0("/BioII/lulab_b/baopengfei/shared_reference/ensembl/biomaRt.ensembl.hsapiens_gene_ensembl.",this.date,".rds"))
        # return(mart)
      },
      error = function(cond){
        asia.ensembl <- "http://asia.ensembl.org"
        message(paste("main site not avail, try mirror site: ",asia.ensembl,", other sites might also be avail"))
        mart <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl", host = asia.ensembl) # try mirror if main site down
        saveRDS(object = mart, file = paste0("/BioII/lulab_b/baopengfei/shared_reference/ensembl/biomaRt.ensembl.hsapiens_gene_ensembl.asia.",this.date,".rds"))
        # return(mart)
      }
    )
  }else{
    message("try offline mart conversion")
    mart <- readRDS("/BioII/lulab_b/baopengfei/shared_reference/ensembl/biomaRt.ensembl.hsapiens_gene_ensembl.2024-04-03.rds")
  }

  ann <- biomaRt::getBM(attributes=c("hgnc_symbol","entrezgene_id","ensembl_gene_id","ensembl_transcript_id","transcript_biotype","description"), # entrezgene --> entrezgene_id
                        "ensembl_transcript_id", unique(query), 
                        mart = mart, useCache = FALSE)
  colnames(ann) <- c("symbol","entrez","ENSG","ENST","biotype","description")
  return(ann)
}

ensg2enstEntrezSymbolDescritptionBiotypeDF <- function(query, online=F){
  if(online){
    message("try online mart conversion")
    this.date <- unlist(sapply(strsplit(x = as.character(Sys.time()), split = " ",fixed = T),"[",1))
    tryCatch(
      {
        message("try main site: www.ensembl.org")
        mart <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl") # need online network
        saveRDS(object = mart, file = paste0("/BioII/lulab_b/baopengfei/shared_reference/ensembl/biomaRt.ensembl.hsapiens_gene_ensembl.",this.date,".rds"))
        # return(mart)
      },
      error = function(cond){
        asia.ensembl <- "http://asia.ensembl.org"
        message(paste("main site not avail, try mirror site: ",asia.ensembl,", other sites might also be avail"))
        mart <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl", host = asia.ensembl) # try mirror if main site down
        saveRDS(object = mart, file = paste0("/BioII/lulab_b/baopengfei/shared_reference/ensembl/biomaRt.ensembl.hsapiens_gene_ensembl.asia.",this.date,".rds"))
        # return(mart)
      }
    )
  }else{
    message("try offline mart conversion")
    mart <- readRDS("/BioII/lulab_b/baopengfei/shared_reference/ensembl/biomaRt.ensembl.hsapiens_gene_ensembl.2024-04-03.rds")
  }
  
  ann <- biomaRt::getBM(attributes=c("hgnc_symbol","entrezgene_id","ensembl_gene_id","ensembl_transcript_id","transcript_biotype","description"), # entrezgene --> entrezgene_id
                        "ensembl_gene_id", unique(query), 
                        mart = mart, useCache = FALSE)
  colnames(ann) <- c("symbol","entrez","ENSG","ENST","biotype","description")
  return(ann)
}
snpId2gnCoord <- function(query){
#snp_ids = c("rs16828074", "rs17232800")
#snp_mart = biomaRt::useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")
#saveRDS(object = snp_mart, file = paste0("/BioII/lulab_b/baopengfei/shared_reference/ensembl/biomaRt.ensembl.hsapiens_snp.asia",".rds"))
snp_mart <- readRDS("/BioII/lulab_b/baopengfei/shared_reference/ensembl/biomaRt.ensembl.hsapiens_snp.asia.2024-06-11.rds")
# return(mart)
ann = biomaRt::getBM(attributes=c("refsnp_id", "chr_name", "chrom_start"), filters="snp_filter", useCache = FALSE,
                      values=query, mart=snp_mart)
return(ann)
}
#



# Bioformat conversion ----
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

## bam
split_cigar <- function(cigarStr) {
  # Regular expression to match numbers followed by letters
  matches <- regmatches(cigarStr, gregexpr("\\d+[MIDNSHP=X]", cigarStr))
  return(matches[[1]])
}
# Function to calculate alignment length from the CIGAR list
aln_length <- function(cigarlist) {
  tlength <- 0
  # Loop through each CIGAR operation
  #Op BAM Description
  #M 0 alignment match (can be a sequence match or mismatch)
  #I 1 insertion to the reference
  #D 2 deletion from the reference
  #N 3 skipped region from the reference
  #S 4 soft clipping (clipped sequences present in SEQ)
  #H 5 hard clipping (clipped sequences NOT present in SEQ)
  #P 6 padding (silent deletion from padded reference)  ???
  #= 7 sequence match
  #X 8 sequence mismatch
  for (cigar in cigarlist) {
    # Extract the operation and length
    operation <- substr(cigar, nchar(cigar), nchar(cigar))
    length <- as.numeric(substr(cigar, 1, nchar(cigar) - 1))
    
    # Check the operation and add length accordingly
    if (operation %in% c("M", "D", "=", "X")) {
      tlength <- tlength + length
    }
  }
  return(tlength)
}
# Example usage
# cigar_string <- "12M66998N9M"
# cigar_list <- split_cigar(cigar_string)
# alignment_length <- aln_length(cigar_list)
# print(alignment_length)


# Coverage plot --------------------------------
get.mat <- function(signal,region,signal.label,region.label,up=50,down=50,bin=10,ratio=0.3,mean_mode = "coverage",bg=0){
  # signal <- "output/GSE71008_NCpool/call_peak_all/tbed_RNA_EM/NCpool.bed.gz"
  # region <- "output/GSE71008_NCpool/call_peak_all/domains_clipper_by_sample/b5_p05/NCpool.bed"
  # up=50,down=50,bin=10,ratio=0.3
  # signal <- "/BioII/lulab_b/baopengfei/projects/RNA-structure/data/cell-2016-icShape/GSE74353_HS_293T_icSHAPE_InVivo_BaseReactivities_filterLen_exPeak_mRNA.bed"
  # region <- "output/GSE71008/call_peak_all/piranha_by_sample/b5_p01/intersect/SAMN03863400.bed6"
  print(signal)
  print(region)
  # signal.gr <- rtracklayer::import(con=as.character(signal)) # import(.bed) got error if input bed not proper formated
  # region.gr <- rtracklayer::import(con=as.character(region))
  signal.gr <- data.table::fread(as.character(signal))
  signal.gr <- GenomicRanges::GRanges(signal.gr$V1, IRanges(signal.gr$V2, signal.gr$V3), name=paste(signal.gr$V4,signal.gr$V1,signal.gr$V2,signal.gr$V3,signal.gr$V6,sep="_"), score=signal.gr$V5, strand=signal.gr$V6)
  region.gr <- data.table::fread(as.character(region))
  region.gr <- GenomicRanges::GRanges(region.gr$V1, IRanges(region.gr$V2, region.gr$V3), name=paste(region.gr$V4,region.gr$V1,region.gr$V2,region.gr$V3,region.gr$V6,sep="_"), score=region.gr$V5, strand=region.gr$V6)
  
  chrs <- intersect(unique(GenomeInfoDb::seqnames(signal.gr)), unique(GenomeInfoDb::seqnames(region.gr)))
  # length(chrs)
  signal.gr <- signal.gr[GenomeInfoDb::seqnames(signal.gr) %in% chrs,]
  region.gr <- region.gr[GenomeInfoDb::seqnames(region.gr) %in% chrs,]
  
  
  #strand info not preserved, not flipped
  mat1 = EnrichedHeatmap::normalizeToMatrix(signal.gr, region.gr, 
                                            #value_column = "cov",   
                                            # extend = c(250, 250), # gn:250,250
                                            extend = c(up, down), # tx:50,50  
                                            target_ratio = ratio, #0.3  
                                            # k = 10, # bins of target regions
                                            mean_mode = mean_mode, # "coverage",  # c("absolute", "weighted", "w0", "coverage")
                                            w = bin, # bins of extended windows 
                                            keep = c(0, 0.99), # Percentiles in the normalized matrix to keep.
                                            background = bg,  # set 0 only for occurrence, not for icSHAPE ?
                                            smooth = TRUE # set TRUE may get negative cov 
  )
  mat1 <- as.data.frame(mat1)
  # will normalize whole matrix !
  # signal <- basename(as.character(signal))
  # signal <- unlist(sapply(strsplit(signal,".",fixed = T),"[",1))
  # region <- basename(as.character(region))
  # region <- unlist(sapply(strsplit(region,".",fixed = T),"[",1))
  
  res.tmp <- as.data.frame(cbind(bed=region.gr$name,region=rep(region.label,nrow(mat1)),signal=rep(signal.label,nrow(mat1)),mat1))
  print(nrow(res.tmp))
  print(ncol(res.tmp))
  return(res.tmp)
}




# IGV plot ------
bw_theme <- theme_void() + 
  theme(plot.title = element_text(size = 20,color="black",hjust = 0.5,face="bold"),
        axis.title = element_text(size = 20,color ="black"), 
        axis.title.y = element_text(size = 20,color ="black"),  #, vjust = 1,angle = 90
        axis.line.y.left = element_blank(),, # colour, size, linetype, lineend, color
        axis.line.x.bottom = element_blank(), 
        axis.ticks = element_line(color = "black"), # ,linewidth = 1 (conficts)
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
  bg <- rtracklayer::import.bw(bw, which=GenomicRanges::GRanges(c(chr), IRanges(start = start+1, width = width))) # warning if no record found, IRanges is 1-based
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
  # bg <- bg[bg$seqnames==chr & bg$start<=start+width & bg$end>=start,] # rtracklayer::import.bw: GenomicRanges::GRanges
  
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

plotBW <- function(single_base_bw,chr,start,width,highlight="Null",ylab,xlab="",title="",color="steelblue4",fill="steelblue4",plotXbreak=TRUE,plotYbreak=TRUE,plotYtext=FALSE,annotate.size=8,maxY=-1,libsize=0){
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
  if(length(highlight)==2){
    tmp <- ggplot(single_base_bw, aes(x=start,y=(score))) +  # , group=sample
      # geom_density(stat = "identity", color="steelblue4", fill="steelblue4")+ # salmon
      geom_bar(stat="identity",color=color,fill=fill) + # position = 'dodge', "steelblue4"
      # geom_rect(aes(xmin=start,xmax=start+1,ymin=0,ymax=score))+
      # ggsci::scale_color_d3()+
      # ggsci::scale_fill_d3()+
      geom_vline(xintercept = highlight, linetype="dashed", color="grey70") +
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
  }else{
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
    }
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




# Genomic interval func. ------------------------
#as.data.frame() seem enough (need convert 0/1-based coord manually)
GRange2bed <- function(gr){
  # gr <- hg38[[1]] # 1-base to 0-base
  options(stringsAsFactors = F)
  bed <- data.frame("chr"=gr@seqnames,"start"=gr@ranges@start-1,"end"=gr@ranges@start,"name"=gr$name,"score"=gr$score,"strand"=gr@strand)
  return(bed)
}
GRange2bed12 <- function(gr){
  # gr <- hg38[[1]] # 1-base to 0-base
  options(stringsAsFactors = F)
  bed <- data.frame("chr"=gr@seqnames,"start"=gr@ranges@start-1,"end"=gr@ranges@start,"name"=gr$name,"score"=gr$score,"strand"=gr@strand,
                    "thickStart"=gr$thick@start-1, "thickEnd"=gr$thick@start+gr$thick@width-1, 
                    "itemRgb"=gr$itemRgb, 
                    "blockCount"=do.call(c,lapply(gr$blocks,length)), "blockSizes"=do.call(c,lapply(X = gr$blocks, FUN = function(x) paste(width(x),collapse=","))),"blockStarts"=do.call(c,lapply(X = gr$blocks, FUN = function(x) paste(start(x)-1,collapse=","))))
  return(bed)
}


readBedDf <- function(x){
  bed <- read.table(x,sep = "\t",header = F)
  if(ncol(bed)==6){
    message("bed6")
    colnames(bed) <- c("tx","start","end","name","score","strand")
  }else if(ncol(bed)==12){
    colnames(bed) <- c("gn","gn_start","gn_end","gn_name","gn_score","gn_strand","gn_thickStart","gn_thickEnd","gn_itemRgb","gn_blockCounts","gn_blockSizes","gn_blockStarts")
  }
  return(bed)
}



getBed12MidPosDf <- function(bed12_record){
  # Extract relevant information from bed12 dataframe
  #chromosome <- bed12_record[,1]
  #strand <- bed12_record[,6]
  block_sizes <- as.numeric(strsplit(bed12_record[,11], ",")[[1]])
  block_sizes_cumsum <- cumsum(block_sizes)
  block_sizes_cumsum <- c(0,block_sizes_cumsum)
  block_starts <- as.numeric(strsplit(bed12_record[,12], ",")[[1]])
  
  # Calculate exon boundaries
  exon_starts <- as.numeric(bed12_record[,2]) + as.numeric(block_starts)
  #exon_ends <- exon_starts + block_sizes
  
  halfSize <- round(sum(block_sizes)/2)
  midIdx <- which(block_sizes_cumsum>halfSize)[1]-1 # max(,1)
  midPos <- halfSize - block_sizes_cumsum[midIdx] + exon_starts[midIdx]
  return(round(midPos))
}

getBed12MidPosGr <- function(gr){
  # Extract relevant information from bed12 GrangeObj
  #chromosome <- gr@seqnames
  #strand <- gr@strand
  block_sizes <- as.numeric(strsplit(do.call(c,lapply(X = gr$blocks, FUN = function(x) paste(width(x),collapse=","))), ",")[[1]])
  block_sizes_cumsum <- cumsum(block_sizes)
  block_sizes_cumsum <- c(0,block_sizes_cumsum)
  block_starts <- as.numeric(strsplit(do.call(c,lapply(X = gr$blocks, FUN = function(x) paste(start(x)-1,collapse=","))), ",")[[1]])
  
  # Calculate exon boundaries
  exon_starts <- as.numeric(gr@ranges@start-1) + as.numeric(block_starts)
  #exon_ends <- exon_starts + block_sizes
  
  halfSize <- round(sum(block_sizes)/2)
  midIdx <- which(block_sizes_cumsum>halfSize)[1]-1 # max(,1)
  midPos <- halfSize - block_sizes_cumsum[midIdx] + exon_starts[midIdx]
  return(round(midPos))
}




# Biol sequence func. -----------------------
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



phastCons100way <- function(x){
  # take aver each bed12 block
  # x <- "chr1:100000-101000"
  tmp <- GenomicScores::gscores(phastCons100way.UCSC.hg38::phastCons100way.UCSC.hg38, GenomicRanges::GRanges(x)) # summaryFunFunction: mean
  return(tmp$default)
}


# Matrix process func. -----------------------
stderror <- function(x) sd(x)/sqrt(length(x))

## scale
maxmin.normalize <- function(x) {  
  #x <- sweep(x, 2, apply(x, 2, min)) 
  #x <- sweep(x, 2, apply(x, 2, max), "/")  # cannot handle zero sd
  x <- apply(x, 2, function(y) (y - min(y))/(max(y)-min(y))^as.logical(sd(y)) ) # can handle zero sd, changed as 230130
  x  # (0,1)
  #2*x - 1  # (-1,1)
}
maxmin.normalize.vec <- function(y) {  
  #x <- sweep(x, 2, apply(x, 2, min)) 
  #x <- sweep(x, 2, apply(x, 2, max), "/")  # cannot handle zero sd
  #y <- as.numeric(mat)
  x <- (y - min(y,na.rm = T))/(max(y,na.rm = T)-min(y,na.rm = T))^as.logical(sd(y,na.rm = T)) # can handle zero sd, changed as 230130
  x  # (0,1)
  #2*x - 1  # (-1,1)
}


reorder_cormat <- function(mat){
  # Use correlation between variables as distance
  #dd <- as.dist((1-mat)/2) # if cor: [-1,1]
  #dd <- as.dist(mat) # if jaccard: [0,1]
  hc <- hclust(dist(mat)) # col as records? 
  mat2 <- mat[hc$order, hc$order]
  return(list("mat"=mat2,"lab"=rownames(mat)[hc$order]))
}
# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

one.hot <- function(mat){
  t <- mat
  t_hot <- NULL
  for (col in colnames(t)) { # seq_along(t[1,])
    # for each unique value in this column (sorted so the resulting
    # columns appear in order)
    #
    for (val in sort(unique(t[, col]))) {
      t_hot <- cbind(t_hot, ifelse(t[, col] == val, 1, 0))
      # make name for this column
      #
      colnames(t_hot)[ncol(t_hot)] <- paste0(col, "_", val)
    }
  }
  rownames(t_hot) <- rownames(t)
  return(t_hot)
}


## smooth
correct.loess <- function(df){
  lo = stats::loess(signal ~ pos, data = df, span = 0.4,degree = 2,family = c("gaussian"))
  df$signal = mean(df$signal, na.rm = TRUE) * df$bc/stats::predict(lo, newdata = df[,c("signal","pos")])
  df$signal = round(df$signal, digits = 2)
  if (any(df$signal < 0, na.rm = TRUE))
    df$signal[df$signal < 0] = 0
  return(df)
}
get.loess <- function(df,span=.4){
  lo = stats::loess(signal ~ pos, data = df, span = span,degree = 2,family = c("gaussian")) # default: span=0.75, higher is smoother
  df$signal = stats::predict(lo, newdata = df$pos)
  df$signal = round(df$signal, digits = 2)
  return(df)
}

get.loess.smooth <- function(df,span=0.4){
  lo = stats::loess.smooth(x = df$pos, y=df$signal, span = span,degree = 2,family = c("gaussian")) # default: span=0.75, higher is smoother
  df <- data.frame(pos=lo$x, signal = lo$y)
  df$signal = round(df$signal, digits = 2)
  return(df)
}
# get.loess <- function(long_refp){
#   lo = stats::loess(signal ~ pos, data = long_refp, span = 0.4) # default: span=0.75, higher is smoother
#   long_refp$signal = stats::predict(lo, newdata = long_refp_tmp$pos)
#   
#   long_refp$signal = round(long_refp$signal, digits = 2)
#   if (any(long_refp$signal < 0, na.rm = TRUE))
#     long_refp$signal[long_refp$signal < 0] = 0
#   return(long_refp)
# }
get.loess.vec <- function(pos,signal,span=0.4){
  lo = stats::loess(signal ~ pos, span = span,degree = 2,family = c("gaussian")) # default: span=0.75, higher is smoother
  signal = stats::predict(lo, newdata = pos)
  signal = round(signal, digits = 2)
  # if (any(df$signal < 0, na.rm = TRUE)){
  #   
  # }
  # df$signal[df$signal < 0] = 0
  return(signal)
}
get.loess.smooth.vec <- function(pos,signal,span=0.4){
  lo = stats::loess.smooth(x = pos, y=signal, span = span,degree = 2,family = c("gaussian")) # default: span=0.75, higher is smoother
  signal = round(lo$y, digits = 2)
  # if (any(df$signal < 0, na.rm = TRUE)){
  #   
  # }
  # df$signal[df$signal < 0] = 0
  return(signal)
}
get.lowess <- function(df,f=0.4){
  #df <- long_refp_tmp[,c("pos","signal")]
  tmpList = lowess(x = df[["pos"]], y = df[["signal"]], f=f) # ,  f=2/3 , iter=3, delta=0
  tmpList <- as.data.frame(do.call(cbind,tmpList))
  df$signal <- tmpList$y
  # df$signal = round(df$signal, digits = 2)
  # if (any(df$signal < 0, na.rm = TRUE))
  #   df$signal[df$signal < 0] = 0
  return(df)
}
#loess is much slower than lowess and sometimes fails when lowess succeeds, so both programs are kept in R.
#lowess seem smoother than loess with same span/f param

#lowess is a little be 'old' in R. Most people use loess now (or in your case loess.smooth to get the same output at lowess). It should be quite a bit faster.
#but loess and loess.smooth has diff default params : family,span,degree
# loess:
#   family = c("gaussian", "symmetric"),span = 0.75, degree = 2,
# loess.smooth:
#   family = c("symmetric", "gaussian"),span = 2/3, degree = 1,


# set.seed(2020617)
# x <- 1:10000
# y <- 5000*rnorm(10000)
# # out.lowess <- lowess(x, y, f=0.3)
# # out.lowess2 <- lowess(x, y, f=0.3, iter=3, delta=0)
# out.loess <- loess(y ~ x,span=0.4,degree = 2,family = c("gaussian"))
# out.loess2 <- loess.smooth(x = x, y = y, span=0.4,degree = 2,family = c("gaussian"))
# # out.loess2 <- predict(out.loess, newdata = out.loess$x)
# plot(x=out.loess$x,y=out.loess$y)
# plot(out.loess2$x,y=out.loess2$y)
# table(out.loess2$x==x)
# # table(out.loess$x==x)
# for (i in c(10000)){
#   set.seed(2020617)
#   x <- 1:i
#   y <- 5000*rnorm(i)
#   time.start <- timeDate::timeDate()
#   out.loess <- loess(y ~ x,span=0.4,degree = 2,family = c("gaussian"))
#   time.end <- timeDate::timeDate()
#   print(paste0(i," loess:",time.end-time.start))
#   time.start <- timeDate::timeDate()
#   out.loess2 <- loess.smooth(x = x, y = y, span=0.4,degree = 2,family = c("gaussian"))
#   time.end <- timeDate::timeDate()
#   print(paste0(i," smoth:",time.end-time.start))
# }




# df <- data.frame(x=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
#                  y=c(1, 4, 7, 13, 19, 24, 20, 15, 13, 11, 15, 18, 22, 27))
# 
# loess50 <- loess(y ~ x, data=df, span=.5)
# smooth50 <- predict(loess50)
# lo2ess50 <- lowess(x=df$x, y=df$y)$y
# 
# loess90 <- loess(y ~ x, data=df, span=.9)
# smooth90 <- predict(loess90)
# lo2ess90 <- lowess(x=df$x, y=df$y, f=.9, iter=3, delta=0)$y
# 
# #create scatterplot with each regression line overlaid
# plot(df$x, df$y, pch=19, main='Loess Regression Models')
# lines(lo2ess50, x=df$x, col='red')
# lines(lo2ess90, x=df$x, col='blue')
# legend('bottomright', legend=c('.5', '.9'),
#        col=c('red', 'blue'), pch=19, title='Smoothing Span')




# DiffTable func. --------------------
## tpm
# countToTpm <- function(counts, effLen)
# {
#   rate <- log(counts) - log(effLen)
#   denom <- log(sum(exp(rate)))
#   exp(rate - denom + log(1e6))
# }
# countToFpkm <- function(counts, effLen)
# {
#   N <- sum(counts)
#   exp( log(counts) + log(1e9) - log(effLen) - log(N) )
# }

# #old version: no TMM and externalLibSize considered !
# tpm <- function(count,gene.len){ 
#   mat <- as.matrix(count)
#   matrix_tpm <- 1000*mat / gene.len
#   matrix_tpm <- t(t(matrix_tpm) * 1e6 / colSums(matrix_tpm))
#   return(matrix_tpm)
# }

#no TMM and externalLibSize considered !
fpkmToTpm <- function(fpkm)
{
  # fpkm[1:3,1:3]
  # TCGA-DM-A288-01A TCGA-QL-A97D-01A TCGA-CM-6164-01A
  # 1:        0.0000000           0.0000           0.0000
  # 2:        0.0051099           0.0000           0.0000
  # 3:        2.7694750           1.9955           3.0041
  # fpkm <- longTPM
  # fpkm <- as.data.frame(fpkm)
  # #str(fpkm)
  # (rowSums(as.matrix(is.na(fpkm))))
  # table(colSums(as.matrix(fpkm))==0)
  # exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
  # table(is.na(1e6*fpkm/colSums(fpkm)))
  # table(is.na(fpkm))
  1e6*fpkm/colSums(fpkm)
}
# 
# countToTpm <- function(counts, effLen)
# {
#   rate <- log(counts) - log(effLen)
#   denom <- log(sum(exp(rate)))
#   exp(rate - denom + log(1e6))
# }
# countToFpkm <- function(counts, effLen)
# {
#   N <- sum(counts)
#   exp( log(counts) + log(1e9) - log(effLen) - log(N) )
# }
# countToEffCounts <- function(counts, len, effLen)
# {
#   counts * (len / effLen)
# }

# countToTpm <- function(count.matrix, useExternalLibSize=FALSE, filterLowExp=FALSE, normMethod="TMM", logValue=FALSE, pseudoCount=1, gene.len){ 
# 
#     
#     return(matrix_tpm)
# }


countToCpm <- function(count.matrix, groups=FALSE, useExternalLibSize=FALSE, filterLowExp=FALSE, normMethod="TMM", logValue=FALSE, pseudoCount=1, outType, gene.len=FALSE){
  #modified run-NormCountMat.R

  # count.matrix = longMat
  # groups = as.numeric(factor(sample.table$group))
  # useExternalLibSize = F
  # filterLowExp = T
  # logValue = F
  # pseudoCount = 1
  # outType = "rpkm"
  # normMethod = "TMM"
  # gene.len = gL$geneLength
  # length(gene.len)
  # dim(count.matrix)
  if(!groups[1]){
    # message("using group info from sample prefix")
    # g <- unlist(lapply(strsplit(s,"-|_|\\."),function(x) x[1]))
    message("not using group info (all as 1)")
    g <- rep(1,ncol(count.matrix))
  } else {
    message(paste0("using group info from vector: ", paste(groups,collapse = " ") ))
    g <- groups # read.delim(groups,header = F,stringsAsFactors = F)$V1
  }
  tissue.types <- factor(g,level=unique(g))
  
  # put normal at reference level
  if (!useExternalLibSize[1]){
    message("use rowSums as lib.size in DGEList")
    y <- edgeR::DGEList(counts=count.matrix,group=tissue.types)  # def: lib.size=colSums()
  }else{
    message("use provided vector as lib.size in DGEList")
    # libSize <- read.delim(ex,header = F,stringsAsFactors = F)$V1
    y <- edgeR::DGEList(counts=count.matrix,group=tissue.types, lib.size=useExternalLibSize)
  }
  print(paste0("dim: ",dim(y$counts)))
  print(paste0("grp: ",table(tissue.types)))
  
  # remove low expressed gene
  if(filterLowExp[1]){
    message("filter low expr by group info")
    keep <- edgeR::filterByExpr(y,group = tissue.types)
    if (!useExternalLibSize[1]){
      message("libSize not kept !")
      y <- y[keep, , keep.lib.sizes=FALSE]   # gene num decrease (may set TRUE； https://support.bioconductor.org/p/9143866/)
    } else{
      message("libSize is kept !")
      y <- y[keep, , keep.lib.sizes=TRUE]   # keep.lib.sizes=FALSE, the lib.size for each sample (cf. the y$samples data.frame) will be recalculated to be the sum of the counts left
    }
    if(gene.len[1]){
      gene.len <- gene.len[keep]
    }
  }else if(!filterLowExp){
    message("do not filter low expr")
  }
  
  # calculate scaling factor for library size
  message(paste0("cal scaling factor using ",normMethod))
  y <- edgeR::calcNormFactors(y, method=normMethod) # method="RLE" etc. is also applicable
  
  if(outType=='cpm' | !gene.len[1]){
  message(paste0("gene length not provided, or outType=cpm setted, will output cpm mat"))
    y.cpm <- as.data.frame(edgeR::cpm(y,normalized.lib.sizes = TRUE,
                                      log = logValue,
                                      prior.count = pseudoCount) )
  }else if(gene.len[1]){
    y.rpkm <- as.data.frame(edgeR::rpkm(y,normalized.lib.sizes = TRUE,
                                       log = logValue,
                                       prior.count = pseudoCount,
                                       gene.length = gene.len) )
    if(outType=='rpkm'){
      message(paste0("output rpkm mat"))
      y.cpm <- y.rpkm
    }else if(outType=='tpm'){
      message(paste0("output tpm mat"))
      y.cpm <- fpkmToTpm(y.rpkm)
    }else{
      message("wrong outType, set to cpm, rpkm, or tpm")
    }
  }
  #note: In the edgeR package, RPKM values are usually calculated using the cpm function. This is because the cpm function (counts per million) is more flexible and can calculate different types of expression standardization, while the rpkm function (reads per kilobase of transcript per million reads) can only calculate RPKM values.
  #normalized.lib.sizes represents the normalized sequencing depth, while lib.size represents the raw sequencing depth.
  #cpmByGroup/rpkmByGroup will use group param and produce merged sample by group, but in cpm/rpkm group param seem not needed
  
  # y.cpm <- cbind(rownames(y.cpm),y.cpm)
  # colnames(y.cpm)[1] <- "gene_id"
    
  return(y.cpm)
}




## def limma diff func.
#with weight
limma.trend <- function(logcpm,group){
  # need convert to log2CPM?
  suppressPackageStartupMessages(library(limma))
  suppressPackageStartupMessages(library(edgeR))
  #y <- DGEList(counts=mat, samples=samples, group=group)
  #y <- calcNormFactors(y, method=norm_method)
  model <- model.matrix(~group)
  #y <- voom(y, model, plot=FALSE)
  #fit <- lmFit(y, model)
  
  ## add weight
  we <- limma::arrayWeights(logcpm, design = model, method = "genebygene")
  fit <- lmFit(logcpm, model, weights = we)
  fit <- eBayes(fit, robust=TRUE, trend=TRUE) # , trend=TRUE (limma.trend), voom比limma-trend更适用于样本库大小不一的情况
  #fit2 <- contrasts.ft(fit)
  #fit2 <- eBayes(fit2, robust=TRUE, trend=TRUE)
  #top_table <- topTable(fit2, sort.by='none', n=Inf)
  top_table <- topTable(fit, coef=2, sort.by='none', n=Inf)
  # rename columns
  mapped_names <- colnames(top_table)
  for(i in 1:ncol(top_table)){
    if(colnames(top_table)[i] == 'logFC'){
      mapped_names[i] <- 'log2FoldChange'
    }else if(colnames(top_table)[i] == 'P.Value'){
      mapped_names[i] <- 'pvalue'
    }else if(colnames(top_table)[i] == 'adj.P.Val') {
      mapped_names[i] <- 'padj'
    }else{
      mapped_names[i] <- colnames(top_table)[i]
    }
  }
  colnames(top_table) <- mapped_names
  res <- top_table
  return(res)
}

# ## limma.trend.paired deprecated !!!
# limma.trend.paired <- function(logcpm,group,patient){
#   #logcpm <- mat.norm.tmp
#   #group <- sample.table.tmp$source
#   #patient <- sample.table.tmp$cell_id
#   suppressPackageStartupMessages(library(limma))
#   suppressPackageStartupMessages(library(edgeR))
#   #y <- DGEList(counts=mat, samples=samples, group=group)
#   #y <- calcNormFactors(y, method=norm_method)
#   model <- model.matrix(~patient+group)
#   #y <- voom(y, model, plot=FALSE)
#   #fit <- lmFit(y, model)
#   fit <- lmFit(logcpm, model)
#   fit <- eBayes(fit, robust=TRUE, trend=TRUE) # , trend=TRUE (limma.trend), voom比limma-trend更适用于样本库大小不一的情况
#   #fit2 <- contrasts.ft(fit)
#   #fit2 <- eBayes(fit2, robust=TRUE, trend=TRUE)
#   #top_table <- topTable(fit2, sort.by='none', n=Inf)
#   top_table <- topTable(fit, coef=ncol(model), sort.by='none', n=Inf)
#   # rename columns
#   mapped_names <- colnames(top_table)
#   for(i in 1:ncol(top_table)){
#     if(colnames(top_table)[i] == 'logFC'){
#       mapped_names[i] <- 'log2FoldChange'
#     }else if(colnames(top_table)[i] == 'P.Value'){
#       mapped_names[i] <- 'pvalue'
#     }else if(colnames(top_table)[i] == 'adj.P.Val') {
#       mapped_names[i] <- 'padj'
#     }else{
#       mapped_names[i] <- colnames(top_table)[i]
#     }
#   }
#   colnames(top_table) <- mapped_names
#   res <- top_table
#   return(res)
# }


## define deseq2 diff func.
diff.deseq2 <- function(mat, group, method, cores){ # patient
  suppressPackageStartupMessages(library(DESeq2))
  suppressPackageStartupMessages(library(BiocParallel))
  register(MulticoreParam(cores))
  dds <- DESeqDataSetFromMatrix(countData = mat,
                                colData = as.matrix(data.frame(group=group)),
                                design = ~group)
  if(method == 'deseq2_wald'){
    dds <- DESeq(dds,test="Wald",parallel=T)
    res <- results(dds, contrast=c('group', 'positive', 'negative'))
  } else if(method == 'deseq2_lrt'){
    dds <- DESeq(dds,test="LRT",reduced = ~ 1)
    res <- results(dds, contrast=c('group', 'positive', 'negative'))
  } else if(method == 'deseq2_sc'){
    dds <- DESeq(dds,test="LRT",reduced = ~ 1,useT=TRUE, minmu=1e-6, minReplicatesForReplace=Inf)
    res <- results(dds, contrast=c('group', 'positive', 'negative'))
  }
  res <- as.data.frame(res)
  return(res)
}


## define edgeR diff func.
#paired patient
diff <- function(mat,samples,group,patient,method,norm_method, filterType="small", featureType){
  suppressPackageStartupMessages(library(edgeR))
  mat <- mat[,samples]
  print(dim(mat))
  y <- DGEList(counts=mat, samples=samples, group=group)
  
  ## filter low expr
  message(filterType)
  counts <- edgeR::getCounts(y)
  
  if (filterType=="small"){
    keep <- rowSums(counts>=1) >= 0.5*length(samples)
  }else if(filterType=="long"){
    if(featureType=="gene"){
      gene.len = as.numeric(unlist(lapply(strsplit(rownames(counts), "|", fixed = T), function(x) x[2])))
    }else if (featureType=="domain"){
      gene.len = as.numeric(unlist(lapply(strsplit(rownames(counts), "|", fixed = T), function(x) x[7]))) - as.numeric(unlist(lapply(strsplit(rownames(counts), "|", fixed = T), function(x) x[6])))
    }
    #tpm_mat <- tpm(count = counts, gene.len = gene.len)
    tpm_mat <- countToCpm(count.matrix =counts, useExternalLibSize = F, filterLowExp = F, normMethod = "none", logValue = F, pseudoCount = 1, outType = "tpm", gene.len = gene.len)
    keep <- rowSums(tpm_mat >= 1) >= 0.5*length(samples)
    
    tpm_mat <- tpm_mat[keep,]
    logtpm <- log2(tpm_mat+1)
    print(dim(logtpm))
    print(logtpm[1:2,1:2])
  }
  #keep <- filterByExpr(y,min.count = 1, min.total.count = 5, min.prop = 0.1) 
  #min.count: Minimum count required for at least some samples.
  #min.total.count: Minimum total count required.
  #min.prop: Minimum proportion of samples in the smallest group that express the gene.
  y <- y[keep,,keep.lib.sizes=TRUE]
  print(dim(y))
  
  y <- calcNormFactors(y, method=norm_method)
  
  logcpm <- edgeR::cpm(y, log=F, normalized.lib.sizes = T,prior.count = 0)
  logcpm <- log2(logcpm+1)
  
  design <- model.matrix(~ factor(patient) + factor(group))  # partial paired mode
  #design <- model.matrix(~  group + patient)  # partial paired mode, design order has great impact
  #design <- model.matrix(~ group) 
  y <- estimateDisp(y, design) 
  if(method == 'edger_glmqlf'){
    fit <- glmQLFit(y, design)
    test <- glmQLFTest(fit, coef=ncol(design))   
    res <- topTags(test, n=nrow(mat), sort.by='none')
  } else if(method == 'edger_glmlrt'){
    fit <- glmFit(y, design)
    test <- glmLRT(fit, coef=ncol(design))  # coef should be the pos col
    res <- topTags(test, n=nrow(mat), sort.by='none')
  } else if(method == 'edger_exact'){
    test <- exactTest(y)
    res <- topTags(test, n=nrow(mat), sort.by='none')
  }
  res <- cbind(res$table, baseMean=2^(res$table$logCPM))
  # rename columns
  mapped_names <- colnames(res)
  for(i in 1:ncol(res)){
    if(colnames(res)[i] == 'logFC'){
      mapped_names[i] <- 'log2FoldChange'
    }else if(colnames(res)[i] == 'PValue'){
      mapped_names[i] <- 'pvalue'
    }else if(colnames(res)[i] == 'FDR') {
      mapped_names[i] <- 'padj'
    }else{
      mapped_names[i] <- colnames(res)[i]
    }
  }
  colnames(res) <- mapped_names
  
  if (filterType=="small"){
    normMat <- logcpm
  }else if(filterType=="long"){
    normMat <- logtpm
  }
  
  outfile <- list()
  outfile[["normMat"]] <- normMat
  outfile[["diffTable"]] <- res
  return(outfile)
}

#not paired
diff.v2 <- function(mat,samples,group,method,norm_method, filterType="small",featureType){
  suppressPackageStartupMessages(library(edgeR))
  mat <- mat[,samples]
  print(dim(mat))
  y <- DGEList(counts=mat, samples=samples, group=group)
  
  ## filter low expr
  #keep <- filterByExpr(y,min.count = 1, min.total.count = 5, min.prop = 0.1)
  message(filterType)
  counts <- edgeR::getCounts(y)
  if (filterType=="small"){
    keep <- rowSums(counts>=1) >= 0.5*length(samples)
  }else if(filterType=="long"){
    if(featureType=="gene"){
      gene.len = as.numeric(unlist(lapply(strsplit(rownames(counts), "|", fixed = T), function(x) x[2])))
    }else if (featureType=="domain"){
      gene.len = as.numeric(unlist(lapply(strsplit(rownames(counts), "|", fixed = T), function(x) x[7]))) - as.numeric(unlist(lapply(strsplit(rownames(counts), "|", fixed = T), function(x) x[6])))
    }
    #tpm_mat <- tpm(count = counts, gene.len = gene.len)
    tpm_mat <- countToCpm(count.matrix =counts, useExternalLibSize = F, filterLowExp = F, normMethod = "none", logValue = F, pseudoCount = 1, outType = "tpm", gene.len = gene.len)
    keep <- rowSums(tpm_mat >= 1) >= 0.5*length(samples)
    
    tpm_mat <- tpm_mat[keep,]
    logtpm <- log2(tpm_mat+1)
    print(dim(logtpm))
    print(logtpm[1:2,1:2])
  }
  #min.count: Minimum count required for at least some samples.
  #min.total.count: Minimum total count required.
  #min.prop: Minimum proportion of samples in the smallest group that express the gene.
  y <- y[keep,,keep.lib.sizes=TRUE] # https://support.bioconductor.org/p/9136787/
  print(dim(y))
  
  y <- calcNormFactors(y, method=norm_method)
  
  logcpm <- edgeR::cpm(y, log=F, normalized.lib.sizes = T,prior.count = 0)
  logcpm <- log2(logcpm+1)
  
  #design <- model.matrix(~ patient + group)  # partial paired mode
  #design <- model.matrix(~  group + patient)  # partial paired mode, design order has great impact
  design <- model.matrix(~ factor(group))
  y <- estimateDisp(y, design)
  if(method == 'edger_glmqlf'){
    fit <- glmQLFit(y, design)
    test <- glmQLFTest(fit, coef=)
    res <- topTags(test, n=nrow(mat), sort.by='none')
  } else if(method == 'edger_glmlrt'){
    fit <- glmFit(y, design)
    test <- glmLRT(fit, coef=2) # coef should be the pos col, ncol(design)=2
    res <- topTags(test, n=nrow(mat), sort.by='none')
  } else if(method == 'edger_exact'){
    test <- exactTest(y)
    res <- topTags(test, n=nrow(mat), sort.by='none')
  }
  res <- cbind(res$table, baseMean=2^(res$table$logCPM))
  # rename columns
  mapped_names <- colnames(res)
  for(i in 1:ncol(res)){
    if(colnames(res)[i] == 'logFC'){
      mapped_names[i] <- 'log2FoldChange'
    }else if(colnames(res)[i] == 'PValue'){
      mapped_names[i] <- 'pvalue'
    }else if(colnames(res)[i] == 'FDR') {
      mapped_names[i] <- 'padj'
    }else{
      mapped_names[i] <- colnames(res)[i]
    }
  }
  colnames(res) <- mapped_names
  
  if (filterType=="small"){
    normMat <- logcpm
  }else if(filterType=="long"){
    normMat <- logtpm
  }
  outfile <- list()
  outfile[["normMat"]] <- normMat
  outfile[["diffTable"]] <- res
  return(outfile)
}



#not paired (optional with batch)
diff.v2.dcb <- function(mat,samples,group, batch="NULL", useExternalLibSize="FALSE", method, norm_method, filterType="NULL",featureType, getCountGt1SmpFrac="Y", getPercExpByGrp="Y",geneLength){
  ## test
  # filterType = "NULL"
  # featureType = "domain"
  # getCountGt1SmpFrac <- "T"
  # getPercExpByGrp = "TRUE"
  # useExternalLibSize <- libSize
  # mat <- mat
  # pos.grp <- "Y"
  # neg.grp <- "N"
  # positive_samples <- sample.table[sample.table$group==pos.grp,"sample"]  # 
  # negative_samples <- sample.table[sample.table$group==neg.grp,"sample"] #
  # samples <- c(positive_samples, negative_samples)
  # group <- c(rep("positive",length(positive_samples)),rep("negative",length(negative_samples)))
  # batch <- "NULL"
  # method <- "edger_glmlrt"
  # norm_method <- "TMM"

  suppressPackageStartupMessages(library(edgeR))
  mat <- mat[,samples]
  print(dim(mat))
  if(useExternalLibSize[1]=="FALSE"){
    message("use default rowsum libsize")
    y <- DGEList(counts=mat, samples=samples, group=group)
  }else if(is.numeric(useExternalLibSize[1])){
    message("use external libsize")
    y <- DGEList(counts=mat, samples=samples, group=group, lib.size=useExternalLibSize[samples])
  }
  
  ## filter low expr
  #keep <- filterByExpr(y,min.count = 1, min.total.count = 5, min.prop = 0.1)
  message(filterType)
  counts <- edgeR::getCounts(y)
  
  if(getCountGt1SmpFrac!="NULL"){
    ## get count>1 sample freq
    counts.neg <- counts[,group=="negative"]
    gt1.neg <- rowSums(counts.neg>=1)/length(samples[group=="negative"])
    counts.neg <- NULL
    counts.pos <- counts[,group=="positive"]
    gt1.pos <- rowSums(counts.pos>=1)/length(samples[group=="positive"])
    counts.pos <- NULL
  }
  
  if (filterType=="small"){
    keep <- rowSums(counts>=1) >= 0.5*length(samples)
  }else if(filterType=="long"){
    if (missing(geneLength)) {
      gene.len = geneLength
    }else{
      if(featureType=="gene"){
        gene.len = as.numeric(unlist(lapply(strsplit(rownames(counts), "|", fixed = T), function(x) x[2])))
      }else if (featureType=="domain"){
        gene.len = as.numeric(unlist(lapply(strsplit(rownames(counts), "|", fixed = T), function(x) x[7]))) - as.numeric(unlist(lapply(strsplit(rownames(counts), "|", fixed = T), function(x) x[6])))
      }
    }
    
    #tpm_mat <- tpm(count = counts, gene.len = gene.len)
    tpm_mat <- countToCpm(count.matrix =counts, useExternalLibSize = useExternalLibSize[samples], filterLowExp = F, normMethod = "none", logValue = F, pseudoCount = 1, outType = "tpm", gene.len = gene.len)
    
    keep <- rowSums(tpm_mat >= 1) >= 0.5*length(samples)
    
    tpm_mat <- tpm_mat[keep,]
    logtpm <- log2(tpm_mat+1)
    # print(dim(logtpm))
    # print(logtpm[1:2,1:2])
  }else{
    message("keep all, not prefilter")
    keep <- rowSums(counts) >= 0
  }
  # if (filterType=="small"){
  #   keep <- rowSums(counts>=1) >= 0.5*length(samples)
  # }else if(filterType=="long"){
  #   tpm_mat <- tpm(count = counts, 
  #                  gene.len = as.numeric(unlist(lapply(strsplit(rownames(counts), "|", fixed = T), function(x) x[2]))))
  #   keep <- rowSums(tpm_mat >= 1) >= 0.5*length(samples)
  # }else{
  #   message("keep all, not prefilter")
  #   keep <- rowSums(counts) >= 0
  # }
  #min.count: Minimum count required for at least some samples.
  #min.total.count: Minimum total count required.
  #min.prop: Minimum proportion of samples in the smallest group that express the gene.
  y <- y[keep,,keep.lib.sizes=TRUE] # https://support.bioconductor.org/p/9136787/
  print(dim(y))
  
  y <- calcNormFactors(y, method=norm_method)
  
  cpm <- edgeR::cpm(y, log=F, normalized.lib.sizes = T,prior.count = 0)
  logcpm <- log2(cpm+1)
  # message(head(row.names(logcpm),2))
  
  if(getPercExpByGrp!="NULL"){
    ## get percentile level by group
    cpm.neg <- cpm[,group=="negative"]
    #dim(cpm)
    p90.neg <- apply(cpm.neg,1, function(x) quantile(x, probs = max(2/ncol(cpm.neg),0.9) )) # 0.9; min 2 samples
    p50.neg <- apply(cpm.neg,1, function(x) quantile(x, probs = max(1/ncol(cpm.neg),0.5) )) # 0.5; min 1 sample
    p10.neg <- apply(cpm.neg,1, function(x) quantile(x, probs = max(1/ncol(cpm.neg),0.1) )) # 0.1; min 1 sample
    mean.neg <- apply(cpm.neg,1, mean ) 
    cpm.neg <- NULL
    cpm.pos <- cpm[,group=="positive"]
    #dim(cpm)
    p90.pos <- apply(cpm.pos,1, function(x) quantile(x, probs = max(2/ncol(cpm.pos),0.9) )) # 0.9; min 2 samples
    p50.pos <- apply(cpm.pos,1, function(x) quantile(x, probs = max(1/ncol(cpm.pos),0.5) )) # 0.5; min 1 sample
    p10.pos <- apply(cpm.pos,1, function(x) quantile(x, probs = max(1/ncol(cpm.pos),0.1) )) # 0.1; min 1 sample
    mean.pos <- apply(cpm.pos,1, mean ) 
    cpm.pos <- NULL
  }
  
  #design <- model.matrix(~ patient + group)  # partial paired mode
  #design <- model.matrix(~  group + patient)  # partial paired mode, design order has great impact
  if(batch[1]=="NULL"){
    design <- model.matrix(~ factor(group))
  }else if(length(batch)==length(group)){
    design <- model.matrix(~ factor(group) + batch)
  }else {
    stop("not same length of batch with group")
  }
  y <- estimateDisp(y, design)
  if(method == 'edger_glmqlf'){
    fit <- glmQLFit(y, design)
    test <- glmQLFTest(fit, coef=2)
    res <- topTags(test, n=nrow(mat), sort.by='none')
  } else if(method == 'edger_glmlrt'){
    fit <- glmFit(y, design)
    test <- glmLRT(fit, coef=2) # coef should be the pos col, ncol(design)=2
    res <- topTags(test, n=nrow(mat), sort.by='none')
  } else if(method == 'edger_exact'){
    test <- exactTest(y) # no design !
    res <- topTags(test, n=nrow(mat), sort.by='none')
  }
  
  res <- cbind(res$table, baseMean=2^(res$table$logCPM))
  if(getCountGt1SmpFrac!="NULL" & getPercExpByGrp!="NULL"){
    res <- cbind(res,negGT1RatioCount=gt1.neg,posGT1RatioCount=gt1.pos,negMeanCPM=mean.neg,negCent90CPM=p90.neg,negCent50CPM=p50.neg,negCent10CPM=p10.neg, posMeanCPM=mean.pos, posCent90CPM=p90.pos,posCent50CPM=p50.pos,posCent10CPM=p10.pos)
  }else if(getCountGt1SmpFrac=="NULL" & getPercExpByGrp!="NULL"){
    res <- cbind(res,negMeanCPM=mean.neg,negCent90CPM=p90.neg,negCent50CPM=p50.neg,negCent10CPM=p10.neg, posMeanCPM=mean.pos, posCent90CPM=p90.pos,posCent50CPM=p50.pos,posCent10CPM=p10.pos)
  }else if(getCountGt1SmpFrac!="NULL" & getPercExpByGrp=="NULL"){
    res <- cbind(res,negGT1RatioCount=gt1.neg,posGT1RatioCount=gt1.pos)
  }
    
  # rename columns
  mapped_names <- colnames(res)
  for(i in 1:ncol(res)){
    if(colnames(res)[i] == 'logFC'){
      mapped_names[i] <- 'log2FoldChange'
    }else if(colnames(res)[i] == 'PValue'){
      mapped_names[i] <- 'pvalue'
    }else if(colnames(res)[i] == 'FDR') {
      mapped_names[i] <- 'padj'
    }else{
      mapped_names[i] <- colnames(res)[i]
    }
  }
  colnames(res) <- mapped_names
  #return(res)
  
  if (filterType!="long"){
    normMat <- logcpm
  }else if(filterType=="long"){
    normMat <- logtpm
  }
  
  # message(head(row.names(logcpm),2))
  # message(head(row.names(res),2))
  
  outfile <- list()
  outfile[["normMat"]] <- normMat
  outfile[["diffTable"]] <- res
  return(outfile)
}

# test for SLE & TCGA (can be removed?)
diff.v3.dcb <- function(mat,sample.table,samples,group,form, useExternalLibSize="FALSE", method, norm_method, filterType="NULL",featureType, getCountGt1SmpFrac="Y", getPercExpByGrp="Y",geneLength){
  # ## test
  # # filterType = "NULL"
  # # featureType = "domain"
  # # getCountGt1SmpFrac <- "T"
  # # getPercExpByGrp = "TRUE"
  # sample.table <- sample.table[sample.table$cellShort=="pDC" & sample.table$group %in% c("HDA","LDA"),] # [!sample.table$patient %in% c("S200","S465","S524"),]
  # # dim(sample.table)
  # useExternalLibSize <- externalLibSize[sample.table$sample]
  # 
  # 
  # mat <- count0[,sample.table$sample]
  # pos.grp <- "HDA"
  # neg.grp <- "LDA"
  # positive_samples <- sample.table[sample.table$group==pos.grp,"sample"]  #
  # negative_samples <- sample.table[sample.table$group==neg.grp,"sample"] #
  # samples <- c(positive_samples, negative_samples)
  # group <- c(rep("positive",length(positive_samples)),rep("negative",length(negative_samples)))
  # #batch <- factor(sample.table$patient)"
  # method <- "edger_glmlrt"
  # norm_method <- "TMM"
  
  # # table(sample.table$group)
  # group <- factor(sample.table$group,levels = c("LDA","HDA") ) # ("HC","LDA","MDA","HDA")
  # Month <- factor(sample.table$Month)
  # # sex <- factor(sample.table$sex)
  # # age <- as.integer(sample.table$age)
  # # # patient <- factor(sample.table$patient)
  # # ExprDate <- factor(sample.table$ExprDate)
  # Operator <- factor(sample.table$Operator)
  

  suppressPackageStartupMessages(library(edgeR))
  # mat <- mat[,samples]
  print(dim(mat))
  if(useExternalLibSize[1]=="FALSE"){
    message("use default rowsum libsize")
    y <- DGEList(counts=mat, samples=samples, group=group)
  }else if(is.numeric(useExternalLibSize[1])){
    message("use external libsize")
    y <- DGEList(counts=mat, samples=samples, group=group, lib.size=useExternalLibSize)
  }

  ## filter low expr
  #keep <- filterByExpr(y,min.count = 1, min.total.count = 5, min.prop = 0.1)
  message(filterType)
  counts <- edgeR::getCounts(y)

  if(getCountGt1SmpFrac!="NULL"){
    ## get count>1 sample freq
    counts.neg <- counts[,group=="negative"]
    gt1.neg <- rowSums(counts.neg>=1)/length(samples[group=="negative"])
    counts.neg <- NULL
    counts.pos <- counts[,group=="positive"]
    gt1.pos <- rowSums(counts.pos>=1)/length(samples[group=="positive"])
    counts.pos <- NULL
  }

  # message(geneLength)
  # message(missing(geneLength))
  if (filterType=="small"){
    keep <- rowSums(counts>=1) >= 0.5*length(samples)
  }else if(filterType=="long"){
    if (!missing(geneLength)) {
      gene.len = geneLength
      message("gL[1]: ",geneLength[1])
      # stop("no gene length provided")
    }else{
      if(featureType=="gene"){
        gene.len = as.numeric(unlist(lapply(strsplit(rownames(counts), "|", fixed = T), function(x) x[2])))
      }else if (featureType=="domain"){
        gene.len = as.numeric(unlist(lapply(strsplit(rownames(counts), "|", fixed = T), function(x) x[7]))) - as.numeric(unlist(lapply(strsplit(rownames(counts), "|", fixed = T), function(x) x[6])))
      }
    }

    
    if(useExternalLibSize[1]=="FALSE"){
      tpm_mat <- countToCpm(count.matrix =counts, filterLowExp = F, normMethod = norm_method, logValue = F, pseudoCount = 1, outType = "tpm", gene.len = gene.len)
    }else if(is.numeric(useExternalLibSize[1])){
      #tpm_mat <- tpm(count = counts, gene.len = gene.len)
      tpm_mat <- countToCpm(count.matrix =counts, filterLowExp = F, normMethod = norm_method, logValue = F, pseudoCount = 1, outType = "tpm", gene.len = gene.len, useExternalLibSize = useExternalLibSize)
    }

    keep <- rowSums(tpm_mat >= 1) >= 0.5*length(samples)

    tpm_mat <- tpm_mat[keep,]
    logtpm <- log2(tpm_mat+1)
    # print(dim(logtpm))
    # print(logtpm[1:2,1:2])
  }else{
    message("keep all, not prefilter")
    keep <- rowSums(counts) >= 0
  }
  # if (filterType=="small"){
  #   keep <- rowSums(counts>=1) >= 0.5*length(samples)
  # }else if(filterType=="long"){
  #   tpm_mat <- tpm(count = counts,
  #                  gene.len = as.numeric(unlist(lapply(strsplit(rownames(counts), "|", fixed = T), function(x) x[2]))))
  #   keep <- rowSums(tpm_mat >= 1) >= 0.5*length(samples)
  # }else{
  #   message("keep all, not prefilter")
  #   keep <- rowSums(counts) >= 0
  # }
  #min.count: Minimum count required for at least some samples.
  #min.total.count: Minimum total count required.
  #min.prop: Minimum proportion of samples in the smallest group that express the gene.
  y <- y[keep,,keep.lib.sizes=TRUE] # https://support.bioconductor.org/p/9136787/
  print(dim(y))

  y <- calcNormFactors(y, method=norm_method)

  plotMD(cpm(y, log=TRUE,normalized.lib.sizes = T))
  
  cpm <- edgeR::cpm(y, log=F, normalized.lib.sizes = T,prior.count = 0) # default: normalized.lib.sizes=T, set to F results to slightly higher FC in high-expr gene
  logcpm <- log2(cpm+1)
  # message(head(row.names(logcpm),2))

  if(getPercExpByGrp!="NULL"){
    ## get percentile level by group
    cpm.neg <- cpm[,group=="negative"]
    #dim(cpm)
    p90.neg <- apply(cpm.neg,1, function(x) quantile(x, probs = max(2/ncol(cpm.neg),0.9) )) # 0.9; min 2 samples
    p50.neg <- apply(cpm.neg,1, function(x) quantile(x, probs = max(1/ncol(cpm.neg),0.5) )) # 0.5; min 1 sample
    p10.neg <- apply(cpm.neg,1, function(x) quantile(x, probs = max(1/ncol(cpm.neg),0.1) )) # 0.1; min 1 sample
    mean.neg <- apply(cpm.neg,1, mean )
    cpm.neg <- NULL
    cpm.pos <- cpm[,group=="positive"]
    #dim(cpm)
    p90.pos <- apply(cpm.pos,1, function(x) quantile(x, probs = max(2/ncol(cpm.pos),0.9) )) # 0.9; min 2 samples
    p50.pos <- apply(cpm.pos,1, function(x) quantile(x, probs = max(1/ncol(cpm.pos),0.5) )) # 0.5; min 1 sample
    p10.pos <- apply(cpm.pos,1, function(x) quantile(x, probs = max(1/ncol(cpm.pos),0.1) )) # 0.1; min 1 sample
    mean.pos <- apply(cpm.pos,1, mean )
    cpm.pos <- NULL
  }

  sample.table <- sample.table[sample.table$sample %in% samples,]
  rownames(sample.table) <- sample.table$sample
  sample.table <- sample.table[samples,]
  design <- model.matrix(form,data = sample.table)
  y <- estimateDisp(y, design)

  # md2 <- plotMD(cpm(y, log=TRUE), column=1) # no diff with md1

  if(method == 'edger_glmqlf'){
    fit <- glmQLFit(y, design)
    test <- glmQLFTest(fit, coef=2)

    # test <- glmQLFTest(fit, contrast=c(0,-1,1) )# 3/2
    res <- topTags(test, n=nrow(mat), sort.by='none')
  } else if(method == 'edger_glmlrt'){
    fit <- glmFit(y, design)
    test <- glmLRT(fit, coef=2) # coef should be the pos col, ncol(design)=2
    res <- topTags(test, n=nrow(mat), sort.by='none')
  } else if(method == 'edger_exact'){
    test <- exactTest(y) # no design !
    res <- topTags(test, n=nrow(mat), sort.by='none')
  }

  res <- cbind(res$table, baseMean=2^(res$table$logCPM))
  if(getCountGt1SmpFrac!="NULL" & getPercExpByGrp!="NULL"){
    res <- cbind(res,negGT1RatioCount=gt1.neg,posGT1RatioCount=gt1.pos,negMeanCPM=mean.neg,negCent90CPM=p90.neg,negCent50CPM=p50.neg,negCent10CPM=p10.neg, posMeanCPM=mean.pos, posCent90CPM=p90.pos,posCent50CPM=p50.pos,posCent10CPM=p10.pos)
  }else if(getCountGt1SmpFrac=="NULL" & getPercExpByGrp!="NULL"){
    res <- cbind(res,negMeanCPM=mean.neg,negCent90CPM=p90.neg,negCent50CPM=p50.neg,negCent10CPM=p10.neg, posMeanCPM=mean.pos, posCent90CPM=p90.pos,posCent50CPM=p50.pos,posCent10CPM=p10.pos)
  }else if(getCountGt1SmpFrac!="NULL" & getPercExpByGrp=="NULL"){
    res <- cbind(res,negGT1RatioCount=gt1.neg,posGT1RatioCount=gt1.pos)
  }

  # rename columns
  mapped_names <- colnames(res)
  for(i in 1:ncol(res)){
    if(colnames(res)[i] == 'logFC'){
      mapped_names[i] <- 'log2FoldChange'
    }else if(colnames(res)[i] == 'PValue'){
      mapped_names[i] <- 'pvalue'
    }else if(colnames(res)[i] == 'FDR') {
      mapped_names[i] <- 'padj'
    }else{
      mapped_names[i] <- colnames(res)[i]
    }
  }
  colnames(res) <- mapped_names
  #return(res)

  if (filterType!="long"){
    normMat <- logcpm
  }else if(filterType=="long"){
    normMat <- logtpm
  }

  # message(head(row.names(logcpm),2))
  # message(head(row.names(res),2))

  outfile <- list()
  outfile[["normMat"]] <- normMat
  outfile[["diffTable"]] <- res
  # outfile[["md1"]] <- (md1) # plot: error
  # outfile[["md2"]] <- draw(md2)
  # outfile[["md2.2"]] <- print(md2)
  return(outfile)
}


## diff wilcox paired 
wilcox.paired <- function(matrix_cpm,group,paired=TRUE){
  # need convert to log2CPM?
  suppressPackageStartupMessages(library(edgeR))
  test_func <- function(x){
    wilcox.test(as.numeric(x[group == "negative"]), as.numeric(x[group == "positive"]), alternative='two.sided', paired = paired)$p.value
  }
  #matrix_cpm <- cpm(mat)
  #table(group == "positive")
  #wilcox.test(as.numeric(matrix_cpm[1,group == "negative"]), as.numeric(matrix_cpm[1,group == "positive"]), alternative='two.sided', paired = paired)$p.value
  #matrix_cpm
  pvalues <- apply(matrix_cpm, 1, test_func)
  treatMeans <- apply(matrix_cpm[,which(group == "positive")], 1, mean)
  ctrlMeans <- apply(matrix_cpm[,which(group == "negative")], 1, mean)
  logFC <- treatMeans - ctrlMeans
  res <- data.frame(log2FoldChange=logFC,
                    pvalue=pvalues,
                    padj=p.adjust(pvalues, method='BH'),
                    baseMean=apply(matrix_cpm, 1, mean),
                    treatMean=treatMeans,
                    ctrlMean=ctrlMeans)
  return(res)
  # outfile <- list()
  # outfile[["normMat"]] <- normMat
  # outfile[["diffTable"]] <- res
  # return(outfile)
}


## gfold
prep.gfold <- function(mat,featureType,outDir){
  # ## test
  # mat <- count0[,samples]
  # featureType <- "domain"
  # outDir <- paste0("./output/",dst,"/gfold")
  
  library(tidyverse)
  expr_count <- mat
  expr_count$feature <- rownames(expr_count)
  if(featureType=="domain"){
    expr_count$Length <- as.numeric(unlist(sapply(strsplit(expr_count$feature,"|",fixed = T), "[", 7))) - as.numeric(unlist(sapply(strsplit(expr_count$feature,"|",fixed = T), "[", 6)))
  }else if(featureType=="gene"){
    expr_count$Length <- as.numeric(unlist(sapply(strsplit(expr_count$feature,"|",fixed = T), "[", 2))) 
  }
  #summary(expr_count$Length)
  #rownames(expr_count) <- expr_count$feature
  
  dir.create(outDir,showWarnings = F,recursive = T)
  for(i in colnames(expr_count)[!colnames(expr_count) %in% c("gene_id","feature","Length")]){
    #i = cn[1]
    outFile <- paste0(outDir,"/",i,".read_cnt")
    if(file.exists(outFile)){
      message(outFile," exist, skip")
      next
    }
    dat <- data.frame(GeneSymbol=rownames(expr_count),
                      GeneName=NA,
                      `Read Count`=expr_count[,i],
                      `Gene exon length`=expr_count[,"Length"],
                      RPKM=NA
    )
    # print(head(dat))
    data.table::fwrite(dat,file = outFile,quote = F,sep = "\t",row.names=F, col.names=F)
  }
}
# # Linux code
# # A_vs_B   注意：s1和s2 为 s2 vs s1
# source activate cfDNA_base
# dst="lulab_oscc_tissue_diff"
# cd /BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/$dst/gfold
# for smp in `cat s2.csv | tr "," " "`
# do 
#   echo "start $smp at `date`"
#   gfold diff -s1 `cat s1.csv` -s2 $smp -sc 0.01 -suf .read_cnt -o ${smp}.diff > ${smp}.log 2>&1
# done
# # gfold diff -s1 BLD,CH -s2 CXS,LSY -sc 0.01 -suf .read_cnt -o test.diff

importGfold <- function(x){
  # x <- "./output/lulab_oscc_tissue_diff/gfold/CH.diff"
  if(!file.exists(x)){
    message(x," not exist")
    return(NULL)
  }
  
  dat <- data.table::fread(cmd=paste0("grep -v '^#' ", x), data.table = F, sep = "\t")
  dat$sample <- x
  return(dat)
}


# significant candidate selection
filterDiffTabByP <- function(diff,pvalue,padj,abs.log2FC){
  return(diff[diff$pvalue <= pvalue & diff$padj <= padj & abs(diff$log2FoldChange)>=abs.log2FC,])
}

filterDiffTabByTopN <- function(diff,N,log2FCcol="log2FoldChange",rankBy="pvalue",countBy="total"){
  # rankBy: padj, pvalue, log2FoldChange
  # countBy: total, trend
  if(rankBy %in% c("padj","pvalue","pval")){
    rank.trend <- F
  }else if(rankBy %in% log2FCcol){
    rank.trend <- T
    rankBy <- paste0("abs.",log2FCcol)
    diff[[rankBy]] <- abs(diff[[log2FCcol]])
  }
  
  diff <- diff[order(diff[[rankBy]],decreasing = rank.trend),]
  
  if(countBy=="total"){
    tmp <- diff[1:min(N,nrow(diff)),]
  }else if(countBy=="trend"){
    diff.up <- diff[diff[[log2FCcol]]>0,]
    diff.dw <- diff[diff[[log2FCcol]]<0,]
    tmp <- as.data.frame(rbind(diff.up[1:min(nrow(diff.up),as.integer(N*0.5)),],diff.dw[1:min(nrow(diff.up),as.integer(N*0.5)),]))
  }
  return(tmp)
}




# Ggplot2 func. --------------------
library(ggplot2)

## theme 
bar_theme <- theme_minimal() +  # base_size=12
  theme(#axis.ticks.x=element_blank(),
    #strip.text.y = element_blank(),
    aspect.ratio = 1,
    strip.text = element_text(size=20),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size = 20,hjust = 1,vjust = 0.5,angle = 90), #
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size=20),
    # strip.text = element_blank(),
    legend.position = "none", #c(0.9,0.8),#,#
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16))
my_theme_point <-   theme(plot.title = element_text(size = 28,color="black",hjust = 0.5),
                    axis.title = element_text(size = 28,color ="black"),
                    axis.text = element_text(size= 24,color = "black"), #,face="bold
                    panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),
                    axis.text.x = element_text( hjust = 0.5 ), # angle = 45,
                    axis.text.y = element_text( hjust = 0 ), # angle = 45,
                    panel.grid=element_blank(),
                    legend.position = "right",#c(0.5,0.3),
                    legend.text = element_text(size= 26,color = "black"),
                    legend.title= element_text(size= 26,color = "black"))
my_theme_box <- theme(aspect.ratio = 1.5,
                      plot.title = element_text(size = 24,color="black",hjust = 0.5),
                      axis.title = element_text(size = 24,color ="black"), 
                      axis.text = element_text(size= 24,color = "black"),
                      axis.text.x = element_text(size= 24,color = "black"),
                      #panel.grid=element_blank(),
                      panel.grid.major.x=element_blank(),
                      panel.grid.minor.x = element_blank(),
                      panel.grid.major.y = element_blank(), #element_line(color = "grey50",linetype = "dashed"), #size= 1,
                      panel.grid.minor.y = element_blank(),
                      #panel.grid.minor.y = element_blank(),
                      panel.border = element_blank(),
                      legend.position = "right",#c(.25,.6),
                      legend.text = element_text(size= 24),
                      legend.title= element_text(size= 24),
                      strip.text.y = element_blank(),
                      strip.text.x = element_text(size=24)
)
my_theme_roc <- theme(aspect.ratio = 1,
                      plot.title = element_text(size = 24,color="black",hjust = 0.5),
                      axis.title = element_text(size = 24,color ="black"), 
                      axis.text = element_text(size= 24,color = "black"),
                      axis.text.x = element_text(size= 24,color = "black"),
                      #panel.grid=element_blank(),
                      # panel.grid.major.x=element_blank(),
                      # panel.grid.minor.x = element_blank(),
                      # panel.grid.major.y = element_blank(), #element_line(color = "grey50",linetype = "dashed"), #size= 1,
                      # panel.grid.minor.y = element_blank(),
                      #panel.grid.minor.y = element_blank(),
                      panel.border = element_blank(),
                      legend.position = "right",#c(.25,.6),
                      legend.text = element_text(size= 16),
                      legend.title= element_text(size= 20),
                      strip.text.y = element_blank(),
                      strip.text.x = element_text(size=24)
)
my_theme_pie <- theme(
  aspect.ratio = 1,
  plot.title = element_text(size = 24,color="black",hjust = 0.5),
  axis.title = element_text(size = 24,color ="black"), 
  axis.ticks = element_blank(),
  panel.grid=element_blank(),
  # panel.grid.major.x=element_blank(),
  # panel.grid.minor.x = element_blank(),
  # panel.grid.major.y = element_line(color = "grey50",linetype = "dashed"), #size= 1,
  #panel.grid.minor.y = element_blank(),
  panel.border = element_blank(),
  axis.text = element_blank(), #element_text(size= 20,color = "black"),
  # axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5), #  , color = c(rep("#003366",length(rna)),rep("darkred",length(dna)))
  legend.position = "right",#c(.25,.6),
  legend.text = element_text(size= 16),
  legend.title= element_text(size= 24)
)
my_theme_volcano <-   theme(aspect.ratio=1,
                            plot.title = element_text(size = 28,color="black",hjust = 0.5),
                            axis.title = element_text(size = 28,color ="black"), 
                            axis.text = element_text(size= 24,color = "black"), #,face="bold
                            panel.grid.minor.y = element_blank(),
                            panel.grid.minor.x = element_blank(),
                            axis.text.x = element_text( hjust = 0.5 ), # angle = 45,
                            axis.text.y = element_text( hjust = 0 ), # angle = 45,
                            panel.grid=element_blank(),
                            legend.position = "right",#c(0.5,0.3),
                            legend.text = element_text(size= 24,color = "black"),
                            legend.title= element_text(size= 26,color = "black")
)

plotViolin <- function(logcpm.sum,sig.size=10,colors=c("grey50","firebrick")){
  p1 <- ggplot(logcpm.sum, aes(x=group,y=value,fill=group))+ # 
    geom_violin(alpha=0.7)+
    geom_boxplot(width=0.1)+
    scale_fill_manual(name="Group",values = colors)+#colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(unique(logcpm.sum$group))))+)
    # facet_grid(.~name,scales = "free")+
    ylab("Scaled peak-index") +
    ggpubr::stat_compare_means(
      label.x.npc = "middle", label.y.npc = "top",
      #size =12,step.increase = 0.08,#paired = TRUE,
      aes(x=group,y=value,
          label = ..p.format..), # p.signif, p.format, group = Group,
      symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01,1),
                       symbols = c( "***", "**", "*", "ns")),
      method.args = list(alternative="greater"),
      ref.group = NC.label,hide.ns=F, size=sig.size,
      method = "wilcox.test"
    ) +
    theme_bw() + 
    my_theme_box
  return(p1)
}
plotROC <- function(logcpm.sum,direction="auto"){
  library(ggplot2)
  library(pROC)
  logcpm.sum <- as.data.frame(logcpm.sum)
  rocobj <- roc(group ~ value, data=logcpm.sum, direction=direction)
  ?roc
  auc <- round(rocobj$auc,4)
  p1 <- ggroc(rocobj, colour= 'steelblue', size = 2) + #  = 'steelblue',
    ggtitle(paste0('', '(AUC = ', auc, ')')) +
    theme_minimal() + 
    my_theme_roc
  return(p1)
}
plotMultiROC <- function(logcpm.sum.list,direction="auto",RNA_colors=RNA_colors){
  library(ggplot2)
  library(pROC)
  library(purrr)
  #library
  library(pROC)
  library(ggplot2)
  library(tidyverse)
  
  logcpm.sum <- as.data.frame(do.call(rbind,logcpm.sum.list))
  #table(logcpm.sum$peak.precursor)
  logcpm.sum$sum.log2cpm <- NULL
  logcpm.sum <- reshape2::dcast(data = logcpm.sum, formula = sample+group~peak.precursor, id.var=c("sample","group"), value.var = "value")
  # logcpm.sum[1:3,]
  # logcpm.sum <- as.data.frame(logcpm.sum)
  
  rocobj.list <- list()
  # auc.list <- list()
  # rocobj <- roc(group ~ sum.log2cpm.scale + mir.sum.log2cpm.scale, data=logcpm.sum)
  for(i in names(logcpm.sum.list)){
    print(i)
    rocobj.list[[i]] <- roc(logcpm.sum[["group"]], logcpm.sum[[i]],direction=direction ) # , data=logcpm.sum
    print(rocobj.list[[i]]$auc)
    # auc.list[[i]] <- round(auc(logcpm.sum[["group"]], logcpm.sum[[i]]),4) # , data=logcpm.sum
  }
  # # extract auc
  # names(rocobj.list)
  # rocobj.list %>% 
  #   map(~tibble(AUC = .x$auc)) %>% 
  #   bind_rows(.id = "name") -> data.auc
  
  
  # generate labels labels
  tibble(name=names(rocobj.list), AUC=as.numeric(do.call(c,lapply(rocobj.list,FUN = function(y) {return(as.numeric(y$auc)) } )) ) ) %>% 
    mutate(label_long=paste0(name,", AUC=",paste(round(AUC,2))),
           label_AUC=paste0("AUC=",paste(round(AUC,2)))) -> data.labels
  
  # # plot a facet plot with AUC within plots
  # ggroc(roc.list) +
  #   facet_wrap(~name) +
  #   
  #   geom_text(data = data.labels,
  #             aes(0.5, 1, 
  #                 label = paste(label_AUC)),
  #             hjust = 1) 
  
  # plot on a single plot with AUC in labels
  p1 <- ggroc(rocobj.list, aes="colour", size = 2, alpha=0.9) + #  = 'steelblue',
    # ggtitle(paste0('', '(AUC = ', auc, ')')) +
    scale_color_manual(values = RNA_colors, labels=data.labels$label_long) +
    guides(fill = guide_legend(title = "Precursor" )) +
    theme_minimal() + 
    my_theme_roc
  
  return(p1)
}
#plotROC(logcpm.sum.list[["all"]])
plotConfusionMat<-function(Actual,Predict,colors=c("white","red4","dodgerblue3"),text.scl=5){
  actual = as.data.frame(table(Actual))
  names(actual) = c("Actual","ActualFreq")
  
  #build confusion matrix
  confusion = as.data.frame(table(Actual, Predict))
  names(confusion) = c("Actual","Predicted","Freq")
  
  #calculate percentage of test cases based on actual frequency
  
  confusion = merge(confusion, actual, by=c('Actual','Actual'))
  confusion$Percent = confusion$Freq/confusion$ActualFreq*100
  confusion$ColorScale<-confusion$Percent*-1
  confusion[which(confusion$Actual==confusion$Predicted),]$ColorScale<-confusion[which(confusion$Actual==confusion$Predicted),]$ColorScale*-1
  # confusion$Label<-paste(round(confusion$Percent,0),"%, n=",confusion$Freq,sep="")
  confusion$Label<-paste(confusion$Freq,sep="")
  confusion$Label[confusion$Label==0] <- ""
  tile <- ggplot() +
    geom_tile(aes(x=Actual, y=Predicted,fill=ColorScale),data=confusion, color="black",size=0.1) +
    labs(x="Actual",y="Predicted")
  tile = tile +
    geom_text(aes(x=Actual,y=Predicted, label=Label),data=confusion, size=text.scl, colour="black") +
    scale_fill_gradient2(low=colors[2],high=colors[3],mid=colors[1],midpoint = 0,guide='none') +
    theme_bw() +
    theme(
      aspect.ratio = 1,
      plot.title = element_text(size = 24,color="black",hjust = 0.5),
      axis.title = element_text(size = 24,color ="black"), 
      axis.text = element_text(size= 24,color = "black"),
      axis.text.x = element_text(size= 24,color = "black"),
      #panel.grid=element_blank(),
      panel.grid.major.x = element_line(color = "grey30",linetype = "dashed"),,
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey30",linetype = "dashed"), #size= 1,
      panel.grid.minor.y = element_blank(),
      #panel.grid.minor.y = element_blank(),
      panel.border = element_blank(),
      legend.position = "right",#c(.25,.6),
      legend.text = element_text(size= 24),
      legend.title= element_text(size= 24),
      strip.text.y = element_blank(),
      strip.text.x = element_text(size=24)
    )
}


plotPie <- function(x){
  #x <- c(1,2,2,3,3,4) 
  # x <- df$RNA
  df2 <- as.data.frame(table(x))
  df2 <- df2[df2$Freq>0,]
  df2 <- dplyr::as_tibble(df2) %>% 
    dplyr::group_by(x) %>% 
    dplyr::summarize(Freq=mean(Freq))
  df2$lab <- round(df2$Freq/sum(df2$Freq),digits = 3)
  df2 <- df2 %>%
    dplyr::mutate(csum = rev(cumsum(rev(Freq))),
                  pos = Freq/2 + lead(csum, 1),
                  pos = if_else(is.na(pos), Freq/2, pos))
  
  RNA.col <- data.frame(RNA=c(rna,dna),col=c(pal_nejm_adaptive()(15)[1:14],"#11838D"))
  RNA.col$RNA <- factor(RNA.col$RNA,levels = c(rna,dna))
  
  # str(df2)
  p <- ggplot(df2, aes(x="", y=Freq, fill=x)) +
    geom_bar(stat="identity", width=1, color="black") +
    coord_polar("y", start=0) +
    # geom_text(aes(label = lab), position = position_stack(vjust=0.5), size = 5, colour="white") +
    ggrepel::geom_label_repel(data = df2, col="white",
                              aes(y = pos, label = paste0(100*lab, "%")),
                              size = 8, nudge_x = 0.75, show.legend = FALSE) +
    labs(x = NULL, y = NULL, title ="") + # paste0("Total repeats peak: ",nrow(peak))
    scale_color_manual("black") +
    scale_fill_manual(name="precursor",values = RNA.col$col[RNA.col$RNA %in% df2$x] ) + 
    # scale_fill_nejm_adaptive(alpha = 0.8) +
    theme_bw() + 
    my_theme_pie
  return(p)
}





# auto color num ggsci
library("ggsci")

#' Adaptive palette (discrete).
#'
#' Create a discrete palette which will use the first n colors from
#' the supplied color values, and interpolate after n.
adaptive_pal <- function(values) {
  force(values)
  function(n = 10) {
    if (n <= length(values)) {
      values[seq_len(n)]
    } else {
      colorRampPalette(values, alpha = TRUE)(n)
    }
  }
}

## npg
pal_npg_adaptive <- function(palette = c("nrc"), alpha = 1) {
  palette <- match.arg(palette)
  
  if (alpha > 1L | alpha <= 0L) stop("alpha must be in (0, 1]")
  
  raw_cols <- ggsci:::ggsci_db$"npg"[[palette]]
  raw_cols_rgb <- col2rgb(raw_cols)
  alpha_cols <- rgb(
    raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], raw_cols_rgb[3L, ],
    alpha = alpha * 255L, names = names(raw_cols),
    maxColorValue = 255L
  )
  
  adaptive_pal(unname(alpha_cols))
}

scale_color_npg_adaptive <- function(palette = c("nrc"), alpha = 1, ...) {
  palette <- match.arg(palette)
  discrete_scale("colour", "npg", pal_npg_adaptive(palette, alpha), ...)
}

scale_fill_npg_adaptive <- function(palette = c("nrc"), alpha = 1, ...) {
  palette <- match.arg(palette)
  discrete_scale("fill", "npg", pal_npg_adaptive(palette, alpha), ...)
}

##nejm
pal_nejm_adaptive <- function(palette = c("default"), alpha = 1) {
  palette <- match.arg(palette)
  
  if (alpha > 1L | alpha <= 0L) stop("alpha must be in (0, 1]")
  
  raw_cols <- ggsci:::ggsci_db$"nejm"[[palette]]
  raw_cols_rgb <- col2rgb(raw_cols)
  alpha_cols <- rgb(
    raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], raw_cols_rgb[3L, ],
    alpha = alpha * 255L, names = names(raw_cols),
    maxColorValue = 255L
  )
  
  adaptive_pal(unname(alpha_cols))
}
scale_color_nejm_adaptive <- function(palette = c("default"), alpha = 1, ...) {
  palette <- match.arg(palette)
  discrete_scale("colour", "nejm", pal_nejm_adaptive(palette, alpha), ...)
}
scale_fill_nejm_adaptive <- function(palette = c("default"), alpha = 1, ...) {
  palette <- match.arg(palette)
  discrete_scale("fill", "nejm", pal_nejm_adaptive(palette, alpha), ...)
}

##d3
#ggsci:::ggsci_db$"d3" # check all palette
pal_d3_adaptive <- function(palette = c("category10"), alpha = 1) {
  palette <- match.arg(palette)
  
  if (alpha > 1L | alpha <= 0L) stop("alpha must be in (0, 1]")
  
  raw_cols <- ggsci:::ggsci_db$"d3"[[palette]]
  raw_cols_rgb <- col2rgb(raw_cols)
  alpha_cols <- rgb(
    raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], raw_cols_rgb[3L, ],
    alpha = alpha * 255L, names = names(raw_cols),
    maxColorValue = 255L
  )
  
  adaptive_pal(unname(alpha_cols))
}
scale_color_d3_adaptive <- function(palette = c("category10"), alpha = 1, ...) {
  palette <- match.arg(palette)
  discrete_scale("colour", "d3", pal_d3_adaptive(palette, alpha), ...)
}
scale_fill_d3_adaptive <- function(palette = c("category10"), alpha = 1, ...) {
  palette <- match.arg(palette)
  discrete_scale("fill", "d3", pal_d3_adaptive(palette, alpha), ...)
}

##lancet
pal_lancet_adaptive <- function(palette = c("lanonc"), alpha = 1) {
  palette <- match.arg(palette)
  
  if (alpha > 1L | alpha <= 0L) stop("alpha must be in (0, 1]")
  
  raw_cols <- ggsci:::ggsci_db$"lancet"[[palette]]
  raw_cols_rgb <- col2rgb(raw_cols)
  alpha_cols <- rgb(
    raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], raw_cols_rgb[3L, ],
    alpha = alpha * 255L, names = names(raw_cols),
    maxColorValue = 255L
  )
  
  adaptive_pal(unname(alpha_cols))
}
scale_color_lancet_adaptive <- function(palette = c("lanonc"), alpha = 1, ...) {
  palette <- match.arg(palette)
  discrete_scale("colour", "lancet", pal_lancet_adaptive(palette, alpha), ...)
}
scale_fill_lancet_adaptive <- function(palette = c("lanonc"), alpha = 1, ...) {
  palette <- match.arg(palette)
  discrete_scale("fill", "lancet", pal_lancet_adaptive(palette, alpha), ...)
}

##igv
#ggsci:::ggsci_db$igv # jama # gsea
pal_igv_adaptive <- function(palette = c("default"), alpha = 1) {
  palette <- match.arg(palette)
  
  if (alpha > 1L | alpha <= 0L) stop("alpha must be in (0, 1]")
  
  raw_cols <- ggsci:::ggsci_db$"igv"[[palette]]
  raw_cols_rgb <- col2rgb(raw_cols)
  alpha_cols <- rgb(
    raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], raw_cols_rgb[3L, ],
    alpha = alpha * 255L, names = names(raw_cols),
    maxColorValue = 255L
  )
  
  adaptive_pal(unname(alpha_cols))
}
scale_color_igv_adaptive <- function(palette = c("default"), alpha = 1, ...) {
  palette <- match.arg(palette)
  discrete_scale("colour", "igv", pal_igv_adaptive(palette, alpha), ...)
}
scale_fill_igv_adaptive <- function(palette = c("default"), alpha = 1, ...) {
  palette <- match.arg(palette)
  discrete_scale("fill", "igv", pal_igv_adaptive(palette, alpha), ...)
}
#pal_igv_adaptive()(4)


##jama
#ggsci:::ggsci_db$jama # jama # gsea
#ggsci:::ggsci_db$jama # jama # gsea
pal_jama_adaptive <- function(palette = c("default"), alpha = 1) {
  palette <- match.arg(palette)
  
  if (alpha > 1L | alpha <= 0L) stop("alpha must be in (0, 1]")
  
  raw_cols <- ggsci:::ggsci_db$"jama"[[palette]]
  raw_cols_rgb <- col2rgb(raw_cols)
  alpha_cols <- rgb(
    raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], raw_cols_rgb[3L, ],
    alpha = alpha * 255L, names = names(raw_cols),
    maxColorValue = 255L
  )
  
  adaptive_pal(unname(alpha_cols))
}
scale_color_jama_adaptive <- function(palette = c("default"), alpha = 1, ...) {
  palette <- match.arg(palette)
  discrete_scale("colour", "jama", pal_jama_adaptive(palette, alpha), ...)
}
scale_fill_jama_adaptive <- function(palette = c("default"), alpha = 1, ...) {
  palette <- match.arg(palette)
  discrete_scale("fill", "jama", pal_jama_adaptive(palette, alpha), ...)
}
pal_jama_adaptive()(4)

##gsea
#ggsci:::ggsci_db$gsea # gsea # gsea
#ggsci:::ggsci_db$gsea # gsea # gsea
pal_gsea_adaptive <- function(palette = c("default"), alpha = 1) {
  palette <- match.arg(palette)
  
  if (alpha > 1L | alpha <= 0L) stop("alpha must be in (0, 1]")
  
  raw_cols <- ggsci:::ggsci_db$"gsea"[[palette]]
  raw_cols_rgb <- col2rgb(raw_cols)
  alpha_cols <- rgb(
    raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], raw_cols_rgb[3L, ],
    alpha = alpha * 255L, names = names(raw_cols),
    maxColorValue = 255L
  )
  
  adaptive_pal(unname(alpha_cols))
}
scale_color_gsea_adaptive <- function(palette = c("default"), alpha = 1, ...) {
  palette <- match.arg(palette)
  discrete_scale("colour", "gsea", pal_gsea_adaptive(palette, alpha), ...)
}
scale_fill_gsea_adaptive <- function(palette = c("default"), alpha = 1, ...) {
  palette <- match.arg(palette)
  discrete_scale("fill", "gsea", pal_gsea_adaptive(palette, alpha), ...)
}
#pal_gsea_adaptive()(4)



# Pathway func. --------------------
suppressPackageStartupMessages(library(argparse))
#remotes::install_github('YuLab-SMU/ggtree')
#BiocManager::install("enrichplot")
#BiocManager::install("clusterProfiler")
suppressPackageStartupMessages(library(enrichplot))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(DOSE))
#BiocManager::install("org.Hs.eg.db")
suppressPackageStartupMessages(library(org.Hs.eg.db))


shortPathDescription <- function(x){
  #KEGG
  x <- sub(" pathway","",x) 
  x <- sub(", "," ",x) 
  x <- sub(" and "," ",x) 
  x <- sub(" - "," ",x) 
  x <- sub(" / "," ",x) 
  x <- gsub("\\s*\\([^\\)]+\\)","",x)  # remove bracket
  #Reactome
  x <- sub("_PATHWAY","",x) 
  return(x)
}



convertGmt <- function(x,source="Reactome"){
  #install.packages("GSA")
  library(GSA)
  gmt_file <-GSA.read.gmt(x)
  
  gmt_file$geneset.names.short <- shortPathDescription(gmt_file$geneset.names)
  
  df.list <- list()
  
  if(sum(grepl("ENSG",gmt_file[[1]][[1]])) >= 1){
    message("treat gene as ENSG...")
    lab <- "ENSG"
  }else{
    message("treat gene as Symbol...")
    lab <- "hgnc_symbol"
  }
  for(i in 1:length(gmt_file[[1]])){
    # i <- 1
    gmt_file_df <- data.frame(n=1:length(gmt_file[[1]][[i]]), ensembl_gene_id="", ENTREZID="", KEGGID=gmt_file$geneset.descriptions[i], DESCRPTION=gmt_file$geneset.names[i], lab=gmt_file[[1]][[i]], short_DESCRPTION=gmt_file$geneset.names.short[i], source=source)
    colnames(gmt_file_df)[colnames(gmt_file_df)=="lab"] <- lab
    df.list[[i]] <- gmt_file_df
  }
  df <- as.data.frame(do.call(rbind,df.list))
  df <- df[df[[lab]]!="",]
  return(df)
}


# 2.Gene Set Enrichment Analysis
gsea.kegg.offline <- function(gene.list,m_t,qValue=0.2, pValue=0.1, adjustPval="BH"){
  # message("start gsea...")
  # gene.list <- rank.list
  # m_t <- m_t2g
  gene.list <- gene.list[!duplicated(names(gene.list))]
  tryCatch(
    {
      set.seed(1234) # important!
      GSEA.KEGG <- GSEA(gene.list, 
                       TERM2GENE = m_t, 
                       nPerm  = 500,
                       minGSSize = 10,
                       maxGSSize = 500,
                       pvalueCutoff = pValue,
                       pAdjustMethod = adjustPval
      )
      GSEA.KEGG <- GSEA.KEGG@result
      if(nrow(GSEA.KEGG)>0){
        GSEA.KEGG$trend <- "all"
        # GSEA.KEGG$minus.log10.pvalue <- -log10(GSEA.KEGG$pvalue)
        # GSEA.KEGG$minus.log10.padj <- -log10(GSEA.KEGG$p.adjust)
        # GSEA.KEGG$minus.log10.qvalues <- -log10(GSEA.KEGG$qvalues)
      }
      return(as.data.frame(GSEA.KEGG))
    },error = function(e) {
      message("error in gsea !")
      return(data.frame(ID="",Description="",setSize="",enrichmentScore="",NES="",pvalue="",p.adjust="",qvalues="",rank="",leading_edge="",core_enrichment="",trend="")) # ,minus.log10.pvalue="",minus.log10.padj="",minus.log10.qvalues=""
    }
    # ,warning = function(w) {
    #   print("warning!") #printing a message when prompted with a warning
    # }
    # , finally ={
    #   print("code executed") #printing a message after code execution
    # }
  )
}

ora.kegg.offline <- function(gene.top.up,gene.top.down,gene.list="Null",m_t,qValue=0.2, pValue=0.1, adjustPval="BH"){
  # 1.Over-Representation Analysis (ncbiid or otherId in TERM2GENE)

  if(gene.list=="Null"){
    gene.list <- unique(m_t[,2])
  }
  
  # 1.1 ORA all gene
  #ORA.KEGG.all=(data.frame(ID="",Description="",setSize="",enrichmentScore="",NES="",pvalue="",p.adjust="",qvalue="",rank="",leading_edge="",core_enrichment="",trend="",minus.log10.pvalue="",minus.log10.padj="",minus.log10.qvalue=""))
  ORA.KEGG.all=data.frame()
  tryCatch(
    {
      ORA.KEGG.all.run <- enricher(unique(c(gene.top.up,gene.top.down)),
                                  universe = gene.list, 
                                  TERM2GENE=m_t,
                                  pvalueCutoff = pValue,
                                  pAdjustMethod = adjustPval,
                                  qvalueCutoff = qValue
                                  #, maxGSSize = 350
      )
      ORA.KEGG.all <- ORA.KEGG.all.run@result
      if(nrow(ORA.KEGG.all)>0){
        ORA.KEGG.all$trend <- "all"
        # ORA.KEGG.all$minus.log10.pvalue <- -log10(ORA.KEGG.all$pvalue)
        # ORA.KEGG.all$minus.log10.padj <- -log10(ORA.KEGG.all$p.adjust)
        # ORA.KEGG.all$minus.log10.qvalue <- -log10(ORA.KEGG.all$qvalue)
      }
    },
    error = function(e) {
      message("not enough valid all.reg gene !")
    }
  )
  
  # 1.2 ORA up gene
  #ORA.KEGG.up=(data.frame(ID="",Description="",setSize="",enrichmentScore="",NES="",pvalue="",p.adjust="",qvalue="",rank="",leading_edge="",core_enrichment="",trend="",minus.log10.pvalue="",minus.log10.padj="",minus.log10.qvalue=""))
  ORA.KEGG.up=data.frame()
  tryCatch(
    {
      ORA.KEGG.up.run <- enricher(gene.top.up,
                                  universe = gene.list, 
                                  TERM2GENE=m_t,
                                  pvalueCutoff = pValue,
                                  pAdjustMethod = adjustPval,
                                  qvalueCutoff = qValue
                                  #, maxGSSize = 350
      )
      ORA.KEGG.up <- ORA.KEGG.up.run@result
      if(nrow(ORA.KEGG.up)>0){
        ORA.KEGG.up$trend <- "up"
        # ORA.KEGG.up$minus.log10.pvalue <- -log10(ORA.KEGG.up$pvalue)
        # ORA.KEGG.up$minus.log10.padj <- -log10(ORA.KEGG.up$p.adjust)
        # ORA.KEGG.up$minus.log10.qvalue <- -log10(ORA.KEGG.up$qvalue)
      }
    },
    error = function(e) {
      message("not enough valid up.reg gene !")
    }
  )
    
  # 1.3 ORA down gene
  #ORA.KEGG.down=(data.frame(ID="",Description="",setSize="",enrichmentScore="",NES="",pvalue="",p.adjust="",qvalue="",rank="",leading_edge="",core_enrichment="",trend="",minus.log10.pvalue="",minus.log10.padj="",minus.log10.qvalue=""))
  ORA.KEGG.down=data.frame()
  tryCatch(
    {
      ORA.KEGG.down.run <- enricher(gene.top.down,
                                  universe = gene.list, 
                                  TERM2GENE=m_t,
                                  pvalueCutoff = pValue,
                                  pAdjustMethod = adjustPval,
                                  qvalueCutoff = qValue
                                  #, maxGSSize = 350
      )
      ORA.KEGG.down <- ORA.KEGG.down.run@result
      if(nrow(ORA.KEGG.down)>0){
        ORA.KEGG.down$trend <- "dw"
        # ORA.KEGG.down$minus.log10.pvalue <- -log10(ORA.KEGG.down$pvalue)
        # ORA.KEGG.down$minus.log10.padj <- -log10(ORA.KEGG.down$p.adjust)
        # ORA.KEGG.down$minus.log10.qvalue <- -log10(ORA.KEGG.down$qvalue)
      }
    },
    error = function(e) {
      message("not enough valid down.reg gene !")
    }
  )
    

  if(nrow(ORA.KEGG.up)>0 & nrow(ORA.KEGG.down)>0){
    return(as.data.frame(rbind(ORA.KEGG.up,ORA.KEGG.down,ORA.KEGG.all)))
  }else if(nrow(ORA.KEGG.up)>0){
    return(as.data.frame(ORA.KEGG.up))
  }else if(nrow(ORA.KEGG.down)>0){
    return(as.data.frame(ORA.KEGG.down))
  }else{
    return(data.frame(ID="",Description="",GeneRatio="",BgRatio="",pvalue="",p.adjust="",qvalue="",geneID="",Count="",trend="")) # ,minus.log10.pvalue="",minus.log10.padj="",minus.log10.qvalue=""
  }
}




# PhastCons Conservation ---------
getPhastCons100 <- function(pos){ # chr1:1000000-1000010
  #mean: pos <- "chr1:1000000-1000010"
  #single base: pos <- c("chr1:1000000-1000001", "chr1:1000001-1000002")
  GenomicScores::gscores(phastCons100way.UCSC.hg38::phastCons100way.UCSC.hg38, GenomicRanges::GRanges(pos))$default
}








# CfPeak utility --------------------
rna <- c("rRNA","pri_miRNA","lncRNA", "mRNA", "piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA')
rna2 <- c("pri_miRNA","piRNA","rRNA","lncRNA","mRNA","snoRNA","snRNA","srpRNA","tRNA","tucpRNA","Y_RNA","ktRNA","Rfam")
dna <- c("intron", "promoter", "enhancer", "repeats")

RNA_colors <- list()
for(i in 1:length(c(rna,dna))){
  j <- c(rna,dna)[i]
  RNA_colors[[j]] <- c(pal_nejm_adaptive()(15)[1:14],"#11838D")[i] #pal_d3_adaptive()(15)[i]
}
RNA_colors[['miRNA']] <- "#5E576FFF"
RNA_colors[['all']] <- "grey50" #"grey50"
RNA_colors[['other']] <- "grey70" #"grey50"
RNA_colors[['RNA']] <- "red2"
RNA_colors[['DNA']] <- "purple" #"#117C8E": similar with repeats
RNA_colors[['spikein']] <- "grey80" 
RNA_colors[['univec']] <- "grey50" 
RNA_colors[['Rfam']] <- "darkred" 
RNA_colors[['ktRNA']] <- "salmon1" 
RNA_colors <- do.call("c",RNA_colors)


# CfPeak cancer Peak index score func. ------------------
#tail gene not correct outlier gene number, just sum up all |zsocre|>3, thus not suited for multi-cancer classify
getRawPeakIndex <- function(feature.lab,logcpm,sample.table,prescale=F,mean.list,sd.list,zscore.cutoff=100,onlyCountZscoreOutlier=F){
  #x <- paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",dst,"/cfPeakCNN_GSE110381_smallDomain_diff_COADvsLAML.cpm")
  #logcpm <- read.table(paste0(pre,"/output/",dst,"/call_peak_all/count_matrix/expeakCNN_b5_d50_p1_GSE110381_CPMrowsum.txt"),check.names = F,header = T) # cpm
  #logcpm <- read.table(x,check.names = F,header = T)  # log2cpm+1
  #rownames(logcpm) <- logcpm$gene_id # feature
  # feature.lab <- "COAD"  #"all"
  cpm <- 2^logcpm-1 # sum before log is more sensitive
  # message(table(rownames(cpm) %in% feature.list[[feature.lab]]))
  cpm <- cpm[rownames(cpm) %in% feature.list[[feature.lab]],] # some may lost some features
  #sample.table <- sample.table.tcga
  #cpm <- cpm[,2:ncol(cpm)]
  sample.table <- sample.table[sample.table$sample %in% colnames(cpm),]
  # table(sample.table$group)
  # positive_samples <- sample.table[sample.table$group==disease,"sample"]
  # negative_samples <- sample.table[sample.table$group!=disease,"sample"]
  # grpNum <- length(unique(sample.table$group))
  # if(grpNum>2){
  #   other_samples <- list()
  #   for(grp in as.character(unique(sample.table$group[!(sample.table$group %in% c(disease,NC.label))] ))   ){
  #     #grp <- "PAAD" 
  #     #grp <- "PRAD"
  #     print(grp)
  #     other_samples[[grp]] <- sample.table[sample.table$group==grp,"sample"]
  #     samples <- c(samples,other_samples[[grp]])
  #     group <- c(group,rep(grp,length(other_samples[[grp]])))
  #   }
  #   sample.table <- sample.table[match(samples,sample.table$sample),]
  # }
  cpm <- cpm[,sample.table$sample]
  message("dim(cpm): ",dim(cpm))
  
  if(prescale){
    message("prescale")
    #standardize by coad&blood group
    cpm <- (cpm-mean.list[[feature.lab]])/sd.list[[feature.lab]]
    #trim to avoid too much outlier, may not need  
    # print(max(cpm))
    # print(min(cpm))
    # zscore.cutoff <- 10
    # tmp <- cpm > zscore.cutoff
    cpm[cpm > zscore.cutoff] <- zscore.cutoff
    cpm[cpm < -zscore.cutoff] <- -zscore.cutoff
    # onlyCountZscoreOutlier <- 3
    if(onlyCountZscoreOutlier){
      message("onlyCountZscoreOutlier: ",onlyCountZscoreOutlier)
      table(cpm>=abs(onlyCountZscoreOutlier))
      cpm[cpm<abs(onlyCountZscoreOutlier)] <- 0 # abs
      cpm[cpm>=abs(onlyCountZscoreOutlier)] <- 1 # abs
      message(cpm[1:2,1:2])
    }
  } else {
    message("not prescale")
  }
  return(cpm)
}

getPeakIndex <- function(feature.lab,logcpm,sample.table,prescale=F,mean.list,sd.list,zscore.cutoff=100,onlyCountZscoreOutlier=F){
  # feature.lab <- "COAD"
  # sample.table <- sample.table.tcga
  # prescale <- T
  # # mean.list <- mean.list[[]]
  # zscore.cutoff <- 100
  cpm <- getRawPeakIndex(feature.lab=feature.lab,logcpm = logcpm,sample.table = sample.table,prescale = prescale,mean.list = mean.list,sd.list = sd.list,zscore.cutoff = zscore.cutoff,onlyCountZscoreOutlier=onlyCountZscoreOutlier)
  library(tidyr)
  library(dplyr)
  conflict_prefer("summarise", "dplyr")
  cpm$feature <- rownames(cpm)
  cpm <- as_tibble(cpm) %>% 
    pivot_longer(cols = 1:(ncol(cpm)-1), names_to = "sample", values_to = "log2cpm")   # cpm acutually !!!
  cpm$group <- sample.table$group[match(cpm$sample,sample.table$sample)]
  #table((cpm$group))
  cpm.sum <- cpm %>% 
    group_by(sample,group) %>% 
    summarise(sum.log2cpm=sum(log2cpm),mean.log2cpm=mean(log2cpm)) 
  cpm.sum$sum.log2cpm.scale <- scale(cpm.sum$sum.log2cpm)[,1] # standize scale
  cpm.sum$sum.log2cpm.scale <- maxmin.normalize.vec(cpm.sum$sum.log2cpm.scale) # min-max scale (better for visualize)
  cpm.sum$mean.log2cpm.scale <- scale(cpm.sum$mean.log2cpm)[,1] # standize scale
  cpm.sum$mean.log2cpm.scale <- maxmin.normalize.vec(cpm.sum$mean.log2cpm.scale) # min-max scale (better for visualize)
  
  # grpNum <- length(unique(sample.table$group))
  # if(grpNum==2){
  #   cpm.sum$group <- factor(cpm.sum$group,levels = c(NC.label,CRC.label))
  # }else if(grpNum>2){
  #   cpm.sum$group <- factor(cpm.sum$group,levels = c(NC.label,CRC.label, unique(sample.table$group[!(sample.table$group %in% c(CRC.label,NC.label))] )))
  # }
  cpm.sum$peak.precursor <- feature.lab
  return(cpm.sum)
}





