#domain intersect with RBPs, G4/i-motif
setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/")
#setwd("/Users/baopengfei/Desktop/lulab/tmp/projects/WCHSU-FTC/exSeek-dev")
options(stringsAsFactors = F)
options(bedtools.path = "/BioII/lulab_b/baopengfei/biosoft") # /BioII/lulab_b/baopengfei/anaconda3/envs/exvariance/bin
suppressPackageStartupMessages(library(bedtoolsr))


# def func. -----------------------------
getMethod <- function(method){
  if (grepl("expeak|cfpeak",method,perl = T)){
    method.name <- "cfPeak"
  }else if  (grepl("localmax",method)){
    method.name <- "LocalMax"
  }else if  (grepl("clam",method)){
    method.name <- "CLAM"
  }else if  (grepl("clipper",method)){
    method.name <- "CLIPper"
  }else if (grepl("piranha",method)){
    method.name <- "Piranha"
  }
  return(method.name)
}

conv.mfe <- function(x){
  # x <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008/call_peak_all/expeak/b5_d50_p1_11RNA.bed.exp.RNAfold"
  print(x)
  rnafold <- readLines(paste0(x))
  rnafold <- as.data.frame(rnafold)
  #head(rnafold,3)
  colnames(rnafold) <- "V1"
  #rnafold$V3 <- unlist(base::sapply(base::strsplit(rnafold$V1,"[",fixed = T),"[",2))
  V3 <- rnafold$V1[!grepl(">",rnafold$V1)]
  
  rnafold <- rnafold[grepl(">",rnafold$V1),,drop=F]
  mir <- rnafold$V1
  mir <- gsub(">","",mir)
  mir <- unlist(sapply(strsplit(mir,"::",fixed = T),"[",1))
  mir <- unlist(sapply(strsplit(mir,"(",fixed = T),"[",1))
  
  V3 <- unlist(sapply(strsplit(as.character(V3),"[",fixed = T),"[",2))
  rnafold$V3 <- as.numeric(gsub("]","",V3 ))
  
  mfe <- rnafold$V3[rnafold$V3!=""]
  mfe <- as.numeric(mfe[!is.na(mfe)])
  rnafold.out <- data.frame(mir=mir)
  rnafold.out$mfe <- mfe
  rnafold.out$path <- x

  # rnafold.out$peak.id <- unlist(sapply(strsplit(rnafold.out$mir,"--"),"[",1))
  rnafold.out$method <- sapply(rnafold.out$path,function(x) getMethod(x))
  # print(summary(rnafold.out$mfe))
  
  # get RNA types
  #x <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008/call_domain_withRepeats_all/domains/b5_p05_11RNA.bed.flank.RNAfold"
  x2 <- gsub(".exp.RNAfold","",x)
  x2 <- gsub(".RNAfold","",x2)
  bed <- read.table(x2)
  
  rnafold.out$txID <- bed$V1[match(rnafold.out$mir,bed$V4)]
  rnafold.out$RNA <- ref$transcript_type[match(rnafold.out$txID,ref$transcript_id)]
  return(rnafold.out)
}


conv.mfe.expeakOnly <- function(x,y){
  # x <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008/call_peak_all/expeak/b5_d50_p1_11RNA.bed.exp.RNAfold"
  print(x)
  rnafold <- readLines(paste0(x))
  rnafold <- as.data.frame(rnafold)
  #head(rnafold,3)
  colnames(rnafold) <- "V1"
  #rnafold$V3 <- unlist(base::sapply(base::strsplit(rnafold$V1,"[",fixed = T),"[",2))
  V3 <- rnafold$V1[!grepl(">",rnafold$V1)]
  
  rnafold <- rnafold[grepl(">",rnafold$V1),,drop=F]
  mir <- rnafold$V1
  mir <- gsub(">","",mir)
  mir <- unlist(sapply(strsplit(mir,"::",fixed = T),"[",1))
  mir <- unlist(sapply(strsplit(mir,"--",fixed = T),"[",1)) # shuffle
  mir <- unlist(sapply(strsplit(mir,"(",fixed = T),"[",1))
  
  V3 <- unlist(sapply(strsplit(as.character(V3),"[",fixed = T),"[",2))
  rnafold$V3 <- as.numeric(gsub("]","",V3 ))
  
  mfe <- rnafold$V3[rnafold$V3!=""]
  mfe <- as.numeric(mfe[!is.na(mfe)])
  rnafold.out <- data.frame(mir=mir)
  rnafold.out$mfe <- mfe
  rnafold.out$path <- x
  
  expeakOnly <- data.table::fread(y,data.table = F,sep = "\t",header = F,check.names = F,stringsAsFactors = F)
  rnafold.out <- rnafold.out[rnafold.out$mir %in% expeakOnly$V4,]
  
  # rnafold.out$peak.id <- unlist(sapply(strsplit(rnafold.out$mir,"--"),"[",1))
  rnafold.out$method <- sapply(rnafold.out$path,function(x) getMethod(x))
  # print(summary(rnafold.out$mfe))
  
  # get RNA types
  #x <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008/call_domain_withRepeats_all/domains/b5_p05_11RNA.bed.flank.RNAfold"
  x2 <- gsub(".exp.RNAfold","",x)
  x2 <- gsub(".RNAfold","",x2)
  bed <- read.table(x2)
  
  rnafold.out$txID <- bed$V1[match(rnafold.out$mir,bed$V4)]
  rnafold.out$RNA <- ref$transcript_type[match(rnafold.out$txID,ref$transcript_id)]
  return(rnafold.out)
}

# Calculate mean ratios, fold change, and confidence intervals for each method
calc_mean_cl_normal <- function(x) {
  mean_val <- mean(x, na.rm = TRUE)  # Calculate mean
  stderr <- sd(x, na.rm = TRUE) / sqrt(length(x))  # Calculate standard error
  error_margin <- 1.96 * stderr  # 95% confidence interval margin (normal distribution)
  c(mean_val, mean_val - error_margin, mean_val + error_margin)  # Return mean, lower, and upper CI
}

calc_mean_sem <- function(x) {
  x <- na.omit(x)  # Remove NA values
  n <- length(x)   # Sample size without NA
  if (n == 0) return(c(NA, NA, NA))  # Handle the case of no valid data
  mean_val <- mean(x)  # Calculate mean
  sem <- sd(x) / sqrt(n)  # Calculate standard error of the mean
  c(mean_val, mean_val - sem, mean_val + sem)  # Return mean, lower SEM, and upper SEM
}

read.intersect <- function(x){
  #x <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008/call_domain_withRepeats_all/domains_by_sample/b5_p10/intersect/SRR2105336.intersect.bed"
  res.df <- data.table::fread(x,data.table = F,sep = "\t",header = T,check.names = F,stringsAsFactors = F)
  res.df$RNA <- ref$transcript_type[match(res.df$chr,ref$transcript_id)]
  #table(is.na(res.df$site_overlap_ratio_product))
  res.df$site_overlap_ratio_product[is.na(res.df$site_overlap_ratio_product)] <- 0
  
  res.df <- dplyr::as_tibble(res.df) %>% 
    dplyr::arrange(desc(site),desc(site_overlap_ratio_product),desc(site_overlap_base)) %>% 
    dplyr::distinct(name, site, .keep_all = TRUE) 
  
  res.df$RNA <- gsub("_rev|_for|\\.for|\\.rev","",res.df$RNA)
  rna <- c("rRNA","pri_miRNA","lncRNA", "mRNA", "piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA')
  #dna <- c("intron","promoter", "enhancer","repeats") # 
  res.df$RNA <- factor(res.df$RNA,levels = c(rna))
  
  res.df$site <- factor(res.df$site,levels = c("AGO2","otherRBPs","RBPhotspot","EV","G4","iM")) # change to before dplyr::distinct?
  # res.df$pos <- paste0(res.df$chr,":",res.df$start,"-",res.df$end) 
  
  res.tbl <- as_tibble(res.df) %>% 
    dplyr::mutate(RBP_structure=dplyr::case_when( (site_overlap_ratio_product>=0.01 ) ~ T,  # | site_overlap_base>=5 may prefer long peak caller
                                                  # site=="structure" & (site_overlap_ratio_product<=0.05) ~ T,
                                                  TRUE ~ F))
  
  mat <- res.tbl %>% 
    dplyr::select(RBP_structure,site,RNA,name) %>% 
    dplyr::group_by(site,RNA) %>% 
    dplyr::mutate(number=n_distinct(name), ratio=sum(RBP_structure)/n_distinct(name)) 
  # tidyr::pivot_wider(id_cols = name, names_from = site, values_from = RBP_structure)
  
  res <- as.data.frame(mat)
  res$path <- x
  res$method <- getMethod(x)
  li[[x]] <- res
}
read.intersect.expeakOnly <- function(x,y){
  # should be both expeak bed input
  # x <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008/call_peak_all/expeak_by_sample/b5_d50_p1/intersect/SAMN03863396.intersect.bed"
  # y <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008/call_peak_all/expeak_by_sample/b5_d50_p1/intersect/SAMN03863396_expeakOnly.bed6"
  print(x)
  res.df <- data.table::fread(x,data.table = F,sep = "\t",header = T,check.names = F,stringsAsFactors = F)
  res.df$peak <- unlist(sapply(strsplit(res.df$name,"--"),"[",1)) # paste0(res.df$chr,"|",res.df$start,"|",res.df$end) #
  #res.df[1:3,]
  
  expeakOnly <- data.table::fread(y,data.table = F,sep = "\t",header = F,check.names = F,stringsAsFactors = F)
  # expeakOnly$peak <- unlist(sapply(strsplit(res.df$name,"--"),"[",1)) # paste0(expeakOnly$V1,"|",expeakOnly$V2,"|",expeakOnly$V3) #
  
  #expeakOnly[1:3,]
  res.df <- res.df[res.df$peak %in% expeakOnly$V4,]
  res.df$peak <- NULL
  res.df$RNA <- ref$transcript_type[match(res.df$chr,ref$transcript_id)]
  #table(is.na(res.df$site_overlap_ratio_product))
  res.df$site_overlap_ratio_product[is.na(res.df$site_overlap_ratio_product)] <- 0
  
  res.df <- dplyr::as_tibble(res.df) %>% 
    dplyr::arrange(desc(site),desc(site_overlap_ratio_product),desc(site_overlap_base)) %>% 
    dplyr::distinct(name, site, .keep_all = TRUE) 
  
  res.df$RNA <- gsub("_rev|_for|\\.for|\\.rev","",res.df$RNA)
  rna <- c("rRNA","pri_miRNA","lncRNA", "mRNA", "piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA')
  #dna <- c("intron","promoter", "enhancer","repeats") # 
  res.df$RNA <- factor(res.df$RNA,levels = c(rna))
  
  res.df$site <- factor(res.df$site,levels = c("AGO2","otherRBPs","RBPhotspot","EV","G4","iM")) # ,"structure"
  # res.df$pos <- paste0(res.df$chr,":",res.df$start,"-",res.df$end) 
  
  res.tbl <- as_tibble(res.df) %>% 
    dplyr::mutate(RBP_structure=dplyr::case_when( (site_overlap_ratio_product>=0.01 ) ~ T,  # | site_overlap_base>=5 may prefer long peak caller
                                                  # site=="structure" & (site_overlap_ratio_product<=0.05) ~ T,
                                                  TRUE ~ F))
  
  mat <- res.tbl %>% 
    dplyr::select(RBP_structure,site,RNA,name) %>% 
    dplyr::group_by(site,RNA) %>% 
    dplyr::mutate(number=n_distinct(name), ratio=sum(RBP_structure)/n_distinct(name)) 
  # tidyr::pivot_wider(id_cols = name, names_from = site, values_from = RBP_structure)
  
  res <- as.data.frame(mat)
  res$path <- x
  res$method <- getMethod(x)
  li[[x]] <- res
}
read.intersect.8DNA <- function(x){
  #x <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008/call_domain_withRepeats_all/domains_by_sample/b5_p10/intersect/SRR2105336.intersect.bed"
  res.df <- data.table::fread(x,data.table = F,sep = "\t",header = T,check.names = F,stringsAsFactors = F)
  wid <- as.numeric(res.df$end)-as.numeric(res.df$start)
  res.df <- res.df[wid>=10 & wid <= 200,] # remove too long (may be exon-spanning)
  
  res.df$RNA <- "8DNA" # ref$transcript_type[match(res.df$chr,ref$transcript_id)]
  #table(is.na(res.df$site_overlap_ratio_product))
  res.df$site_overlap_ratio_product[is.na(res.df$site_overlap_ratio_product)] <- 0
  
  res.df <- dplyr::as_tibble(res.df) %>% 
    dplyr::arrange(desc(site),desc(site_overlap_ratio_product),desc(site_overlap_base)) %>% 
    dplyr::distinct(name, site, .keep_all = TRUE) 
  
  # res.df$RNA <- gsub("_rev|_for|\\.for|\\.rev","",res.df$RNA)
  # rna <- c("rRNA","pri_miRNA","lncRNA", "mRNA", "piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA')
  #dna <- c("intron","promoter", "enhancer","repeats") # 
  # res.df$RNA <- factor(res.df$RNA,levels = c(rna))
  
  res.df$site <- factor(res.df$site,levels = c("AGO2","otherRBPs","RBPhotspot","EV","G4","iM")) # ,"structure"
  # res.df$pos <- paste0(res.df$chr,":",res.df$start,"-",res.df$end) 
  
  res.tbl <- as_tibble(res.df) %>% 
    dplyr::mutate(RBP_structure=dplyr::case_when( (site_overlap_ratio_product>=0.01 ) ~ T,  # | site_overlap_base>=5 may prefer long peak caller
                                                  # site=="structure" & (site_overlap_ratio_product<=0.05) ~ T,
                                                  TRUE ~ F))
  
  mat <- res.tbl %>% 
    dplyr::select(RBP_structure,site,name,RNA) %>% 
    dplyr::group_by(site,RNA) %>% 
    dplyr::mutate(number=n_distinct(name), ratio=sum(RBP_structure)/n_distinct(name)) 
  # tidyr::pivot_wider(id_cols = name, names_from = site, values_from = RBP_structure)
  
  res <- as.data.frame(mat)
  res$path <- x
  res$method <- getMethod(x)
  li[[x]] <- res
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







# plot intersect enrichment: sample-wise boxplot/barplot (11RNA, dot is sample) -----------------------------------
#need run domain_intersect_RBPG4.R before
ref <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)
table(ref$transcript_type)

li <- list()
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
dst="GSE71008" # "GSE71008","GSE110381","GSE123972" # test use multi samplses
bg.type <- "Background"  # "Background" "Flank"
region.type <- c("",".flank",".shuffle") # .shuffle2: shufSameRNATxInPriority, seem little effect
methods <- c("piranha_by_sample/b5_p01","clipper_by_sample/b5_p05","clam_by_sample/b5_p005","expeakCNN_by_sample/b5_d50_p1") # ,"_localmax/b5_d05_p05"
smps <- read.table(paste0(pre,"/exSeek-dev/data/",dst,"/sample_ids_NCpool_test15.txt"))$V1
fs <- as.data.frame(expand_grid(methods=methods,smps=smps,region=region.type))
fs$path <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_all/",fs$methods,"/intersect/",fs$smps,".intersect.bed",fs$region) # 11RNA
fs <- unique(fs$path)

li2 <- parallel::mclapply(fs,read.intersect,mc.cores = 5) # 11RNA
res.tbl <- as_tibble(do.call(rbind,li2))
#table(res.tbl$method)
res.tbl$sample <- unlist(sapply(strsplit(res.tbl$path,"/",fixed=T),"[",13)) # 14
#table(res.tbl$sample )
res.tbl$type <- "Peak"
res.tbl$type[grepl("shuffle",res.tbl$sample)] <- "Background"
res.tbl$type[grepl("flank",res.tbl$sample)] <- "Flank"
res.tbl$sample <- unlist(sapply(strsplit(res.tbl$sample,".",fixed=T),"[",1))
res.tbl$type <- factor(res.tbl$type, levels = c("Peak","Background","Flank"))
res.tbl <- res.tbl %>% 
  dplyr::filter(type %in% c(bg.type,"Peak")) # Background, Flank
table(res.tbl$type)


## plot
library(ggplot2)
library(ggpubr) # need lib explicit
res.df2.perSample <- res.tbl %>% 
  dplyr::select(type,site,sample,RNA,number,ratio,method) %>% 
  dplyr::distinct( .keep_all = TRUE) %>% 
  dplyr::group_by(type,site,sample,method) %>% 
  dplyr::summarise(total.num=sum(number,na.rm=T),sig.number=sum(number*ratio,na.rm=T),ratio=sig.number/total.num) 
res.df2.perSample2 <- res.tbl %>% 
  dplyr::mutate(site2=dplyr::case_when(site %in% c("AGO2","otherRBPs") ~ "RBP",site=="G4" ~ "G4",site=="iM" ~ "iM",site=="RBPhotspot" ~ "RBPhotspot",site=="EV" ~ "EV",TRUE ~ "NA") ) %>% 
  dplyr::filter(site2 %in% c("RBP","G4","EV")) %>%  # rm hotspot,iM
  dplyr::group_by(type,sample,name,site2,method) %>% # add sample,type
  dplyr::mutate(RBP_structure2=as.logical(sum(RBP_structure))) %>% 
  dplyr::distinct(type,sample,name,site2,method,.keep_all = TRUE) %>% 
  dplyr::group_by(type,sample,site2,method) %>% 
  dplyr::summarise(number=n_distinct(name), ratio=sum(RBP_structure2)/n_distinct(name)) %>% # re-cal combined num/ratio
  dplyr::distinct( .keep_all = TRUE) # %>% 


#  ggpubr::stat_compare_means(ref.group="Peak", label = "p.signif",  hide.ns = F, paired = T)
res.df2.perSample <- as.data.frame(res.df2.perSample)
res.df2.perSample$type <- as.factor(res.df2.perSample$type)
res.df2.perSample$site <- factor(res.df2.perSample$site, levels = c("AGO2","otherRBPs","RBPhotspot","EV","G4","iM"))
res.df2.perSample$method <- factor(res.df2.perSample$method, levels = c("Piranha","CLIPper","CLAM","cfPeak"))
#table(res.df2.perSample$method)
#res.df2.perSample

p11 <- ggbarplot(data = res.df2.perSample, x = "method", y = "ratio", fill = "type", facet.by=c("site"), add="mean_se", short.panel.labs = T, position = position_dodge()) + # facet.by=c("site"), mean_ci mean_se
  stat_compare_means(aes(x=method, y=ratio, group=type), # ref.group="Flank",  # 
                     # comparisons = list(c("Peak",bg.type)),
                     label="p.signif",  # ..p.signif../..p.format..
                     # symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01,1),symbols = c( "***", "**", "*", "ns")),
                     method = "wilcox.test", method.args = list(alternative = "less"),  # greater means ref.group less
                     label.x.npc = 0.8, label.y.npc = 0.8, vjust = 0.6, #step.increase = 0.08,
                     hide.ns=T,size =10, paired = F
  ) +
  # facet_grid(site~.,scales = "free")+ # facet_grid seem not work well with ggbarplot
  scale_fill_manual(values = c('firebrick','salmon'))+ #c("red","grey90")
  # scale_x_discrete(label=c("Peak Ratio","MIR Ratio"))+
  labs(title="",x="", y = "")+
  # xlim(c(0,1))+   
  ylim(c(0,1))+
  theme_minimal() +  # base_size=12
  theme(#axis.ticks.x=element_blank(),
    #strip.text.y = element_blank(),
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
p11
#

res.df2.perSample2 <- as.data.frame(res.df2.perSample2)
res.df2.perSample2$type <- as.factor(res.df2.perSample2$type)
#str(res.df2.perSample2)
#table(res.df2.perSample2$site2)
res.df2.perSample2 <- res.df2.perSample2[res.df2.perSample2$site2!="iM",]
res.df2.perSample2$site2 <- factor(res.df2.perSample2$site2, levels = c("RBP","EV","G4","iM"))
res.df2.perSample2$method <- factor(res.df2.perSample2$method, levels = c("Piranha","CLIPper","CLAM","cfPeak"))

res.df2.perSample2.2 <- res.df2.perSample2[res.df2.perSample2$method=="cfPeak",]

# res.df2.perSample2.2
# ggplot(res.df2.perSample2.2)+
#   geom_boxplot(aes(x=site2,y=ratio,fill=type)) +
#   facet_grid(site2~.)

# tmp1 <- res.df2.perSample2.2$ratio[res.df2.perSample2.2$site2=="G4" & res.df2.perSample2.2$type=="Peak"]
# tmp2 <- res.df2.perSample2.2$ratio[res.df2.perSample2.2$site2=="G4" & res.df2.perSample2.2$type=="Background"]
# wilcox.test(x=tmp1,y=tmp2)
p12 <- ggbarplot(data = res.df2.perSample2.2, x = "method", y = "ratio", fill = "type", facet.by=c("site2"), add="mean_se", short.panel.labs = T, position = position_dodge()) + # add.params=list(size=4,color="black"), mean_ci mean_se
  stat_compare_means(aes(x=method, y=ratio, group=type), # ref.group="Flank",  #
                     # comparisons = list(c("Peak",bg.type)),
                     label="p.signif",  # ..p.signif../..p.format..
                     # symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01,1),symbols = c( "***", "**", "*", "ns")),
                     method = "wilcox.test", method.args = list(alternative = "less"),  # greater means ref.group less
                     label.x.npc = 0.8, label.y.npc = 0.8, vjust = 0.6, #step.increase = 0.08,
                     hide.ns=T,size =10, paired = F
  ) +
  # facet_grid(site~.,scales = "free")+ # facet_grid seem not work well with ggbarplot
  scale_fill_manual(values = c('firebrick','salmon'))+ #c("red","grey90")
  # scale_x_discrete(label=c("Peak Ratio","MIR Ratio"))+
  labs(title="",x="", y = "")+
  # xlim(c(0,1))+
  ylim(c(0,1))+
  theme_minimal() +  # base_size=12
  theme(#axis.ticks.x=element_blank(),
    #strip.text.y = element_blank(),
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
p12
# ggsave(filename = "domain_enrich_sum.pdf", width = 10, height = 7) # 11RNA
ggsave(filename = "figure/domain_enrich_sum_expeak.pdf", width = 6, height = 7) # expeak


fc_data <- res.df2.perSample2.2 %>%
  group_by(method, site2) %>%
  summarise(mean_firebrick = mean(ratio[type == 'Peak']),
            ci_firebrick = calc_mean_cl_normal(ratio[type == 'Peak']),
            mean_salmon = mean(ratio[type == 'Background']),
            ci_salmon = calc_mean_cl_normal(ratio[type == 'Background']),
            FC = mean_firebrick / mean_salmon,
            delta = mean_firebrick-mean_salmon) %>%
  mutate(ci_firebrick_text = sprintf("%.4f ± %.4f", ci_firebrick[1], ci_firebrick[3] - ci_firebrick[1]),
         ci_salmon_text = sprintf("%.4f ± %.4f", ci_salmon[1], ci_salmon[3] - ci_salmon[1]),
         FC_label = sprintf("FC=%.2f\nDelta=%.2f\nPeak: %s\nBackground: %s", FC, delta, ci_firebrick_text, ci_salmon_text)) %>%
  ungroup()

# Merge FC data back with the original dataframe
res.df2.perSample2.2 <- left_join(res.df2.perSample2.2, fc_data, by = c("method", "site2"))

# Create the bar plot with ggplot2
p12 <- ggplot(res.df2.perSample2.2, aes(x = method, y = ratio, fill = type)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge()) + 
  geom_errorbar(stat = "summary", fun.data = "mean_cl_normal", position = position_dodge(0.9), width = 0.2) + # calc_mean_sem mean_cl_normal
  stat_compare_means(aes(group = type), method = "wilcox.test", method.args = list(alternative = "less"),
                     label = "p.signif", label.x.npc = "right", label.y.npc = "top", hide.ns = TRUE, size = 10) +
  facet_wrap(~ site2, scales = "free") + 
  scale_fill_manual(values = c('firebrick', 'salmon')) +
  
  # Add FC and 95% CI labels for each method
  geom_text(data = fc_data, aes(x = method, y = 0.9, label = FC_label), 
            inherit.aes = FALSE, size = 5) + 
  
  labs(title = "", x = "", y = "") +
  ylim(c(0, 1)) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 20, hjust = 1, vjust = 0.5, angle = 90),
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size = 20),
    legend.position = "none",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16)
  )
p12
# Save the plot
# ggsave(filename = "figure/domain_enrich_sum_expeak.pdf", plot = p12, width = 6, height = 8)
ggsave(filename = "figure/domain_enrich_sum_expeak.pdf", plot = p12, width = 10, height = 8)



library(ggplot2)
# res.df <- as.data.frame(res.tbl[,c("ratio","std","site","RNA","type")])
#res.tbl <- res.tbl[sample(1:nrow(res.tbl),size = 0.1*nrow(res.tbl),replace = F),]
res.df3 <- res.tbl %>% 
  dplyr::select(type,site,sample,RNA,number,ratio,method) %>% 
  dplyr::distinct( .keep_all = TRUE) %>% 
  dplyr::group_by(type,site,sample,method) %>% 
  dplyr::mutate(total.num=sum(number,na.rm=T),sig.number=sum(number*ratio,na.rm=T)) %>% 
  dplyr::group_by(type,site,RNA,method) %>% 
  dplyr::mutate(std=sd(ratio,na.rm=T), mean.ratio=mean(ratio,na.rm=T), num.std=sd(number,na.rm=T), mean.number=mean(number,na.rm=T))
res.df3 <- as.data.frame(res.df3)
res.df3$type <- factor(res.df3$type, levels = c("Peak","Background","Flank"))
res.df3$method <- factor(res.df3$method, levels = c("Piranha","CLIPper","CLAM","exPeak"))
rna <- c("rRNA","pri_miRNA","lncRNA", "mRNA", "piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA')
#dna <- c("intron","promoter", "enhancer","repeats") # 
res.df3$RNA <- factor(res.df3$RNA,levels = c(rna))
res.df3 <- res.df3 %>% 
  dplyr::filter(type %in% c(bg.type,"Peak")) #  Flank

for(tmp in unique(res.df3$site)){
  # tmp <- "RBPhotspot"
  print(tmp)
  res.df3.tmp <- res.df3[res.df3$site==tmp,]
  p22 <- ggplot(res.df3.tmp, aes(x=method, y=mean.ratio, fill=type)) +
    geom_bar(position=position_dodge(), stat = "summary",color="black")+ # 
    geom_errorbar(position=position_dodge(.9), aes(ymin=mean.ratio-std, ymax=mean.ratio+std, group=type), width=.2)+
    ggpubr::stat_compare_means(aes(x=method, y=mean.ratio, group=type), 
                               label = "p.signif", # ..p.signif../..p.format..
                               # ref.group = "Peak", comparisons,
                               # symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01,1),symbols = c( "***", "**", "*", "ns")),
                               method = "wilcox.test", paired = T, method.args = list(alternative = "less"),
                               label.x.npc = 0.4, label.y.npc = 0.85, step.increase = 0.08, 
                               hide.ns=T,size =10
    )+
    scale_fill_manual(values = c('firebrick','salmon'))+ #c("red","grey90")
    # scale_x_discrete(label=c("Peak Ratio","MIR Ratio"))+
    labs(title="",x="", y = "")+
    # facet_grid(method~.)+
    ylim(c(0,1))+
    # xlim(c(0,1))+
    theme_minimal() +  # base_size=12
    theme(#axis.ticks.x=element_blank(),  
      #strip.text.y = element_blank(),
      #strip.text.x = element_text(face="bold",family="arial",size=20),
      axis.title.x = element_text(size=20),
      axis.title.y = element_text(size=20),
      axis.text.x = element_text(size = 20,hjust = 1,vjust = 0.5,angle = 90), # 
      axis.text.y = element_text(size = 20),
      plot.title = element_text(size=20),
      strip.text = element_text(size = 20),
      legend.position = "right", #c(0.9,0.8),#,# 
      legend.text = element_text(size= 16),
      legend.title= element_text(size= 16)) +
    facet_grid(RNA~.,scales = "free")
  p22
  ggsave(filename = paste0("domain_enrich_",tmp,".pdf"), width = 7, height = 22) # 11RNA only
}
  
  
  
  


# plot intersect enrichment: sample-wise boxplot/barplot (11RNA, dot is sample, expeak only OR expeak gold) -----------------------------------
#need run domain_intersect_RBPG4.R before
ref <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)
table(ref$transcript_type)

li <- list()

pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
dst="GSE71008" # "GSE71008","GSE110381","GSE123972" # test use multi samplses
bg.type <- "Background"  # "Background" "Flank"
region.type <- c("",".flank",".shuffle") # .shuffle, .shuffle2
methods <- c("expeakCNN_by_sample/b5_d50_p1") # "piranha_by_sample/b5_p01","clipper_by_sample/b5_p05","clam_by_sample/b5_p005",
smps <- read.table(paste0(pre,"/exSeek-dev/data/",dst,"/sample_ids_NCpool_test15.txt"))$V1
fs <- as.data.frame(expand_grid(methods=methods,smps=smps,region=region.type))
fs$path <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_all/",fs$methods,"/intersect/",fs$smps,".intersect.bed",fs$region) # 11RNA
fs <- unique(fs$path)
#fs <- Sys.glob(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_all/",methods,"/intersect/",smps,".intersect.bed*")) # 11RNA
fs2 <- paste0(unlist(sapply(strsplit(fs,".intersect.bed",fixed=T),"[",1)), "_expeakCNNOnly.bed6") # _expeakCNNOnly.bed6, _expeakCNNGold.bed6


# li2 <- parallel::mclapply(1:length(fs),FUN = function(i) {read.intersect.expeakOnly(x=fs[i],y=fs2[i])},mc.cores = 1) # 11RNA
li2 <- lapply(1:length(fs),FUN = function(i) {read.intersect.expeakOnly(x=fs[i],y=fs2[i])} )
# lapply(1:nrow(cmb), function(j) {get.mat(signal=as.character(cmb$signal[j]), region=as.character(cmb$region[j]),
#                                          signal.label=as.character(cmb$signal.label[j]),region.label=as.character(cmb$region.label[j]),
#                                          up = up, down = down, bin = bin, ratio = ratio)} )
#names(li2)
# li2 <- parallel::mclapply(fs,read.intersect.8DNA,mc.cores = 10) # 8DNA
res.tbl <- as_tibble(do.call(rbind,li2))
#table(res.tbl$method)
res.tbl$sample <- unlist(sapply(strsplit(res.tbl$path,"/",fixed=T),"[",13)) # 14
res.tbl$type <- "Peak"
#table(res.tbl$sample)
res.tbl$type[grepl("shuffle",res.tbl$sample)] <- "Background"
res.tbl$type[grepl("flank",res.tbl$sample)] <- "Flank"
res.tbl$sample <- unlist(sapply(strsplit(res.tbl$sample,".",fixed=T),"[",1))
res.tbl$type <- factor(res.tbl$type, levels = c("Peak","Background","Flank"))
res.tbl <- res.tbl %>% 
  dplyr::filter(type %in% c(bg.type,"Peak")) # Background, Flank
table(res.tbl$type)


## plot
library(ggplot2)
library(ggpubr) # need lib explicit
res.df2.perSample <- res.tbl %>% 
  dplyr::select(type,site,sample,RNA,number,ratio,method) %>% 
  dplyr::distinct( .keep_all = TRUE) %>% 
  dplyr::group_by(type,site,sample,method) %>% 
  dplyr::summarise(total.num=sum(number,na.rm=T),sig.number=sum(number*ratio,na.rm=T),ratio=sig.number/total.num) 
res.df2.perSample2 <- res.tbl %>% 
  dplyr::mutate(site2=dplyr::case_when(site %in% c("AGO2","otherRBPs") ~ "RBP",site=="G4" ~ "G4",site=="iM" ~ "iM",site=="RBPhotspot" ~ "RBPhotspot",site=="EV" ~ "EV",TRUE ~ "NA") ) %>% 
  dplyr::filter(site2 %in% c("RBP","G4","iM","EV")) %>%  # rm hotspot
  dplyr::group_by(type,sample,name,site2,method) %>% # add sample,type
  dplyr::mutate(RBP_structure2=as.logical(sum(RBP_structure))) %>% 
  dplyr::distinct(type,sample,name,site2,method,.keep_all = TRUE) %>% 
  dplyr::group_by(type,sample,site2,method) %>% 
  dplyr::summarise(number=n_distinct(name), ratio=sum(RBP_structure2)/n_distinct(name)) %>% # re-cal combined num/ratio
  dplyr::distinct( .keep_all = TRUE) # %>% 


#  ggpubr::stat_compare_means(ref.group="Peak", label = "p.signif",  hide.ns = F, paired = T)
res.df2.perSample <- as.data.frame(res.df2.perSample)
res.df2.perSample$type <- as.factor(res.df2.perSample$type)
res.df2.perSample$site <- factor(res.df2.perSample$site, levels = c("AGO2","otherRBPs","RBPhotspot","EV","G4","iM"))
res.df2.perSample$method <- factor(res.df2.perSample$method, levels = c("Piranha","CLIPper","CLAM","exPeak"))
# table(res.df2.perSample$method)
#res.df2.perSample

p11 <- ggbarplot(data = res.df2.perSample, x = "method", y = "ratio", fill = "type", facet.by=c("site"), add="mean_se", short.panel.labs = T, position = position_dodge()) + # facet.by=c("site"),
  stat_compare_means(aes(x=method, y=ratio, group=type), # ref.group="Flank",  # 
                     # comparisons = list(c("Peak",bg.type)),
                     label="p.signif",  # ..p.signif../..p.format..
                     # symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01,1),symbols = c( "***", "**", "*", "ns")),
                     method = "wilcox.test", method.args = list(alternative = "less"),  # greater means ref.group less
                     label.x.npc = 0.8, label.y.npc = 0.8, vjust = 0.6, #step.increase = 0.08,
                     hide.ns=T,size =10, paired = F
  ) +
  # facet_grid(site~.,scales = "free")+ # facet_grid seem not work well with ggbarplot
  scale_fill_manual(values = c('firebrick','salmon'))+ #c("red","grey90")
  # scale_x_discrete(label=c("Peak Ratio","MIR Ratio"))+
  labs(title="",x="", y = "")+
  # xlim(c(0,1))+   
  ylim(c(0,1))+
  theme_minimal() +  # base_size=12
  theme(#axis.ticks.x=element_blank(),
    #strip.text.y = element_blank(),
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
p11
res.df2.perSample2 <- as.data.frame(res.df2.perSample2)
res.df2.perSample2$type <- as.factor(res.df2.perSample2$type)
#str(res.df2.perSample2)
#table(res.df2.perSample2$site2)
res.df2.perSample2 <- res.df2.perSample2[res.df2.perSample2$site2!="iM",]
res.df2.perSample2$site2 <- factor(res.df2.perSample2$site2, levels = c("RBP","EV","G4","iM"))
res.df2.perSample2$method <- factor(res.df2.perSample2$method, levels = c("Piranha","CLIPper","CLAM","exPeak"))
res.df2.perSample2[1:3,]

# res.df2.perSample2.tmp <- res.df2.perSample2 %>% 
#   dplyr::group_by(type,site2,method) %>% 
#   dplyr::summarise(mean.ratio=mean(ratio),std=sd(ratio))
# ggplot(res.df2.perSample2.tmp, aes(x=site2, y=mean.ratio, fill=type)) +
#       geom_bar(position=position_dodge(), stat = "summary",color="black")+ #
#       geom_errorbar(position=position_dodge(.9), aes(ymin=mean.ratio-std, ymax=mean.ratio+std, group=type), width=.2)
#   
p12 <- ggbarplot(data = res.df2.perSample2, x = "method", y = "ratio", fill = "type", facet.by=c("site2"), add="mean_se", short.panel.labs = T, position = position_dodge()) +
  stat_compare_means(aes(x=method, y=ratio, group=type), # ref.group="Flank",  # 
                     # comparisons = list(c("Peak",bg.type)),
                     label="p.signif",  # ..p.signif../..p.format..
                     # symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01,1),symbols = c( "***", "**", "*", "ns")),
                     method = "wilcox.test", method.args = list(alternative = "less"),  # greater means ref.group less
                     label.x.npc = 0.8, label.y.npc = 0.8, vjust = 0.6, #step.increase = 0.08,
                     hide.ns=T,size =10, paired = F
  ) +
  # facet_grid(site~.,scales = "free")+ # facet_grid seem not work well with ggbarplot
  scale_fill_manual(values = c('firebrick','salmon'))+ #c("red","grey90")
  # scale_x_discrete(label=c("Peak Ratio","MIR Ratio"))+
  labs(title="",x="", y = "")+
  # xlim(c(0,1))+   
  ylim(c(0,1))+
  theme_minimal() +  # base_size=12
  theme(#axis.ticks.x=element_blank(),
    #strip.text.y = element_blank(),
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
p12
# cowplot::plot_grid(plotlist = list(p11,p12),rel_heights = c(8,4), nrow = 2, align = "hv") # 11RNA
# cowplot::plot_grid(plotlist = list(p11,p12),rel_heights = c(4,4), nrow = 2, align = "hv") # 8DNA
ggsave(filename = "domain_enrich_sum.pdf", width = 6, height = 7) # 11RNA
# ggsave(filename = "domain_enrich_sum.pdf", width = 8, height = 10) # 8DNA

# library(ggplot2)
# # res.df <- as.data.frame(res.tbl[,c("ratio","std","site","RNA","type")])
# #res.tbl <- res.tbl[sample(1:nrow(res.tbl),size = 0.1*nrow(res.tbl),replace = F),]
# res.df3 <- res.tbl %>% 
#   dplyr::select(type,site,sample,RNA,number,ratio,method) %>% 
#   dplyr::distinct( .keep_all = TRUE) %>% 
#   dplyr::group_by(type,site,sample,method) %>% 
#   dplyr::mutate(total.num=sum(number,na.rm=T),sig.number=sum(number*ratio,na.rm=T)) %>% 
#   dplyr::group_by(type,site,RNA,method) %>% 
#   dplyr::mutate(std=sd(ratio,na.rm=T), mean.ratio=mean(ratio,na.rm=T), num.std=sd(number,na.rm=T), mean.number=mean(number,na.rm=T))
# res.df3 <- as.data.frame(res.df3)
# res.df3$type <- factor(res.df3$type, levels = c("Peak","Background","Flank"))
# res.df3$method <- factor(res.df3$method, levels = c("Piranha","CLIPper","CLAM","exPeak"))
# rna <- c("rRNA","pri_miRNA","lncRNA", "mRNA", "piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA')
# #dna <- c("intron","promoter", "enhancer","repeats") # 
# res.df3$RNA <- factor(res.df3$RNA,levels = c(rna))
# res.df3 <- res.df3 %>% 
#   dplyr::filter(type %in% c(bg.type,"Peak")) #  Flank
# 
# for(tmp in unique(res.df3$site)){
#   # tmp <- "RBPhotspot"
#   print(tmp)
#   res.df3.tmp <- res.df3[res.df3$site==tmp,]
#   p22 <- ggplot(res.df3.tmp, aes(x=method, y=mean.ratio, fill=type)) +
#     geom_bar(position=position_dodge(), stat = "summary",color="black")+ # 
#     geom_errorbar(position=position_dodge(.9), aes(ymin=mean.ratio-std, ymax=mean.ratio+std, group=type), width=.2)+
#     ggpubr::stat_compare_means(aes(x=method, y=mean.ratio, group=type), 
#                                label = "p.signif", # ..p.signif../..p.format..
#                                # ref.group = "Peak", comparisons,
#                                # symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01,1),symbols = c( "***", "**", "*", "ns")),
#                                method = "wilcox.test", paired = T, method.args = list(alternative = "less"),
#                                label.x.npc = 0.4, label.y.npc = 0.85, step.increase = 0.08, 
#                                hide.ns=T,size =10
#     )+
#     scale_fill_manual(values = c('firebrick','salmon'))+ #c("red","grey90")
#     # scale_x_discrete(label=c("Peak Ratio","MIR Ratio"))+
#     labs(title="",x="", y = "")+
#     # facet_grid(method~.)+
#     ylim(c(0,1))+
#     # xlim(c(0,1))+
#     theme_minimal() +  # base_size=12
#     theme(#axis.ticks.x=element_blank(),  
#       #strip.text.y = element_blank(),
#       #strip.text.x = element_text(face="bold",family="arial",size=20),
#       axis.title.x = element_text(size=20),
#       axis.title.y = element_text(size=20),
#       axis.text.x = element_text(size = 20,hjust = 1,vjust = 0.5,angle = 90), # 
#       axis.text.y = element_text(size = 20),
#       plot.title = element_text(size=20),
#       strip.text = element_text(size = 20),
#       legend.position = "right", #c(0.9,0.8),#,# 
#       legend.text = element_text(size= 16),
#       legend.title= element_text(size= 16)) +
#     facet_grid(RNA~.,scales = "free")
#   p22
#   ggsave(filename = paste0("domain_enrich_",tmp,".pdf"), width = 7, height = 22) # 11RNA only
# }








# plot intersect enrichment: sample-wise boxplot/barplot (8DNA, dot is sample) -----------------------------------
#need run domain_intersect_RBPG4.R before
ref <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)
table(ref$transcript_type)

li <- list()

dst="GSE71008" # "GSE71008","GSE110381"
bg.type <- "Background"
methods <- c("piranha_by_sample/b5_p01","clipper_by_sample/b5_p05","clam_by_sample/b5_p005","expeakCNN_by_sample/b5_d50_p1") # ,"_localmax/b5_d05_p05"
# fs <- Sys.glob(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_all/domains",methods,"/intersect/*.intersect.bed*")) # 11RNA
fs <- Sys.glob(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_all/",methods,"/intersect_8DNA/*.intersect.bed*")) # 8DNA


# li2 <- parallel::mclapply(fs,read.intersect,mc.cores = 10) # 11RNA
li2 <- parallel::mclapply(fs,read.intersect.8DNA,mc.cores = 1) # 8DNA
res.tbl <- as_tibble(do.call(rbind,li2))
#table(res.tbl$method)
res.tbl$sample <- unlist(sapply(strsplit(res.tbl$path,"/",fixed=T),"[",13)) # 14
res.tbl$type <- "Peak"
res.tbl$type[grepl("shuffle",res.tbl$sample)] <- "Background"
res.tbl$type[grepl("flank",res.tbl$sample)] <- "Flank"
res.tbl$sample <- unlist(sapply(strsplit(res.tbl$sample,".",fixed=T),"[",1))
res.tbl$type <- factor(res.tbl$type, levels = c("Peak","Background","Flank"))
res.tbl <- res.tbl %>% 
  dplyr::filter(type %in% c(bg.type,"Peak")) # Background, Flank


## plot
library(ggplot2)
library(ggpubr) # need lib explicit
res.df2.perSample <- res.tbl %>% 
  dplyr::select(type,site,sample,RNA,number,ratio,method) %>% 
  dplyr::distinct( .keep_all = TRUE) %>% 
  dplyr::group_by(type,site,sample,method) %>% 
  dplyr::summarise(total.num=sum(number,na.rm=T),sig.number=sum(number*ratio,na.rm=T),ratio=sig.number/total.num) 
res.df2.perSample2 <- res.tbl %>% 
  dplyr::mutate(site2=dplyr::case_when(site %in% c("AGO2","otherRBPs") ~ "RBP",site=="G4" ~ "G4",site=="iM" ~ "iM",site=="RBPhotspot" ~ "RBPhotspot",site=="EV" ~ "EV",TRUE ~ "NA") ) %>% 
  dplyr::filter(site2 %in% c("RBP","G4","iM","EV")) %>%  # rm hotspot
  dplyr::group_by(type,sample,name,site2,method) %>% # add sample,type
  dplyr::mutate(RBP_structure2=as.logical(sum(RBP_structure))) %>% 
  dplyr::distinct(type,sample,name,site2,method,.keep_all = TRUE) %>% 
  dplyr::group_by(type,sample,site2,method) %>% 
  dplyr::summarise(number=n_distinct(name), ratio=sum(RBP_structure2)/n_distinct(name)) %>% # re-cal combined num/ratio
  dplyr::distinct( .keep_all = TRUE) # %>% 


#  ggpubr::stat_compare_means(ref.group="Peak", label = "p.signif",  hide.ns = F, paired = T)
res.df2.perSample <- as.data.frame(res.df2.perSample)
res.df2.perSample$type <- as.factor(res.df2.perSample$type)
res.df2.perSample$site <- factor(res.df2.perSample$site, levels = c("AGO2","otherRBPs","RBPhotspot","EV","G4","iM"))
res.df2.perSample$method <- factor(res.df2.perSample$method, levels = c("Piranha","CLIPper","CLAM","exPeak"))
# table(res.df2.perSample$method)
#res.df2.perSample

p11 <- ggbarplot(data = res.df2.perSample, x = "method", y = "ratio", fill = "type", facet.by=c("site"), add="mean_se", short.panel.labs = T, position = position_dodge()) + # facet.by=c("site"),
  stat_compare_means(aes(x=method, y=ratio, group=type), # ref.group="Flank",  # 
                    # comparisons = list(c("Peak",bg.type)),
                     label="p.signif",  # ..p.signif../..p.format..
                     # symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01,1),symbols = c( "***", "**", "*", "ns")),
                     method = "wilcox.test", method.args = list(alternative = "less"),  # greater means ref.group less
                     label.x.npc = 0.8, label.y.npc = 0.8, vjust = 0.6, #step.increase = 0.08,
                     hide.ns=T,size =10, paired = F
                     ) +
  # facet_grid(site~.,scales = "free")+ # facet_grid seem not work well with ggbarplot
  scale_fill_manual(values = c('firebrick','salmon'))+ #c("red","grey90")
  # scale_x_discrete(label=c("Peak Ratio","MIR Ratio"))+
  labs(title="",x="", y = "")+
  # xlim(c(0,1))+   
  ylim(c(0,1))+
  theme_minimal() +  # base_size=12
  theme(#axis.ticks.x=element_blank(),
    #strip.text.y = element_blank(),
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
p11
res.df2.perSample2 <- as.data.frame(res.df2.perSample2)
res.df2.perSample2$type <- as.factor(res.df2.perSample2$type)
res.df2.perSample2 <- res.df2.perSample2[res.df2.perSample2$site2!="iM",]
res.df2.perSample2$site2 <- factor(res.df2.perSample2$site2, levels = c("RBP","EV","G4"))
res.df2.perSample2$method <- factor(res.df2.perSample2$method, levels = c("Piranha","CLIPper","CLAM","exPeak"))
p12 <- ggbarplot(data = res.df2.perSample2, x = "method", y = "ratio", fill = "type", facet.by=c("site2"), add="mean_se", short.panel.labs = T, position = position_dodge()) +
    stat_compare_means(aes(x=method, y=ratio, group=type), # ref.group="Flank",  # 
      # comparisons = list(c("Peak",bg.type)),
      label="p.signif",  # ..p.signif../..p.format..
      # symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01,1),symbols = c( "***", "**", "*", "ns")),
      method = "wilcox.test", method.args = list(alternative = "less"),  # greater means ref.group less
      label.x.npc = 0.8, label.y.npc = 0.8, vjust = 0.6, #step.increase = 0.08,
      hide.ns=T,size =10, paired = F
    ) +
  # facet_grid(site~.,scales = "free")+ # facet_grid seem not work well with ggbarplot
  scale_fill_manual(values = c('firebrick','salmon'))+ #c("red","grey90")
  # scale_x_discrete(label=c("Peak Ratio","MIR Ratio"))+
  labs(title="",x="", y = "")+
  # xlim(c(0,1))+   
  ylim(c(0,1))+
  theme_minimal() +  # base_size=12
  theme(#axis.ticks.x=element_blank(),
    #strip.text.y = element_blank(),
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
p12
# cowplot::plot_grid(plotlist = list(p11,p12),rel_heights = c(8,4), nrow = 2, align = "hv") # 11RNA
cowplot::plot_grid(plotlist = list(p11,p12),rel_heights = c(4,4), nrow = 2, align = "hv") # 8DNA
# ggsave(filename = "domain_enrich_sum.pdf", width = 8, height = 18) # 11RNA
ggsave(filename = "domain_enrich_sum.pdf", width = 8, height = 10) # 8DNA

res.df2.perSample3 <- res.df2.perSample2[res.df2.perSample2$method=="exPeak",]
p13 <- ggbarplot(data = res.df2.perSample3, x = "type", y = "ratio", fill = "type", add="mean_se", short.panel.labs = T, facet.by=c("site2"), position = position_dodge()) + # 
  stat_compare_means(aes(x=type, y=ratio, group=type), # ref.group="Flank",  # 
                     # comparisons = list(c("Peak",bg.type)),
                     label="p.signif",  # ..p.signif../..p.format..
                     # symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01,1),symbols = c( "***", "**", "*", "ns")),
                     method = "wilcox.test", method.args = list(alternative = "less"),  # greater means ref.group less
                     label.x.npc = 0.8, label.y.npc = 0.8, vjust = 0.6,  hjust = 1, #step.increase = 0.08,
                     hide.ns=T,size =10, paired = F
  ) +
  # facet_grid(site~.,scales = "free")+ # facet_grid seem not work well with ggbarplot
  scale_fill_manual(values = c('firebrick','salmon'))+ #c("red","grey90")
  # scale_x_discrete(label=c("Peak Ratio","MIR Ratio"))+
  labs(title="",x="", y = "")+
  # xlim(c(0,1))+   
  ylim(c(0,1))+
  theme_minimal() +  # base_size=12
  theme(#axis.ticks.x=element_blank(),
    #strip.text.y = element_blank(),
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
p13
ggsave(filename = "domain_enrich_sum2.pdf", width = 6, height = 6) # 8DNA
#




# plot MFE       enrichment: group-wise  boxplot/barplot (11RNA, dot is peak)  ----------------------------------------------------------------
ref <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)
table(ref$transcript_type)

options(stringsAsFactors = F)
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
dst <- "GSE71008"
bg.type <- "Background"
methods <- c("piranha/b5_p01","clipper/b5_p05","clam/b5_p005","expeakCNN/b5_d50_p1") # 11RNA

inFile <- paste0(pre,"/output/",dst,"/call_peak_all/",methods,"_11RNA.bed.exp.RNAfold") # 11RNA
# methods <- c("/b5_p05","_localmax_EM2/b5_d05_p05") # 8DNA
# inFile <- paste0(pre,"/output/",dst,"/call_peak_all/domains",methods,"_8DNA.bed.exp.RNAfold") # 8DNA

mfe <- lapply(inFile,conv.mfe)
mfe <- do.call("rbind",mfe)
mfe$type <- "Peak" # "Peak"

inFile2 <- paste0(pre,"/output/",dst,"/call_peak_all/",methods,"_11RNA.bed.exp.shuffle.RNAfold") # 11RNA: flank, shuffle
# inFile2 <- paste0(pre,"/output/",dst,"/call_peak_all/domains",methods,"_8DNA.bed.exp.shuffle.RNAfold") # 8DNA: flank, shuffle
mfe2 <- lapply(inFile2,conv.mfe)
mfe2 <- do.call("rbind",mfe2)
mfe2$type <- bg.type # "Flank" # "Peak"


mfe <- rbind(mfe,mfe2)
mfe$method <- factor(mfe$method,levels = c("Piranha","CLIPper","CLAM","cfPeak")) # 11RNA
rna <- c("rRNA","pri_miRNA","lncRNA", "mRNA", "piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA') # 11RNA
mfe$RNA <- factor(mfe$RNA,levels = c(rna)) # 11RNA
# mfe$method <- factor(mfe$method,levels = c("Piranha","exPeak")) # 8DNA
# dna <- c("intron","promoter", "enhancer","repeats") # 8DNA
# mfe$RNA <- factor(mfe$RNA,levels = c(dna)) # 8DNA

mfe$type <- factor(mfe$type,levels = c("Peak",bg.type))
# table(mfe$mfe<=-20)
table(mfe2$method)
mfe.2 <- mfe[mfe$method=="cfPeak",]
p3 <- 
  # ggplot(mfe,aes(x = method, y = mfe, fill = type)) +
  # geom_boxplot()+
  # facet_grid(method~.,scales = "free_y")
  ggbarplot(data = mfe.2, x = "method", y = "mfe", fill = "type", add="mean_se",position = position_dodge(), short.panel.labs = T) + #  facet.by=c("method"),
  stat_compare_means(
    aes(x=method, y=mfe, group=type), # ref.group="Flank",  #
    # comparisons = list(c("Peak","Flank")),
    label="p.signif",  # p.signif / p.format
    # symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01,1),symbols = c( "***", "**", "*", "ns")),
    method = "wilcox.test", method.args = list(alternative = "greater"),  # greater means ref.group less, 
    label.x.npc = 0.8, label.y.npc = 0.8, # 11RNA
    step.increase = 0.08, vjust = -0.6, # 11RNA
    # label.x.npc = 0.8, label.y.npc = 0.8, # 8DNA
    # bracket.size = 0.5, step.increase = 0.08, vjust = 9.6, # 8DNA
    hide.ns=T,size =0, paired = F # F?
  ) + 
  # geom_bracket(
  #   xmin = "0.5", xmax = "1", y.position = 30,
  #   label = "t-test, p < 0.05"
  # ) + 
  # facet_grid(site~.,scales = "free")+ # facet_grid seem not work well with ggbarplot
  scale_fill_manual(values = c('firebrick','salmon'))+ #c("red","grey90")
  # scale_x_discrete(label=c("Peak Ratio","MIR Ratio"))+
  labs(title="",x="", y = "")+
  # xlim(c(0,1))+   
  # ylim(c(-20,0))+
  theme_minimal() +  # base_size=12
  theme(#aspect.ratio = 1.5,
    #axis.ticks.x=element_blank(),
    #strip.text.y = element_blank(),
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
p3
#cowplot::plot_grid(plotlist = list(p11,p12),rel_heights = c(8,4), nrow = 2, align = "hv") # , axis = "b"
# ggsave(filename = "domain_mfe_sum.pdf", width = 7, height = 7) # 11RNA
ggsave(filename = "domain_mfe_sum_expeak.pdf", width = 3, height = 7) # expeak


fc_data <- mfe.2 %>%
  group_by(method) %>%
  summarise(mean_firebrick = mean(mfe[type == 'Peak']),
            ci_firebrick = calc_mean_cl_normal(mfe[type == 'Peak']),
            mean_salmon = mean(mfe[type == 'Background']),
            ci_salmon = calc_mean_cl_normal(mfe[type == 'Background']),
            FC = mean_firebrick / mean_salmon,
            delta = mean_firebrick-mean_salmon) %>%
  mutate(ci_firebrick_text = sprintf("%.4f ± %.4f", ci_firebrick[1], ci_firebrick[3] - ci_firebrick[1]),
         ci_salmon_text = sprintf("%.4f ± %.4f", ci_salmon[1], ci_salmon[3] - ci_salmon[1]),
         FC_label = sprintf("FC=%.2f\nDelta=%.2f\nPeak: %s\nBackground: %s", FC, delta, ci_firebrick_text, ci_salmon_text)) %>%
  ungroup()
mfe.2 <- left_join(mfe.2, fc_data, by = c("method"))

p12 <- ggplot(mfe.2, aes(x = method, y = mfe, fill = type)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge()) + 
  geom_errorbar(stat = "summary", fun.data = "mean_cl_normal", position = position_dodge(0.9), width = 0.2) + # calc_mean_sem mean_cl_normal
  stat_compare_means(aes(group = type), method = "wilcox.test", method.args = list(alternative = "less"),
                     label = "p.signif", label.x.npc = "right", label.y.npc = "top", hide.ns = TRUE, size = 10) +
  # facet_wrap(~ site2, scales = "free") + 
  scale_fill_manual(values = c('firebrick', 'salmon')) +
  
  # Add FC and 95% CI labels for each method
  geom_text(data = fc_data, aes(x = method, y = 0.9, label = FC_label), 
            inherit.aes = FALSE, size = 5) + 
  
  labs(title = "", x = "", y = "") +
  # ylim(c(0, 1)) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 20, hjust = 1, vjust = 0.5, angle = 90),
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size = 20),
    legend.position = "none",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16)
  )
p12
# Save the plot
# ggsave(filename = "figure/domain_enrich_sum_expeak.pdf", plot = p12, width = 6, height = 8)
ggsave(filename = "domain_mfe_sum_expeak2.pdf", plot = p12, width = 6, height = 13)
#


mfe.df <- mfe %>% 
  # dplyr::select(type,site,sample,RNA,number,ratio,method) %>% 
  # dplyr::distinct( .keep_all = TRUE) %>% 
  dplyr::group_by(type,RNA,method) %>% 
  dplyr::mutate(std=sd(mfe,na.rm=T), mean.mfe=mean(mfe,na.rm=T))

# p33 <- 
#   # ggplot(mfe,aes(x = method, y = mfe, fill = type)) +
#   # geom_boxplot()+
#   # facet_grid(method~.,scales = "free_y")
#   ggbarplot(data = mfe, x = "method", y = "mfe", fill = "type", add="mean_se",position = position_dodge(),facet.by=c("RNA"), short.panel.labs = T) + #  
#   stat_compare_means(
#     aes(x=method, y=mfe, group=type), # ref.group="Flank",  #
#     # comparisons = list(c("Peak","Flank")),
#     label="p.signif",  # ..p.signif../..p.format..
#     # symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01,1),symbols = c( "***", "**", "*", "ns")),
#     method = "wilcox.test", method.args = list(alternative = "greater"),  # greater means ref.group less, 
#     label.x.npc = 0.8, label.y.npc = 0.8, vjust = 0.6, #step.increase = 0.08,
#     hide.ns=T,size =14, paired = F # F?
#   ) +
#   # facet_grid(site~.,scales = "free")+ # facet_grid seem not work well with ggbarplot
#   scale_fill_manual(values = c('firebrick','salmon'))+ #c("red","grey90")
#   # scale_x_discrete(label=c("Peak Ratio","MIR Ratio"))+
#   labs(title="",x="", y = "")+
#   # xlim(c(0,1))+   
#   # ylim(c(-8,10))+
#   theme_minimal() +  # base_size=12
#   theme(#axis.ticks.x=element_blank(),
#     #strip.text.y = element_blank(),
#     strip.text = element_text(size=20),
#     axis.title.x = element_text(size=20),
#     axis.title.y = element_text(size=20),
#     axis.text.x = element_text(size = 20,hjust = 1,vjust = 0.5,angle = 90), #
#     axis.text.y = element_text(size = 20),
#     plot.title = element_text(size=20),
#     # strip.text = element_blank(),
#     legend.position = "none", #c(0.9,0.8),#,#
#     legend.text = element_text(size= 16),
#     legend.title= element_text(size= 16))
#mfe.df[1:3,]
table(mfe.df$RNA)
p33 <- ggplot(mfe.df, aes(x=method, y=mean.mfe, fill=type)) +
  geom_bar(position=position_dodge(), stat = "summary",color="black")+ #
  geom_errorbar(position=position_dodge(.9), aes(ymin=mean.mfe-std, ymax=mean.mfe+std, group=type), width=.2)+
  ggpubr::stat_compare_means(aes(x=method, y=mean.mfe, group=type),
                             label = "p.signif", # ..p.signif../..p.format..
                             # ref.group = "Peak", comparisons,
                             # symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01,1),symbols = c( "***", "**", "*", "ns")),
                             method = "wilcox.test", paired = T, method.args = list(alternative = "less"),
                             label.x.npc = 0.4, label.y.npc = 0.85, step.increase = 0.08,
                             hide.ns=T,size =10
  )+
  scale_fill_manual(values = c('firebrick','salmon'))+ #c("red","grey90")
  # scale_x_discrete(label=c("Peak Ratio","MIR Ratio"))+
  labs(title="",x="", y = "")+
  # facet_grid(method~.)+
  # ylim(c(0,1))+
  # xlim(c(0,1))+
  theme_minimal() +  # base_size=12
  theme(#axis.ticks.x=element_blank(),
    #strip.text.y = element_blank(),
    #strip.text.x = element_text(face="bold",family="arial",size=20),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size = 20,hjust = 1,vjust = 0.5,angle = 90), #
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size=20),
    strip.text = element_text(size = 20),
    legend.position = "right", #c(0.9,0.8),#,#
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16)) +
  facet_grid(RNA~.,scales = "free")
p33
ggsave(filename = "domain_mfe_sum_by_RNA_type.pdf", width = 7, height = 22)#p1+p2+patchwork::guide_area+patchwork::plot_layout()
#no need for 8DNA !!!!






# plot MFE       enrichment: group-wise  boxplot/barplot (11RNA, dot is peak, expeak only OR expeak gold)  ----------------------------------------------------------------
ref <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)
table(ref$transcript_type)

options(stringsAsFactors = F)
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
dst <- "GSE71008"
bg.type <- "Background"
methods <- c("expeakCNN/b5_d50_p1") # "piranha/b5_p01","clipper/b5_p05","clam/b5_p005",

inFile <- paste0(pre,"/output/",dst,"/call_peak_all/",methods,"_11RNA.bed.exp.RNAfold") # 11RNA
# methods <- c("/b5_p05","_localmax_EM2/b5_d05_p05") # 8DNA
# inFile <- paste0(pre,"/output/",dst,"/call_peak_all/domains",methods,"_8DNA.bed.exp.RNAfold") # 8DNA 
inFile.2 <- paste0(unlist(sapply(strsplit(inFile,".bed",fixed=T),"[",1)), "_expeakCNNOnly.bed") # _expeakCNNOnly.bed, _expeakCNNGold.bed
#mfe <- lapply(inFile,conv.mfe)
# mfe <- lapply(1:length(inFile),conv.mfe)
mfe <- lapply(1:length(inFile),FUN = function(i) {conv.mfe.expeakOnly(x=inFile[i],y=inFile.2[i])} )
mfe <- do.call("rbind",mfe)
mfe$type <- "Peak" # "Peak"

inFile2 <- paste0(pre,"/output/",dst,"/call_peak_all/",methods,"_11RNA.bed.exp.shuffle.RNAfold") # 11RNA: flank, shuffle
inFile2.2 <- inFile.2
# mfe2 <- lapply(inFile2,conv.mfe)
mfe2 <- lapply(1:length(inFile2),FUN = function(i) {conv.mfe.expeakOnly(x=inFile2[i],y=inFile2.2[i])} )

mfe2 <- do.call("rbind",mfe2)
mfe2$type <- bg.type # "Flank" # "Peak"


mfe <- rbind(mfe,mfe2)
# mfe$method <- factor(mfe$method,levels = c("Piranha","CLIPper","CLAM","exPeak")) # 11RNA
rna <- c("rRNA","pri_miRNA","lncRNA", "mRNA", "piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA') # 11RNA
mfe$RNA <- factor(mfe$RNA,levels = c(rna)) # 11RNA
# mfe$method <- factor(mfe$method,levels = c("Piranha","exPeak")) # 8DNA
# dna <- c("intron","promoter", "enhancer","repeats") # 8DNA
# mfe$RNA <- factor(mfe$RNA,levels = c(dna)) # 8DNA

mfe$type <- factor(mfe$type,levels = c("Peak",bg.type))
# table(mfe$mfe<=-20)
table(mfe$RNA)
p3 <- 
  # ggplot(mfe,aes(x = method, y = mfe, fill = type)) +
  # geom_boxplot()+
  # facet_grid(method~.,scales = "free_y")
  ggbarplot(data = mfe, x = "method", y = "mfe", fill = "type", add="mean_se",position = position_dodge(), short.panel.labs = T) + #  facet.by=c("method"),
  stat_compare_means(
    aes(x=method, y=mfe, group=type), # ref.group="Flank",  #
    # comparisons = list(c("Peak","Flank")),
    label="p.signif",  # p.signif / p.format
    # symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01,1),symbols = c( "***", "**", "*", "ns")),
    method = "wilcox.test", method.args = list(alternative = "greater"),  # greater means ref.group less, 
    label.x.npc = 0.8, label.y.npc = 0.8, # 11RNA
    step.increase = 0.08, vjust = -0.6, # 11RNA
    # label.x.npc = 0.8, label.y.npc = 0.8, # 8DNA
    # bracket.size = 0.5, step.increase = 0.08, vjust = 9.6, # 8DNA
    hide.ns=T,size =0, paired = F # F?
  ) + 
  # geom_bracket(
  #   xmin = "0.5", xmax = "1", y.position = 30,
  #   label = "t-test, p < 0.05"
  # ) + 
  # facet_grid(site~.,scales = "free")+ # facet_grid seem not work well with ggbarplot
  scale_fill_manual(values = c('firebrick','salmon'))+ #c("red","grey90")
  # scale_x_discrete(label=c("Peak Ratio","MIR Ratio"))+
  labs(title="",x="", y = "")+
  # xlim(c(0,1))+   
  # ylim(c(-20,0))+
  theme_minimal() +  # base_size=12
  theme(aspect.ratio = 3,
        #axis.ticks.x=element_blank(),
        #strip.text.y = element_blank(),
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
p3
#cowplot::plot_grid(plotlist = list(p11,p12),rel_heights = c(8,4), nrow = 2, align = "hv") # , axis = "b"
ggsave(filename = "domain_mfe_sum.pdf", width = 3, height = 7) # 11RNA
# ggsave(filename = "domain_mfe_sum.pdf", width = 4, height = 7) # 8DNA


mfe.df <- mfe %>% 
  # dplyr::select(type,site,sample,RNA,number,ratio,method) %>% 
  # dplyr::distinct( .keep_all = TRUE) %>% 
  dplyr::group_by(type,RNA,method) %>% 
  dplyr::mutate(std=sd(mfe,na.rm=T), mean.mfe=mean(mfe,na.rm=T))

# p33 <- 
#   # ggplot(mfe,aes(x = method, y = mfe, fill = type)) +
#   # geom_boxplot()+
#   # facet_grid(method~.,scales = "free_y")
#   ggbarplot(data = mfe, x = "method", y = "mfe", fill = "type", add="mean_se",position = position_dodge(),facet.by=c("RNA"), short.panel.labs = T) + #  
#   stat_compare_means(
#     aes(x=method, y=mfe, group=type), # ref.group="Flank",  #
#     # comparisons = list(c("Peak","Flank")),
#     label="p.signif",  # ..p.signif../..p.format..
#     # symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01,1),symbols = c( "***", "**", "*", "ns")),
#     method = "wilcox.test", method.args = list(alternative = "greater"),  # greater means ref.group less, 
#     label.x.npc = 0.8, label.y.npc = 0.8, vjust = 0.6, #step.increase = 0.08,
#     hide.ns=T,size =14, paired = F # F?
#   ) +
#   # facet_grid(site~.,scales = "free")+ # facet_grid seem not work well with ggbarplot
#   scale_fill_manual(values = c('firebrick','salmon'))+ #c("red","grey90")
#   # scale_x_discrete(label=c("Peak Ratio","MIR Ratio"))+
#   labs(title="",x="", y = "")+
#   # xlim(c(0,1))+   
#   # ylim(c(-8,10))+
#   theme_minimal() +  # base_size=12
#   theme(#axis.ticks.x=element_blank(),
#     #strip.text.y = element_blank(),
#     strip.text = element_text(size=20),
#     axis.title.x = element_text(size=20),
#     axis.title.y = element_text(size=20),
#     axis.text.x = element_text(size = 20,hjust = 1,vjust = 0.5,angle = 90), #
#     axis.text.y = element_text(size = 20),
#     plot.title = element_text(size=20),
#     # strip.text = element_blank(),
#     legend.position = "none", #c(0.9,0.8),#,#
#     legend.text = element_text(size= 16),
#     legend.title= element_text(size= 16))
#mfe.df[1:3,]
table(mfe.df$RNA)
p33 <- ggplot(mfe.df, aes(x=method, y=mean.mfe, fill=type)) +
  geom_bar(position=position_dodge(), stat = "summary",color="black")+ #
  geom_errorbar(position=position_dodge(.9), aes(ymin=mean.mfe-std, ymax=mean.mfe+std, group=type), width=.2)+
  ggpubr::stat_compare_means(aes(x=method, y=mean.mfe, group=type),
                             label = "p.signif", # ..p.signif../..p.format..
                             # ref.group = "Peak", comparisons,
                             # symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01,1),symbols = c( "***", "**", "*", "ns")),
                             method = "wilcox.test", paired = T, method.args = list(alternative = "less"),
                             label.x.npc = 0.4, label.y.npc = 0.85, step.increase = 0.08,
                             hide.ns=T,size =10
  )+
  scale_fill_manual(values = c('firebrick','salmon'))+ #c("red","grey90")
  # scale_x_discrete(label=c("Peak Ratio","MIR Ratio"))+
  labs(title="",x="", y = "")+
  # facet_grid(method~.)+
  # ylim(c(0,1))+
  # xlim(c(0,1))+
  theme_minimal() +  # base_size=12
  theme(#axis.ticks.x=element_blank(),
    #strip.text.y = element_blank(),
    #strip.text.x = element_text(face="bold",family="arial",size=20),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size = 20,hjust = 1,vjust = 0.5,angle = 90), #
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size=20),
    strip.text = element_text(size = 20),
    legend.position = "right", #c(0.9,0.8),#,#
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16)) +
  facet_grid(RNA~.,scales = "free")
p33
ggsave(filename = "domain_mfe_sum_by_RNA_type.pdf", width = 7, height = 22)#p1+p2+patchwork::guide_area+patchwork::plot_layout()
#no need for 8DNA !!!!




# plot MFE       enrichment: group-wise  boxplot/barplot (8DNA, dot is peak)  ----------------------------------------------------------------
ref <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)
table(ref$transcript_type)

options(stringsAsFactors = F)
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
dst <- "GSE71008"
bg.type <- "Background"
# methods <- c("/b5_p05","_clipper/b5_p05","_clam/b5_p05","_localmax_EM2/b5_d05_p05") # 11RNA
# inFile <- paste0(pre,"/output/",dst,"/call_peak_all/domains",methods,"_11RNA.bed.RNAfold") # 11RNA
methods <- c("expeakCNN/b5_d50_p1") # 8DNA ("piranha/b5_p05",)
inFile <- paste0(pre,"/output/",dst,"/call_peak_all/",methods,"_8DNA_gn.bed.exp.RNAfold") # 8DNA

mfe <- lapply(inFile,conv.mfe)
mfe <- do.call("rbind",mfe)
mfe$type <- "Peak" # "Peak"

# inFile2 <- paste0(pre,"/output/",dst,"/call_peak_all/domains",methods,"_11RNA.bed.shuffle.RNAfold") # 11RNA: flank, shuffle
inFile2 <- paste0(pre,"/output/",dst,"/call_peak_all/",methods,"_8DNA_gn.bed.exp.shuffle.RNAfold") # 8DNA: flank, shuffle
mfe2 <- lapply(inFile2,conv.mfe)
mfe2 <- do.call("rbind",mfe2)
mfe2$type <- bg.type # "Flank" # "Peak"


mfe <- rbind(mfe,mfe2)
# mfe$method <- factor(mfe$method,levels = c("Piranha","exPeak")) # 8DNA
dna <- c("intron","promoter", "enhancer","repeats") # 8DNA
mfe$RNA <- factor(mfe$RNA,levels = c(dna)) # 8DNA

mfe$type <- factor(mfe$type,levels = c("Peak",bg.type))
# table(mfe$mfe<=-20)
# table(mfe2$RNA)
mfe2 <- mfe[mfe$method=="exPeak",]
p3 <- 
  # ggplot(mfe,aes(x = method, y = mfe, fill = type)) +
  # geom_boxplot()+
  # facet_grid(method~.,scales = "free_y")
  ggbarplot(data = mfe2, x = "type", y = "mfe", fill = "type", add="mean_se",position = position_dodge(), short.panel.labs = T) + #  facet.by=c("method"),
  stat_compare_means(
    aes(x=type, y=mfe, group=type), # ref.group="Flank",  #
    # comparisons = list(c("Peak","Flank")),
    label="p.format",  # p.signif / p.format
    # symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01,1),symbols = c( "***", "**", "*", "ns")),
    method = "wilcox.test", method.args = list(alternative = "greater"),  # greater means ref.group less, 
    # label.x.npc = 0.8, label.y.npc = 0.8, # 11RNA 
    # step.increase = 0.08, vjust = -0.6, # 11RNA 
    label.x.npc = 0.8, label.y.npc = 0.8, # 8DNA
    bracket.size = 0.5, step.increase = 0.08, vjust = 6.6, hjust = 1, # 8DNA
    hide.ns=T,size =6, paired = F # F?
  ) + 
  # geom_bracket(
  #   xmin = "0.5", xmax = "1", y.position = 30,
  #   label = "t-test, p < 0.05"
  # ) + 
  # facet_grid(site~.,scales = "free")+ # facet_grid seem not work well with ggbarplot
  scale_fill_manual(values = c('firebrick','salmon'))+ #c("red","grey90")
  # scale_x_discrete(label=c("Peak Ratio","MIR Ratio"))+
  labs(title="",x="", y = "")+
  # xlim(c(0,1))+   
  # ylim(c(-20,0))+
  theme_minimal() +  # base_size=12
  theme(aspect.ratio = 2,
    #axis.ticks.x=element_blank(),
    #strip.text.y = element_blank(),
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
p3
#cowplot::plot_grid(plotlist = list(p11,p12),rel_heights = c(8,4), nrow = 2, align = "hv") # , axis = "b"
# ggsave(filename = "domain_mfe_sum.pdf", width = 6, height = 8) # 11RNA
ggsave(filename = "domain_mfe_sum.pdf", width = 3, height = 8) # 8DNA





# plot intersect enrichment: overview freq hist + upset (11RNA, Fig5) -------------------------------------------------
### get tx type
ref <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)
#ref <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length.txt",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)
table(ref$transcript_type)

### read (shuffle) domain bedtools intersect
dst="GSE71008" # "GSE71008","GSE110381", "GSE94533"
methods <- "expeakCNN" #
x <- (paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_all/",methods,"/b5_d50_p1_11RNA.intersect.bed")) # 11RNA
res.df <- data.table::fread(x,data.table = F,sep = "\t",header = T,check.names = F,stringsAsFactors = F)
# res.df <- res.df[order(res.df$site_overlap_ratio_product,res.df$site_overlap_base,decreasing = T),]  # dedup ?
# res.df <- res.df[!duplicated(paste0(res.df$name,res.df$site)),]
res.df$RNA <- ref$transcript_type[match(res.df$chr,ref$transcript_id)]
table(is.na(res.df$RNA))
table(is.na(res.df$site_overlap_ratio_product))
res.df$site_overlap_ratio_product[is.na(res.df$site_overlap_ratio_product)] <- 0
res.df$peakID <- unlist(sapply(strsplit(res.df$name,"--"),"[",1))
#table(res.df$site)
res.df[1:3,]


x2 <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_all/expeakCNN/b5_d50_p1_11RNA.bed.exp.fa.csv")
#x2 <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008/call_peak_all/domains_localmax_EM/b5_d05_p01.bed.shuffle.fa.csv"
#x2 <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE133684/call_peak_dedupByPos/domains_localmax_EM/b5_d05_p01.bed.fa.csv"
# x2 <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE133684/call_peak_dedupByPos/domains_localmax_EM/b5_d05_p01.bed.shuffle.fa.csv"
res.df2 <- data.table::fread(x2,data.table = F,sep = ",",header = F,check.names = F,stringsAsFactors = F)
res.df2$V1 <- gsub(">","",res.df2$V1)
colnames(res.df2) <- c("name2","site_overlap_ratio_product","mfe") # site_overlap_ratio_product: structure_p
res.df2$site <- "structure"
res.df2[1:3,]
res.df2$peakID <- unlist(sapply(strsplit(res.df2$name2,"\\(|:",perl = T,),"[",1)) # (),::; 
res.df2$chr <- res.df$chr[match(res.df2$peakID,res.df$peakID)]
res.df2$RNA <- ref$transcript_type[match(res.df2$chr,ref$transcript_id)]
table((res.df2$RNA))
res.df2$site_overlap_ratio_product[is.na(res.df2$site_overlap_ratio_product)] <- 1
#table(res.df2$site_overlap_ratio_product<=0.1)

num <- nrow(res.df) # total peak num
#dplyr::full_join(x = , y = )
res.df3 <- rbind(res.df[,c("chr","site_overlap_ratio_product","site","RNA","peakID")], res.df2[,c("chr","site_overlap_ratio_product","site","RNA","peakID")])
#table(res.df3$RNA)

library(ggplot2)
res.df3$RNA <- gsub("_rev|_for|\\.for|\\.rev","",res.df3$RNA)
rna <- c("rRNA","pri_miRNA","lncRNA", "mRNA", "piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA')
#dna <- c("intron","promoter", "enhancer","repeats") # 
res.df3$RNA <- factor(res.df3$RNA,levels = c(rna)) #dna
table(is.na(res.df3$RNA))

res.df3 <- res.df3[res.df3$site!="iM",]
res.df3$site <- factor(res.df3$site,levels = c("AGO2","otherRBPs","RBPhotspot","EV","G4","structure"))
# res.df3$pos <- paste0(res.df3$chr,":",res.df3$start,"-",res.df3$end) 
res.df3$site2 <- res.df3$site
res.df3$site2 <- gsub("AGO2|otherRBPs|RBPhotspot","RBP",res.df3$site2,perl = T)
#res.df3$site2[grepl("AGO2|otherRBPs|RBPhotspot",as.character(res.df3$site),perl = T)] <- "RBP"
table(res.df3$site2)

res.df3$site2 <- factor(res.df3$site2, levels = c("RBP","EV","G4","structure"))


# peak.gold <- data.table::fread(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_all/",methods,"/b5_d50_p1_11RNA_expeakGold.bed"),
#                             data.table = F,sep = "\t",header = F,check.names = F,stringsAsFactors = F)
# table(res.df3$peakID %in% peak.gold$V4) # ~50%
# res.df3 <- res.df3[res.df3$peakID %in% peak.gold$V4,] # filter gold overlapped peak


ggplot(res.df3, aes(x=site_overlap_ratio_product, color=site2, fill=site2)) +
  geom_histogram(alpha=0.5)+
  # geom_histogram(aes(y=..density..), position="identity", alpha=0.5)+
  # geom_density(alpha=0.1)+
  # geom_vline(data=mu, aes(xintercept=grp.mean, color=method),linetype="dashed")+
  # scale_color_manual(values=c("firebrick", "steelblue4", "grey30"))+
  # scale_fill_manual(values=c("firebrick", "steelblue4", "grey30"))+
  # geom_vline(xintercept = 0,color="grey30",linetype="dashed")+
  # geom_label(data=mu, aes(label=ratio, color=method, fill=method )) +
  #annotate(x = 0.5, y=0.5)+
  # geom_text()+
  ggsci::scale_fill_nejm()+
  ggsci::scale_color_nejm()+
  labs(title="",x="site_overlap_ratio_product", y = "Domain Number")+
  facet_grid(RNA~site2, scales = "free_y")+
  # ylim(c(0,1))+xlim(c(-10,10))+
  theme_classic(base_size=12) + 
  theme(#axis.ticks.x=element_blank(),  
    #strip.text.y = element_blank(),
    #strip.text.x = element_text(face="bold",family="arial",size=20),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size = 20,hjust = 0.5,vjust = 0.5), # ,angle = 90
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size=20),
    strip.text = element_text(size = 20),
    legend.position = "right",# 
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16))
ggsave("intersect.pdf", width = 22, height = 20)

res.tbl <- as_tibble(res.df3) %>% 
  dplyr::mutate(RBP_structure=dplyr::case_when( site2=="RBP" & (site_overlap_ratio_product>=0.01) ~ T,   #site2
                                                # site2=="otherRBPs" & (site_overlap_ratio_product>=0.01) ~ T,
                                                # site2=="RBPhotspot" & (site_overlap_ratio_product>=0.01) ~ T,
                                                site2=="EV" & (site_overlap_ratio_product>=0.01) ~ T,
                                                site2=="G4" & (site_overlap_ratio_product>=0.01) ~ T,
                                                site2=="structure" & (site_overlap_ratio_product<=0.1) ~ T,
                                                TRUE ~ F
                                                )
                )


mat <- res.tbl %>% 
  dplyr::select(RBP_structure,site2,peakID) %>% 
  dplyr::group_by(site2,peakID) %>% 
  dplyr::mutate(RBP_structure=sum(RBP_structure)) %>% 
  dplyr::distinct(RBP_structure,site2,peakID,.keep_all = TRUE) %>% 
  tidyr::pivot_wider(id_cols = peakID, names_from = site2, values_from = RBP_structure)

#library("ggvenn") # not suport more than 5
#create Venn diagram and display all sets
# ggvenn(as.data.frame(mat),set_name_size=4,text_size = 5,show_percentage=F, text_color = "white") + ggsci::scale_fill_nejm() + 
#   theme(aspect.ratio = 0.8)
mat <- as.data.frame(mat)
rownames(mat) <- mat$peakID
mat$peakID <- NULL
mat[1:3,]
# table(is.na(mat)) # some structure is na (seem filtered before rnafold)
mat[is.na(mat)] <- 0
mat <- as.matrix(mat)
mat[mat>=1] <- 1
mat[mat==0] <- 0
mat <- as.data.frame(mat)
# tmp <- mat[mat$G4iM==1 & rowSums(mat)==1,]
# tmp <- mat[is.na(mat[,5]),]
# table(is.na(mat[,5]))
colnames(mat) <- c("RBP bind","G4 Struc.","EV sort","2nd Struc.")
pdf("upset.pdf",width = 6,height = 3)
UpSetR::upset(data = mat, sets = c("RBP bind","EV sort","G4 Struc.","2nd Struc."), nsets = 4, keep.order = T,
              matrix.color = "salmon",set_size.numbers_size = T,
              sets.bar.color = pal_nejm_adaptive()(4)#colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(4)
              )
dev.off()

#get G4 or EV only eg.
# tmp <- mat[mat$G4iM==1 & rowSums(mat)==1,]  # G4iM, EV
# x <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008/call_peak_all/domains_localmax_EM2/b5_d05_p05_11RNA.bed"
# res.df3 <- data.table::fread(x,data.table = F,sep = "\t",header = F,check.names = F,stringsAsFactors = F)
# res.df3 <- res.df3[res.df3$V4 %in% rownames(tmp),]
# res.df3$enst <- unlist(sapply(strsplit(res.df3$V1,"_"),"[",1))
# library(clusterProfiler)
# library(org.Hs.eg.db)
# convert.df <- bitr(res.df3$enst , fromType="ENSEMBLTRANS", toType=c("SYMBOL","ENTREZID"), OrgDb="org.Hs.eg.db")
# #seem no RP mRNA related
# res.df3$gene.symbol<- convert.df$SYMBOL[match(res.df3$enst,convert.df$ENSEMBLTRANS)]

# #colnames(mat)
# nrow(mat)
# colSums(mat)
#small,domain
#AGO2  otherRBPs RBPhotspot       G4iM  structure 
#396        563        342         44        118 
#small,bg
#AGO2  otherRBPs RBPhotspot       G4iM  structure 
#350        498        342         22        188 
#long,domain,4659
#AGO2  otherRBPs RBPhotspot       G4iM  structure 
#99        152        171         53        274 
#AGO2  otherRBPs RBPhotspot       G4iM  structure 
#99        152        171         53        367 
#long,bg
#AGO2  otherRBPs RBPhotspot       G4iM  structure 
#90        154        169         60        525 
#AGO2  otherRBPs RBPhotspot       G4iM  structure 
#90        154        169         60        644

# ## hyper geometric test
# domain.T <- 239
# bg.T <- 302
# domain.F <- nrow(mat)-domain.T
# bg.F <- nrow(mat)-bg.T
# df.hyper <- data.frame(row.names = c("Domain","Background"),withRBP=c(domain.T,bg.T),noRBP=c(domain.F,bg.F))
# test <- fisher.test(df.hyper,)
# test$p.value
# #df.hyper$Type <- rownames(df.hyper)



mat.list <- list(rownames(mat)[mat$`RBP bind`==1],rownames(mat)[mat$`EV sort`==1],rownames(mat)[mat$`G4 Struc.`==1],rownames(mat)[mat$`2nd Struc.`==1])
names(mat.list) <- colnames(mat)
library("ggvenn")
#names(mat.list) 
ggvenn(data = mat.list, columns = c("RBP bind","EV sort","G4 Struc.","2nd Struc.") , set_name_size=6, text_size = 8,show_percentage=F, stroke_color = "white", text_color = "white") + ggsci::scale_fill_nejm() + theme(aspect.ratio = 0.8)
ggsave("tmpVenn.pdf",width = 10,height = 5, limitsize=F)

table(rowSums(mat)>0)
sum(rowSums(mat)>0)/nrow(mat)
#GSE110381.raw: 0.2723504
#GSE110381.gold: 0.5998914
#GSE94533.raw:0.5515424
#GSE94533.gold:0.5923182


#RBP:
sum(mat$`RBP bind`)/nrow(mat)
#GSE71008.raw: 38.1%

# 
# ## select 4 overlap peak id
# peak.id <- rownames(mat)[rowSums(mat)==4] 
# #peak_1837 (MIR532-(3p)-201)
# #tx: ENST00000385025_____1 56 76 peak_1837 50  +
# #gn: chrX    50003203        50003223        peak_1837       50      +       0       0       0       1       20      0
# #>peak_1837::ENST00000385025_____1:36-91(+)   exp 15bp
# #ACCGTTGGCATCTTAATTACCCTCCCACACCCAAGGCTTGCAGAAGAGCGAGCCT
# #                   CCCTCCCACACCCAAGGCTTG
# #15,15
# #CGACTTGCTTTCTCTCCTCCATGCCTTGAGTGTAGGACCGTTGGCATCTTAATTACCCTCCCACACCCAAGGCTTGCAGAAGAGCGAGCCT
# #full-len: 91nt
# 
# 
# ## select all RBP peak id 
# #peak.id <- rownames(mat)[mat$`RBP bind`==1 & mat$`G4 Struc.`==1 & mat$`EV sort`==0 & mat$`2nd Struc.`==1]
# peak.id <- rownames(mat)[mat$`RBP bind`==1]
# 
# y <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008/call_peak_all/expeakCNN/b5_d50_p1_11RNA.bed"
# peak <- read.table(y)
# peak <- peak[peak$V4 %in% peak.id,]
# write.table(peak,"tmp/RBP_peak.bed",quote = F,sep = "\t",row.names = F,col.names = F)
# 
# #bash
# pre="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
# /usr/bin/Rscript  $pre/exSeek-dev/bin/peak_intersect_G4iM_RBP_EV.R --coord tx --RNAtype RNA \
#   -i tmp/RBP_peak.bed \
#   -o tmp/RBP_peak.bed.intersect


# plot RBP distribution each RNA types (11RNA , Fig5, adapted from Taiwei) ----
suppressPackageStartupMessages(library(dplyr))
conflict_prefer("desc", "dplyr")

# get G4 RBP ref
g4.ref <- data.table::fread("/lulabdata/baopengfei/shared_reference/structure/quadratlas/quadratlas_RG4BPs.txt",skip = 2,sep = "\t")
#warning is normal for this is a mix tab
g4.ref[1:3,1:3]
length(unique(g4.ref$genename))
#table(g4.ref$isKnownRBP)
#table(unique(as.character(df$RBP)) %in% unique(as.character(g4.ref$genename)))

# get tx ref
tx_anno <- ref


# read RBP intersect
#bash
x <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/tmp/RBP_peak.bed.intersect"
RBP_enrich <- data.table::fread(x,data.table = F,sep = "\t",header = T,check.names = F,stringsAsFactors = F)
RBP_enrich$RNA <- ref$transcript_type[match(RBP_enrich$chr,ref$transcript_id)]
table(is.na(RBP_enrich$RNA))
table(is.na(RBP_enrich$site_overlap_ratio_product))
RBP_enrich$site_overlap_ratio_product[is.na(RBP_enrich$site_overlap_ratio_product)] <- 0
RBP_enrich$peakID <- unlist(sapply(strsplit(RBP_enrich$name,"--"),"[",1))
table(RBP_enrich$site)
RBP_enrich <- RBP_enrich[RBP_enrich$site %in% c("AGO2","otherRBPs"),]
RBP_enrich <- RBP_enrich[RBP_enrich$site_overlap_ratio_product>=0.01,]
RBP_enrich$RBP <- unlist(sapply(strsplit(RBP_enrich$site_name,"-"),"[",1))
table(is.na(RBP_enrich$RBP))
RBP_enrich <- RBP_enrich[,c("chr","site_name","peakID","site_overlap_base","site_overlap_ratio_product","RBP")]
colnames(RBP_enrich) <- paste0("V",1:6)
#not need keep max overlap_ratio_product each peak (already done in intersect_RBP_G4 script)

#RBP_enrich <- data.table::fread("/BioII/lulab_b/wangtaiwei/peak_calling/RBP_enrich/final.bed", sep = "\t")
RBP_enrich[1:2,]
#V1                                               V2       V3 V4      V5   V6
#1: ENST00000384997_____3 AGO2-chr1:1167120-1167140|chr1|1167120|1167140|+ peak_928 16 0.64000 AGO2
#2: ENST00000384997_____3 AGO2-chr1:1167120-1167149|chr1|1167120|1167149|+ peak_928 20 0.68966 AGO2


RBP_enrich <- merge(RBP_enrich, tx_anno, by.x = "V1", by.y = "transcript_id") %>%
  transmute(transcript_id=V1,position=V2,peak_id=V3,overlap_ratio=V5,RBP_type=V6,transcript_type=transcript_type)


RNApeak_RBP_enrich <- data.frame(unclass(table(RBP_enrich$RBP_type,RBP_enrich$transcript_type))) #unclass() to remove class object 
RNApeak_RBP_enrich <- RNApeak_RBP_enrich[,c("rRNA","pri_miRNA","lncRNA","mRNA","piRNA","snoRNA","snRNA","srpRNA","tRNA","tucpRNA","Y_RNA")]

RNApeak_RBP_percentage <- apply(RNApeak_RBP_enrich,2 , function(x) {x*100/sum(x,na.rm = T)} )
RBP_rank_desc <- apply(RNApeak_RBP_enrich,2,function(x) row.names(RNApeak_RBP_enrich)[order(desc(x))])

#prepare percentage barplot matrix
barplot_mat <- RNApeak_RBP_percentage[unique(c(RBP_rank_desc[1:5,])),]
for (i in colnames(barplot_mat)) barplot_mat[-which(row.names(barplot_mat) %in% RBP_rank_desc[1:5,i]),i] = 0 #other RNA type's top RBP set to 0
barplot_mat <- rbind(barplot_mat, 100 - as.vector(apply(barplot_mat, 2 , sum)))
row.names(barplot_mat)[nrow(barplot_mat)] <- "Others"

# library(RColorBrewer)
# coul <- brewer.pal(8, "Pastel2")
library(paletteer)
coul <- paletteer_d("ggthemes::Tableau_20")
#barplot(barplot_mat, col = coul,border="white")+geom_text(aes(label = row.names(barplot_mat)))

# ggplot
library(ggplot2)
df <- melt(barplot_mat, varnames = c("RBP","RNA"), value.name = c("percentage"))
df <- df[which(df$percentage != 0),]
conflict_prefer("arrange", "dplyr")
df2 <- df %>% 
  dplyr::filter(RBP!="Others") %>% 
  dplyr::group_by(RBP) %>% 
  dplyr::summarise(mean.ratio=mean(percentage)) %>% 
  arrange(desc(mean.ratio))
df$RBP <- factor(df$RBP,levels = c(unique(as.character(df2$RBP)),"Others")) # rank RBP by all RNA mean.perc

df3 <- df # rank RNA by all RBP sequentially
#df[df$RBP %in% unique(as.character(df2$RBP))[1:5],] # rank RNA by top5 RBP sequentially
#df3$RBP <- factor(df3$RBP,levels = unique(as.character(df2$RBP))[1:5])
table(df3$RBP)
df3 <- df3[order(df3$RBP,-df3$percentage,decreasing=F),]
df$RNA <- factor(df$RNA,levels = c(unique(as.character(df3$RNA))))

unique(as.character(df$RBP))[unique(as.character(df$RBP)) %in% unique(as.character(g4.ref$genename))] # 9
df$G4RBP <- "N"
df$G4RBP[(as.character(df$RBP)) %in% unique(as.character(g4.ref$genename))] <- "Y"
df$G4RBP <- factor(df$G4RBP,levels = c("N","Y"))
table(df$G4RBP)

ggplot(df, aes(x = RNA, y = percentage, RBP = RBP))+
  geom_col(aes(fill = RBP), position = position_stack(reverse = T))+ # ,color=G4RBP
  geom_text(aes(label = RBP), position = position_stack(vjust = 0.5,reverse = T), size = 4)+ # ,angle = -90
  scale_fill_d3_adaptive()+
  scale_color_manual(values = c("grey90","black"))+
  ylab("Percentage") + 
  # coord_flip() +
  # scale_fill_manual(values = ggthemes::tableau_color_pal(palette = "Tableau 20",type="ordered-sequential")(18)) +
  theme_minimal() +  # base_size=12
  theme(
    aspect.ratio = 0.5,
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    #axis.ticks.x=element_blank(),
    #strip.text.y = element_blank(),
    strip.text = element_text(size=26),
    axis.title.x = element_text(size=26),
    axis.title.y = element_text(size=26),
    axis.text.x = element_text(size = 26,hjust = 1,vjust = 0.5,angle = 90), #hjust = 1,vjust = 0.5
    axis.text.y = element_text(size = 26),
    plot.title = element_text(size=26),
    # strip.text = element_blank(),
    legend.position = "right", #c(0.9,0.8),#,#
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 26))
ggsave("RBP_dist.pdf",width = 13,height = 10)

#

# plot intersect enrichment: overview freq hist + upset (8DNA) -------------------------------------------------
### get tx type
ref <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)
#ref <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length.txt",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)
table(ref$transcript_type)

### read (shuffle) domain bedtools intersect
dst="GSE71008" # "GSE71008","GSE110381"
methods <- "expeak" #
#x <- (paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_all/",methods,"/b5_d50_p1_11RNA.intersect.bed")) # 11RNA
x <- (paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_all/",methods,"/b5_d50_p1_8DNA_gn.intersect.bed")) # 8DNA

res.df <- data.table::fread(x,data.table = F,sep = "\t",header = T,check.names = F,stringsAsFactors = F)
# res.df <- res.df[order(res.df$site_overlap_ratio_product,res.df$site_overlap_base,decreasing = T),]  # dedup ?
# res.df <- res.df[!duplicated(paste0(res.df$name,res.df$site)),]
res.df$RNA <- ref$transcript_type[match(res.df$chr,ref$transcript_id)]
table(is.na(res.df$RNA))
table(is.na(res.df$site_overlap_ratio_product))
res.df$site_overlap_ratio_product[is.na(res.df$site_overlap_ratio_product)] <- 0
res.df$peakID <- unlist(sapply(strsplit(res.df$name,"--"),"[",1))
#table(res.df$site)
res.df[1:3,]


x2 <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008/call_peak_all/expeak/b5_d50_p1_8DNA_gn.bed.exp.fa.csv"
res.df2 <- data.table::fread(x2,data.table = F,sep = ",",header = F,check.names = F,stringsAsFactors = F)
res.df2$V1 <- gsub(">","",res.df2$V1)
colnames(res.df2) <- c("name2","site_overlap_ratio_product","mfe") # site_overlap_ratio_product: structure_p
res.df2$site <- "structure"
res.df2[1:3,]
res.df2$peakID <- unlist(sapply(strsplit(res.df2$name,"::",fixed = T),"[",1)) # ()
res.df2$chr <- res.df$chr[match(res.df2$peakID,res.df$peakID)]
res.df2$RNA <- ref$transcript_type[match(res.df2$chr,ref$transcript_id)]
table((res.df2$RNA))
res.df2$site_overlap_ratio_product[is.na(res.df2$site_overlap_ratio_product)] <- 1
#table(res.df2$site_overlap_ratio_product<=0.1)

num <- nrow(res.df) # total peak num
#dplyr::full_join(x = , y = )
res.df3 <- rbind(res.df[,c("chr","site_overlap_ratio_product","site","RNA","peakID")], res.df2[,c("chr","site_overlap_ratio_product","site","RNA","peakID")])
#table(res.df3$RNA)

library(ggplot2)
# res.df3$RNA <- gsub("_rev|_for|\\.for|\\.rev","",res.df3$RNA)
# rna <- c("rRNA","pri_miRNA","lncRNA", "mRNA", "piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA')
# #dna <- c("intron","promoter", "enhancer","repeats") # 
# res.df3$RNA <- factor(res.df3$RNA,levels = c(rna)) #dna
# table(is.na(res.df3$RNA))

res.df3 <- res.df3[res.df3$site!="iM",]
res.df3$site <- factor(res.df3$site,levels = c("AGO2","otherRBPs","RBPhotspot","EV","G4","structure"))
# res.df3$pos <- paste0(res.df3$chr,":",res.df3$start,"-",res.df3$end) 
res.df3$site2 <- res.df3$site
res.df3$site2 <- gsub("AGO2|otherRBPs|RBPhotspot","RBP",res.df3$site2,perl = T)
#res.df3$site2[grepl("AGO2|otherRBPs|RBPhotspot",as.character(res.df3$site),perl = T)] <- "RBP"
table(res.df3$site2)

res.df3$site2 <- factor(res.df3$site2, levels = c("RBP","EV","G4","structure"))



res.tbl <- as_tibble(res.df3) %>% 
  dplyr::mutate(RBP_structure=dplyr::case_when( site2=="RBP" & (site_overlap_ratio_product>=0.01) ~ T,   #site2
                                                # site2=="otherRBPs" & (site_overlap_ratio_product>=0.01) ~ T,
                                                # site2=="RBPhotspot" & (site_overlap_ratio_product>=0.01) ~ T,
                                                site2=="EV" & (site_overlap_ratio_product>=0.01) ~ T,
                                                site2=="G4" & (site_overlap_ratio_product>=0.01) ~ T,
                                                site2=="structure" & (site_overlap_ratio_product<=0.1) ~ T,
                                                TRUE ~ F
  )
  )


mat <- res.tbl %>% 
  dplyr::select(RBP_structure,site2,peakID) %>% 
  dplyr::group_by(site2,peakID) %>% 
  dplyr::mutate(RBP_structure=sum(RBP_structure)) %>% 
  dplyr::distinct(RBP_structure,site2,peakID,.keep_all = TRUE) %>% 
  tidyr::pivot_wider(id_cols = peakID, names_from = site2, values_from = RBP_structure)

#library("ggvenn") # not suport more than 5
#create Venn diagram and display all sets
# ggvenn(as.data.frame(mat),set_name_size=4,text_size = 5,show_percentage=F, text_color = "white") + ggsci::scale_fill_nejm() + 
#   theme(aspect.ratio = 0.8)
mat <- as.data.frame(mat)
rownames(mat) <- mat$peakID
mat$peakID <- NULL
mat[1:3,]
# table(is.na(mat)) # some structure is na (seem filtered before rnafold)
mat[is.na(mat)] <- 0
mat <- as.matrix(mat)
mat[mat>=1] <- 1
mat[mat==0] <- 0
mat <- as.data.frame(mat)
# tmp <- mat[mat$G4iM==1 & rowSums(mat)==1,]
# tmp <- mat[is.na(mat[,5]),]
# table(is.na(mat[,5]))
colnames(mat) <- c("RBP bind","G4 Struc.","EV sort","2nd Struc.")
pdf("upset.pdf",width = 6,height = 3)
UpSetR::upset(data = mat, sets = c("RBP bind","EV sort","G4 Struc.","2nd Struc."), nsets = 4, keep.order = T,
              matrix.color = "salmon",set_size.numbers_size = T,
              sets.bar.color = pal_nejm_adaptive()(4)#colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(4)
)
dev.off()




mat.list <- list(rownames(mat)[mat$`RBP bind`==1],rownames(mat)[mat$`EV sort`==1],rownames(mat)[mat$`G4 Struc.`==1],rownames(mat)[mat$`2nd Struc.`==1])
names(mat.list) <- colnames(mat)
library("ggvenn")
#names(mat.list) 
ggvenn(data = mat.list, columns = c("RBP bind","EV sort","G4 Struc.","2nd Struc.") , set_name_size=6, text_size = 8,show_percentage=F, stroke_color = "white", text_color = "white") + ggsci::scale_fill_nejm() + theme(aspect.ratio = 0.8)
ggsave("tmpVenn.pdf",width = 10,height = 5, limitsize=F)

table(rowSums(mat)>0)
sum(rowSums(mat)>0)/nrow(mat)


## select overlap peak id
# peak.id1 <- rownames(mat)[mat$`RBP bind`==0 & mat$`G4 Struc.`==1 & mat$`EV sort`==1 & mat$`2nd Struc.`==0]
# peak.id2 <- rownames(mat)[mat$`RBP bind`==1 & mat$`G4 Struc.`==1 & mat$`EV sort`==0 & mat$`2nd Struc.`==1]
# peak.id3 <- rownames(mat)[mat$`RBP bind`==1 & mat$`G4 Struc.`==0 & mat$`EV sort`==1 & mat$`2nd Struc.`==1]
peak.id <- rownames(mat)[rowSums(mat)>=2]

y <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008/call_peak_all/expeak/b5_d50_p1_8DNA.bed"
peak <- read.table(y)
peak <- peak[peak$V4 %in% c(peak.id),] # c(peak.id1,peak.id2,peak.id3)
write.table(peak,"tmp/RBP_peak_8DNA.bed",quote = F,sep = "\t",row.names = F,col.names = F)

peak$RNA <- ref$transcript_type[match(peak$V1,ref$transcript_id)]
peak$RNA <- gsub("_rev|_for|\\.for|\\.rev","",peak$RNA)
#rna <- c("rRNA","pri_miRNA","lncRNA", "mRNA", "piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA')
dna <- c("intron","promoter", "enhancer","repeats") #
peak$RNA <- factor(peak$RNA,levels = c(dna)) #dna,rna
table(peak$RNA)

#filter peak of interest manually
AluJb__chr15___68202323____68202601_pos,135,155 # 278 (loose: >=1)
G041621__chr19___53051683____53057401_neg,550,564 # 5718
tRNA____Cys____TGY__chr6___149027805____149027840_neg,0,23 # 35
G____rich__chr6___7146362____7146415_pos,29,49 # 53
promoter__chr12___121806128____121809528_neg,914,934 # 3400
promoter__chr16___30486277____30489677_pos,2926,2944 # 3400
enhancer__chr9___2824800____2829600_neg,1639,1652 # 4800
enhancer__chr9___128906479____128909279_pos,568,586 # 2800

#bash
# pre="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
# /usr/bin/Rscript  $pre/exSeek-dev/bin/peak_intersect_G4iM_RBP_EV.R --coord tx --RNAtype RNA \
# -i tmp/RBP_peak.bed \
# -o tmp/RBP_peak.bed.intersect
#RBP_peak.bed.intersect used for RBP distribution plot
#may filter from old intersect ? skip this step ?





# sum isSHAPE structure barplot (icshape, exPeak, Fig5) -----------------------------------------------------------------
setwd("/BioII/lulab_b/baopengfei/projects/RNA-structure/")
options(stringsAsFactors = F)
suppressPackageStartupMessages(library(bedtoolsr))
options(bedtools.path = "/BioII/lulab_b/baopengfei/biosoft/")


## read ref
ref <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)
table(ref$transcript_type)
ref <- ref[ref$transcript_type=="mRNA",]
ref$transcript_id2 <- unlist(sapply(strsplit(ref$transcript_id,"_____"),"[",1))

## icSHaPE
#"./data/cell-2016-icShape/GSE74353_HS_293T_icSHAPE_InVitro_BaseReactivities_filterLen.txt"
#icshape <- rio::import(file_path,fill=TRUE)
#icshape <- data.table::fread(file_path,col.names = rep("sss",6000),sep = "\t", header = FALSE,fill=TRUE)
#icshape <- read.table("./data/cell-2016-icShape/GSE74353_HS_293T_icSHAPE_InVitro_BaseReactivities_filterLen.txt", header = FALSE,fill=T)
file_path <- "./data/cell-2016-icShape/GSE74353_HS_293T_icSHAPE_InVivo_BaseReactivities_filterLen.txt"
#file_path <- "./data/cell-2016-icShape/GSE74353_HS_293T_icSHAPE_InVitro_BaseReactivities_filterLen.txt"
icshape <- readLines(file_path)
icshape[1:2]
length(icshape)
icshape.l <- sapply(strsplit(icshape,"\t"),"[")
length(icshape.l[[4]])
icshape.df <- do.call("rbind",icshape.l)
icshape.df <- as.data.frame(icshape.df)
# dim(icshape.df)
# table(duplicated(icshape.df$V1))
g <- icshape.df$V1
gl <- as.numeric(icshape.df$V2)

icshape.df <- icshape.df[,4:ncol(icshape.df)]
icshape.df[1:4,1:10]
# table(icshape.df=="NULL")
# FALSE     TRUE 
# 39936730 59468710
icshape.df[icshape.df=="NULL" | icshape.df=="NA"] <- NaN
r <- nrow(icshape.df)
icshape.df <- matrix(as.numeric(as.matrix(icshape.df)),nrow=r)
# table(is.nan(icshape.df))
# FALSE     TRUE 
# 39936730 59468710
# table(icshape.df>2)
# FALSE     TRUE 
# 39567102   246401
# hist(icshape.df,xlim = c(0,10),breaks = 100000)
# str(icshape.df[1:4,1:100])
table(icshape.df>1)
hist(icshape.df,breaks = 100000,xlim = c(0,2)) # usually in [0,1]: 0 indicating a highly structured or inaccessible region and 1 representing an unstructured or accessible region.
icshape.df[icshape.df>1] <- 1

# hist(icshape.df,xlim = c(0,10),breaks = 100000)
# dim(icshape.df)
rownames(icshape.df) <- g

#only keep the same mRNA tx version/length
table(g %in% ref$transcript_id2)
# FALSE  TRUE 
# 1359  8801 
table(paste(g,gl) %in% paste(ref$transcript_id2,ref$tx.length))
# FALSE  TRUE 
# 9124  1036
icshape.df <- icshape.df[ paste(g,gl) %in% paste(ref$transcript_id2,ref$tx.length) , ]
g2 <- g[paste(g,gl) %in% paste(ref$transcript_id2,ref$tx.length)]
gl2 <- gl[paste(g,gl) %in% paste(ref$transcript_id2,ref$tx.length)]
dim( icshape.df )



#"FTC_small_localmax","FTC_long_piranha" 
#ds <- c("GSE104251","GSE123972")
# for (ds in c("GSE104251","GSE123972")){
#   # ds <- "GSE123972"
ds <- "GSE71008"
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
dst <- "GSE71008"
smps <- read.table(paste0(pre,"/exSeek-dev/data/",dst,"/sample_ids_NC_test15.txt"))$V1

for(smp in smps){
  #smp <- "SAMN03863400"
  print(smp)
  
  # ref.bed.df <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/transcript_table/all_newTxID.txt", data.table = F, header = T)
  # ref.bed.df <- ref.bed.df[,c("transcript_id","transcript_type","strand")]
  # #ref.bed.df$transcript_id2 <- paste0(ref.bed.df$transcript_id,".",ref.bed.df$strand)
  # table(duplicated(ref.bed.df$transcript_id))  # ENST00000607781.1
  # ref.bed.df <- ref.bed.df[!duplicated(ref.bed.df$transcript_id),]
  # # table(ref.bed.df$transcript_type)
  # rownames(ref.bed.df) <- ref.bed.df$transcript_id
  # ref.bed.df[1:4,]
  # tail(ref.bed.df)
  
  rna.type <- "mRNA" # c("lncRNA","mRNA","snoRNA","snRNA","srpRNA","tRNA","tucpRNA","Y_RNA")
  #for (rna.type in c("Y_RNA","lncRNA","mRNA","snoRNA","snRNA","srpRNA","tRNA","tucpRNA")){
  message(rna.type)
  
  ## read peak.bed
  #peak.bed <- read.table(paste0("/BioII/lulab_b/baopengfei/shared_reference/RBP/data/",ds),sep = "\t",header = F)
  # peak.bed <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",ds,"/call_peak_all/expeakCNN_by_sample/b5_d50_p1/intersect/",smp,".bed6"),sep = "\t",header = F)
  peak.bed <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",ds,"/call_peak_all/expeakCNN/b5_d50_p1_11RNA.bed.exp"),sep = "\t",header = F)
  peak.bed <- bedtoolsr::bt.slop(s=T,l=20,r=20,g = paste0(pre,"/exSeek-dev/genome/hg38/chrom_sizes/transcriptome_genome_sort_uniq_newTxID"),i=peak.bed)# slop/expand
  
  # peak.bed$type <- unlist(sapply(strsplit(peak.bed$V1,"|",fixed=T),tail,1)) # option1
  # table(peak.bed$type)
  #ref.bed.df <- NULL
  peak.bed$type <- ref$transcript_type[match(peak.bed$V1,ref$transcript_id)] # option2
  table((peak.bed$type))
  # table(is.na(peak.bed$type))
  peak.bed <- peak.bed[!is.na(peak.bed$type),]
  peak.bed <- peak.bed[peak.bed$type %in% rna.type,]
  # peak.bed$V1 <- gsub(paste0("|",rna.type,"|"), "", peak.bed$V1, fixed = T )
  # peak.bed[1:3,]
  
  peak.bed <- peak.bed[peak.bed$V3 < 10000,]
  peak.bed$V1 <- unlist(sapply(strsplit(peak.bed$V1,"_",fixed=T),"[",1))
  #peak.bed <- peak.bed[!duplicated(peak.bed$V1),]
  #peak.bed <- peak.bed[grepl("ENST",peak.bed$V1),]
  peak.bed <- peak.bed[peak.bed$V1 %in% rownames(icshape.df),]
  #table(peak.bed$type)
  peak.bed$mean.reactivity <- NA
  peak.bed$median.reactivity <- NA
  #summary(peak.bed$V3-peak.bed$V2)
  
  for(i in 1:nrow(peak.bed)){
    #i <- 2
    # print(i)
    enst <- peak.bed$V1[i]
    start <- peak.bed$V2[i]
    end <- peak.bed$V3[i]
    peak.bed$mean.reactivity[i] <- mean(as.numeric(icshape.df[enst,start:end]),na.rm = T)
    peak.bed$median.reactivity[i] <- median(as.numeric(icshape.df[enst,start:end]),na.rm = T)
  }
  
  
  ## shuffle (bg only)
  #rna.fa <- read.table(paste0("/BioII/lulab_b/baopengfei/tmp/index/mRNA.transcripts.fa"),sep = "\t",header = F)
  rna.fa <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/fasta/mRNA.fa"),sep = "\t",header = F)
  l1 <- rna.fa$V1[grepl(">",rna.fa$V1)]
  l1 <- gsub(">","",l1,fixed = T)
  rna.fa <- rna.fa[!grepl(">",rna.fa$V1),,drop=F]
  rna.fa$V2 <- l1
  rna.fa$len <- nchar(rna.fa$V1)
  rna.fa <- rna.fa[rna.fa$len < 10000,]
  rna.fa <- rna.fa[,2:3]
  g <- unlist(sapply(strsplit(rna.fa$V2,".",fixed=T),"[",1))
  rna.fa <- rna.fa[g %in% rownames(icshape.df),]
  #head(rna.fa)
  peak.bed.shuffle <- bedtoolsr::bt.shuffle(seed = 1234, noOverlapping=T, i = peak.bed[,1:3], g = rna.fa) #, chrom = T  # excl = args[4], seems not working
  
  # all(peak.bed.shuffle$V3==peak.bed$V3)
  # all((peak.bed.shuffle$V3-peak.bed.shuffle$V2)==(peak.bed$V3-peak.bed$V2))
  #peak.bed[,1:3] <- peak.bed.shuffle[,1:3]
  
  peak.bed.shuffle <- peak.bed.shuffle[peak.bed.shuffle$V3 < 10000,]
  peak.bed.shuffle$V1 <- unlist(sapply(strsplit(peak.bed.shuffle$V1,".",fixed=T),"[",1))
  #peak.bed.shuffle <- peak.bed.shuffle[!duplicated(peak.bed.shuffle$V1),]
  peak.bed.shuffle <- peak.bed.shuffle[peak.bed.shuffle$V1 %in% rownames(icshape.df),]
  peak.bed.shuffle$mean.reactivity <- NA
  peak.bed.shuffle$median.reactivity <- NA
  
  for(i in 1:nrow(peak.bed.shuffle)){
    #i <- 2
    # print(i)
    enst <- peak.bed.shuffle$V1[i]
    start <- peak.bed.shuffle$V2[i]
    end <- peak.bed.shuffle$V3[i]
    peak.bed.shuffle$mean.reactivity[i] <- mean(as.numeric(icshape.df[enst,start:end]),na.rm = T)
    peak.bed.shuffle$median.reactivity[i] <- median(as.numeric(icshape.df[enst,start:end]),na.rm = T)
  }
  
  
  # summary(peak.bed$mean.reactivity,breaks = 10000,xlim = c(0,1))
  # summary(peak.bed$median.reactivity,breaks = 10000,xlim = c(0,1))
  # 
  # summary(peak.bed.shuffle$mean.reactivity,breaks = 10000,xlim = c(0,1))
  # summary(peak.bed.shuffle$median.reactivity,breaks = 10000,xlim = c(0,1))
  # 
  # wilcox.test(peak.bed.shuffle$median.reactivity,peak.bed$median.reactivity)
  # wilcox.test(peak.bed.shuffle$mean.reactivity,peak.bed$mean.reactivity)
  
  dat <- as.data.frame(cbind(peak.bed$mean.reactivity,peak.bed.shuffle$mean.reactivity))
  colnames(dat) <- c("Domain","Background")
  #dat[1:3,]
  
  dat.m <- reshape2::melt(dat)
  colnames(dat.m) <- c("variable","value")
  dat.m$sample <- smp
  # # head(dat.m)
  # # str(dat.m)
  # dat.m1 <- na.omit(dat.m[dat.m$variable=="Domain",])
  # Rmisc::CI(dat.m1$value,ci = 0.95)
  # # Domain:
  # # upper    mean   lower 
  # # 0.19834 0.18361 0.16888
  # dat.m2 <- na.omit(dat.m[dat.m$variable=="Background",])
  # Rmisc::CI(dat.m2$value,ci = 0.95)
  # # Background:
  # # upper    mean   lower 
  # # 0.23040 0.21108 0.19175
  data.table::fwrite(dat.m,paste0("tmp/icshapeInVivo_",dst,"_",smp,".txt2"),quote = F,sep = "\t",row.names = F,col.names = T)
}


dat.list <- list()
for(smp in smps){
  #smp <- "SAMN03863400"
  print(smp)
  dat.list[[smp]] <- data.table::fread(paste0("tmp/icshapeInVivo_",dst,"_",smp,".txt2"),sep = "\t",header = T,check.names = F,stringsAsFactors = F)
}
dat.m <- as.data.frame(do.call(rbind,dat.list))


## plot  
library(ggplot2)
table(dat.m$sample)
table(dat.m$variable)
dat.m1 <- na.omit(dat.m[dat.m$variable=="Domain",])
Rmisc::CI(dat.m1$value,ci = 0.95)
dat.m2 <- na.omit(dat.m[dat.m$variable=="Background",])
Rmisc::CI(dat.m2$value,ci = 0.95)
dat.m$variable <- factor(dat.m$variable,levels = c("Domain","Background"))

ggbarplot(data = dat.m, x = "variable", y = "value", fill = "variable", add="mean_se", short.panel.labs = T, position = position_dodge()) + # add.params=list(size=4,color="black"),
  stat_compare_means(aes(x=variable, y=value, group=variable), # ref.group="Flank",  # 
                     # comparisons = list(c("Peak",bg.type)),
                     label="p.signif",  # ..p.signif../..p.format..
                     # symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01,1),symbols = c( "***", "**", "*", "ns")),
                     method = "wilcox.test", method.args = list(alternative = "greater"),  # greater means ref.group less
                     label.x.npc = 0.8, label.y.npc = 0.8, vjust = 0.6, #step.increase = 0.08,
                     hide.ns=T,size =0, paired = F
  ) +
  # facet_grid(site~.,scales = "free")+ # facet_grid seem not work well with ggbarplot
  scale_fill_manual(values = c('firebrick','salmon')) + #c("red","grey90")
  # scale_x_discrete(label=c("Peak Ratio","MIR Ratio")) +
  labs(title="",x="", y = "")+
  # xlim(c(0,1))+   
  ylim(c(0,1))+
  theme_minimal() +  # base_size=12
  theme(#axis.ticks.x=element_blank(),
    #strip.text.y = element_blank(),
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
# ggplot(dat.m,aes(x = variable, y = value, fill = variable)) +
#   geom_boxplot() +
#   # facet_wrap( ~ sample) +
#   ylim(c(0,1))+
#   labs(title = dst ) + 
#   ylab("icSHAPE reactivity (trim mean)")+
#   scale_fill_manual(name="Type",values = c('firebrick','salmon')) +
#   theme_bw() + 
#   ggpubr::stat_compare_means(    
#     label.x.npc = 0.3,
#     label.y.npc = 0.8,
#     hide.ns=F,
#     size = 8, #paired = TRUE,  
#     #aes(group = Group,label = p.format), # p.signif,p.format
#     label = "p.format",
#     method = "wilcox.test",
#     method.args = list(alternative = "greater")
#   )+
#   theme_minimal() +  # base_size=12
#   theme(#axis.ticks.x=element_blank(),
#     #strip.text.y = element_blank(),
#     #strip.text.x = element_text(face="bold",family="arial",size=20),
#     axis.title.x = element_text(size=20),
#     axis.title.y = element_text(size=20),
#     axis.text.x = element_text(size = 20,hjust = 1,vjust = 0.5,angle = 90), #
#     axis.text.y = element_text(size = 20),
#     plot.title = element_text(size=20),
#     strip.text = element_text(size = 20),
#     legend.position = "none", #c(0.9,0.8),#,#
#     legend.text = element_text(size= 16),
#     legend.title= element_text(size= 16))
dir.create(paste0("./output/",dst,"/"),showWarnings = F,recursive = T)
#ggsave(paste0("./output/",dst,"/exPeak_icSHPAPE_box.pdf"))
#ggsave(filename = paste0("./exPeak_icSHPAPE_box.pdf"), width = 3, height = 7) #p1+p2+patchwork::guide_area+patchwork::plot_layout()
ggsave(filename = paste0("./exPeak_icSHPAPE_bar.pdf"), width = 3, height = 7)
# }


mfe.2 <- na.omit(dat.m)
mfe.2$method <- "cfPeak"
mfe.2$mfe <- mfe.2$value
mfe.2$value <- NULL
mfe.2$type <- mfe.2$variable
mfe.2$variable <- NULL
mfe.2$type <- gsub("Domain","Peak",mfe.2$type)
mfe.2$type <- factor(mfe.2$type,levels = c("Peak","Background"))
table(mfe.2$type)
fc_data <- mfe.2 %>%
  group_by(method) %>%
  summarise(mean_firebrick = mean(mfe[type == 'Peak']),
            ci_firebrick = calc_mean_cl_normal(mfe[type == 'Peak']),
            mean_salmon = mean(mfe[type == 'Background']),
            ci_salmon = calc_mean_cl_normal(mfe[type == 'Background']),
            FC = mean_firebrick / mean_salmon,
            delta = mean_firebrick-mean_salmon) %>%
  mutate(ci_firebrick_text = sprintf("%.4f ± %.4f", ci_firebrick[1], ci_firebrick[3] - ci_firebrick[1]),
         ci_salmon_text = sprintf("%.4f ± %.4f", ci_salmon[1], ci_salmon[3] - ci_salmon[1]),
         FC_label = sprintf("FC=%.2f\nDelta=%.2f\nPeak: %s\nBackground: %s", FC, delta, ci_firebrick_text, ci_salmon_text)) %>%
  ungroup()
mfe.2 <- left_join(mfe.2, fc_data, by = c("type","method"))

p12 <- ggplot(mfe.2, aes(x = method, y = mfe, fill = type)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge()) + 
  geom_errorbar(stat = "summary", fun.data = "mean_cl_normal", position = position_dodge(0.9), width = 0.2) + # calc_mean_sem mean_cl_normal
  stat_compare_means(aes(group = type), method = "wilcox.test", method.args = list(alternative = "less"),
                     label = "p.signif", label.x.npc = "right", label.y.npc = "top", hide.ns = TRUE, size = 10) +
  # facet_wrap(~ site2, scales = "free") + 
  scale_fill_manual(values = c('firebrick', 'salmon')) +
  
  # Add FC and 95% CI labels for each method
  geom_text(data = fc_data, aes(x = method, y = 0.9, label = FC_label), 
            inherit.aes = FALSE, size = 5) + 
  
  labs(title = "", x = "", y = "") +
  # ylim(c(0, 1)) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 20, hjust = 1, vjust = 0.5, angle = 90),
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size = 20),
    legend.position = "none",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16)
  )
p12
# Save the plot
# ggsave(filename = "figure/domain_enrich_sum_expeak.pdf", plot = p12, width = 6, height = 8)
ggsave(filename = "domain_icSHAPE_sum_expeak2.pdf", plot = p12, width = 6, height = 13)
#






