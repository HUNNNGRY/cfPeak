# func. for cfPeak
# last 2305 by bpf 
# b.p.f@qq.com

library(tidyverse)
library(dplyr)
library(ggpubr)

tocapital <- function(my_string){
  capital_string <- paste(toupper(substring(my_string, 1, 1)), substring(my_string, 2), sep = "")
  return(capital_string)
}

# test filtering too many standard peak number (deprecated) ------
test <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/tmp/4dst_filter2_11RNA.bed")
colnames(test) <- c("seqnames","start","end","name","score","strand")
test$width <- test$end-test$start
table(test$width<=10)
hist(test$width)
#test <- test[test$width>=10 & test$width<=200,] # filterLen
test <- test[order(test$V5,decreasing = T),]
for (top in c(500, 1000, 2000, 5000, 10000)) { # 
  test.tmp <- test[1:min(top,nrow(test)),]
  write.table(test.tmp[,1:6],paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/tmp/4dst_filter2_11RNA_top",top,".bed"),quote = F,sep = "\t",row.names = F,col.names = F)
}





#keep >=2 dsts detect reads
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
#bwType <- "tbigwig_11RNA_primary"
bwType <- "tbigwig_RNA_EM"
cmb <- data.frame(dst=c("GSE71008_NCpool","GSE94533_NCpool","GSE123972_NCpool",    "GSE94582_NCpool","GSE94582_NCpool","GSE94582_NCpool","GSE94582_NCpool"),
                  smp=c("NCpool","NCpool","NCpool",     "NCpool_CleanTag","NCpool_N4","NCpool_NEBNext","NCpool_TruSeq"))

res.list <- list()
for (i in 1:nrow(cmb)){## consensus peak
  #i <- 1
  dst <- cmb$dst[i]
  smp <- cmb$smp[i]
  # if(i %% 1000 ==0){
  #   print(i)
  # }
  count.tmp <- read.table(paste0(pre,"/exSeek-dev/tmp/4dst_filter2_11RNA.bed.",dst,".",smp,".count"),header = F )
  print(paste0(dst," ",smp,": ",(nrow(count.tmp))))
  colnames(count.tmp)[5] <- paste0(dst,"-",smp)
  res.list[[paste0(dst,smp)]] <- count.tmp[,paste0(dst,"-",smp),drop=F]
}
#table(res.list[[1]]$V1==res.list[[5]]$V1)
df <- as.data.frame(do.call(cbind,res.list))
# tmp <- (res.list[[2]])

test <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/tmp/4dst_filter2_11RNA.bed")
colnames(test) <- c("seqnames","start","end","name","score","strand")

cutoff.count <- 2
recur <- rowSums(df>=cutoff.count)
#length(recur)==nrow(test)
cutoff.recur <- 1
table(recur>=cutoff.recur)

for (cutoff.recur in c(1:4)) { # 
  test.tmp <- test[recur>=cutoff.recur,]
  print(nrow(test.tmp))
  write.table(test.tmp[,1:6],paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/tmp/4dst_filter2_11RNA_exist",nrow(test.tmp),".bed"),quote = F,sep = "\t",row.names = F,col.names = F)
}
#


# calculate confusion matrix (v1: DIY in R, deprecated)  -----------------------------------------------------------------------
bedtoolsPath <- "/BioII/lulab_b/baopengfei/biosoft"
options(bedtools.path = bedtoolsPath) # /BioII/lulab_b/baopengfei/anaconda3/envs/exvariance/bin
chrSize <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/transcriptome_genome_sort_uniq_newTxID"
standardFile <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/tmp/3dst_filter2_11RNA.bed" # "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/domains_clipper_b5_p10_recur_2.bed"
seed <- 1234

## read gold standard
gold.standard.peak <- data.table::fread(standardFile,data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)
gold.standard.peak <- gold.standard.peak[,1:6]
colnames(gold.standard.peak) <- c("chr","start","end","name","score","strand")
# k <- list(gold.standard.peak)
# l <- c("gold.standard.peak")


res <- list()
for (sample_id in sample_ids){
#sample_id <- "SRR2105127"
print(sample_id)
for (method in methods){
#method <- "domains_by_sample/b5_p10"
inputFile <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_domain_withRepeats_all/",method,"/intersect/",sample_id,".bed") # sample-wise with pval

if (grepl("EM2",method)){
  method.name <- "exPeak"
}else if  (grepl("localmax",method)){
  method.name <- "LocalMax"
}else if  (grepl("clam",method)){
  method.name <- "CLAM"
}else if  (grepl("clipper",method)){
  method.name <- "CLIPper"
}else {
  method.name <- "Piranha"
}
print(method.name)

## read domain bed
domain0 <- as.data.frame(data.table::fread(inputFile,data.table = F,sep = '\t',check.names = F,stringsAsFactors = F))
domain0[1:3,]
if (method.name=="exPeak" | method.name=="LocalMax"){
  domain0 <- domain0[,c(1:4,9,6)]
} else if (method.name=="CLAM"){
  domain0$V5 <- unlist(sapply(strsplit(domain0$V4,","),"[",1))
  domain0$V5 <- as.numeric(unlist(sapply(strsplit(domain0$V5,":"),"[",3)))
  domain0 <- domain0[,1:6]
} else if (method.name=="CLIPper"){
  domain0 <- domain0[,1:6]
} else if (method.name=="Piranha"){
  domain0 <- domain0[,c(1:4,7,6)]
}

colnames(domain0) <- c("chr","start","end","name","pval","strand")

# use count, not pval
if (useCountReplacePval) {
  count <- as.data.frame(data.table::fread(paste0(inputFile,".count"),data.table = F,sep = '\t',check.names = F,stringsAsFactors = F))
  count <- count[match(domain0$name,count$V4),]
  count <- na.omit(count)
  domain0$pval <- count$V5
}


domain0$name <- paste0(domain0$name,"--",domain0$chr,"--",domain0$start,"--",domain0$end)
domain0$pval[is.nan(domain0$pval) | is.na(domain0$pval)] <- 0
# domain0$pval.group <- cut(domain0$pval, breaks=pval_seq, include.lowest = TRUE, right = TRUE )

## get real peak number (all peaks overlap with gold standard)
domain0.intersect <- bedtoolsr::bt.intersect(a = domain0[,1:6],wao = T,s = T,b = gold.standard.peak)
x <- "tmp"
domain0.intersect <- as_tibble(domain0.intersect) %>% 
  dplyr::mutate(overlap_ratio_product = V13/(V3-V2) * (V13/(V9-V8)), site = x) %>% 
  dplyr::arrange(desc(overlap_ratio_product)) %>% 
  dplyr::distinct(V4, .keep_all = TRUE)   # only keep max overlap_ratio_product each peak
colnames(domain0.intersect) <- c( c("chr","start","end","name","score","strand"),paste0("site","_",c("chr","start","end","name","score","strand","overlap_base","overlap_ratio_product") ) ,"site" ) 
domain0.intersect$type <- "real.neg"
domain0.intersect$type[domain0.intersect$site_overlap_ratio_product>=0.01 ] <- "real.pos" # | domain0.intersect$site_overlap_base>=5
table(domain0.intersect$type)
real.neg <- sum(domain0.intersect$type == "real.neg")
real.pos <- sum(domain0.intersect$type == "real.pos")
if (real.pos*real.neg==0){
  print("error: real.pos*real.neg==0 !")
}


pval_seq2 <- unique(c(domain0$pval)) #c(-1,pval_seq)
# if(length(pval_seq2)>=50){
#   pval_seq2 <- pval_seq2[sample(1:length(pval_seq2),50,replace = F)]
# }
# pval_seq2 <- unique(c(-1,pval_seq2,max(domain0$pval)+1))
for (i in 1:length(pval_seq2)){
  # cutoff <- 0
  cutoff <- pval_seq2[i]
  print(cutoff)
  if (useCountReplacePval) {
    pred.pos <- sum(domain0$pval>=cutoff)
    pred.neg <- sum(domain0$pval<cutoff)
    TP <- sum(domain0$pval>=cutoff & domain0.intersect$type == "real.pos")
    TN <- sum(domain0$pval<cutoff & domain0.intersect$type == "real.neg")
  }else{
    pred.pos <- sum(domain0$pval<=cutoff)
    pred.neg <- sum(domain0$pval>cutoff)
    TP <- sum(domain0$pval<=cutoff & domain0.intersect$type == "real.pos")
    TN <- sum(domain0$pval>cutoff & domain0.intersect$type == "real.neg")
  }
  FP <- pred.pos-TP
  FN <- pred.neg-TN
  if (pred.pos*pred.neg==0){
    print("error: pred.pos*pred.neg==0 !")
  }
  if (TP==0 | TN==0 | FP==0 | FN==0){
    print("error: TP*TN*FP*FN==0 !")
  }

  id <- paste0(method.name,sample_id,i)
  res[[id]] <- data.frame("sample"=sample_id, "method"=method.name, "cutoff"=cutoff, "TP"=TP, "FP"=FP, "TN"=TN, "FN"=FN)
}
}
}
res.df <- do.call(rbind,res)
res.df <- dplyr::as_tibble(res.df) %>% 
  dplyr::mutate("precision" = (TP)/(TP+FP),
                "fpr" = (FP)/(FP+TN),
                "recall" = (TP)/(TP+FN))
  # dplyr::mutate("precision" = (TP+1)/(TP+FP+1),
  #               "fpr" = (FP+1)/(FP+TN+1),
  #               "recall" = (TP+1)/(TP+FN+1))
# res.df$precision[is.nan(res.df$precision)] <- 0
# res.df$fpr[is.nan(res.df$fpr)] <- 0
# res.df$precision[is.nan(res.df$precision)] <- 0




# calculate confusion matrix (test params: 2022,GenomeBiol,GoPeaks) -----------------------------------------------------------------------
## run in bash (need run RBP/G4 intersect before get matrix in R)
setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/")
dst <- "GSE71008" # "GSE71008,WSQ_SMARTer_NEB"
sample_ids <- Sys.glob("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/*multi3*_*_decay*_pval*_t*p*_SRR2105340.bed") # SRR2105340,NC_PKU-2392860_smart_PNK_1
sample_ids <- sapply(sample_ids, function(x) basename(x))
sample_ids <- gsub(".bed","",sample_ids)

#dst="WSQ_SMARTer_NEB"
#smp="NC_PKU-2392860_smart_PNK_1"


## get matrix in R
res <- list()
for (sample_id in sample_ids){
  #sample_id <- "notrecursive_multi3local_background_decay0.1_pval1.0_t12p1_SRR2105340"
  print(sample_id)
  # for (method in methods){
    #method <- "domains_by_sample/b5_p10"
  #inputFile <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_domain_withRepeats_all/",method,"/intersect/",sample_id,".bed.count.txt") # sample-wise count
  inputFile <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/tmp/intersect/",sample_id,".bed.count.txt") # test params
  # if (grepl("EM2",method)){
  #   method.name <- "exPeak"
  # }else if  (grepl("localmax",method)){
  #   method.name <- "LocalMax"
  # }else if  (grepl("clam",method)){
  #   method.name <- "CLAM"
  # }else if  (grepl("clipper",method)){
  #   method.name <- "CLIPper"
  # }else {
  #   method.name <- "Piranha"
  # }
  # print(method.name)
  # 
  tmp <- read.table(inputFile,header = T,sep = "")
  colnames(tmp)[colnames(tmp)=="replicate"] <- "sample"
  colnames(tmp)[colnames(tmp)=="counts"] <- "cutoff"
  tmp$method <- sample_id
  
  id <- paste0(sample_id)
  res[[id]] <- tmp
# }
}
res.df <- do.call(rbind,res)



# calculate confusion matrix (fig3: adapted from 2022,GenomeBiol,GoPeaks) -----------------------------------------------------------------------
setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/") 

# dsts <- c("GSE71008_NCpool") # "GSE71008","GSE94533", "GSE110381", "WSQ_SMARTer_NEB","GSE94582" ,"GSE50676"
dsts <- c("GSE71008_NCpool","GSE94533_NCpool","GSE123972_NCpool", "GSE94582_NCpool","WSQ_SMARTer_NEB_NCpool","GSE148861_GSE148862_NCpool", "Phospho-RNA-seq_NCpool", "GSE110381_NCpool", "AGO2_IP_NCpool")
#"GSE50676_NCpool"
#sample_ids <- c("NCpool") #c("SRR2105340","SRR2105341","SRR2105342","SRR2105346","SRR2105335")  #c("SRR5230634","SRR5230661","SRR5230662","SRR5230667") # c("SRR2105336","SRR2105128","SRR2105129","SRR2105335","SRR2105336"), c("NC_PKU-2392860_smart_1","NC_PKU-2392860_smart_PNK_1","NC_PKU-2392860_1"),c(CleanTag_Lab5,NEBNext_Lab1,NEBNext_Lab2,TruSeq_Lab1,TruSeq_Lab2,N4_A_Lab5,N4_B_Lab2),c("GM12878","GM12891")
# methods <- c("piranha_by_sample/b5_p01","clipper_by_sample/b5_p05","clam_by_sample/b5_p005","expeakCNN_by_sample/b5_d50_p1") # ,"domains_localmax_significant/b5_d05_p05"
methods <- c("piranha/b5_p01","clipper/b5_p05","clam/b5_p005","expeakCNN/b5_d50_p1") # ,"domains_localmax_significant/b5_d05_p05"


## get matrix in R
res <- list()
for (frac in c("")){ # "0.1","0.5"
  for (top in c("")){ # "500","1000","2000","5000","10564"; '9426', '7077', '4893','3419'
    for (dst in dsts){
      print(dst)
      # dst <- "AGO2_IP_NCpool"
      # tmp <- tmp_df[tmp_df$dst==dst,]
      sample_ids <- read.table(paste0("data/",dst,"/sample_ids.txt"))$V1
      for (sample_id in sample_ids){
        #sample_id <- "SRR2105127"
        # sample_id <- "NCpool"
        print(sample_id)
        for (method in methods){
          # method <- "domains_by_sample/b5_p05"
          
          if(grepl("AGO2_IP",dst)){
            dedup <- "dedup"
          }else{
            dedup <- "all"
          }
          # inputFile <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_",dedup,"/",method,"/intersect/",sample_id,".bed.count.txt") # sample-wise count
          # inputFile <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_",dedup,"/",method,"_11RNA.bed.count.txt") # sample-consensus count
          inputFile <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_",dedup,"/",method,"_11RNA.bed.count",".filterReads6dst",frac,top,".txt") # sample-consensus count
          
          # inputFile <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/tmp/intersect/",sample_id,".bed.count") # test params
          if (grepl("expeak",method)){
            method.name <- "exPeak"
          }else if  (grepl("localmax",method)){
            method.name <- "LocalMax"
          }else if  (grepl("clam",method)){
            method.name <- "CLAM"
          }else if  (grepl("clipper",method)){
            method.name <- "CLIPper"
          }else {
            method.name <- "Piranha"
          }
          print(method.name)
          
          tmp <- read.table(inputFile,header = T,sep = "")
          colnames(tmp)[colnames(tmp)=="replicate"] <- "sample"
          colnames(tmp)[colnames(tmp)=="counts"] <- "cutoff"
          tmp$method <- method.name
          tmp$dst <- dst
          tmp$top <- top
          tmp$frac <- frac
          # tmp$smp <- sample_id
          
          res[[paste0(dst,sample_id,method.name,frac,top)]] <- tmp
        }
      }
    }
  }
}
res.df <- do.call(rbind,res)
table(res.df$dst, res.df$method)
#res.df[1:3]
# summary(res.df$fpr)
# omit any record that not in [0,1]
res.df <- res.df[res.df$precision>=0 & res.df$precision<=1 & res.df$recall>=0 & res.df$recall<=1 & res.df$fpr>=0 & res.df$fpr<=1 & res.df$f1>=0 & res.df$f1<=1 ,]
# res.df.bak <- res.df
# res.df.bak$method <- gsub("exPeak","cfPeak",res.df.bak$method)
# table(res.df.bak$method)
# data.table::fwrite(res.df.bak,"/BioII/lulab_b/baopengfei/smb/cfpeak/GSE71008/benchmark.txt",quote = F,row.names = F,col.names = T,sep = "\t")



## plot standard-overlap/recall/precison pie/bar plot 
# pval_seq2 <- c(0,10^seq(-15,-2,5),0.05,0.1)
# res.df2 <- res.df[res.df$cutoff %in% pval_seq2,]
res.df2 <- res.df
#res.df2$sample <- ""
#res.df2$recurrence <- res.df2$TP/(res.df2$TP+res.df2$FP+res.df2$TN+res.df2$FN)
res.df2 <- dplyr::as_tibble(res.df2) %>% 
  dplyr::group_by(dst,method,sample) %>%
  dplyr::arrange(cutoff) %>% # only select min cutoff (all peaks)
  dplyr::distinct(method, sample, dst, .keep_all = TRUE) %>% 
  tidyr::pivot_longer(cols = precision:f1, names_to = "group",values_to = "number") %>%
  dplyr::mutate(total.peak.num=TP+FP+TN+FN, no.number=1-number) %>% # TN+FN=0
  tidyr::pivot_longer(cols = c("number","no.number"), names_to = "group.fill",values_to = "number") %>%
  dplyr::mutate(method2=paste0(method,".",group.fill))
table(res.df2$group.fill)
res.df2$method2 <- factor(res.df2$method2, levels = c("Piranha.no.number","Piranha.number","CLIPper.no.number","CLIPper.number","CLAM.no.number","CLAM.number","exPeak.no.number","exPeak.number"))
res.df2$method <- gsub("exPeak","cfPeak",res.df2$method)
res.df2$method <- factor(res.df2$method, levels = c("Piranha","CLIPper","CLAM","cfPeak"))
#table(res.df2$method2 )
#str(res.df2)
#unique(res.df2$method2)
mypal <- c("#1F77B4E5","#FF7F0EE5","#2CA02CE5","#D62728E5") #pal_d3_adaptive("category10", alpha = 0.9)(4)
black.col <- alpha("grey90",0.3)
#scales::show_col(mypal)
#table(res.df2$group)
#tmp <- res.df2[res.df2$group=="TP",]
res.df2 <- res.df2[res.df2$group!="fpr",]
res.df2$label <- base::format(res.df2$number,digits = 1)
res.df2$label[res.df2$group.fill=="no.number"] <- ""




### plot single dst (GSE71008_NCpool) (fig3)
this.dst <- "GSE71008_NCpool"  # GSE71008_NCpool
for (g in unique(res.df2$group)){
  # g <- "precision"
  print(g)
  tmp <- res.df2[res.df2$group==g & res.df2$dst==this.dst,]
  tmp <- tmp[tmp$group.fill=="number",]
  ggpubr::ggbarplot(data = tmp, x = "method", y = "number", fill = "method", position = position_stack(), width = 0.5) + # , add="mean_se", , facet.by=c("group")
    # stat_compare_means(#aes(x=type, y=ratio, group=site), # ref.group="Flank",  # 
    #   comparisons = list(c("Domain","Flank")),
    #   label="p.signif",  # ..p.signif../..p.format..
    #   # symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01,1),symbols = c( "***", "**", "*", "ns")),
    #   method = "wilcox.test", method.args = list(alternative = "greater"),  # greater means ref.group less
    #   label.x.npc = 0.8, label.y.npc = 0.8, vjust = 0.9, #step.increase = 0.08,
    #   hide.ns=T,size =14, paired = T
    # ) +
    # facet_grid(site~.,scales = "free")+ # facet_grid seem not work well with ggbarplot
    # scale_fill_manual(values = c('firebrick','salmon'))+ #c("red","grey90")
    # scale_x_discrete(label=c("Domain Ratio","MIR Ratio"))+
  # scale_fill_d3_adaptive()+
  # scale_fill_manual(values = c(black.col,mypal[1],black.col,mypal[2],black.col,mypal[3],black.col,mypal[4])) +
  scale_fill_manual(values = c(mypal[1],mypal[2],mypal[3],mypal[4])) +
    geom_text(aes(x=method,y=0.15*max(tmp$number),label=label),size=6,color="black",angle=90)+
    labs(title="",x="", y = tocapital(g))+
    # xlim(c(0,1))+   
    # ylim(c(0,1))+
    theme_minimal() +  # base_size=12
    theme(#axis.ticks.x=element_blank(),
      #strip.text.y = element_blank(),
      aspect.ratio = 1,
      strip.text = element_text(size=20),
      axis.title.x = element_text(size=20),
      axis.title.y = element_text(size=24),
      axis.text.x = element_text(size = 20,vjust = 1,hjust = 1, angle = 45, color="black"), #
      axis.text.y = element_text(size = 20),
      plot.title = element_text(size=20),
      # strip.text = element_blank(),
      legend.position = "none", #c(0.9,0.8),#,#
      legend.text = element_text(size= 16),
      legend.title= element_text(size= 16))
  ggsave(paste0(this.dst,"_",g,".png"),width = 5,height = 5)
}





### plot multi dsts (GSE71008_NCpool...) (suppl fig3 ?)
#dsts <- c("GSE71008_NCpool", "GSE94533_NCpool")  # , "GSE110381_NCpool", "GSE94582_NCpool"
#"GSE123972_NCpool", "Phospho-RNA-seq_NCpool", "GSE110381_NCpool","GSE94582_NCpool","WSQ_SMARTer_NEB_NCpool","GSE148861_GSE148862_NCpool","AGO2_IP_NCpool","GSE50676_NCpool"
res.df2$X <- gsub("_NCpool","",res.df2$dst )
res.df2$X <- gsub("_NEBNext","",res.df2$X )
res.df2$X <- gsub("WSQ_SMARTer_NEB","in-house",res.df2$X )
res.df2$X <- gsub("Phospho-RNA-seq","GSE126051",res.df2$X )
res.df2$X <- gsub("GSE148861_GSE148862","GSE148861",res.df2$X )
res.df2$X <- factor(res.df2$X, levels = c("GSE71008", "GSE94533","GSE123972", "GSE110381","GSE94582", "GSE126051", "in-house","AGO2_IP","GSE148861"))

#str(res.df2)
for (g in unique(res.df2$group)){
  # g <- "precision"
  print(g)
  tmp <- res.df2[res.df2$group==g,]
  tmp <- tmp[tmp$group.fill=="number",]
  
  ggpubr::ggbarplot(data = tmp, x = "X", y = "number", fill = "method", position = position_dodge(width = 0.8)) + #, add="mean_se", facet.by=c("dataset")
    ggsci::scale_fill_d3()+
    labs(title="",x="", y = stringr::str_to_title(g))+
    # xlim(c(0,1))+   
    # ylim(c(0,1))+
    theme_minimal() +  # base_size=12
    theme(#axis.ticks.x=element_blank(),
      #strip.text.y = element_blank(),
      strip.text = element_text(size=20),
      axis.title.x = element_text(size=20),
      axis.title.y = element_text(size=20),
      axis.text.x = element_text(size = 20,hjust = 1,vjust = 0.5,angle = 90), #
      axis.text.y = element_text(size = 24),
      plot.title = element_text(size=20),
      # strip.text = element_blank(),
      legend.position = "right", #c(0.9,0.8),#,#
      legend.text = element_text(size= 16),
      legend.title= element_text(size= 16)) #+
    # facet_grid(frac~top)
  #ggsave("multi_AUPR.pdf",width = 12,height = 8)
  ggsave(paste0("multi_",g,".png"),width = 12,height = 5)
}


## pie chart (deprecated)
# res.df3 <- res.df2 
# # %>% 
# #   tidyr::pivot_longer(cols = "recurrence":"unoverlap.ratio", names_to = "Overlap", values_to = "Value")
# # res.df3$Overlap <- factor(res.df3$Overlap, levels = c("unoverlap.ratio","recurrence"))
# 
# # clr <- c("#1B61A5","#FC6910","#269321","#C9121E")
# methods <- c("Piranha","CLIPper","CLAM","exPeak")#as.character(unique(res.df3$method)) 
# p.list <- list()
# for (i in 1:length(methods)){
#   method <- methods[i]
#   # method <- "CLIPper"
#   print(method)
#   res.df3.tmp <- res.df3[res.df3$method==method,]
#   
#   res.df3.tmp$lab <- ifelse(res.df3.tmp$group=="TP",paste0(round(res.df3.tmp$overlap.ratio, digits = 2)*100, "%"),"")
#   res.df3.tmp$ratio <- ifelse(res.df3.tmp$group=="FP",res.df3.tmp$unoverlap.ratio,res.df3.tmp$overlap.ratio)
#   res.df3.tmp$group <- factor(res.df3.tmp$group)
#   #res.df3.tmp
#   tmp.p <- ggplot(res.df3.tmp, aes(x="", y=ratio, fill=group)) +
#     geom_bar(stat="identity", width=1, color="black") +
#     coord_polar("y", start=0) +
#     # geom_text(aes(label = lab), position = position_stack(vjust=0.5), size = 15, colour="white") +
#     labs(x = NULL, y = NULL) +
#     scale_color_manual("black") +
#     scale_fill_manual(values = c('grey60',mypal[i])) + 
#     theme_bw() + 
#     theme(plot.title = element_text(size = 24,color="black",hjust = 0.5,face="bold"),
#           axis.title = element_text(size = 24,color ="black"), 
#           axis.ticks = element_blank(),
#           panel.grid=element_blank(),
#           # panel.grid.major.x=element_blank(),
#           # panel.grid.minor.x = element_blank(),
#           # panel.grid.major.y = element_line(color = "grey50",linetype = "dashed"), #size= 1,
#           #panel.grid.minor.y = element_blank(),
#           panel.border = element_blank(),
#           axis.text = element_blank(), #element_text(size= 20,color = "black"),
#           # axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5), #  , color = c(rep("#003366",length(rna)),rep("darkred",length(dna)))
#           legend.position = "none",#c(.25,.6),
#           legend.text = element_text(size= 16),
#           legend.title= element_text(size= 16))
#   p.list[[i]] <- tmp.p
#   # theme_classic() +
#   # theme(axis.line = element_blank(),
#   #       axis.text = element_blank(),
#   #       axis.ticks = element_blank()) +
#   # scale_fill_brewer(palette="Blues")
#   # ggsave(filename = paste0(method,"_recur_pie.pdf"), plot = tmp.p)
# }
# p <- gridExtra::grid.arrange(p.list[[1]],p.list[[2]],p.list[[3]],p.list[[4]], ncol=4)
# ggsave(filename = paste0("recur_pie_precision.pdf"), plot = p,height = 3, width = 13)
# #0.987,0.788,0.701,0.813
# 
# standard.num <- 10564
# res.df3$recall.real <- res.df3$total.peak.num * res.df3$precision / standard.num
# p.list <- list()
# for (i in 1:length(methods)){
#   method <- methods[i]
#   # method <- "CLIPper"
#   print(method)
#   res.df3.tmp <- res.df3[res.df3$method==method,]
#   
#   res.df3.tmp$lab <- ifelse(res.df3.tmp$group=="TP",paste0(round(res.df3.tmp$recall.real, digits = 2)*100, "%"),"")
#   res.df3.tmp$ratio <- ifelse(res.df3.tmp$group=="FP",1-res.df3.tmp$recall.real,res.df3.tmp$recall.real)
#   res.df3.tmp$group <- factor(res.df3.tmp$group)
#   #res.df3.tmp
#   tmp.p <- ggplot(res.df3.tmp, aes(x="", y=ratio, fill=group)) +
#     geom_bar(stat="identity", width=1, color="black") +
#     coord_polar("y", start=0) +
#     # geom_text(aes(label = lab), position = position_stack(vjust=0.5), size = 15, colour="white") +
#     labs(x = NULL, y = NULL) +
#     scale_color_manual("black") +
#     scale_fill_manual(values = c('grey60',mypal[i])) + 
#     theme_bw() + 
#     theme(plot.title = element_text(size = 24,color="black",hjust = 0.5,face="bold"),
#           axis.title = element_text(size = 24,color ="black"), 
#           axis.ticks = element_blank(),
#           panel.grid=element_blank(),
#           # panel.grid.major.x=element_blank(),
#           # panel.grid.minor.x = element_blank(),
#           # panel.grid.major.y = element_line(color = "grey50",linetype = "dashed"), #size= 1,
#           #panel.grid.minor.y = element_blank(),
#           panel.border = element_blank(),
#           axis.text = element_blank(), #element_text(size= 20,color = "black"),
#           # axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5), #  , color = c(rep("#003366",length(rna)),rep("darkred",length(dna)))
#           legend.position = "none",#c(.25,.6),
#           legend.text = element_text(size= 16),
#           legend.title= element_text(size= 16))
#   p.list[[i]] <- tmp.p
#   # theme_classic() +
#   # theme(axis.line = element_blank(),
#   #       axis.text = element_blank(),
#   #       axis.ticks = element_blank()) +
#   # scale_fill_brewer(palette="Blues")
#   # ggsave(filename = paste0(method,"_recur_pie.pdf"), plot = tmp.p)
# }
# p <- gridExtra::grid.arrange(p.list[[1]],p.list[[2]],p.list[[3]],p.list[[4]], ncol=4)
# ggsave(filename = paste0("recur_pie_recall.pdf"), plot = p,height = 3, width = 13)
# #0.0357,0.244,0.128,0.321
# #
# 




# ## plot ROC,PR (deprecated)
# {
# library(dplyr)
# library(stringr)
# #devtools::install_github("r-lib/conflicted")
# library(conflicted)
# library(ggplot2)
# #library(rio)
# #library(data.table)
# library(tidyr)
# conflict_prefer("filter", "dplyr")
# conflict_prefer("xlim", "ggplot2")
# 
# # model_df$replicate <- as.factor(model_df$replicate)
# tmp_df <- as.data.frame(res.df) #%>% 
# tmp_df <- tmp_df[order(tmp_df$sample,tmp_df$method,tmp_df$cutoff),]
# #head(tmp_df)
# 
# # plot ROC and PR
# my_theme <- theme_minimal() +  # base_size=12
#   theme(#axis.ticks.x=element_blank(),
#     #strip.text.y = element_blank(),
#     aspect.ratio = 1,
#     strip.text = element_text(size=20),
#     axis.title.x = element_text(size=20),
#     axis.title.y = element_text(size=20),
#     axis.text.x = element_text(size = 20), # angle = 90,vjust = 0.5,hjust = 1
#     axis.text.y = element_text(size = 20),
#     plot.title = element_text(size=20),
#     # strip.text = element_blank(),
#     legend.position = "none", #c(0.9,0.8),#,#
#     legend.text = element_text(size= 16),
#     legend.title= element_text(size= 16))
# }
# 
# {
# tmp_df$method <- factor(tmp_df$method,levels = c("Piranha","CLIPper","CLAM","exPeak")) # ,"LocalMax"
# #str(tmp_df)
# tmp_roc <- ggplot(tmp_df, aes(x = fpr, y = recall, color = method)) +
#   # geom_point() +
#   geom_path(size = 1) +
#   xlim(c(0,1)) +
#   # ylim(c(0,1)) +
#   xlab("False Positive Rate") +
#   ylab("Recall") +
#   scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1.00), limits = c(0,1),labels = c("0.00","0.25","0.50","0.75","1.00")) +
#   geom_abline(slope = 1, linetype = 2) +
#   # facet_wrap(sample ~ ., nrow = 1,) +
#   ggsci::scale_color_d3()+
#   # scale_color_brewer(palette = "RdYlBu", direction = -1) +
#   # theme_minimal() +
#   my_theme 
#   # theme(legend.position = "")
# # theme(plot.background = element_rect(fill = "white"))
# # tmp_roc
# tmp_pr <- ggplot(tmp_df, aes(x = recall, y = precision, color = method)) +
#   # geom_point()+
#   geom_path(size = 1) +
#   xlim(c(0,1)) +
#   xlab("Recall") +
#   ylab("Precision") +
#   scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1.00), limits = c(0,1),labels = c("0.00","0.25","0.50","0.75","1.00")) +
#   # facet_wrap(sample ~ ., nrow = 1, scales = "free_y") +
#   ggsci::scale_color_d3()+
#   # scale_color_brewer(palette = "RdYlBu", direction = -1) +
#   # theme_minimal() +
#   my_theme
# # theme(plot.background = element_rect(fill = "white"))
# # tmp_pr
# ggpubr::ggarrange(plotlist = list(tmp_roc,tmp_pr), heights=c(1,1), ncol = 1, nrow = 2, align = "hv")
# }
# ggsave("./roc_pr.pdf", width = 4.5, height = 10, dpi = 600)
# # dir.create("data/figures/roc", recursive = TRUE)
# # dir.create("data/figures/pr", recursive = TRUE)
# # ggsave("./roc.pdf", tmp_roc, width = 10, height = 5, dpi = 600)
# # ggsave("./pr.pdf", tmp_pr, width = 10, height = 5, dpi = 600)
# # }
# # 








# plot (test expeak params, deprecated) -----------------------------------------------------------------------
library(dplyr)
library(stringr)
#devtools::install_github("r-lib/conflicted")
library(conflicted)
library(ggplot2)
#library(rio)
#library(data.table)
library(tidyr)
conflict_prefer("filter", "dplyr")
conflict_prefer("xlim", "ggplot2")

# model_df$replicate <- as.factor(model_df$replicate)
tmp_df <- as.data.frame(res.df) #%>% 
# tmp_df <- tmp_df[order(tmp_df$sample,tmp_df$method,tmp_df$cutoff),]

convertParamLab <- function(expeak_path){
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
  
  tmp <- paste0(recur,"_",mode,"_",boundary,"_",d,"_",p)
  return(tmp)
}

tmp_df$method <- sapply(tmp_df$method,function(x) convertParamLab(x))

tmp_df$method <- factor(tmp_df$method,levels = unique(tmp_df$method)) # 
tmp_df$sample <- "X"
tmp_df <- tmp_df[order(tmp_df$sample,tmp_df$method,tmp_df$cutoff),]


# plot ROC and PR
my_theme <- theme_minimal() + 
  theme(plot.title = element_text(size = 20,color="black",hjust = 0.5,face="bold"),
        axis.title = element_text(size = 20,color ="black"), 
        axis.title.y = element_text(size = 20,color ="black"),  #, vjust = 1,angle = 90
        axis.text = element_text(size= 16,color = "black"),
        #panel.grid=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(), #size= 1,
        panel.grid.minor.y = element_blank(),
        # panel.border = element_blank(),
        axis.text.x = element_text(), # angle = 45, hjust = 1 
        axis.text.y = element_text(size= 16,color = "black"),#, angle=90
        legend.position = "right",
        legend.text = element_text(size= 16),
        legend.title= element_text(size= 16),
        strip.text = element_text(size= 16) ) #+
# tmp_df$method <- factor(tmp_df$method,levels = c("Piranha","CLIPper","CLAM","LocalMax","exPeak")) # comment for test expeak params

#str(tmp_df)
# roc1 <- pROC::roc(tmp_df$,)
# pROC::ci.auc(roc1)
tmp_roc <- ggplot(tmp_df, aes_string(x = "fpr", y = "recall", color = "method")) +
  # geom_point()
  geom_path(size = 1) +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  xlab("False Positive Rate") +
  ylab("Recall") +
  geom_abline(slope = 1, linetype = 2) +
  facet_wrap(sample ~ ., nrow = 1,) +
  scale_color_d3_adaptive()+
  # scale_color_brewer(palette = "RdYlBu", direction = -1) +
  # theme_minimal() +
  my_theme
  # theme(plot.background = element_rect(fill = "white"))
tmp_roc
tmp_pr <- ggplot(tmp_df, aes_string(x = "recall", y = "precision", color = "method")) +
  # geom_point()+
  geom_path(size = 1) +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  xlab("Recall") +
  ylab("Precision") +
  facet_wrap(sample ~ ., nrow = 1, scales = "free_y") +
  scale_color_d3_adaptive()+
  # scale_color_brewer(palette = "RdYlBu", direction = -1) +
  # theme_minimal() +
  my_theme
  # theme(plot.background = element_rect(fill = "white"))
tmp_pr
ggpubr::ggarrange(plotlist = list(tmp_roc,tmp_pr), widths = c(2,2), ncol = 2, nrow = 1, align = "hv")
ggsave("./roc_pr.pdf", width = 16, height = 5, dpi = 600)
# dir.create("data/figures/roc", recursive = TRUE)
# dir.create("data/figures/pr", recursive = TRUE)
# ggsave("./roc.pdf", tmp_roc, width = 10, height = 5, dpi = 600)
# ggsave("./pr.pdf", tmp_pr, width = 10, height = 5, dpi = 600)
# }













# plot mult-dataset AUROC/AUPR barplot (fig3, calculate area under curve manually, deprecated) ---------------------------------------
# devtools::install_github('smin95/smplot2')
library(smplot2)
setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/")
dsts <- c("GSE71008_NCpool", "GSE94533_NCpool")  # , "GSE110381_NCpool", "GSE94582_NCpool"
#"WSQ_SMARTer_NEB_NCpool", "Phospho-RNA-seq_NCpool", "GSE110381_NCpool","GSE94582_NCpool","GSE123972_NCpool","GSE148861_GSE148862_NCpool" 
#"AGO2_IP_NCpool","GSE50676_NCpool"
# methods <- c("piranha_by_sample/b5_p01","clipper_by_sample/b5_p05","clam_by_sample/b5_p005","expeakCNN_by_sample/b5_d50_p1") # ,"domains_localmax_significant/b5_d05_p05"
methods <- c("piranha/b5_p01","clipper/b5_p05","clam/b5_p005","expeakCNN/b5_d50_p1") # ,"domains_localmax_significant/b5_d05_p05"

res2 <- list()
for (frac in c("")){ # "0.1","0.5"
  for (top in c("")){ # "500","1000","2000","5000","10564"; '9426', '7077', '4893','3419'
    for (dst in dsts){
      print(dst)
      # dst <- "AGO2_IP_NCpool"
      # tmp <- tmp_df[tmp_df$dst==dst,]
      sample_ids <- read.table(paste0("data/",dst,"/sample_ids.txt"))$V1
      for (sample_id in sample_ids){
        #sample_id <- "SRR2105127"
        # sample_id <- "NCpool"
        print(sample_id)
        for (method in methods){
          # method <- "domains_by_sample/b5_p05"
          
          if(grepl("AGO2_IP",dst)){
            dedup <- "dedup"
          }else{
            dedup <- "all"
          }
          # inputFile <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_",dedup,"/",method,"/intersect/",sample_id,".bed.count.txt") # sample-wise count
          # inputFile <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_",dedup,"/",method,"_11RNA.bed.count.txt") # sample-consensus count
          inputFile <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_",dedup,"/",method,"_11RNA.bed.count",".filterReads",frac,top,".txt") # sample-consensus count
          
          # inputFile <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/tmp/intersect/",sample_id,".bed.count") # test params
          if (grepl("expeak",method)){
            method.name <- "exPeak"
          }else if  (grepl("localmax",method)){
            method.name <- "LocalMax"
          }else if  (grepl("clam",method)){
            method.name <- "CLAM"
          }else if  (grepl("clipper",method)){
            method.name <- "CLIPper"
          }else {
            method.name <- "Piranha"
          }
          print(method.name)
          
          tmp <- read.table(inputFile,header = T,sep = "")
          colnames(tmp)[colnames(tmp)=="replicate"] <- "sample"
          colnames(tmp)[colnames(tmp)=="counts"] <- "cutoff"
          tmp$method <- method.name
          tmp$dst <- dst
          tmp$top <- top
          tmp$frac <- frac
          # tmp$smp <- sample_id

          res2[[paste0(dst,sample_id,method.name,frac,top)]] <- tmp
        }
      }
    }
  }
}
res2.df <- do.call("rbind",res2) 
res2.df <- res2.df[res2.df$precision>=0 & res2.df$precision<=1 & res2.df$recall>=0 & res2.df$recall<=1 & res2.df$fpr>=0 & res2.df$fpr<=1 & res2.df$f1>=0 & res2.df$f1<=1 ,]
#table(res2.df$sample, res2.df$dst)
#table(res2.df$top,res2.df$frac)
#table(res2.df$dst)

res2.list <- list()
for (frac in c("")){ # "0.1","0.5"
  for (top in c("")){ # "500","1000","2000","5000","10564"; '9426', '7077', '4893','3419'
    for (dst in unique(res2.df$dst)){
      # dst <- "WSQ_SMARTer_NEB_NCpool"
      print(dst)
      tmp <- res2.df[res2.df$dst==dst & res2.df$frac==frac & res2.df$top==top,]
      for (smp in unique(tmp$sample)){
        # smp <- "NCpool_smart"
        print(smp)
        tmp2 <- tmp[tmp$sample==smp,]
        for (method in unique(tmp2$method)){
          # print(method)
          # if (grepl("EM2",method)){
          #   method.name <- "exPeak"
          # }else if  (grepl("localmax",method)){
          #   method.name <- "LocalMax"
          # }else if  (grepl("clam",method)){
          #   method.name <- "CLAM"
          # }else if  (grepl("clipper",method)){
          #   method.name <- "CLIPper"
          # }else {
          #   method.name <- "Piranha"
          # }
          method.name <- method
          print(method.name)
          tmp3 <- tmp2[tmp2$method==method,]
          res2.list[[paste0(frac,top,dst,smp,method.name)]] <- data.frame("dataset"=dst,"method"=method.name,"sample"=smp,"frac"=frac,"top"=top,
                                                        "AUROC"=smplot2::sm_auc(tmp3$fpr, tmp3$recall),
                                                        "AUPR"=smplot2::sm_auc(tmp3$recall, tmp3$precision))
        }
      }
    }
  }
}
res2.list.df <- do.call("rbind",res2.list)
#table(res2.list.df$sample, res2.list.df$dataset)
#table(res2.list.df$dataset)

# ## plot
# # res2.list.df$dataset <- 
# table(res2.list.df$sample )
res2.list.df$method <- factor(res2.list.df$method, levels = c("Piranha","CLIPper","CLAM","exPeak"))
# res2.list.df$X <- paste0(res2.list.df$dataset) # ,"_",res2.list.df$sample
# res2.list.df$X <- gsub("_NCpool","",res2.list.df$X )
# res2.list.df$X <- gsub("_NEBNext","",res2.list.df$X )
# res2.list.df$X <- gsub("WSQ_SMARTer_NEB","inhouse",res2.list.df$X )
# res2.list.df$X <- gsub("Phospho-RNA-seq","GSE126051",res2.list.df$X )
# res2.list.df$X <- gsub("GSE148861_GSE148862","GSE148861",res2.list.df$X )
# 
# table(res2.list.df$X)
# 
# res2.list.df <- res2.list.df[!(res2.list.df$X %in% c("GSE71008","GSE110381","AGO2_IP","inhouse_smart_PNK","inhouse_smart","inhouse","GSE126051_PNK","GSE126051_none","GSE94582_TruSeq","GSE94582_CleanTag","GSE94582_N4") ),]
# res2.list.df$X <- factor(res2.list.df$X,levels = c("GSE94533","GSE94582","GSE123972","GSE148861","GSE50676"))
# # res2.list.df <- res2.list.df[res2.list.df$sample %in% c("NCpool_NEBNext","NCpool"),]
# # res2.list.df <- res2.list.df[!(res2.list.df$dataset %in% c("GSE110381_NCpool") ),]
# # table(res2.list.df$dataset )
# #res2.list.df[1:3,]
# # ggplot(data = res2.list.df, aes(x = "X", y = "AUROC", fill = "method")) + 
# #   geom_bar(stat = "summary", position="dodge") + 
#   # facet_grid(dataset~.,scale="free") +
# table(res2.list.df$dataset)
ggpubr::ggbarplot(data = res2.list.df, x = "dataset", y = "AUROC", fill = "method", position = position_dodge(width = 0.8)) + #, add="mean_se", facet.by=c("dataset")
  scale_fill_d3()+
  labs(title="",x="", y = "AUROC")+
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
    legend.position = "right", #c(0.9,0.8),#,#
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16)) + 
  facet_grid(frac~top)
#ggsave("multi_AUROC.pdf",width = 12,height = 8)
ggsave("multi_AUROC.pdf",width = 24,height = 10)
#table(res2.list.df$X)
ggpubr::ggbarplot(data = res2.list.df, x = "dataset", y = "AUPR", fill = "method", position = position_dodge(width = 0.8)) + #, add="mean_se", facet.by=c("dataset")
  scale_fill_d3()+
  labs(title="",x="", y = "AUPR")+
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
    legend.position = "right", #c(0.9,0.8),#,#
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16)) +
  facet_grid(frac~top)
#ggsave("multi_AUPR.pdf",width = 12,height = 8)
ggsave("multi_AUPR.pdf",width = 24,height = 10)



