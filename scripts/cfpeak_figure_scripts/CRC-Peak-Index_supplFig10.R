# test localmax

# last 2205 by bpf 
# b.p.f@qq.com
# transmit form local to bioii at 220604

setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/")
#setwd("/Users/baopengfei/Desktop/lulab/tmp/projects/WCHSU-FTC/exSeek-dev")
options(stringsAsFactors = F)
library(ggplot2)
library(conflicted)

# load all func.
source("./util.R")





# CRC-peak-index TCGA_small+GSE110381+PRJNA (GSE110381 peak) (cfpeak suppl Fig10)  ------------------------------------

## TCGA 
dst <- "TCGA_small_diff3"
#p <- 0.01
#feature <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_all/count_matrix/expeakCNN_b5_d50_p1_GSE110381.txt.diff.tissueHighP",p))$V1 # use CRC diff-high seem poor
feature <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_all/count_matrix/cfpeakCNN_b5_d50_p1_GSE110381.txt.diff.dcb"))$V1
#head(feature)

df <- data.frame(peak.id = unlist(sapply(strsplit(feature,"|",fixed=T),"[",4)),
                 RNA = unlist(sapply(strsplit(feature,"|",fixed=T),"[",2)))
df$RNA <- gsub("_rev|_for","",df$RNA,perl=T)
df$RNA <- factor(df$RNA,levels = c(rna,dna))
table(df$RNA)



## plot diff peak RNA type (pie plot)
#df <- df[df$peak.id %in% peak.ids[['top1000']],] # 50
df2 <- as.data.frame(table(df$RNA))
df2 <- df2[df2$Freq>0,]
df2 <-dplyr:: as_tibble(df2) %>% 
  dplyr::group_by(Var1) %>% 
  dplyr::summarize(Freq=mean(Freq))
df2$lab <- round(df2$Freq/sum(df2$Freq),digits = 3)
df2 <- df2 %>%
  dplyr::mutate(csum = rev(cumsum(rev(Freq))),
                pos = Freq/2 + lead(csum, 1),
                pos = if_else(is.na(pos), Freq/2, pos))

RNA.col <- data.frame(RNA=c(rna,dna),col=c(pal_nejm_adaptive()(15)[1:14],"#11838D"))
RNA.col$RNA <- factor(RNA.col$RNA,levels = c(rna,dna))

# str(df2)
ggplot(df2, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width=1, color="black") +
  coord_polar("y", start=0) +
  # geom_text(aes(label = lab), position = position_stack(vjust=0.5), size = 5, colour="white") +
  ggrepel::geom_label_repel(data = df2, col="white",
                            aes(y = pos, label = paste0(100*lab, "%")),
                            size = 8, nudge_x = 0.75, show.legend = FALSE) +
  labs(x = NULL, y = NULL, title ="") + # paste0("Total repeats peak: ",nrow(peak))
  scale_color_manual("black") +
  scale_fill_manual(name="precursor",values = RNA.col$col[RNA.col$RNA %in% df2$Var1] ) + 
  # scale_fill_nejm_adaptive(alpha = 0.8) +
  theme_bw() + 
  my_theme_pie
ggsave("TDCPs_pie.pdf",width = 7,height = 5)





## diff volcano plot of TDCP
dst <- "TCGA_small_diff3"
disease.label <- "COAD"
normal.label <- "LAML"

diff.plot <- list()
# tmp <- read.table(paste0("./output/",dst,"/exPeak_smallDomain_diff_",disease.label,"vs",normal.label,".diff"),sep = "\t",header = T,check.names = F)
tmp <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",dst,"/cfPeakCNN_GSE110381_smallDomain_diff_",disease.label,"vs",normal.label,".diff"),sep = "\t",header = T,check.names = F,stringsAsFactors = F)

tmp$minus.log10.pvalue <- -log10(tmp$pvalue)
tmp$log10.baseMean <- log10(tmp$baseMean)
table(feature %in% rownames(tmp)) # all True
tmp$TDCP <- "No"
tmp$TDCP[rownames(tmp) %in% feature] <- "Yes"
table(tmp$TDCP)
tmp <- tmp[order(tmp$TDCP),]
tmp$alpha <- maxmin.normalize.vec(tmp$minus.log10.pvalue)
diff.plot[[disease.label]] <- ggplot(tmp, aes(x=log10.baseMean,y=log2FoldChange)) + 
  geom_point(aes(fill=TDCP),shape=21,color="white",alpha=0.9) + 
  labs(title = paste0(disease.label," vs. ",normal.label)) +
  ylab("log2.FoldChange") +
  xlab("log10.Mean.CPM") +
  # ggraph::scale_fill_viridis(direction = -1) +
  scale_fill_manual(values = c("grey80","firebrick")) +
  geom_hline(yintercept = 0,linetype="dashed",color="salmon") + 
  theme_bw()+
  my_theme_volcano

# p <- ggpubr::ggarrange(plotlist = diff.plot, rel_widths = c(1,1,1), nrow = 1,  common.legend = T,align = "hv", axis = "b", legend = "right") #  labels = c("RBPs", "", "EV", "", "G4"),
ggsave(plot = diff.plot[[disease.label]], paste0("./TDCPs_log2FC_basemean.pdf"), width=6, height=6) # "_",sample






# TCGA_small_diff CRCvsNC AUROC
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
dst <- "TCGA_small_diff3"
NC.label <- "Blood"
CRC.label <- "COAD"

### read smp table
sample.table.tcga <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/archive/",dst,"/sample_table.txt"),check.names = F, sep = "\t", header = T)
#sample.table.tcga <- sample.table.tcga[,1:2]
sample.table.tcga <- sample.table.tcga[,c(2,5)]
colnames(sample.table.tcga) <- c("sample","group")
sample.table.tcga$sample <- gsub(".bam","",sample.table.tcga$sample)
rownames(sample.table.tcga) <- sample.table.tcga$sample
#table(sample.table.tcga$sample %in% colnames(cpm.tcga))
#sample.table.tcga <- sample.table.tcga[sample.table.tcga$sample %in% colnames(cpm.tcga),]
sample.table.tcga$group <- gsub("TCGA-","",sample.table.tcga$group)
table(sample.table.tcga$group)
sample.table.tcga[1:3,]
# sample.table.tcga$group <- gsub("COAD","CRC",sample.table.tcga$group)
sample.table.tcga$group <- gsub("LAML","Blood",sample.table.tcga$group)
sample.table.tcga <- sample.table.tcga[sample.table.tcga$group %in% c("COAD", "Blood"),] # this peak index can not classify other cancer types !
sample.table.tcga$group <- factor(sample.table.tcga$group,levels = c("Blood","COAD")) # , "PAAD", "PRAD"
sample.table.tcga <- sample.table.tcga[order(sample.table.tcga$group),]
sample.table.tcga <- sample.table.tcga[order(sample.table.tcga$group),]


#logcpm <- read.table("count_matrix",check.names = F,header = T) # cpm
#logcpm <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",dst,"/cfPeakCNN_GSE110381_smallDomain_diff_COADvsLAML.cpm"),check.names = F,header = T)  # log2cpm+1
logcpm <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_all/count_matrix/cfpeakCNN_b5_d50_p1_GSE110381_CPMrowsum.txt"),check.names = F,header = T)  # log2cpm+1
rownames(logcpm) <- logcpm$gene_id # feature
logcpm <- logcpm[,2:ncol(logcpm)] ## cpm, not log2cpm !!!
max(logcpm[,1])
cpm <- logcpm

# logcpm <- logcpm[feature,]
sample.table.tcga <- sample.table.tcga[sample.table.tcga$sample %in% colnames(logcpm),]
table(sample.table.tcga$group)
# Blood  COAD 
# 60    60 
#sample.table.tcga$group <- gsub("adenoma","CRC",sample.table.tcga$group)
#sample.table.tcga <- sample.table.tcga[order(sample.table.tcga$group),]
positive_samples <- sample.table.tcga[sample.table.tcga$group==CRC.label,"sample"]
negative_samples <- sample.table.tcga[sample.table.tcga$group==NC.label,"sample"]
samples <- c(positive_samples, negative_samples)
# sample.table.tcga <- sample.table.tcga[match(samples,sample.table.tcga$sample),]
# group <- c(rep("positive",length(positive_samples)),rep("negative",length(negative_samples)))
# logcpm <- logcpm[,sample.table.tcga$sample]



## build ref db
nrow(df)==length(feature) # TRUE !!!
rownames(df) <- feature
feature.lab.sort <- c(rna,dna)
feature.lab.sort <- feature.lab.sort[feature.lab.sort %in% as.character(unique(df$RNA))]
feature.lab.sort <- c(feature.lab.sort, "RNA", "DNA","all","highROC")

highROC <- c("pri_miRNA","lncRNA","mRNA","tucpRNA","intron","promoter","enhancer","repeats") # #only keep RNA species with roc >0.9 in TCGA

feature.list <- list()
mean.list <- list()
sd.list <- list()
for(i in feature.lab.sort){
  print(i)
  if(i=="all"){
    feature.list[[i]] <- rownames(df)
  }else if(i=="RNA"){
    feature.list[[i]] <- rownames(df)[df$RNA %in% rna]
  }else if(i=="DNA"){
    feature.list[[i]] <- rownames(df)[df$RNA %in% dna]
  }else if(i=="highROC"){
    feature.list[[i]] <- rownames(df)[df$RNA %in% highROC]
  }else{
    feature.list[[i]] <- rownames(df)[df$RNA==i]
  }
  
  #extract fixed coad&blood group
  blood.mat <- cpm[feature.list[[i]], negative_samples] # cpm, not log2cpm !!!
  coad.mat <- cpm[feature.list[[i]], positive_samples]
  mean.list[[i]] <- apply(blood.mat,1,mean) # length(mean.vec)
  sd.list[[i]] <- apply(cbind(coad.mat,blood.mat),1,sd) # blood mat has many 0 : Zhu - 2021 - NC - Tissue-specific cell-free DNA degradation...
  # sd.list <- (rep(sd(as.matrix(blood.mat)),length(mean.vec)))
  #print(table(sd.list==0))
  print( table(sd.list[[i]]==0) ) # need all FALSE !!! or else enlarge cohort when calculating sd (eg: add other groups)
}




### read count mat
#x <- paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",dst,"/cfPeakCNN_GSE110381_smallDomain_diff_COADvsLAML.cpm")
#logcpm <- read.table(x,check.names = F,header = T)  # log2cpm+1
#dim(logcpm)
logcpm <- log2(logcpm+1)
#max(logcpm[,1])


logcpm.sum.list <- list()
for(i in feature.lab.sort ){
  print(i)
  if(i=="all"){
    feature.list[[i]] <- rownames(df)
  }else if(i=="RNA"){
    feature.list[[i]] <- rownames(df)[df$RNA %in% rna]
  }else if(i=="DNA"){
    feature.list[[i]] <- rownames(df)[df$RNA %in% dna]
  }else if(i=="highROC"){
    feature.list[[i]] <- rownames(df)[df$RNA %in% highROC]
  }else{
    feature.list[[i]] <- rownames(df)[df$RNA==i]
  }
  logcpm.sum.list[[i]] <- getPeakIndex(feature.lab = i,logcpm = logcpm, sample.table = sample.table.tcga, prescale = T, mean.list = mean.list, sd.list = sd.list, zscore.cutoff = 100, onlyCountZscoreOutlier = F)
}


### plot boxplot
for(y in names(logcpm.sum.list)) {
  logcpm.sum.list[[y]]$feature.lab <- y
  logcpm.sum.list[[y]]$value <- logcpm.sum.list[[y]]$mean.log2cpm.scale # sum.log2cpm.scale, mean.log2cpm
  logcpm.sum.list[[y]]$group <- factor(logcpm.sum.list[[y]]$group,levels = c(NC.label,CRC.label))
} 
logcpm.sum.df <- as.data.frame(do.call(rbind,logcpm.sum.list))
logcpm.sum.df$feature.lab <- factor(logcpm.sum.df$feature.lab,levels = feature.lab.sort)

plotViolin(logcpm.sum.df,sig.size = 4) + # 2e-16
  facet_grid(~feature.lab)
ggsave(filename = paste0("./",dst,"_CRC-idx_box_all19.pdf"),width = 48,height = 8)
plotViolin(logcpm.sum.df[logcpm.sum.df$feature.lab=="all",],sig.size = 4) # 2e-16
ggsave(filename = paste0("./",dst,"_CRC-idx_box_all.pdf"),width = 8,height = 8)

plotMultiROC(logcpm.sum.list, direction = "<", RNA_colors = RNA_colors)
ggsave(filename = paste0("./",dst,"_CRC-idx_multiroc_all19.pdf"),width = 9,height = 9)

# #only keep RNA species with roc >0.9 in TCGA
# highROC <- c("pri_miRNA","lncRNA","mRNA","tucpRNA","intron","promoter","enhancer","repeats")






## GSE110381 CRCvsNC 
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
dst <- "GSE110381_diff"
NC.label <- "NC"
CRC.label <- "CRC"
### read smp table
sample.table <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/GSE110381_diff/sample_table.txt",check.names = F,header = T) 
# passQC <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/GSE110381_diff/sample_ids_diff_passQC2_rmBatch.txt",check.names = F,header = F)$V1
# sample.table <- sample.table[sample.table$sample %in% passQC,]
rownames(sample.table) <- sample.table$sample
sample.table$group <- gsub("adenoma","CRC",sample.table$group)

x <- paste0(pre,"/output/",dst,"/call_peak_all/count_matrix/expeakCNN_b5_d50_p1_CPMrowsum.txt") # use all samples: better auroc for this dst
#x <- paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",dst,"/filerSmp_exPeakCNN_smallDomain_diff_CRCvsNC.cpm")
logcpm <- read.table(x,check.names = F,header = T)  # not log2cpm+1
max(logcpm[,2])
rownames(logcpm) <- logcpm$gene_id
logcpm$gene_id <- NULL
max(as.matrix(logcpm))
logcpm <- log2(logcpm+1)


logcpm.sum.list <- list()
for(i in feature.lab.sort ){
  print(i)
  if(i=="all"){
    feature.list[[i]] <- rownames(df)
  }else if(i=="RNA"){
    feature.list[[i]] <- rownames(df)[df$RNA %in% rna]
  }else if(i=="DNA"){
    feature.list[[i]] <- rownames(df)[df$RNA %in% dna]
  }else if(i=="highROC"){
    feature.list[[i]] <- rownames(df)[df$RNA %in% highROC]
  }else{
    feature.list[[i]] <- rownames(df)[df$RNA==i]
  }
  logcpm.sum.list[[i]] <- getPeakIndex(feature.lab = i,logcpm = logcpm, sample.table = sample.table, prescale = T, mean.list = mean.list, sd.list = sd.list, zscore.cutoff = 100, onlyCountZscoreOutlier = F)
}


### plot boxplot
for(y in names(logcpm.sum.list)) {
  logcpm.sum.list[[y]]$feature.lab <- y
  logcpm.sum.list[[y]]$value <- logcpm.sum.list[[y]]$mean.log2cpm.scale # sum.log2cpm.scale, mean.log2cpm
  logcpm.sum.list[[y]]$group <- factor(logcpm.sum.list[[y]]$group,levels = c(NC.label,CRC.label))
} 
logcpm.sum.df <- as.data.frame(do.call(rbind,logcpm.sum.list))
logcpm.sum.df$feature.lab <- factor(logcpm.sum.df$feature.lab,levels = feature.lab.sort)

plotViolin(logcpm.sum.df,sig.size = 4) + # 2e-16
  facet_grid(~feature.lab)
ggsave(filename = paste0("./",dst,"_CRC-idx_box_all19.pdf"),width = 48,height = 8)
plotViolin(logcpm.sum.df[logcpm.sum.df$feature.lab=="all",],sig.size = 4) # 2e-16
ggsave(filename = paste0("./",dst,"_CRC-idx_box_all.pdf"),width = 8,height = 8)

plotMultiROC(logcpm.sum.list, direction = "<", RNA_colors = RNA_colors)
ggsave(filename = paste0("./",dst,"_CRC-idx_multiroc_all19.pdf"),width = 9,height = 9)





## PRJNA540919_diff CRCvsNC
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
dst <- "PRJNA540919_diff"
NC.label <- "NC"
CRC.label <- "CRC"
### read smp table
sample.table <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/PRJNA540919_diff/sample_table.txt",check.names = F,header = T)
colnames(sample.table) <- c("sample","group")
table(sample.table$group)

### read cpm
logcpm <- read.table(paste0(pre,"/output/",dst,"/call_peak_all/count_matrix/cfpeakCNN_b5_d50_p1_GSE110381_CPMrowsum.txt"),check.names = F,header = T)
# x <- paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",dst,"/filterSmp_exPeakCNN_smallDomain_diff_CRCvsNC.cpm")
# logcpm <- read.table(x,check.names = F,header = T)  # log2cpm+1
rownames(logcpm) <- logcpm$gene_id
logcpm$gene_id <- NULL
max(as.matrix(logcpm))
logcpm <- log2(logcpm+1)

logcpm.sum.list <- list()
for(i in feature.lab.sort ){
  print(i)
  if(i=="all"){
    feature.list[[i]] <- rownames(df)
  }else if(i=="RNA"){
    feature.list[[i]] <- rownames(df)[df$RNA %in% rna]
  }else if(i=="DNA"){
    feature.list[[i]] <- rownames(df)[df$RNA %in% dna]
  }else if(i=="highROC"){
    feature.list[[i]] <- rownames(df)[df$RNA %in% highROC]
  }else{
    feature.list[[i]] <- rownames(df)[df$RNA==i]
  }
  logcpm.sum.list[[i]] <- getPeakIndex(feature.lab = i,logcpm = logcpm, sample.table = sample.table, prescale = T, mean.list = mean.list, sd.list = sd.list, zscore.cutoff = 100, onlyCountZscoreOutlier = F)
}


### plot boxplot
for(y in names(logcpm.sum.list)) {
  logcpm.sum.list[[y]]$feature.lab <- y
  logcpm.sum.list[[y]]$value <- logcpm.sum.list[[y]]$mean.log2cpm.scale # sum.log2cpm.scale, mean.log2cpm
  logcpm.sum.list[[y]]$group <- factor(logcpm.sum.list[[y]]$group,levels = c(NC.label,CRC.label))
} 
logcpm.sum.df <- as.data.frame(do.call(rbind,logcpm.sum.list))
logcpm.sum.df$feature.lab <- factor(logcpm.sum.df$feature.lab,levels = feature.lab.sort)

plotViolin(logcpm.sum.df,sig.size = 4) + # 2e-16
  facet_grid(~feature.lab)
ggsave(filename = paste0("./",dst,"_CRC-idx_box_all19.pdf"),width = 48,height = 8)
plotViolin(logcpm.sum.df[logcpm.sum.df$feature.lab=="all",],sig.size = 4) # 2e-16
ggsave(filename = paste0("./",dst,"_CRC-idx_box_all.pdf"),width = 8,height = 8)
plotViolin(logcpm.sum.df[logcpm.sum.df$feature.lab %in% c("all","RNA","DNA"),],sig.size = 0) + # 2e-16
  facet_grid(~feature.lab)
ggsave(filename = paste0("./",dst,"_CRC-idx_box_all3.pdf"),width = 8,height = 8)

plotMultiROC(logcpm.sum.list, direction = "<", RNA_colors = RNA_colors) # poor if use onlyCountZscoreOutlier=3 !!!!
ggsave(filename = paste0("./",dst,"_CRC-idx_multiroc_all19.pdf"),width = 9,height = 9)
plotMultiROC(list("all"=logcpm.sum.list[["all"]],"RNA"=logcpm.sum.list[["RNA"]],"DNA"=logcpm.sum.list[["DNA"]]), direction = "<", RNA_colors = RNA_colors) # poor if use onlyCountZscoreOutlier=3 !!!!
ggsave(filename = paste0("./",dst,"_CRC-idx_multiroc_all3.pdf"),width = 9,height = 9)
# plotMultiROC(list("all"=logcpm.sum.list[["all"]],"pri_miRNA"=logcpm.sum.list[["pri_miRNA"]])) 
# ggsave(filename = paste0("./",dst,"_CRC-idx_multi2roc_2.pdf"),width = 9,height = 9)

#



