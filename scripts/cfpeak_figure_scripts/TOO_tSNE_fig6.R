# test t-SNE
# last 220322 by bpf
# b.p.f@qq.com

setwd("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/")
options(scipen = 4,stringsAsFactors = F,digits = 5)


# load all func.
source("./util.R")





#################
# GSE71008
#################
setwd("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/")
options(scipen = 4,stringsAsFactors = F,digits = 5)


dst <- "GSE71008_diff"
# sample.table <- read.table("./meta/lulab/FTC/sample_table.txt",sep = "\t",header = T)
sample.table0 <- read.table("../WCHSU-FTC/exSeek-dev/data/GSE71008/sample_table.txt",sep = "\t",header = T)
#sample.table <- reshape2::melt(sample.table,id.var=c("patient_id","group","name"),variable.name="source",value.name = "sample")
#sample.table$group <- substr(sample.table$patient_id,1,3)
colnames(sample.table0)[c(1,3)] <- c("sample","group")
#sample.table <- sample.table[sample.table$sample!="",]
table(sample.table0$group)

## read mat
count <- read.table(paste0("../WCHSU-FTC/output/",dst,"/call_peak_all/count_matrix/cfpeakCNN_b5_d50_p1.txt"),header = T,sep = "\t",stringsAsFactors = F,check.names = F)

rownames(count) <- count$feature
mat <- count[,2:ncol(count)]
# mat <- mat[,!(colnames(mat) %in% fail)]


## diff
peak.list <- list()
nums <- c(50,100,200,500,1000,5000)
for (top in nums){
  peak.list[[paste0("top",top)]] <- list()
}
for (disease in c("Colorectal Cancer","Pancreatic Cancer","Prostate Cancer","Healthy Control")){
  #disease <- "Colorectal Cancer" # "Colorectal Cancer" "Healthy Control"   "Pancreatic Cancer" "Prostate Cancer"
    if(disease=="Colorectal Cancer"){
      disease.label <- "CRC"
    }else if(disease=="Pancreatic Cancer"){
      disease.label <- "PACA"
    }else if(disease=="Prostate Cancer"){
      disease.label <- "PRCA"
    }else if(disease=="Healthy Control"){
      disease.label <- "NC"
    }
    print(disease)
    #normal <- "" #"Healthy Control"
    normal.label <- "R" #"NC" 
  
  tmp <- read.table(paste0("./output/",dst,"/allSmp_rmBatch_OvR_cfPeakCNN_smallDomain_diff_",disease.label,"vs",normal.label,".diff"),sep = "\t",header = T,check.names = F,stringsAsFactors = F)
    
    tmp <- tmp[order(as.numeric(tmp$pvalue),decreasing = F),]
  rownames(tmp) <- unlist(sapply(strsplit(rownames(tmp),"|",fixed=T),"[",4))

  # op1: fixed cutoff
  #tmp.sig <- tmp[tmp$pvalue<0.00001,] 

  # op2: fixed top
  for (i in nums){
    peak.id <- rownames(tmp)[1:i]
    # print(length(peak.id))
    peak.list[[paste0("top",i)]][[disease]] <- peak.id
  }
}

peak.ids <- list()
for (top in nums){
  peak.ids[[paste0("top",top)]] <- unique(do.call("c",peak.list[[paste0("top",top)]])) 
  print(length( peak.ids[[paste0("top",top)]] ))
}


## mds/PCA/tSNE plot (must TMM-CPM after filtering, not read previous *.cpm or *.rowsum.cpm !!!!)
# y <- DGEList(counts=mat) # cpm of all samples
# counts <- edgeR::getCounts(y)
# y <- calcNormFactors(y, method="TMM")
# logcpm <- edgeR::cpm(y, log=TRUE, normalized.lib.sizes = T,prior.count = 1)
# logcpm <- read.table(paste0("./output/",dst,"/filter112Smp_OvR_exPeakCNN_smallDomain_diff_",disease.label,"vs",normal.label,".cpm"), header = T, sep="\t",check.names = F)
disease.label <- "NC" 
normal.label <- "R"
logcpm <- read.table(paste0("./output/",dst,"/allSmp_rmBatch_OvR_cfPeakCNN_smallDomain_diff_",disease.label,"vs",normal.label,".cpm"), header = T, sep="\t",check.names = F)
max(logcpm[,3]) # log2CPM+1
# logcpm <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_all/count_matrix/expeakCNN_b5_d50_p1_CPMrowsum.txt"), header = T, sep="\t",check.names = F, row.names = 1)
rownames(logcpm) <- unlist(sapply(strsplit(rownames(logcpm),"|",fixed=T),"[",4))
sample.table0 <- sample.table0[sample.table0$sample %in% colnames(logcpm),]
sample.table0$group <- factor(sample.table0$group, levels = c("Healthy Control","Colorectal Cancer","Prostate Cancer","Pancreatic Cancer"))
sample.table0 <- sample.table0[order(sample.table0$group),] # draw point in group order !!!!
table(sample.table0$group)
logcpm <- logcpm[,sample.table0$sample]

library(Rtsne)
library(scatterplot3d)
library(plyr)
library(dplyr)

#CategColors <- c("blue", "green", "gray", "purple", "orangered", "black", "brown", "yellow3")
#PerplexityVec = seq(5, 15, by=5) # 70 too large, 20
PerplexityVec = seq(5, 60, by=10) # (80-1)/3
PerplexityVec
Iter1000_Rtsne_Results = vector(mode="list", length=length(PerplexityVec))

Phenotype <- sample.table0$group[match(colnames(logcpm),sample.table0$sample)]
table(Phenotype)
# Healthy Control Colorectal Cancer   Prostate Cancer Pancreatic Cancer 
# 50               100                36                 6
Colours4Classes_TCGA=c("grey", "salmon", "purple", "steelblue") # "orangered", "brown", "black", "gray",salmon
AllSamplesColors_TCGA=mapvalues(Phenotype, unique(Phenotype), Colours4Classes_TCGA)
AllSamplesColors <- as.character(AllSamplesColors_TCGA)

for (q in 1:length(PerplexityVec)) {
    #q <- 1
    # pdf(paste0("rm112_OvR_tSNE_GSE71008_iter100_perplexity", PerplexityVec[q], ".pdf"), width = 13,height = 10) # 
    pdf(paste0("all_rmBatch_OvR_tSNE_GSE71008_iter100_perplexity", PerplexityVec[q], ".pdf"), width = 13,height = 10) # 
    par(mfrow=c(3,4)) # length(nums)
    #par(mfrow=c(2,2))
    
    for (top in nums){
      #top <- 50
      print(top)
      logcpm2 <- logcpm[peak.ids[[paste0("top",top)]],sample.table0$sample]
      print(dim(logcpm2))
      
      j <- q
      # for (j in 1:length(PerplexityVec)) {
      #   print(j)
        set.seed(42)  # Sets seed for reproducibility
        Iter1000_Rtsne_Results[[j]] <- Rtsne(t(logcpm2), dims=3, max_iter=1000, perplexity=PerplexityVec[j])
        
        # set.seed(42)  # Sets seed for reproducibility
        # Iter3000_Rtsne_Results[[j]] <- Rtsne(t(logcpm2), dims=3, max_iter=3000, perplexity=PerplexityVec[j])
        
        # set.seed(42)  # Sets seed for reproducibility
        # Iter5000_Rtsne_Results[[j]] <- Rtsne(t(logcpm2), dims=3, max_iter=5000, perplexity=PerplexityVec[j])
      # }
      
        scatterplot3d(Iter1000_Rtsne_Results[[q]]$Y,color=AllSamplesColors, bg=AllSamplesColors, pch=20, cex.symbols=1, main=paste0("Perplexity ", PerplexityVec[q], " & Iterations 1000"),
                      xlab="", ylab="", zlab="") # pch=21,
        
        plot(Iter1000_Rtsne_Results[[q]]$Y[,1], 
             Iter1000_Rtsne_Results[[q]]$Y[,2], 
             col=AllSamplesColors, bg=AllSamplesColors, pch=20, cex=1.5, cex.axis = 1.5, main=paste0("Perplexity ", PerplexityVec[q], " & Iterations 1000"), # pch=21, col="white",
             xlab="", ylab="") #, xlim=c(-10, 10), ylim=c(-10,15))  
        
        plot(Iter1000_Rtsne_Results[[q]]$Y[,1],
             Iter1000_Rtsne_Results[[q]]$Y[,3], 
             col=AllSamplesColors, bg=AllSamplesColors, pch=20, cex=1.5, cex.axis = 1.5, main=paste0("Perplexity ", PerplexityVec[q], " & Iterations 1000"), #  pch=21, col="white", 
             xlab="", ylab="") #, xlim=c(-10, 10), ylim=c(-10,15))    
        
        plot(Iter1000_Rtsne_Results[[q]]$Y[,2], 
             Iter1000_Rtsne_Results[[q]]$Y[,3], 
             col=AllSamplesColors, bg=AllSamplesColors, pch=20, cex=1.5, cex.axis = 1.5, main=paste0("Perplexity ", PerplexityVec[q], " & Iterations 1000"), # pch=20, col="white", 
             xlab="", ylab="") #, xlim=c(-10, 10), ylim=c(-10,15))    
    }
  dev.off()
  
}
#



## plot diff peak RNA type (pie plot)
df <- data.frame(peak.id = unlist(sapply(strsplit(rownames(mat),"|",fixed=T),"[",4)),
                 RNA = unlist(sapply(strsplit(rownames(mat),"|",fixed=T),"[",2)))
df$RNA <- gsub("_rev|_for","",df$RNA,perl=T)
rna <- c("rRNA","pri_miRNA","lncRNA", "mRNA", "piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA')
dna <- c("intron", "promoter", "enhancer", "repeats")
df$RNA <- factor(df$RNA,levels = c(rna,dna))
table(df$RNA)

df <- df[df$peak.id %in% peak.ids[['top500']],] # top200:556; top500:1810

df2 <- as.data.frame(table(df$RNA))
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
  # scale_fill_manual(name="precursor",values = c(pal_nejm_adaptive()(15)[1:14],"#11838D") ) + 
  # scale_fill_nejm_adaptive(alpha = 0.8) +
  theme_bw() + 
  theme(
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
    legend.title= element_text(size= 24))
ggsave(paste0(dst,"_diff_pie.pdf"),width = 12,height = 7)





#diff volcano plot
dst <- "GSE71008_diff"
my_theme <-   theme(aspect.ratio=1,
                    plot.title = element_text(size = 28,color="black",hjust = 0.5),
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

#diff.plot <- list()
for (disease in c("Colorectal Cancer","Pancreatic Cancer","Prostate Cancer","Healthy Control")){
  #disease <- "Colorectal Cancer" # "Colorectal Cancer" "Healthy Control"   "Pancreatic Cancer" "Prostate Cancer"
  if(disease=="Colorectal Cancer"){
    disease.label <- "CRC"
  }else if(disease=="Pancreatic Cancer"){
    disease.label <- "PACA"
  }else if(disease=="Prostate Cancer"){
    disease.label <- "PRCA"
  }else if(disease=="Healthy Control"){
    disease.label <- "NC"
  }
  print(disease)
  # normal <- "Healthy Control"
  normal.label <- "R" #"NC"

  # tmp <- read.table(paste0("./output/",dst,"/exPeak_smallDomain_diff_",disease.label,"vs",normal.label,".diff"),sep = "\t",header = T,check.names = F)
  # tmp <- read.table(paste0("./output/",dst,"/filter112Smp_exPeakCNN_smallDomain_diff_",disease.label,"vs",normal.label,".diff"),sep = "\t",header = T,check.names = F,stringsAsFactors = F)
  tmp <- read.table(paste0("./output/",dst,"/allSmp_rmBatch_OvR_cfPeakCNN_smallDomain_diff_",disease.label,"vs",normal.label,".diff"),sep = "\t",header = T,check.names = F,stringsAsFactors = F)
  
  tmp$minus.log10.pvalue <- -log10(tmp$pvalue)
  tmp$log10.baseMean <- log10(tmp$baseMean)
  
  tmp$alpha <- maxmin.normalize.vec(tmp$minus.log10.pvalue)
  # ggplot(tmp, aes(x=log2FoldChange,y=minus.log10.pvalue)) + 
  #   geom_point() 
  # ggplot(tmp, aes(x=baseMean.pos,y=baseMean.neg)) + 
  #   geom_point() 
  p <- ggplot(tmp, aes(x=log10.baseMean,y=log2FoldChange, fill=minus.log10.pvalue)) + 
    geom_point(aes(alpha=alpha),shape=21,color="white") + 
    labs(title = paste0(disease.label," vs. ",normal.label)) +
    ggraph::scale_fill_viridis(direction = -1) +
    geom_hline(yintercept = 0,linetype="dashed",color="salmon") + 
    theme_bw()+
    my_theme
  ggsave(plot = p, paste0(dst,"_diff_basemean_",disease.label,".pdf"), width=9, height=9) # "_",sample
}






#################
# TCGA-small
#################
setwd("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/")
options(scipen = 4,stringsAsFactors = F,digits = 5)
dst <- "TCGA_small_diff3"

# sample.table <- read.table("./meta/lulab/FTC/sample_table.txt",sep = "\t",header = T)
sample.table0 <- read.table(paste0("../WCHSU-FTC/exSeek-dev/data/",dst,"/sample_table.txt"),sep = "\t",header = T)
#sample.table <- reshape2::melt(sample.table,id.var=c("patient_id","group","name"),variable.name="source",value.name = "data_id")
#sample.table$group <- substr(sample.table$patient_id,1,3)
#table(sample.table0$Project.ID) # TCGA-LAML, TCGA-PRAD, TCGA-PAAD, TCGA-COAD
sample.table0 <- sample.table0[,c(2,5)]
colnames(sample.table0)[c(1,2)] <- c("sample","group")
sample.table0$sample <- gsub(".bam","",sample.table0$sample)
sample.table0$group <- gsub("TCGA-","",sample.table0$group)
#sample.table <- sample.table[sample.table$sample!="",]

## read mat
count <- read.table(paste0("../WCHSU-FTC/output/",dst,"/call_peak_all/count_matrix/cfpeakCNN_b5_d50_p1_GSE71008.txt"),header = T,sep = "\t",stringsAsFactors = F,check.names = F)

rownames(count) <- count$feature
mat <- count[,2:ncol(count)]
#mat[1:3,1:3]
#colnames(mat) # 72

## diff

## mds/PCA plot
y <- DGEList(counts=mat)
# counts <- edgeR::getCounts(y)
y <- calcNormFactors(y, method="TMM")
# logcpm <- edgeR::cpm(y, log=TRUE, normalized.lib.sizes = T,prior.count = 1)
logcpm <- edgeR::cpm(y, log=F, normalized.lib.sizes = T)
logcpm <- log2(logcpm+1)
# logcpm <- read.table(paste0("./output/",dst,"/cfPeakCNN_GSE110381_smallDomain_diff_",disease.label,"vs",normal.label,".cpm"), header = T, sep="\t",check.names = F)
dim(logcpm)
rownames(logcpm) <- unlist(sapply(strsplit(rownames(logcpm),"|",fixed=T),"[",4))

sample.table0 <- sample.table0[sample.table0$sample %in% colnames(logcpm),]
table(sample.table0$group)
sample.table0$group <- factor(sample.table0$group, levels = c("LAML","COAD","PAAD","PRAD"))
sample.table0 <- sample.table0[order(sample.table0$group),]
logcpm <- logcpm[,sample.table0$sample]

#CategColors <- c("blue", "green", "gray", "purple", "orangered", "black", "brown", "yellow3")
PerplexityVec=seq(5, 70, by=20) # 70 too large, 20
PerplexityVec
Iter1000_Rtsne_Results=vector(mode="list", length=length(PerplexityVec))
Phenotype <- sample.table0$group[match(colnames(mat),sample.table0$sample)]
table(Phenotype)
# COAD LAML PAAD PRAD 
# 8   30    4   30
# LAML COAD PAAD PRAD 
# 60   60   60   60 
Colours4Classes_TCGA=c("grey", "salmon",  "steelblue", "purple") # "orangered", "brown", "black", "gray",
AllSamplesColors_TCGA=mapvalues(Phenotype, c("COAD","LAML","PAAD","PRAD"), Colours4Classes_TCGA)
AllSamplesColors <- as.character(AllSamplesColors_TCGA)

for (q in 1:length(PerplexityVec)) {
  #q <- 1
  pdf(paste0("tSNE_TCGA3_iter1000_perplexity", PerplexityVec[q], ".pdf"), width = 13,height = 3.5) # 
  par(mfrow=c(1,4))
  #par(mfrow=c(2,2))
  
  
  for (top in nums){
    #top <- 50
    print(top)
    logcpm2 <- logcpm[peak.ids[[paste0("top",top)]],sample.table0$sample]
    print(dim(logcpm2))
    
    j <- q
    # for (j in 1:length(PerplexityVec)) {
    #   print(j)
    set.seed(42)  # Sets seed for reproducibility
    Iter1000_Rtsne_Results[[j]] <- Rtsne(t(logcpm2), dims=3, max_iter=1000, perplexity=PerplexityVec[j])
    
    # set.seed(42)  # Sets seed for reproducibility
    # Iter3000_Rtsne_Results[[j]] <- Rtsne(t(logcpm2), dims=3, max_iter=3000, perplexity=PerplexityVec[j])
    
    # set.seed(42)  # Sets seed for reproducibility
    # Iter5000_Rtsne_Results[[j]] <- Rtsne(t(logcpm2), dims=3, max_iter=5000, perplexity=PerplexityVec[j])
    # }
    
    scatterplot3d(Iter1000_Rtsne_Results[[q]]$Y,color=AllSamplesColors, bg=AllSamplesColors, pch=20, cex.symbols=1, main=paste0("Perplexity ", PerplexityVec[q], " & Iterations 1000"),
                  xlab="", ylab="", zlab="") # pch=21,
    
    plot(Iter1000_Rtsne_Results[[q]]$Y[,1], 
         Iter1000_Rtsne_Results[[q]]$Y[,2], 
         col=AllSamplesColors, bg=AllSamplesColors, pch=20, cex=1.5, cex.axis = 1.5, main=paste0("Perplexity ", PerplexityVec[q], " & Iterations 1000"), # pch=21, col="white",
         xlab="", ylab="") #, xlim=c(-10, 10), ylim=c(-10,15))  
    
    plot(Iter1000_Rtsne_Results[[q]]$Y[,1],
         Iter1000_Rtsne_Results[[q]]$Y[,3], 
         col=AllSamplesColors, bg=AllSamplesColors, pch=20, cex=1.5, cex.axis = 1.5, main=paste0("Perplexity ", PerplexityVec[q], " & Iterations 1000"), #  pch=21, col="white", 
         xlab="", ylab="") #, xlim=c(-10, 10), ylim=c(-10,15))    
    
    plot(Iter1000_Rtsne_Results[[q]]$Y[,2], 
         Iter1000_Rtsne_Results[[q]]$Y[,3], 
         col=AllSamplesColors, bg=AllSamplesColors, pch=20, cex=1.5, cex.axis = 1.5, main=paste0("Perplexity ", PerplexityVec[q], " & Iterations 1000"), # pch=20, col="white", 
         xlab="", ylab="") #, xlim=c(-10, 10), ylim=c(-10,15))     
  }
  
  dev.off()
}

