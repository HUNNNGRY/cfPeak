setwd("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/")
# .libPaths(c("/home/baopengfei/R/x86_64-redhat-linux-gnu-library/3.6","/usr/lib64/R/library","/usr/share/R/library")) # ,"/BioII/lulab_b/baopengfei/R/x86_64-redhat-linux-gnu-library/3.6"
# .libPaths()
options(scipen = 4,stringsAsFactors = F,digits = 5)
options(bedtools.path = "/BioII/lulab_b/baopengfei/biosoft")
library(bedtoolsr)
library(dplyr)
#library(magrittr)
library(conflicted)
#unloadNamespace("conflicted")
conflict_prefer("colnames", "base")
conflict_prefer("count", "dplyr")
conflict_prefer("mutate", "dplyr")

# load all func.
source("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/scripts/util.R")
source("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/scripts/lulab/sncRNA_utils.R")
source("/BioII/lulab_b/baopengfei/shared_utils/stats.R")

#Sys.setlocale(category = "LC_ALL", locale = "en_US.UTF-8")
#Sys.getlocale(category = "LC_ALL")


#bedtools intersect -split -u -f 0.5 -s -wa -a ../WCHSU-FTC/output/SLE/call_peak_dedup/cfpeakCNN/b5_d50_p1_gn.bed -b <(cat output/lulab/*NAsmallDomain_EVenrich_gn.bed | sort -k5,5n  | tail -n 100 | LC_COLLATE=C sort -k1,1 -k2,2n) | LC_COLLATE=C sort -k1,1 -k2,2n > ../WCHSU-FTC/output/SLE/call_peak_dedup/cfpeakCNN/b5_d50_p1_gn_intersectEV_top100EV.bed
#bedtools intersect -split -u -f 0.5 -s -wa -a ../WCHSU-FTC/output/SLE/call_peak_dedup/cfpeakCNN/b5_d50_p1_gn.bed -b <(cat output/lulab/*NAsmallDomain_EVenrich_gn.bed | LC_COLLATE=C sort -k1,1 -k2,2n) | LC_COLLATE=C sort -k1,1 -k2,2n > ../WCHSU-FTC/output/SLE/call_peak_dedup/cfpeakCNN/b5_d50_p1_gn_intersectEV.bed
#bedtools intersect -split -f 0.5 -s -wa -wb -a ../WCHSU-FTC/output/SLE/call_peak_dedup/cfpeakCNN/b5_d50_p1_gn.bed -b <(cat output/lulab/*NAsmallDomain_EVenrich_gn.bed | LC_COLLATE=C sort -k1,1 -k2,2n) | LC_COLLATE=C sort -k1,1 -k2,2n > ../WCHSU-FTC/output/SLE/call_peak_dedup/cfpeakCNN/b5_d50_p1_gn_intersectEV.txt
plasmaEV <- read.table("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/lulab/smallDomain_EVenrich_gn.bed",header = F,sep = '\t')

SLE.EV.top100 <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/SLE/call_peak_dedup/cfpeakCNN/b5_d50_p1_gn_intersectEV_top100EV.bed",header = F,sep = '\t')
SLE.EV <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/SLE/call_peak_dedup/cfpeakCNN/b5_d50_p1_gn_intersectEV.bed",header = F,sep = '\t')
#table(SLE.EV$V4 %in% peak.id)
#SLE.EV <- SLE.EV[SLE.EV$V4 %in% peak.id,]
plasma.SLE.EV <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/SLE/call_peak_dedup/cfpeakCNN/b5_d50_p1_gn_intersectEV.txt",header = F,sep = '\t')
plasma.SLE.EV <- plasma.SLE.EV[,c("V4","V16")]
plasma.SLE.EV <- plasma.SLE.EV[!grepl(",",plasma.SLE.EV$V16),]
colnames(plasma.SLE.EV) <- c("SLE","plasma")


# read sncRNA/peak tab/mat ------------
dst <- "SLE"
dedup <- "_dedup"
diffpre <- "/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA"
pre="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
outpre="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
dir.create(paste0(diffpre,"/output/",dst),recursive = T,showWarnings = F) # ,"/",disease.label,"vs",normal.label

sample.table0 <- read.csv( paste0(pre,"/exSeek-dev/data/SLE/240316/Anno_0315.csv.bak"),check.names = F,sep = ",",header = T,encoding="UTF-8")
colnames(sample.table0)[1:2] <- c("sample","patient")
rownames(sample.table0) <- sample.table0$sample
sample.table0 <- sample.table0[sample.table0$pass=="TRUE",]
sample.table0 <- sample.table0[!grepl("J2",sample.table0$sample),] # 1691 pass
table(sample.table0$group)
sample.table0 <- sample.table0[sample.table0$group %in% c("HDA","Healthy_control","Inactive","LDA","MDA"),]
# Disease_control             HDA Healthy_control        Inactive             LDA             MDA         NO_INFO        NoDNaseI 
# 49             371             496              80             270             415               1               9 
sample.table0$group <- gsub("Healthy_control","HC",sample.table0$group)
table(sample.table0$cell)
#EV, PBMC, Platelet, Plasma
table(sample.table0$lineage)

sample.table0 <- sample.table0[sample.table0$lineage %in% c("Plasma" ,"EV"),] #
# B        CD4        CD8         DC         EV   Monocyte Neutrophil         NK       PBMC     Plasma   Platelet 
# 305        556        250        125         55        182        118         61          1         19         18
# tmp <- sample.table0$sample


## append operator et al
tmp <- read.csv( paste0(pre,"/exSeek-dev/data/SLE/240316/Anno_0315_utf8.csv"),check.names = F,sep = ",",header = T,encoding="UTF-8")
#colnames(tmp)
colnames(tmp)[1:2] <- c("sample","patient")
update.fields <- c("Age","Gender","SLEDAI","Anti-SSA","Organ","Threatment","Threatment Response","treatmentnaive","Prednisolone (mg/day), median (IQR)","Hydroxychloroquine","Mycofenolate mofetil","Tacrolimus","Methotrexate","Experiment Date","Operator")
tmp <- tmp[,c("sample",update.fields)]
#colnames(tmp) <- c("sample",)
sample.table0 <- dplyr::left_join(sample.table0[,1:6],tmp,by="sample")
rownames(sample.table0) <- sample.table0$sample


# table(sample.table0$Gender)
sample.table0$Gender <- gsub("female","Female",sample.table0$Gender,perl = T)
sample.table0$Gender <- gsub("male","Male",sample.table0$Gender,perl = T)
table(sample.table0$Gender)

sample.table0$Month <- substr(sample.table0$`Experiment Date`,1,6)
sample.table0$Month2[sample.table0$Month %in% c("202306","202307")] <- "202306"
sample.table0$Month2[sample.table0$Month %in% c("202308","202309")] <- "202308"
sample.table0$Month2[sample.table0$Month %in% c("202310","202311","202312")] <- "202310"
table(sample.table0$Month2)
table(sample.table0$Month2,sample.table0$Operator)
sample.table0$Month2.Operator <- paste0(sample.table0$Month2,"_",sample.table0$Operator)
table(sample.table0$Month2.Operator )
colnames(sample.table0)[colnames(sample.table0)=="Prednisolone (mg/day), median (IQR)"] <- "Prednisolone"


table(duplicated(sample.table0$patient))
sample.table0 <- sample.table0[sample.table0$patient %in% unique(sample.table0$patient[duplicated(sample.table0$patient)]),]
sample.table0 <- sample.table0[sample.table0$group=="HC",]

table(sample.table0$lineage)
table(duplicated(sample.table0$patient))

## table2: read count mat
count0 <- data.table::fread( paste0(outpre,"/output/",dst,"/call_peak",dedup,"/count_matrix/cfpeakCNN_b5_d50_p1.txt"), header = T,sep = "\t",stringsAsFactors = F,check.names = F, data.table = F)
rownames(count0) <- count0$feature
count0$feature <- NULL
count0 <- count0[,sample.table0$sample]




peak.table <- data.table::fread(paste0(outpre,"/output/",dst,"/call_peak",dedup,"/cfpeakCNN/b5_d50_p1.summary"),header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",data.table = F)
rownames(peak.table) <- peak.table$name

libSizePath <- paste0(outpre,"/output/",dst,"/call_peak",dedup,"/count_matrix/EM.libSize")
tmp <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/SLE/240316/QC_0315.txt",header = T,sep = "\t",check.names = F,row.names = 1)
# tmp <- tmp[colnames(count),]
tmp <- tmp[sample.table0$sample,]
externalLibSize <- tmp[,"star_hg38_v38_dedup"]
#names(externalLibSize) <- colnames(count) # [2:ncol(count)]
names(externalLibSize) <- sample.table0$sample # [2:ncol(count)
#



# diff of sncRNA ----------------------------
table(sample.table0$patient)
this.patients <- c("HD18","HD19","HD20") #c("HD29","HD30","HD31") # "HD18","HD19","HD20"
sample.table.pair <- sample.table0[(sample.table0$patient %in% this.patients),] 
sample.table.pair$patient_id <- sample.table.pair$patient
sample.table.pair$data_id <- sample.table.pair$sample
sample.table.pair$source <- sample.table.pair$cell
length(unique(sample.table.pair$patient_id))
positive_samples <- sample.table.pair[sample.table.pair$source=="EV","data_id"]
negative_samples <- sample.table.pair[sample.table.pair$source=="Plasma","data_id"]
samples <- c(positive_samples, negative_samples)
sample.table.pair <- sample.table.pair[match(samples,sample.table.pair$data_id),]
table(sample.table.pair$cell,sample.table.pair$Month2.Operator)
sample.table.pair$group <- sample.table.pair$cell

mat <- count0[,samples]
dim(mat)

#colnames(mat)
group <- c(rep("positive",length(positive_samples)),rep("negative",length(negative_samples)))
sample.table.pair$patient_id <- factor(sample.table.pair$patient_id,levels = unique(as.character(sample.table.pair$patient_id)))
patient <- sample.table.pair$patient_id[(sample.table.pair$source=="EV" | sample.table.pair$source=="Plasma")]
table(patient)
length(patient)
method <- "edger_glmlrt"
norm_method <- "TMM"
#res <- diff(mat = mat, samples = samples, group = group, patient = patient, method = method, norm_method = norm_method, filterType = "small")
form <- ~ group + patient #+ Month + Operator
res <- diff.v3.dcb(mat = mat, sample.table = sample.table.pair, form = form, useExternalLibSize = externalLibSize[samples], samples = samples, group = group, method = method, norm_method = norm_method, filterType = "NULL",featureType = "domain")
#write.table(res[["diffTable"]],paste0("./output/SLE/","diff_EVvsCF_pair.txt"),quote = F,sep = "\t",row.names = T,col.names = T)
write.table(res[["diffTable"]],paste0("./output/SLE/","diff_EVvsCF_pair3.txt"),quote = F,sep = "\t",row.names = T,col.names = T)

# tmp <- res[["diffTable"]]
tmp <- read.table(paste0("./output/SLE/","diff_EVvsCF_pair3.txt"))
#filter <- tmp$pvalue<0.1
filter <- tmp$padj<0.01 & tmp$log2FoldChange>0.8
peak.id <- unlist(sapply(strsplit(rownames(tmp[filter,]),"|",fixed=T),"[",4))
length(peak.id)
length(SLE.EV$V4)
num.real <- length(unique(intersect(peak.id,SLE.EV$V4)))
num.real

perm.id.list <- list()
inter.list <- list()
for (i in 1:500){
  set.seed(i)
  perm.id <- unlist(sapply(strsplit(rownames(tmp[sample(1:nrow(tmp),sum(filter),replace = F),]),"|",fixed=T),"[",4))
  print(length(unique(intersect(perm.id,SLE.EV$V4))))
  perm.id.list[[i]] <- perm.id
  inter.list[[i]] <- length(unique(intersect(perm.id,SLE.EV$V4)))
}
num.distribution <- do.call(c,inter.list)
hist(num.distribution)

#Calculate FDR
num.greater_or_equal <- round(sum(num.distribution >= num.real),digits = 6)
fdr <- round((num.greater_or_equal) / length(num.distribution),digits = 6)
print(sprintf("FDR: %.18f", fdr))
#****

# Create a pretty ggplot histogram
df <- data.frame(num_distribution = num.distribution)
ggplot(df, aes(x = num_distribution)) +
  geom_histogram(binwidth = 1, fill = "steelblue4", color = "black", alpha = 0.7) +
  geom_vline(aes(xintercept = num.real), color = "red", linetype = "dashed", size = 1) +
  labs(title = "Histogram of Permuted Intersections",
       x = "Number of Intersections",
       y = "Frequency",
       subtitle = sprintf("****")) +
  theme_minimal() +
  theme(plot.title = element_text(size = 24,color="black",hjust = 0.5,face="bold"),
        axis.title = element_text(size = 24,color ="black"),
        axis.text = element_text(size= 20,color = "black"),
        plot.subtitle = element_text(hjust = 0.5,size=30),
        #panel.grid=element_blank(),
        # panel.grid.major.x=element_blank(),
        # panel.grid.minor.x = element_blank(),
        # panel.grid.major.y = element_line(color = "grey50",linetype = "dashed"), #size= 1,
        #panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(), #   , color = c(rep("steelblue",length(rna)),rep("firebrick",length(dna)))
        # legend.position = "right",#c(.25,.6),
        legend.text = element_text(size= 16),
        legend.title= element_text(size= 16))
ggsave("./figure/EV_inter.png",height = 5, width = 7)

# # comp log2FC cor (depre, no plasma FC for merged EV)
# rownames(tmp) <- unlist(sapply(strsplit(rownames(tmp),"|",fixed=T),"[",4))
# plasma.SLE.EV <- plasma.SLE.EV[plasma.SLE.EV$SLE %in% rownames(tmp),]
# tmp <- tmp[plasma.SLE.EV$SLE,]
# 
# plasma.diff <- read.table(paste0("./output/lulab/","expeakCNN","_smallDomain_diff_EVvsCF_pair.txt"))
# rownames(plasma.diff) <- unlist(sapply(strsplit(rownames(plasma.diff),"|",fixed=T),"[",4))
# plasma.diff <- plasma.diff[plasma.SLE.EV$plasma,]
# 
# cor.test(plasma.diff$log2FoldChange,tmp$log2FoldChange)
# #



# plot IGV: multi sample bigwig ------------------------
sample.table.pair <- sample.table0[(sample.table0$patient %in% this.patients),]
table(sample.table.pair$cell)
sample.table.pair$group <- sample.table.pair$cell
sample.table.pair$cell <- "plasmaEV"


library(IRanges)
library(GenomicRanges)
# read candidate summary table
#dst <- "SLE"
#SLE.EV.top100$V4
#"peak_18076", "peak_8701",  "peak_10252", "peak_11103", "peak_27366", "peak_11114", "peak_11107"
candidate.ids <- c("peak_27782", "peak_9935",  "peak_28110", "peak_16095", "peak_11082", "peak_24114", "peak_9923",  "peak_11104", "peak_6346",  "peak_27755", "peak_26334","peak_18076", "peak_8701",  "peak_10252", "peak_11103", "peak_27366", "peak_11114", "peak_11107")  #c("peak_9935", "peak_16095", "peak_11082") # 
this.cells <- c("plasmaEV")

#sncCandidate <- read.table( paste0(diffpre,"/output/",dst,"/candidate/",candidate.class,".txt"), header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "")
sncCandidate <- peak.table[peak.table$name %in% candidate.ids, ]
# peakMaxY <- list()
# peakMaxY[["peak_9935"]] <- 12
# peakMaxY[["peak_16095"]] <- 16
# peakMaxY[["peak_11082"]] <- 80
for(cellID in c(this.cells)){
  # cellID <- "Plasma"
  # this.cells
  #cellID <- "Naive.CD8.T" #cellShort
  # sncCandidate <- read.table( paste0(diffpre, "/output/",dst,"/",cellID,"_candidate_diff_summary.txt"),header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "")
  # sncCandidate <- sncCandidate[!(is.na(sncCandidate$tx)),]
  # sncCandidate <- sncCandidate[sample(1:nrow(sncCandidate),10,replace = F),]
  # #outFile <- paste0(outpre,"/output/",dst,"/call_peak",dedup,"/cfpeakCNN/b5_d50_p1.summary")
  
  sample.table <- sample.table.pair[sample.table.pair$cell==cellID,] # 1042
  #table(sample.table0$group,sample.table0$treatmentnaive) # seems all inactive was treated
  
  # prepare meta.list
  meta.col <- list()
  genes.list <- list()
  #annotation_row_merge <- read.table(file = "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/tmp/PRNNA-TCGA3_annotation_row_merge.txt",header = T,sep = "\t",check.names = F,stringsAsFactors = F)
  annotation_row_merge <- sncCandidate[,c("tx","start","end","name","tx_len")]
  colnames(annotation_row_merge)[c(1,4)] <- c("chr","peak")
  annotation_row_merge$width <- annotation_row_merge$end-annotation_row_merge$start
  expFold <- 2
  
  annotation_row_merge$exp_start <- 0
  annotation_row_merge$exp_end <- annotation_row_merge$tx_len
  annotation_row_merge$exp_width <- 0
  for(i in 1:nrow(annotation_row_merge)){
    annotation_row_merge$exp_start[i] <- max(0,annotation_row_merge$start[i]-expFold*annotation_row_merge$width[i])
    annotation_row_merge$exp_end[i] <- min(annotation_row_merge$tx_len[i],annotation_row_merge$end[i]+expFold*annotation_row_merge$width[i])
    annotation_row_merge$exp_width[i] <- annotation_row_merge$exp_end[i]-annotation_row_merge$exp_start[i]
    genes.list[[i]] <- data.frame(chr=annotation_row_merge[i,"chr"],start=annotation_row_merge[i,"exp_start"],width=annotation_row_merge[i,"exp_width"],peak=annotation_row_merge[i,"peak"],peak.start=annotation_row_merge[i,"start"],peak.width=annotation_row_merge[i,"width"]) #mir
  }
  
  
  # prepare IGV meta table
  randNum <- 3
  library(tidyverse)
  library(tidydr)
  #conflict_prefer("select", "dplyr")
  
  # dataset: >18 cell types
  dsts.list <- list()
  
  #smp.tbl <- read.table(paste0(pre,"/exSeek-dev/data/",dst,"/sample_table_202311.txt"),header = T, sep = "\t",quote = "",check.names = F,stringsAsFactors = F)
  smp.tbl <- sample.table
  smp.tbl <- smp.tbl[!(is.na(smp.tbl$group)),]
  #colnames(smp.tbl)
  #smp.tbl[1:3,]
  #colnames(smp.tbl)[c(1,4,5)] <- c("sample","group")
  smp.tbl <- smp.tbl %>% 
    dplyr::select(c("sample","group")) %>% 
    # dplyr::mutate(group = case_when(group=="No" ~ 'Localized',
    #                                 group=="Yes" ~ 'Metastasis')) %>% 
    dplyr::group_by(group) %>% 
    dplyr::sample_n(randNum)
  smp.tbl$bwPath <- paste0(outpre,"/output/",dst,"/call_peak",dedup,"/tbigwig_RNA_EM/",smp.tbl$sample,".transcriptome.bigWig")
  smp.tbl$dataset <- cellID
  #table(smp.tbl$origin) # Buccal Mucosa Floor of mouth    Oral Cavity    Oral Tongue
  dsts.list[[dst]] <- smp.tbl
  
  
  dsts.df <- do.call(rbind,dsts.list)
  #dsts.df$dataset <- factor(dsts.df$dataset,levels = cellIDshort)
  groupShort <- c("Plasma","EV") #c("HC","Inactive","LDA","MDA","HDA")
  groupColor <- c("Plasma"="steelblue","EV"="salmon") 
  dsts.df$group <- factor(dsts.df$group,levels = groupShort)
  dsts.df <- dsts.df[order(dsts.df$dataset,dsts.df$group),] # bw order !
  #table(dsts.df$dataset)
  #table(dsts.df$group)
  
  
  # add color each dst group
  dsts <- unique(as.character(dsts.df$dataset)) # pDC
  #cols <- c("seagreen4","steelblue","salmon","salmon3","salmon4") # groupShortColors
  cols <- groupColor
  #names(cols) <- groupShort
  dsts.df$color <- "grey"#FA80721A
  #j2 <- 1
  for (i in 1:length(dsts)){
    tmp.df <- dsts.df[dsts.df$dataset==dsts[i],]
    #tmp.df <- tmp.df[order(tmp.df$group),]
    grps <- unique(as.character(tmp.df$group))
    #print(grps)
    for (j in 1:length(grps)){
      dsts.df$color[dsts.df$dataset==dsts[i] & dsts.df$group==grps[j]] <- cols[grps[j]]
      # j2 <- j2+1
    }
  }
  #peak_5859
  
  
  
  for (i in 1:length(genes.list)){
    # i <- 1 #"ENST00000385214_____1_30_49_+|pri_miRNA|ENST00000385214_____1|peak_6181|ENST00000385214_____1|30|49"
    chr <- genes.list[[i]]$chr
    start <- as.numeric(genes.list[[i]]$start)
    width <- as.numeric(genes.list[[i]]$width)
    peak.start <- as.numeric(genes.list[[i]]$peak.start)
    peak.width <- as.numeric(genes.list[[i]]$peak.width)
    peakID <- genes.list[[i]]$peak
    print(paste0(chr,".",start,".",width))
    
    res.list <- list()
    for (dst.tmp in dsts){ # names(meta.list)
      smps <- dsts.df$sample[dsts.df$dataset==dst.tmp]
      grps <- dsts.df$group[dsts.df$dataset==dst.tmp]    
      for (j in 1:length(smps)){
        #j <- 6
        smp <- smps[j]
        grp <- grps[j]
        # print(smp)
        # print(grp)
        bw.color <- dsts.df$color[dsts.df$dataset==dst.tmp & dsts.df$sample==smp]
        plotXbreak <- ifelse(dst.tmp==dsts[length(dsts)] & smp==smps[length(smps)], TRUE, FALSE)
        dst.lab <- dst.tmp
        bw.res.df <- readBW(bw=dsts.df$bwPath[dsts.df$dataset==dst.tmp & dsts.df$sample==smp],chr=chr,start=start,width=width)
        # title <- ifelse(dst.tmp==names(meta.list)[1] & smp==smps[1], 
        #                 paste0(chr,":",start,"-",start+width),
        #                 "")
        title <- "" #paste0(chr,":",start,"-",start+width)#""
        res.list[[paste0(dst.lab,"_",smp)]] <- plotBW(single_base_bw=bw.res.df,chr=chr,start=start,width=width,highlight=c(peak.start,peak.start+peak.width),ylab="",title=title,plotXbreak = plotXbreak, plotYbreak = T, plotYtext = F, color = bw.color, fill = bw.color , libsize = externalLibSize[smp]) # maxY=20, ,  # ,maxY=peakMaxY[[peakID]]
      }
    }
    ## collect
    tmp <- ggpubr::ggarrange(plotlist = res.list, ncol = 1,  align = "hv") + #heights = c(7,1,1,1,1,1,1,1,1,1,1), nrow = length(res.list),
      ggtitle(label = paste0(c(dsts,unique(grps))))
    # ggsave(plot = tmp,filename = paste0(diffpre,"/output/",dst,"/candidate/",peakID,"_",cellID,"_bw","-",chr,".",start,".",width,".pdf"),width = 15, height = 1.3*length(res.list)) # 15
    ggsave(plot = tmp,filename = paste0("./figure/",peakID,"_",cellID,"_bw","-",chr,".",start,".",width,".png"),width = 15, height = 1.3*length(res.list)) # 15
  }
}
#
#
