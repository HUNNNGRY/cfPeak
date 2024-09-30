setwd("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/")
options(scipen = 4,stringsAsFactors = F,digits = 5)

# load all func.
source("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/scripts/util.R")
source("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/scripts/lulab/sncRNA_utils.R")
setwd("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/")
options(scipen = 4,stringsAsFactors = F,digits = 5)

## Cell free RNA peaks
# Figure 6E
dst <- "GSE71008_diff"
# sample.table <- read.table("./meta/lulab/FTC/sample_table.txt",sep = "\t",header = T)
sample.table0 <- read.table("../WCHSU-FTC/exSeek-dev/data/GSE71008/sample_table.txt",sep = "\t",header = T)
#sample.table <- reshape2::melt(sample.table,id.var=c("patient_id","group","name"),variable.name="source",value.name = "sample")
#sample.table$group <- substr(sample.table$patient_id,1,3)
colnames(sample.table0)[c(1,3)] <- c("sample","group")
#sample.table <- sample.table[sample.table$sample!="",]
table(sample.table0$group)

# # filter sample by batch/QC
# #SAMN union of (bottom 30% qc rank &) bottom 30% piRNA ratio & top 30% miRNA ratio
# cutoff <- 0.7 # 0, 0.4, 0.5, 0.6, 0.7, 0.82
# qc <- read.table("../WCHSU-FTC/exSeek-dev/data/GSE71008_diff/qc.txt",header = T)
# fail.mir <- qc[order(qc$pri_miRNA,decreasing = T),]
# fail.mir <- fail.mir$sample[1:(nrow(fail.mir)*cutoff)]
# fail.pir <- qc[order(qc$piRNA,decreasing = F),]
# fail.pir <- fail.pir$sample[1:(nrow(fail.pir)*cutoff)]
# fail <- intersect(fail.pir,fail.mir) #0, 49, 66, 88, 112, 136 smps

## read mat
count <- read.table(paste0("../WCHSU-FTC/output/",dst,"/call_peak_all/count_matrix/expeakCNN_b5_d50_p1.txt"),header = T,sep = "\t",stringsAsFactors = F,check.names = F)

rownames(count) <- count$feature
mat <- count[,2:ncol(count)]
# mat <- mat[,!(colnames(mat) %in% fail)]


## diff
#sample.table <- sample.table[!(sample.table$patient_id %in% c("FTA-11","FTA-12","FTA-13","FTA-14","FTA-15","FTA-16","FTA-18","FTA-19","FTA-20","FTA-23","FTA-24","FTA-9",     "FTA-13","FTA-17","FTA-4","FTA-6")),]
#sample.table$sample <- unlist(sapply(strsplit(sample.table$sample,"_",fixed=T),"[",1)) # trim samples surfix !
#colnames(mat) <- unlist(sapply(strsplit(colnames(mat),"_",fixed=T),"[",1)) # trim samples surfix !
#write.table(sample.table0, file = "../WCHSU-FTC/exSeek-dev/data/GSE71008/sample_table_filterBatch.txt", sep = "\t", quote = F, row.names = F, col.names = T)
#sample.table0[1:3,1:3]

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
  
  # write.table(res[["normMat"]],paste0("./output/",dst,"/filterSmp_exPeak_smallDomain_diff_",disease.label,"vs",normal.label,".cpm"),quote = F,sep = "\t",row.names = T,col.names = T)
  # write.table(res[["diffTable"]],paste0("./output/",dst,"/filterSmp_exPeak_smallDomain_diff_",disease.label,"vs",normal.label,".diff"),quote = F,sep = "\t",row.names = T,col.names = T)
  # #piranha,clipper,clam,localmax,exPeak
  #tmp <- read.table(paste0("./output/",dst,"/filterSmp_exPeakCNN_smallDomain_diff_",disease.label,"vs",normal.label,".diff"),sep = "\t",header = T,check.names = F,stringsAsFactors = F)
  # tmp <- read.table(paste0("./output/",dst,"/filter49Smp_exPeakCNN_smallDomain_diff_",disease.label,"vs",normal.label,".diff"),sep = "\t",header = T,check.names = F,stringsAsFactors = F)
  # tmp <- read.table(paste0("./output/",dst,"/filter0Smp_exPeakCNN_smallDomain_diff_",disease.label,"vs",normal.label,".diff"),sep = "\t",header = T,check.names = F,stringsAsFactors = F)
  # tmp <- read.table(paste0("./output/",dst,"/filter66Smp_exPeakCNN_smallDomain_diff_",disease.label,"vs",normal.label,".diff"),sep = "\t",header = T,check.names = F,stringsAsFactors = F)
  #tmp <- read.table(paste0("./output/",dst,"/filter88Smp_exPeakCNN_smallDomain_diff_",disease.label,"vs",normal.label,".diff"),sep = "\t",header = T,check.names = F,stringsAsFactors = F)
  # tmp <- read.table(paste0("./output/",dst,"/filter136Smp_exPeakCNN_smallDomain_diff_",disease.label,"vs",normal.label,".diff"),sep = "\t",header = T,check.names = F,stringsAsFactors = F)
  # tmp <- read.table(paste0("./output/",dst,"/filter112Smp_exPeakCNN_smallDomain_diff_",disease.label,"vs",normal.label,".diff"),sep = "\t",header = T,check.names = F,stringsAsFactors = F)
  # tmp <- read.table(paste0("./output/",dst,"/filter112Smp_OvR_exPeakCNN_smallDomain_diff_",disease.label,"vs",normal.label,".diff"),sep = "\t",header = T,check.names = F,stringsAsFactors = F)
  tmp <- read.table(paste0("./output/",dst,"/allSmp_rmBatch_OvR_exPeakCNN_smallDomain_diff_",disease.label,"vs",normal.label,".diff"),sep = "\t",header = T,check.names = F,stringsAsFactors = F)
    
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
logcpm <- read.table(paste0("./output/",dst,"/allSmp_rmBatch_OvR_exPeakCNN_smallDomain_diff_",disease.label,"vs",normal.label,".cpm"), header = T, sep="\t",check.names = F)
max(logcpm[,3]) # log2CPM+1
# logcpm <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_all/count_matrix/expeakCNN_b5_d50_p1_CPMrowsum.txt"), header = T, sep="\t",check.names = F, row.names = 1)
rownames(logcpm) <- unlist(sapply(strsplit(rownames(logcpm),"|",fixed=T),"[",4))
sample.table0 <- sample.table0[sample.table0$sample %in% colnames(logcpm),]
sample.table0$group <- factor(sample.table0$group, levels = c("Healthy Control","Colorectal Cancer","Prostate Cancer","Pancreatic Cancer"))
sample.table0 <- sample.table0[order(sample.table0$group),] # draw point in group order !!!!
table(sample.table0$group)
logcpm <- logcpm[,sample.table0$sample]

## run PCA&UMAP (2408)
logcpm0 <- logcpm
rownames(sample.table0) <- sample.table0$sample
for (top in nums){
  #top <- 50
  print(top)
  logcpm2 <- logcpm[peak.ids[[paste0("top",top)]],sample.table0$sample]
  print(dim(logcpm2))

  sample.table <- sample.table0
  table(sample.table$group)
  sample.table$group <- gsub("Healthy Control","NC",sample.table$group)
  sample.table$group <- gsub("Colorectal Cancer","CRC",sample.table$group)
  sample.table$group <- gsub("Prostate Cancer","PRAD",sample.table$group)
  sample.table$group <- gsub("Pancreatic Cancer","PAAD",sample.table$group)
  groupColor <- c("NC"="grey", "CRC"="salmon", "PRAD"="purple", "PAAD"="steelblue")
  sample.table$Color <- sample.table$group
  sample.table$Fill <- sample.table$group
  sample.table$Shape <- factor(21,levels = 1:21) #sample.table$group
  # plotPCA(sample.table = sample.table, logcpm = logcpm2, outFile = paste0("./output/",dst,"/PCA/allSmp_cfPeakCNN_smallDomain_diff.cpm.PCA.","top",top,".pdf"), 
  #        plotGrpColor = groupColor, topVarNum = top, plotUMAP = T, plotDotSize = 2, plotDotAlpha = 0.9, figWid = 8, figHigh = 5)
  # Redirct output
  plotPCA(sample.table = sample.table, logcpm = logcpm2, outFile = paste0("/BioII/lulab_b/wangtaiwei/peak_calling/review/output/",dst,"/PCA/allSmp_cfPeakCNN_smallDomain_diff.cpm.PCA.","top",top,".pdf"), 
          plotGrpColor = groupColor, topVarNum = top, plotUMAP = T, plotDotSize = 2, plotDotAlpha = 0.9, figWid = 8, figHigh = 5)
  #
}
#

## compare within vs. inter- gorup distance from UMAP distance
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)

calculate_distances <- function(umap_data, metadata) {
  # Merge the umap_data and metadata
  merged_data <- merge(umap_data, metadata, by = "sample")
  
  # Calculate pairwise Euclidean distances
  dist_matrix <- as.matrix(dist(merged_data[, c("UMAP_Dim1", "UMAP_Dim2")]))
  rownames(dist_matrix) <- merged_data$sample
  colnames(dist_matrix) <- merged_data$sample
  
  # Initialize result lists
  within_group_distances <- list()
  inter_group_distances <- list()
  
  # Get unique groups
  groups <- unique(merged_data$group)
  groups <- as.character(groups)
  
  # Calculate within-group distances
  for (group in groups) {
    group_samples <- merged_data %>% filter(group == !!group) %>% dplyr::select(sample)
    group_distances <- dist_matrix[group_samples$sample, group_samples$sample]
    within_group_distances[[group]] <- mean(group_distances[upper.tri(group_distances)])
  }
  
  # Calculate inter-group distances
  for (i in 1:(length(groups) - 1)) {
    for (j in (i + 1):length(groups)) {
      group1_samples <- merged_data %>% filter(group == !!groups[i]) %>% dplyr::select(sample)
      group2_samples <- merged_data %>% filter(group == !!groups[j]) %>% dplyr::select(sample)
      inter_distances <- dist_matrix[group1_samples$sample, group2_samples$sample]
      # Ensure the order of groups is consistent
      if (groups[i] < groups[j]) {
        inter_group_distances[[paste(groups[i], groups[j], sep = "_vs_")]] <- mean(inter_distances)
      }
      else {
        inter_group_distances[[paste(groups[j], groups[i], sep = "_vs_")]] <- mean(inter_distances)
      }    
    }
  }
  
  # Combine results into data frames
  within_group_df <- data.frame(group = names(within_group_distances), avg_within_distance = unlist(within_group_distances))
  inter_group_df <- data.frame(group_pair = names(inter_group_distances), avg_inter_distance = unlist(inter_group_distances))
  
  return(list(within_group = within_group_df, inter_group = inter_group_df))
}

top <- 500
# umap_data <- read.csv(paste0("./output/",dst,"/PCA/allSmp_cfPeakCNN_smallDomain_diff.cpm.PCA.","top",top,".pdf.UMAP.txt"),sep = "\t")
umap_data <- read.csv(paste0("/BioII/lulab_b/wangtaiwei/peak_calling/review/output/",dst,"/PCA/allSmp_cfPeakCNN_smallDomain_diff.cpm.PCA.","top",top,".pdf.UMAP.txt"),sep = "\t")
rownames(umap_data) <- umap_data$sample
distances <- calculate_distances(umap_data=umap_data, metadata = sample.table)
print(distances$within_group)
print(distances$inter_group)
#write.csv(PCA_distance,"./Figure 3/Expression/PCA/PCA_distance.csv")

# plot bar
within_group_df <- distances$within_group
inter_group_df <- distances$inter_group
within_group_df$group <- as.character(within_group_df$group)
inter_group_df$group_pair <- as.character(inter_group_df$group_pair)
combined_df <- inter_group_df %>%
  separate(group_pair, into = c("group1", "group2"), sep = "_vs_") %>%
  left_join(within_group_df, by = c("group1" = "group")) %>%
  left_join(within_group_df, by = c("group2" = "group"), suffix = c("_group1", "_group2")) %>%
  mutate(inter_intra_fraction = avg_inter_distance / ((avg_within_distance_group1 + avg_within_distance_group2) / 2))
# ggplot(combined_df, aes(x = paste(group1, group2, sep = "_vs_"), y = inter_intra_fraction)) +
#   geom_bar(stat = "identity", fill = "steelblue") +
#   geom_hline(yintercept = 1,size=1, linetype = "dashed", color = "grey30") +
#   theme_minimal(base_size = 14) +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     axis.title.x = element_text(margin = margin(t = 10)),
#     axis.title.y = element_text(margin = margin(r = 10)),
#     plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
#   ) +
#   labs(
#     x = "Group Pair",
#     y = "Inter vs Intra Group Distance Fraction",
#     title = "Group Separation based on UMAP Distances"
#   ) +
#   ylim(0, max(combined_df$inter_intra_fraction) * 1.1)

# Perform permutation test
set.seed(123)
num_permutations <- 100

# Initialize data frame to store permuted fractions
permuted_fractions <- data.frame(matrix(ncol = nrow(combined_df), nrow = num_permutations))
colnames(permuted_fractions) <- combined_df$group_pair

for (i in 1:num_permutations) {
  message(i)
  permuted_metadata <- sample.table
  permuted_metadata$group <- sample(permuted_metadata$group)
  permuted_distances <- calculate_distances(umap_data, permuted_metadata)
  # Ensure the order of groups is consistent
  permuted_distances$within_group <- permuted_distances$within_group %>% arrange(group)
  permuted_distances$inter_group <- arrange(permuted_distances$inter_group, group_pair)

  permuted_combined_df <- permuted_distances$inter_group %>%
    separate(group_pair, into = c("group1", "group2"), sep = "_vs_") %>%
    left_join(permuted_distances$within_group, by = c("group1" = "group")) %>%
    left_join(permuted_distances$within_group, by = c("group2" = "group"), suffix = c("_group1", "_group2")) %>%
    mutate(inter_intra_fraction = avg_inter_distance / ((avg_within_distance_group1 + avg_within_distance_group2) / 2))

  # Ensure permutation results are in the same format as combined_df
  permuted_combined_df <- permuted_combined_df %>%
    mutate(group_pair = paste(group1, group2, sep = "vs")) %>%
    dplyr::select(group_pair, inter_intra_fraction) %>%
    mutate(perm_id = i)
  
  combined_df <- combined_df %>% mutate(group_pair = paste(group1, group2, sep = "vs"))
  # Initialize permuted_fractions
  if (i ==1){
    colnames(permuted_fractions) <- permuted_combined_df$group_pair
  }
  for (pair in combined_df$group_pair) {
    permuted_fractions[i, pair] <- mean(permuted_combined_df$inter_intra_fraction[permuted_combined_df$group_pair == pair], na.rm = TRUE)
  }
}

# Calculate p-values
# Not work
# p_values <- apply(permuted_fractions, 2, function(x) mean(x >= combined_df$inter_intra_fraction, na.rm = TRUE))
p_values <- c()
for (pair in combined_df$group_pair) {
  p_values <- c(p_values, mean(permuted_fractions[, pair] >= combined_df$inter_intra_fraction[combined_df$group_pair == pair], na.rm = TRUE))
}

combined_df$p_value <- p_values

# Add significance symbols
combined_df$significance <- cut(combined_df$p_value,
                                breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                                labels = c("***", "**", "*", ""))

# Plot bar plot with significance and add * on the bar
pdf(paste0("/BioII/lulab_b/wangtaiwei/peak_calling/review/output/",dst,"/allSmp_cfPeakCNN_smallDomain_diff_group_sep_","top",top,".pdf"), width = 8, height = 5)
ggplot(combined_df, aes(x = group_pair, y = inter_intra_fraction, fill = significance)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 1, size = 1, linetype = "dashed", color = "grey30") +
  geom_text(aes(label = significance), vjust = -0.5, size = 8) +
  scale_fill_manual(values = c("***" = "#953B3E", "**" = "#BA6A68", "*" = "#E98677"), labels = c("***" = "*** p < 0.001", "**" = "**  p < 0.01", "*" = "*   p < 0.05")) +
  theme_minimal(base_size = 14) +
  theme(
    # panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    # panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
    axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  ) +
  labs(
    x = "Group Pair",
    y = "Inter vs Intra Group Distance Fraction",
    title = "Group Separation based on UMAP Distances (Plasma)",
    fill = "Significance"
  ) +
  ylim(0, max(combined_df$inter_intra_fraction) * 1.1)
dev.off()




## TCGA Tissue
# Figure 6F
setwd("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/")
options(scipen = 4,stringsAsFactors = F,digits = 5)
dst <- "TCGA_small_diff3"

# sample.table <- read.table("./meta/lulab/FTC/sample_table.txt",sep = "\t",header = T)
sample.table0 <- read.table(paste0("../WCHSU-FTC/exSeek-dev/data/archive/",dst,"/sample_table.txt"),sep = "\t",header = T)
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
dim(count)

rownames(count) <- count$feature
mat <- count[,2:ncol(count)]
#mat[1:3,1:3]
#colnames(mat) # 72

## diff
#sample.table <- sample.table[!(sample.table$patient_id %in% c("FTA-11","FTA-12","FTA-13","FTA-14","FTA-15","FTA-16","FTA-18","FTA-19","FTA-20","FTA-23","FTA-24","FTA-9",     "FTA-13","FTA-17","FTA-4","FTA-6")),]
#sample.table$sample <- unlist(sapply(strsplit(sample.table$sample,"_",fixed=T),"[",1)) # trim samples surfix !
#colnames(mat) <- unlist(sapply(strsplit(colnames(mat),"_",fixed=T),"[",1)) # trim samples surfix !
#write.table(sample.table0, file = "../WCHSU-FTC/exSeek-dev/data/GSE71008/sample_table_filterBatch.txt", sep = "\t", quote = F, row.names = F, col.names = T)
#sample.table0[1:3,1:3]

peak.list <- list()
#nums <- c(50,100,200,500,1000)
nums <- c(500)
for (top in nums){
  peak.list[[paste0("top",top)]] <- list()
}
for (disease in c("Colorectal Cancer","Pancreatic Cancer","Prostate Cancer","Healthy Control")){
  #disease <- "Colorectal Cancer" # "Colorectal Cancer" "Healthy Control"   "Pancreatic Cancer" "Prostate Cancer"
  message(disease)
  if(disease=="Colorectal Cancer"){
    disease.label <- "COAD"
  }else if(disease=="Pancreatic Cancer"){
    disease.label <- "PAAD"
  }else if(disease=="Prostate Cancer"){
    disease.label <- "PRAD"
  }else if(disease=="Healthy Control"){
    disease.label <- "LAML"
  }
  print(disease)
  #normal <- "" #"Healthy Control"
  normal.label <- "R" #"NC" 
  
  tmp <- read.table(paste0("./output/",dst,"/OvR_cfPeakCNN_GSE71008_smallDomain_diff_",disease.label,"vs",normal.label,".diff"),sep = "\t",header = T,check.names = F,stringsAsFactors = F)
  
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


library(edgeR)
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
dim(logcpm)



## run PCA&UMAP (2408)
logcpm0 <- logcpm
rownames(sample.table0) <- sample.table0$sample
table(sample.table0$group)
for (top in 500){
  #top <- 500
  print(top)
  logcpm2 <- logcpm[peak.ids[[paste0("top",top)]],sample.table0$sample]
  print(dim(logcpm2))

  sample.table <- sample.table0
  table(sample.table$group)
  # sample.table$group <- gsub("Healthy Control","NC",sample.table$group)
  # sample.table$group <- gsub("Colorectal Cancer","CRC",sample.table$group)
  # sample.table$group <- gsub("Prostate Cancer","PRAD",sample.table$group)
  # sample.table$group <- gsub("Pancreatic Cancer","PAAD",sample.table$group)
  groupColor <- c("LAML"="grey", "COAD"="salmon", "PRAD"="purple", "PAAD"="steelblue")
  sample.table$Color <- sample.table$group
  sample.table$Fill <- sample.table$group
  sample.table$Shape <- factor(21,levels = 1:21) #sample.table$group
  # plotPCA(sample.table = sample.table, logcpm = logcpm2, outFile = paste0("./output/",dst,"/PCA/allSmp_cfPeakCNN_smallDomain_diff.cpm.PCA.","top",top,".pdf"), 
  #         plotGrpColor = groupColor, topVarNum = top, plotUMAP = T, plotDotSize = 2, plotDotAlpha = 0.9, figWid = 8, figHigh = 5)
  plotPCA(sample.table = sample.table, logcpm = logcpm2, outFile = paste0("/BioII/lulab_b/wangtaiwei/peak_calling/review/output/",dst,"/PCA/allSmp_cfPeakCNN_smallDomain_diff.cpm.PCA.","top",top,".pdf"), 
          plotGrpColor = groupColor, topVarNum = top, plotUMAP = T, plotDotSize = 2, plotDotAlpha = 0.9, figWid = 8, figHigh = 5)
  #
}

# compare within vs. inter- gorup distance from UMAP distance
top <- 500
# umap_data <- read.csv(paste0("./output/",dst,"/PCA/allSmp_cfPeakCNN_smallDomain_diff.cpm.PCA.","top",top,".pdf.UMAP.txt"),sep = "\t")
umap_data <- read.csv(paste0("/BioII/lulab_b/wangtaiwei/peak_calling/review/output/",dst,"/PCA/allSmp_cfPeakCNN_smallDomain_diff.cpm.PCA.","top",top,".pdf.UMAP.txt"),sep = "\t")
rownames(umap_data) <- umap_data$sample
distances <- calculate_distances(umap_data=umap_data, metadata = sample.table)
print(distances$within_group)
print(distances$inter_group)

within_group_df <- distances$within_group
inter_group_df <- distances$inter_group
within_group_df$group <- as.character(within_group_df$group)
inter_group_df$group_pair <- as.character(inter_group_df$group_pair)
combined_df <- inter_group_df %>%
  separate(group_pair, into = c("group1", "group2"), sep = "_vs_") %>%
  left_join(within_group_df, by = c("group1" = "group")) %>%
  left_join(within_group_df, by = c("group2" = "group"), suffix = c("_group1", "_group2")) %>%
  mutate(inter_intra_fraction = avg_inter_distance / ((avg_within_distance_group1 + avg_within_distance_group2) / 2))

# Perform permutation test
set.seed(123)
num_permutations <- 100

# Initialize data frame to store permuted fractions
permuted_fractions <- data.frame(matrix(ncol = nrow(combined_df), nrow = num_permutations))
colnames(permuted_fractions) <- combined_df$group_pair

for (i in 1:num_permutations) {
  message(i)
  permuted_metadata <- sample.table
  permuted_metadata$group <- sample(permuted_metadata$group)
  permuted_distances <- calculate_distances(umap_data, permuted_metadata)
  # Ensure the order of groups is consistent
  permuted_distances$within_group <- permuted_distances$within_group %>% arrange(group)
  permuted_distances$inter_group <- arrange(permuted_distances$inter_group, group_pair)

  permuted_combined_df <- permuted_distances$inter_group %>%
    separate(group_pair, into = c("group1", "group2"), sep = "_vs_") %>%
    left_join(permuted_distances$within_group, by = c("group1" = "group")) %>%
    left_join(permuted_distances$within_group, by = c("group2" = "group"), suffix = c("_group1", "_group2")) %>%
    mutate(inter_intra_fraction = avg_inter_distance / ((avg_within_distance_group1 + avg_within_distance_group2) / 2))

  # Ensure permutation results are in the same format as combined_df
  permuted_combined_df <- permuted_combined_df %>%
    mutate(group_pair = paste(group1, group2, sep = "vs")) %>%
    dplyr::select(group_pair, inter_intra_fraction) %>%
    mutate(perm_id = i)
  
  combined_df <- combined_df %>% mutate(group_pair = paste(group1, group2, sep = "vs"))
  # Initialize permuted_fractions
  if (i ==1){
    colnames(permuted_fractions) <- permuted_combined_df$group_pair
  }
  for (pair in combined_df$group_pair) {
    permuted_fractions[i, pair] <- mean(permuted_combined_df$inter_intra_fraction[permuted_combined_df$group_pair == pair], na.rm = TRUE)
  }
}

# Calculate p-values
# Not work
# p_values <- apply(permuted_fractions, 2, function(x) mean(x >= combined_df$inter_intra_fraction, na.rm = TRUE))
p_values <- c()
for (pair in combined_df$group_pair) {
  p_values <- c(p_values, mean(permuted_fractions[, pair] >= combined_df$inter_intra_fraction[combined_df$group_pair == pair], na.rm = TRUE))
}

combined_df$p_value <- p_values

# Add significance symbols
combined_df$significance <- cut(combined_df$p_value,
                                breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                                labels = c("***", "**", "*", ""))

# Plot bar plot with significance and add * on the bar
pdf(paste0("/BioII/lulab_b/wangtaiwei/peak_calling/review/output/",dst,"/allSmp_cfPeakCNN_smallDomain_diff_group_sep_","top",top,".pdf"), width = 4, height = 2.5)
ggplot(combined_df, aes(x = group_pair, y = inter_intra_fraction, fill = significance)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 1, size = 1, linetype = "dashed", color = "grey30") +
  geom_text(aes(label = significance), vjust = -0.5, size = 8) +
  scale_fill_manual(values = c("***" = "#953B3E", "**" = "#BA6A68", "*" = "#E98677"), labels = c("***" = "*** p < 0.001", "**" = "**  p < 0.01", "*" = "*   p < 0.05")) +
  theme_minimal(base_size = 14) +
  theme(
    # panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    # panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
    axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  ) +
  labs(
    x = "Group Pair",
    y = "Inter vs Intra Group Distance Fraction",
    title = "Group Separation based on UMAP Distances (Cancer tissue)",
    fill = "Significance"
  ) +
  ylim(0, max(combined_df$inter_intra_fraction) * 1.1)
dev.off()