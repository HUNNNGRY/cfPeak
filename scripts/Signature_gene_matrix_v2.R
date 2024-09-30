#! /usr/bin/env Rscript
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(biomaRt))
# UMAP in R (source: wechat 降维你还在用多种R包吗？快来试试这个宝藏R包吧！)
#remotes::install_github("jlmelville/uwot")
#install.packages("tidydr")
suppressPackageStartupMessages(library(tidydr))#降维
suppressPackageStartupMessages(library(ggplot2))#可视化
#install.packages("mlr3")
#library(mlr3)

parser <- ArgumentParser(description='Signature Matrix v2')
parser$add_argument('-m', '--matrix', type='character', required=TRUE,
                    help='Input count matrix. Rows are genes(TPM). Columns are samples.')
parser$add_argument('-s', '--sample_info', type='character', required=TRUE,
                    help='Sample information for each sample. CSV file with 2 columns. Column 1: sample, Column 2: group')
parser$add_argument('-e', '--TPM_cutoff', type='double', required=FALSE,default=1,
                    help='TPM cutoff, maximum expressionn larger than this cutoff. Default: 1.')
parser$add_argument('-ts','--TSS_cutoff', type='double', required=FALSE,default=0.5,
                    help='TSS cutoff, minimum tissue specific score larger than this cutoff. Default: 0.5.')
parser$add_argument('-n','--Number', type='integer', required=FALSE,default=20,
                    help='Signature gene number. Default: 20.')
parser$add_argument('-o', '--outdir', type='character', required=TRUE,
                    help='output directory')
args <- parser$parse_args()

#cl <- makeCluster(1) #not to overload your computer
#registerDoParallel(cl)

message('read count matrix: ', args$matrix)
message('output: ', args$outdir)
message('Start read in files at: ',Sys.time())
dir.create(args$outdir,recursive=T,showWarnings=F)

# def function
mean_TSS <- function(mat,sample_info,TPM_cutoff=1,TSS_cutoff=0.5,Num=20,outdir){
  sample_info_ex <- as.data.frame(list(c("a","b","c","d","e","f","g","h","i"),c("A","A","A","B","B","B","C","C","C")))
  colnames(sample_info_ex) <- c("sample","group")
  mat_ex <- as.data.frame(list(c(100,10,1),c(110,8,2),c(120,6,3),c(10,110,1),c(8,120,2),c(6,100,3),c(1,10,110),c(2,8,120),c(3,6,100)))
  colnames(mat_ex) <- c("a","b","c","d","e","f","g","h","i")
  rownames(mat_ex) <- c("gene1","gene2","gene3")
  
  message("Example mat:")
  print(mat_ex)
  message("Example sample_info:")
  print(sample_info_ex)
  message("Example TPM_cutoff: 1")
  message("Example TSS_cutoff: 0.5")
  message("Example Number: 20")
  
  group <- unique(sort(sample_info$group))
  i=1
  while(i<=length(group)){
    sample <- sample_info[which(sample_info$group==group[i]),]$sample
    median <- as.data.frame(rowMeans(as.matrix(mat[,as.character(sample)])))
    colnames(median) <- group[i]
    rownames(median) <- rownames(mat)
    if(i==1){
      median_matrix <- median
    } else {
      median_matrix <- cbind(median_matrix,median)
    }
    i=i+1
  }
  
  median_matrix <- median_matrix[which(rowMaxs(as.matrix(median_matrix)) > TPM_cutoff),]
  
  j=1
  sum <- as.data.frame(rowSums(median_matrix))
  while(j<=ncol(median_matrix)){
    TSS <- median_matrix[,j]/sum
    colnames(TSS) <- colnames(median_matrix)[j]
    if(j==1){
      TSS_matrix <- TSS
    } else {
      TSS_matrix <- cbind(TSS_matrix,TSS)
    }
    j=j+1
  }

  n=1
  TSG_plot={}
  while(n<=ncol(TSS_matrix)){
    test_sorted <- TSS_matrix[order(TSS_matrix[,n],decreasing = TRUE),]
    test_cutoff <- test_sorted[test_sorted[,n]>TSS_cutoff,]
    test_cutoff$gene <- rownames(test_cutoff)
    TSG_plot <- rbind(TSG_plot,test_cutoff[1:min(Num,sum(test_sorted[,n]>TSS_cutoff)),])
    #     TSG_plot <- rbind(TSG_plot,test_cutoff[1:max(2,min(Num,sum(test_sorted[,n]>TSS_cutoff))),])
    n=n+1
  }  

  TSG_plot <- TSG_plot[!duplicated(TSG_plot$gene),]
  rownames(TSG_plot) <- TSG_plot$gene
  TSG_plot$gene <- NULL
  # print(rowMaxs(as.matrix(TSG_plot)))
  TSG_plot <- TSG_plot[which(rowMaxs(as.matrix(TSG_plot)) > TSS_cutoff),]
  
  k=1
  TSG_plot$group <- NA 
  while(k<=length(group)){
    if(length(grep("TRUE",TSG_plot[,as.character(group[k])] > TSS_cutoff))>0){
      TSG_plot[(TSG_plot[,as.character(group[k])] > TSS_cutoff),]$group <- as.character(group[k])
      k=k+1
    } else {
      k=k+1
    }
  }
  
  #TSG_plot <- TSG_plot[order(TSG_plot$group),]
  signature_matrix <- median_matrix[rownames(TSG_plot),]
  signature_matrix <- signature_matrix[order(TSG_plot$group),]
  rownames(signature_matrix) <- as.character(lapply(strsplit(rownames(signature_matrix),".",fixed = TRUE),function(x) x[1]))
  write.table(signature_matrix,paste0(outdir,"/","signature_genes_",TPM_cutoff,"_",TSS_cutoff,"_",Num,".txt"),quote = FALSE,sep = "\t",row.names = TRUE,col.names = TRUE)
  row_annotation <- data.frame(rownames(TSG_plot),TSG_plot[,which(colnames(TSG_plot)=="group")])
  colnames(row_annotation) <- c("gene","group")
  rownames(row_annotation) <- row_annotation$gene
  row_annotation$gene <- NULL
  write.table(row_annotation,paste0(outdir,"/","gene_list_",TPM_cutoff,"_",TSS_cutoff,"_",Num,".txt"),quote=FALSE,sep="\t",row.names = TRUE, col.names = TRUE)
  return(TSG_plot)
}

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}



# get params
outdir <- args$outdir
count_matrix <- read.table(args$matrix,sep = "\t",header = TRUE, row.names = 1, check.names = FALSE)
sample_info <- read.csv(args$sample_info,header = TRUE)
count_matrix <- count_matrix[,colnames(count_matrix) %in% sample_info$sample]
if(min(count_matrix[,1])<0){
  message("minus value exist in mat, treat mat as logtpm and convert to tpm")
  count_matrix <- 2^count_matrix
}else{
  message("treat mat as tpm")
}

TPM_cutoff <- args$TPM_cutoff
TSS_cutoff <- args$TSS_cutoff
Num <- args$Number



# ## test SLE
# outdir <- "/data2/lulab1/bpf/projects/WCHSU-FTC/output/SLE/call_peak_dedup/TOO"
# count_matrix <- read.table("/data2/lulab1/bpf/projects/WCHSU-FTC/output/SLE/call_peak_dedup/count_matrix/cfpeakCNN_b5_d50_p1_log2TPM_rmMonthOperator.txt",sep = "\t",header = TRUE, row.names = 1, check.names = FALSE)
# sample_info <- read.csv("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/SLE/240316/TOO_HC.csv",header = TRUE)
# count_matrix <- count_matrix[,colnames(count_matrix) %in% sample_info$sample]
# count_matrix <- 2^count_matrix
# 
# TPM_cutoff <- 0.1
# TSS_cutoff <- 0.5
# Num <- 200


## test old
# setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/")
# outdir <- "../output/GSE163534/call_domain_withRepeats_all/too"
# count_matrix <- read.table("../output/GSE163534/call_domain_withRepeats_all/count_matrix/domains_long_localmax2_b5_d05_p05_CPM-TMM_filter.txt",sep = "\t",header = TRUE, row.names = 1, check.names = FALSE)
# sample_info <- read.csv("./data/GSE163534/sample_table_filter.csv",header = TRUE)
# all_sample <- as.data.frame(colnames(count_matrix))
# colnames(all_sample) <- "sample"
# all_sample$exist <- "Y"
# sample_info <- dplyr::left_join(sample_info,all_sample, by=c("sample"="sample"))
# sample_info <- sample_info[which(sample_info$exist=="Y"),]
# sample_info$exist <- NULL
# mat <- count_matrix[,as.character(sample_info$sample)]
# TPM_cutoff <- 1
# TSS_cutoff <- 0.5
# Num <- 50


all_sample <- as.data.frame(colnames(count_matrix))
colnames(all_sample) <- "sample"
all_sample$exist <- "Y"
sample_info <- sample_info[order(sample_info$group),]
sample_info <- left_join(sample_info,all_sample, by=c("sample"="sample"))
sample_info <- sample_info[which(sample_info$exist=="Y"),]
sample_info$exist <- NULL

mat <- count_matrix[,as.character(sample_info$sample)]
#table(duplicated(rownames(mat)))

message('TPM cutoff: ', TPM_cutoff)
message('TSS cutoff: ', TSS_cutoff)
message('Signature gene number: ', Num)
message('Start calculating tissue specific genes at: ',Sys.time())


# run
TSG_plot <- mean_TSS(mat,sample_info,TPM_cutoff,TSS_cutoff,Num,outdir)

message('End calculating tissue specific genes at: ',Sys.time())
message('Start ploting signature genes at: ',Sys.time())

col_annotation <- sample_info
rownames(col_annotation) <- col_annotation$sample
col_annotation <- col_annotation[order(sample_info$group),]
col_annotation$sample <- NULL
#col_annotation <- col_annotation[,c("group")]
row_annotation <- data.frame(row.names = rownames(TSG_plot), 
                             species=unlist(sapply(strsplit(rownames(TSG_plot),"|",fixed = T),"[",2)) )
row_annotation$species=unlist(sapply(strsplit(rownames(TSG_plot),"|",fixed = T),"[",2)) 
row_annotation$species <- gsub("_for|_rev|\\.for|\\.rev","",row_annotation$species,perl = T)
rna <- c("rRNA","pri_miRNA","lncRNA", "mRNA","piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA')
dna <- c("intron", "promoter", "enhancer","repeats") 
row_annotation$species <- factor(row_annotation$species,levels = c(rna,dna))
#row_annotation <- row_annotation$[,c("species","group")]

mat_plot <- mat[rownames(TSG_plot),order(sample_info$group)]
bk = unique(c(seq(0,1, length=100)))
plot_sample <- pheatmap(
  mat_plot[order(TSG_plot$group),],
  breaks = bk,
  annotation_col = col_annotation,
  annotation_row = row_annotation,
  scale = "row",
  cluster_cols = FALSE,cluster_rows = FALSE,
  show_colnames=FALSE, show_rownames=FALSE,
  #colorRampPalette(c("blue","white","red"))(100),
  # color = (colorRampPalette(RColorBrewer::brewer.pal(7, "Reds"))(100)), #RdBu, colorRampPalette(c("steelblue","white","salmon"))(100),
  color = viridis::viridis_pal()(100),
  # annotation_colors = (colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(unique(col_annotation$sample)))),
  fontsize_row = 5,
  clustering_distance_cols = "euclidean")
save_pheatmap_pdf(plot_sample, paste0(outdir,"/","signature_genes_samples_",TPM_cutoff,"_",TSS_cutoff,"_",Num,".pdf"))
dev.off()

#plot signature gene matrix
#mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
#id <- as.character(lapply(strsplit(rownames(TSG_plot),"\\."),function(x) x[1]))
#transversion <- data.frame("raw"=rownames(TSG_plot),"id"=id)
#gene_names <- getBM(attributes=c("ensembl_gene_id", "external_gene_name","gene_biotype"),
#                    filters = "ensembl_gene_id",
#                    values=transversion$id, mart= mart,useCache = FALSE)
#annotated <- left_join(transversion,gene_names, by= c("id"="ensembl_gene_id"))
#rownames(TSG_plot) <- paste(annotated$gene_biotype,annotated$raw,annotated$external_gene_name,sep = "|")

row_annotation <- data.frame(rownames(TSG_plot),TSG_plot[,which(colnames(TSG_plot)=="group")])
colnames(row_annotation) <- c("gene","group")
rownames(row_annotation) <- row_annotation$gene
row_annotation$species=unlist(sapply(strsplit(rownames(TSG_plot),"|",fixed = T),"[",2)) 
row_annotation$species <- gsub("_for|_rev|\\.for|\\.rev","",row_annotation$species,perl = T)
rna <- c("rRNA","pri_miRNA","lncRNA", "mRNA","piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA')
dna <- c("intron", "promoter", "enhancer","repeats") 
row_annotation$species <- factor(row_annotation$species,levels = c(rna,dna))
row_annotation$gene <- NULL
#row_annotation <- row_annotation[order(TSG_plot$group),]
#head(row_annotation[order(TSG_plot$group),])

bk = unique(c(seq(0,1, length=100)))
plot_tissue <- pheatmap(
  TSG_plot[order(TSG_plot$group),-which(colnames(TSG_plot)=="group")],
  breaks = bk,
  #annotation_col = col_annotation,
  annotation_row = row_annotation,
  scale = "row",
  cluster_cols = FALSE,cluster_rows = FALSE,
  show_colnames=TRUE, show_rownames=FALSE,  
  # color = (colorRampPalette(RColorBrewer::brewer.pal(7, "Reds"))(100)), #colorRampPalette(c("steelblue","white","salmon"))(100),
  color = viridis::viridis_pal()(100),
  # annotation_colors = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(unique(row_annotation$gene))),
  fontsize_row = 5,
  clustering_distance_cols = "euclidean")

save_pheatmap_pdf(plot_tissue, paste0(outdir,"/","signature_genes_tissues_",TPM_cutoff,"_",TSS_cutoff,"_",Num,".pdf"))
dev.off()



# plot sig_gene/row_anno seperately
#row_annotation[1:3,]
ggplot(row_annotation,aes(x=group,fill=species))+
  geom_bar(stat="count",position="stack")+ #,position="fill"
  geom_vline(xintercept = 0,linetype="dashed",color="grey10")+
  #ggsci::scale_fill_d3()+
  scale_fill_manual(name="Species",values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(unique(row_annotation$species))))+ # c(RColorBrewer::brewer.pal(8, "Accent"),RColorBrewer::brewer.pal(9, "Set1"))
  # paletteer::scale_fill_paletteer_d("RColorBrewer::Dark2")+
  #scale_y_continuous(trans = "log10")+
  #xlim(c(0,500))+
  #geom_hline(yintercept = c(0))+
  #ggraph::scale_fill_viridis(option = "C") + # name = "value",
  theme_bw() + 
  theme(
    plot.title = element_text(size = 24,color="black",hjust = 0.5),
    axis.title = element_text(size = 24,color ="black"), 
    axis.text = element_text(size= 24,color = "black"),
    axis.text.x = element_text(size= 24,color = "black",angle = 90,vjust = 0.5,hjust = 1),
    #panel.grid=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey50",linetype = "dashed"), #size= 1,
    #panel.grid.minor.y = element_blank(),
    panel.border = element_blank(),
    legend.position = "right",#c(.25,.6),
    legend.text = element_text(size= 24),
    legend.title= element_text(size= 24))
ggsave(filename = paste0(outdir,"/","signature_stat_bar_",TPM_cutoff,"_",TSS_cutoff,"_",Num,".pdf") )



# plot tsne
#tidydr::available_methods()
#x <- tidydr::dr(data = t(mat), fun = Rtsne::Rtsne,perplexity = 10)# uwot::tumap, uwot::umap, uwot::lvish, Rtsne::Rtsne, stats::prcomp
for(i in seq(5,50,10)){
tryCatch(
expr={
x <- tidydr::dr(data = t(mat), fun = Rtsne::Rtsne,perplexity = i,max_iter=1000,dims=2)# uwot::tumap, uwot::umap, uwot::lvish, Rtsne::Rtsne, stats::prcomp
autoplot(x, aes(fill=group),alpha=0.9, metadata = sample_info[, "group", drop=FALSE] ) +
  geom_point(size=5,alpha=0.99,shape=21,color="white")+
  theme_dr()+
  scale_fill_manual(name="Group",values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(unique(sample_info$group))))+ 
  theme(
    plot.title = element_text(size = 24,color="black",hjust = 0.5),
    axis.title = element_text(size = 24,color ="black"), 
    axis.text = element_blank(),
    panel.grid=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(), #element_line(color = "grey50",linetype = "dashed"), #size= 1,
    panel.grid.minor.y = element_blank(),
    panel.border = element_blank(),
    legend.position = "right",#c(.25,.6),
    legend.text = element_text(size= 24),
    legend.title= element_text(size= 24),
    strip.text.y = element_blank(),
    strip.text.x = element_text(size=24)
  )
ggsave(filename = paste0(outdir,"/","tSNE_all_",TPM_cutoff,"_",TSS_cutoff,"_",Num,"_perplexity",i,".pdf"), width =  18, height = 12)
},
error = function(e){
message("seem err in tSNE")
}
)
}

sig <- read.table(paste0(outdir,"/","gene_list_",TPM_cutoff,"_",TSS_cutoff,"_",Num,".txt"),sep="\t",row.names = 1, header = T,stringsAsFactors = F,check.names = F)
mat <- mat[rownames(mat) %in% rownames(sig),]
for(i in seq(5,50,10)){
tryCatch(
expr={x <- tidydr::dr(data = t(mat), fun = Rtsne::Rtsne,perplexity = i,max_iter=1000,dims=2)# uwot::tumap, uwot::umap, uwot::lvish, Rtsne::Rtsne, stats::prcomp
autoplot(x, aes(fill=group), metadata = sample_info[, "group", drop=FALSE] ) +
  geom_point(size=5,alpha=0.99,shape=21,color="white")+
  theme_dr()+
  scale_fill_manual(name="Group",values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(unique(sample_info$group))))+ 
  theme(
    plot.title = element_text(size = 24,color="black",hjust = 0.5),
    axis.title = element_text(size = 24,color ="black"), 
    axis.text = element_blank(),
    panel.grid=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(), #element_line(color = "grey50",linetype = "dashed"), #size= 1,
    panel.grid.minor.y = element_blank(),
    panel.border = element_blank(),
    legend.position = "right",#c(.25,.6),
    legend.text = element_text(size= 24),
    legend.title= element_text(size= 24),
    strip.text.y = element_blank(),
    strip.text.x = element_text(size=24)
  )
ggsave(filename = paste0(outdir,"/","tSNE_sig_",TPM_cutoff,"_",TSS_cutoff,"_",Num,"_perplexity",i,".pdf"), width =  18, height = 12)
},
error = function(e){
message("seem err in tSNE")
}
)
}
