#!/usr/bin/env Rscript

#2017, BIB, A benchmark of gene expression tissue-specificity metrics

options( stringsAsFactors = F)
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(tidydr))#降维
suppressPackageStartupMessages(library(ggplot2))#可视化
suppressPackageStartupMessages(library(ggsci))#可视化

parser <- ArgumentParser(description=paste(
  "Signature Matrix v3",
  "benchmark: 2017, BIB, A benchmark of gene expression tissue-specificity metrics",
  "This script demonstrates the creation of a sample information dataframe and a matrix dataframe.",
  "",
  "Example:",
  "Sample Info:",
  "sample group",
  "0      a     A",
  "1      b     A",
  "2      c     A",
  "3      d     B",
  "4      e     B",
  "5      f     B",
  "6      g     C",
  "7      h     C",
  "8      i     C",
  "",
  "Matrix:",
  "      a    b    c    d    e    f    g    h    i",
  "gene1  100  110  120   10    8    6    1    2    3",
  "gene2   10    8    6  110  120  100   10    8    6",
  "gene3    1    2    3    1    2    3  110  120  100",
  collapse="\n"))
parser$add_argument('-m', '--matrix', type='character', required=TRUE,
                    help='Input count matrix. Rows are genes (TPM). Columns are samples.')
parser$add_argument('-s', '--sample_info', type='character', required=TRUE,
                    help='Sample information for each sample. CSV file with 2 columns. Column 1: sample, Column 2: group.')
parser$add_argument('-e', '--TPM_cutoff', type='double', required=FALSE, default=1,
                    help='TPM cutoff, maximum expression larger than this cutoff. Default: 1.')
parser$add_argument('-ts', '--TSS_cutoff', type='double', required=FALSE, default=0.5,
                    help='TSS cutoff, minimum tissue-specific score larger than this cutoff. Default: 0.5.')
parser$add_argument('-n', '--Number', type='integer', required=FALSE, default=20,
                    help='Signature gene number. Default: 20.')
parser$add_argument('-o', '--outdir', type='character', required=TRUE,
                    help='Output directory.')
parser$add_argument('-method', '--method', type='character', required=FALSE, default='mean',
                    help="Method to select signature genes: 'mean' for mean expression norm by rowSums (TSI). TODO: add other options, like zscore. Default: 'mean'.")
parser$add_argument('--color_list', type='character', required=FALSE, default='NULL',
                    help="Rdata path that include color_list, which contains 'group' and 'species'. Default: 'NULL'. e.g.: /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/TCGA_small16/color.rds")
parser$add_argument('--rnaColumn', type='integer', required=FALSE, default=8,
                    help='| seperated Column index for RNA species label for Signature gene. Default: 8. Opt: 2')

# Parse command line arguments
args <- parser$parse_args()
for(i in 1:length(args)){
  v.name <- names(args)[i]
  assign(v.name,args[[i]])
  message(paste0(v.name,": ",get(v.name)))
}

# auto color num ggsci
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

## cal each gene single value
# calculate_tau <- function(expression_matrix) {
#   if (is.vector(expression_matrix)) {
#     expression_matrix <- matrix(expression_matrix, nrow = 1)
#   } else if (!is.matrix(expression_matrix)) {
#     expression_matrix <- as.matrix(expression_matrix)
#   }
#   max_expression <- rowMaxs(expression_matrix)
#   relative_expression <- expression_matrix / max_expression
#   tau_values <- apply(relative_expression, 1, function(x) sum(1 - x) / (length(x) - 1))
#   return(tau_values)
# }

#calculate TSI 
mean_TSS <- function(mat, sample_info, TPM_cutoff=1, TSS_cutoff=0.5, Num=20, outdir, method='mean') {
  group <- unique(sort(sample_info$group))
  median_matrix <- sapply(group, function(g) rowMeans(mat[, sample_info$sample[sample_info$group == g], drop=FALSE]))
  median_matrix <- median_matrix[which(rowMaxs(median_matrix) > TPM_cutoff),]
  #dim(median_matrix)
  
  # if (method == 'tau') {
  #   tau_mat <- data.frame(row.names = rownames(median_matrix), score=calculate_tau(median_matrix) )
  #   # lapply(group, function(g) {
  #   #   tau_mat[[g]] <- s
  #   # }
  #   # #hist(TSS_matrix)
  #   TSS_matrix <- TSS_matrix
  #   colnames(TSS_matrix) <- group
  # } else {
    sum_matrix <- rowSums(median_matrix)
    TSS_matrix <- sweep(median_matrix, 1, sum_matrix, '/')
  # }
    #table(rowSums((TSS_matrix))) 
   # dim(TSS_matrix)

  #head(sig.tab)
  sig.tab <- data.frame( gene=rownames(TSS_matrix), group=colnames(TSS_matrix)[apply(TSS_matrix, 1, which.max)] ) # each gene assgin one group with max expr
  data.table::fwrite(sig.tab, paste0(outdir, "/all_gene_list_", TPM_cutoff, "_", TSS_cutoff, "_", Num,"_",method, ".txt"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
  
  TSG_plot <- (do.call(rbind, lapply(group, function(g) {
    # g <- "Resting.Memory.CD8.T"
    # message(g)
    selected_genes <- as.data.frame(TSS_matrix[order(TSS_matrix[, g], decreasing = TRUE), ])
    keepNum <- sum(selected_genes[, g] > TSS_cutoff)
    if(keepNum>=2){
      # keepNum <- 5
      keepNum <- min(Num,keepNum)
      message(g,": use ",keepNum)
      selected_genes <- as.data.frame(selected_genes[selected_genes[, g] > TSS_cutoff, ])
      # dim(selected_genes)
      selected_genes <- selected_genes[1:min(Num,keepNum), , drop=FALSE]
    }else if(keepNum<2){
      pseudoNum <- 2
      message(g,": not enough gene left, use ",pseudoNum)
      selected_genes <- selected_genes[1:min(Num,pseudoNum), , drop=FALSE]
    }

    selected_genes$gene <- rownames(selected_genes)
    # print(dim(selected_genes))
    return(selected_genes)
  })))
  colnames(TSG_plot)[1:(ncol(TSG_plot)-1)] <- as.character(group)
  TSG_plot <- TSG_plot[!duplicated(TSG_plot$gene),]
  rownames(TSG_plot) <- TSG_plot$gene
  TSG_plot$gene <- NULL

  # TSG_plot <- (TSG_plot[which(rowMaxs(as.matrix(TSG_plot)) > TSS_cutoff), ])
  #TSG_plot$group <- apply(TSG_plot, 1, function(row) group[which.max(row)])
  TSG_plot$group <- sig.tab$group[match(rownames(TSG_plot),sig.tab$gene)]
  #table(TSG_plot$group)
  tmp <- as.data.frame(median_matrix)[rownames(TSG_plot), ]
  tmp <- cbind(gene=rownames(tmp),tmp)
  data.table::fwrite(tmp, paste0(outdir, "/signature_genes_", TPM_cutoff, "_", TSS_cutoff, "_", Num,"_",method, ".txt"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
  tmp <- data.frame(gene=rownames(TSG_plot), group=TSG_plot$group)
  data.table::fwrite(tmp, paste0(outdir, "/gene_list_", TPM_cutoff, "_", TSS_cutoff, "_", Num,"_",method, ".txt"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
  return(TSG_plot)
}

save_pheatmap_pdf <- function(x, filename, device="pdf", width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  if(device=="pdf"){
    pdf(filename, width=width, height=height)
  }else if(device=="png"){
    png(filename, width=width, height=height, units="in", res=800)
  }
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}



## read input
#count_matrix <- read.table(args$matrix, sep="\t", header=TRUE, row.names=1, check.names=FALSE)
count_matrix <- data.table::fread(args$matrix,sep = "\t",header = TRUE, stringsAsFactors = F, check.names = FALSE, data.table = F)
rownames(count_matrix) <- count_matrix[,1]
count_matrix[,1] <- NULL
sample_info <- read.csv(args$sample_info, header=TRUE, stringsAsFactors = F)

TPM_cutoff <- args$TPM_cutoff
TSS_cutoff <- args$TSS_cutoff
Num <- args$Number
method <- args$method
outdir <- args$outdir

# ## test SLE
# outdir <- "/data2/lulab1/bpf/projects/WCHSU-FTC/output/SLE/call_peak_dedup/TOO_CPM"
# dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
# #count_matrix <- read.table("/data2/lulab1/bpf/projects/WCHSU-FTC/output/SLE/call_peak_dedup/count_matrix/cfpeakCNN_b5_d50_p1_log2TPM_rmMonthOperator.txt",sep = "\t",header = TRUE, row.names = 1, check.names = FALSE)
# #count_matrix <- data.table::fread("/data2/lulab1/bpf/projects/WCHSU-FTC/output/SLE/call_peak_dedup/count_matrix/cfpeakCNN_b5_d50_p1_CPM_rmMonthOperator.txt",sep = "\t",header = TRUE, stringsAsFactors = F, check.names = FALSE, data.table = F)
# count_matrix <- data.table::fread("/data2/lulab1/bpf/projects/WCHSU-FTC/output/SLE/call_peak_dedup/count_matrix/CPM/cfpeakCNN_b5_d50_p1_CPM.txt",sep = "\t",header = TRUE, stringsAsFactors = F, check.names = FALSE, data.table = F)
# rownames(count_matrix) <- count_matrix[,1]
# count_matrix[,1] <- NULL
# #sum(count_matrix$`20230606-SLE-HD10_BC12`)
# #summary(count_matrix$`20230606-SLE-HD10_BC12`)
# sample_info <- read.csv("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/SLE/240316/TOO_HC.csv",header = TRUE)
# count_matrix <- count_matrix[,colnames(count_matrix) %in% sample_info$sample]
# #count_matrix <- 2^count_matrix
# TPM_cutoff <- 0.1
# TSS_cutoff <- 0.1
# Num <- 200
# method <- "mean" #mean

# ## test TCGA
# outdir <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/TCGA_small7/call_peak_all/TOO_all_tmm"
# #count_matrix <- read.table("/data2/lulab1/bpf/projects/WCHSU-FTC/output/SLE/call_peak_dedup/count_matrix/cfpeakCNN_b5_d50_p1_log2TPM_rmMonthOperator.txt",sep = "\t",header = TRUE, row.names = 1, check.names = FALSE)
# count_matrix <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/TCGA_small7/call_peak_all/count_matrix/cfpeakCNN_b5_d50_p1_CPMtmm.txt",sep = "\t",header = TRUE, stringsAsFactors = F, check.names = FALSE, data.table = F)
# rownames(count_matrix) <- count_matrix[,1]
# count_matrix[,1] <- NULL
# sample_info <- read.csv("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/TCGA_small7/TOO_All.csv",header = TRUE)
# count_matrix <- count_matrix[,colnames(count_matrix) %in% sample_info$sample]
# TPM_cutoff <- 0.1
# TSS_cutoff <- 0.1
# Num <- 200
# method <- "mean" #mean
# color_list <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/TCGA_small16/color.rds"

dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

count_matrix <- count_matrix[, colnames(count_matrix) %in% sample_info$sample]
if (min(count_matrix) < 0) {
  message("Negative values found in matrix, treating matrix as logTPM and converting to TPM.")
  count_matrix <- 2^count_matrix
} else {
  message("Treating matrix as TPM.")
}

sample_info <- sample_info[sample_info$sample %in% colnames(count_matrix),]
sample_info <- sample_info[order(sample_info$group), ]
mat <- count_matrix[, as.character(sample_info$sample)]

message('TPM cutoff: ', TPM_cutoff)
message('TSS cutoff: ', TSS_cutoff)
message('Signature gene number: ', Num)
message('Method: ', method)
message('Start calculating tissue-specific genes at: ', Sys.time())


# run
TSG_plot <- mean_TSS(mat,sample_info,TPM_cutoff,TSS_cutoff,Num,outdir)
#TSG_plot <- TSG_plot[order(TSG_plot$group),] # not sort this


message('End calculating tissue specific genes at: ',Sys.time())
message('Start ploting signature genes at: ',Sys.time())

rna <- c("rRNA","pri_miRNA","lncRNA", "mRNA","piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA')
rna2 <- c("pri_miRNA","piRNA","rRNA","lncRNA","mRNA","snoRNA","snRNA","srpRNA","tRNA","tucpRNA","Y_RNA","ktRNA","Rfam")
dna <- c("intron", "promoter", "enhancer","repeats") 

col_annotation <- sample_info
rownames(col_annotation) <- col_annotation$sample
col_annotation$sample <- NULL

row_annotation <- data.frame(row.names = rownames(TSG_plot), species=unlist(sapply(strsplit(rownames(TSG_plot),"|",fixed = T),"[",2)) )
if(rnaColumn==8){
row_annotation$species=unlist(sapply(strsplit(rownames(TSG_plot),"|",fixed = T),"[",8))  # 2, 8
row_annotation$species <- gsub("_for|_rev|\\.for|\\.rev","",row_annotation$species,perl = T)
row_annotation$species <- factor(row_annotation$species,levels = c(rna2,dna)) # rna
}else if(rnaColumn==2){
row_annotation$species=unlist(sapply(strsplit(rownames(TSG_plot),"|",fixed = T),"[",2))  # 2, 8
row_annotation$species <- gsub("_for|_rev|\\.for|\\.rev","",row_annotation$species,perl = T)
row_annotation$species <- factor(row_annotation$species,levels = c(rna,dna)) # rna
}
row_annotation$group <- TSG_plot$group #[match(rownames(row_annotation),rownames(TSG_plot))]
row_annotation <- row_annotation[order(row_annotation$group,row_annotation$species),] # not sort this ?

#table(row_annotation$group)
if(color_list=="NULL"){
  message("color_list not given")
if(rnaColumn==8){
  color_list <- list(group=pal_npg_adaptive()( length((unique(col_annotation$group)) )),
                     species=c(pal_nejm_adaptive()(15)[1:14],"#11838D","grey70","grey30") #pal_d3_adaptive()(15)[i]
  )
  names(color_list$group) <- unique(col_annotation$group)
names(color_list$species) <- c(rna2,dna)
}else if(rnaColumn==2){
  color_list <- list(group=pal_npg_adaptive()( length((unique(col_annotation$group)) )),
                     species=c(pal_nejm_adaptive()(15)[1:14],"#11838D") #pal_d3_adaptive()(15)[i]
  )
  names(color_list$group) <- unique(col_annotation$group)
names(color_list$species) <- c(rna,dna)
}

}else if(file.exists(color_list)){
  message("color_list from ",color_list)
  color_list <- readRDS(color_list)
}else{
  stop("no color rda provided")
}


mat_plot <- mat[rownames(row_annotation),rownames(col_annotation)]
bk = unique(c(seq(0,1, length=100)))
plot_sample <- pheatmap(
  mat_plot, # [order(TSG_plot$group),]
  breaks = bk,border_color = "NA",
  annotation_colors = color_list,
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
save_pheatmap_pdf(plot_sample, paste0(outdir,"/","signature_genes_samples_",TPM_cutoff,"_",TSS_cutoff,"_",Num,"_",method,".png") ,device = "png", width=7, height=7)
dev.off()


#plot signature gene matrix
col_annotation <- data.frame(row.names = colnames(TSG_plot)[1:(ncol(TSG_plot)-1)],group=colnames(TSG_plot)[1:(ncol(TSG_plot)-1)]) #sample_info[!duplicated(sample_info$group),]
rownames(col_annotation) <- col_annotation$group
col_annotation$sample <- NULL

bk = unique(c(seq(0,1, length=100)))
plot_tissue <- pheatmap(
  TSG_plot[rownames(row_annotation),-which(colnames(TSG_plot)=="group")],
  breaks = bk, border_color = "NA",
  annotation_col = col_annotation,
  annotation_row = row_annotation,
  annotation_colors = color_list,
  scale = "row",
  cluster_cols = FALSE,cluster_rows = FALSE,
  show_colnames=F, show_rownames=FALSE,  
  # color = (colorRampPalette(RColorBrewer::brewer.pal(7, "Reds"))(100)), #colorRampPalette(c("steelblue","white","salmon"))(100),
  color = viridis::viridis_pal()(100),
  # annotation_colors = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(unique(row_annotation$gene))),
  fontsize_row = 5,
  clustering_distance_cols = "euclidean")

save_pheatmap_pdf(plot_tissue, paste0(outdir,"/","signature_genes_tissues_",TPM_cutoff,"_",TSS_cutoff,"_",Num,"_",method,".pdf"))
dev.off()



# plot sig_gene/row_anno seperately
#row_annotation[1:3,]
#dim(row_annotation)
ggplot(row_annotation,aes(x=group,fill=species))+
  geom_bar(stat="count",position="stack")+ #,position="fill"
  geom_vline(xintercept = 0,linetype="dashed",color="grey10")+
  #ggsci::scale_fill_d3()+
  scale_fill_manual(name="Species",values = color_list[['species']] )+ # c(RColorBrewer::brewer.pal(8, "Accent"),RColorBrewer::brewer.pal(9, "Set1"))
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
ggsave(filename = paste0(outdir,"/","signature_stat_bar_",TPM_cutoff,"_",TSS_cutoff,"_",Num,"_",method,".pdf") , width = 12, height = 12)



# # plot venn/upset
# sig.tab <- data.table::fread(paste0(outdir, "/all_gene_list_", TPM_cutoff, "_", TSS_cutoff, "_", Num, ".txt"), sep="\t", check.names = F, header = T)
# sig.tab <- reshape2::acast(data = sig.tab, formula = gene ~ group, value.var = "group", id.var="gene")


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
        scale_fill_manual(name="Group",values = color_list[['group']])+ 
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
      ggsave(filename = paste0(outdir,"/","tSNE_all_",TPM_cutoff,"_",TSS_cutoff,"_",Num,"_",method,"_perplexity",i,".pdf"), width =  18, height = 12)
    },
    error = function(e){
      message("seem err in tSNE")
    }
  )
}

sig <- read.table(paste0(outdir,"/","gene_list_",TPM_cutoff,"_",TSS_cutoff,"_",Num,"_",method,".txt"),sep="\t",row.names = 1, header = T,stringsAsFactors = F,check.names = F)
mat <- mat[rownames(mat) %in% rownames(sig),]
for(i in seq(5,50,10)){
  tryCatch(
    expr={x <- tidydr::dr(data = t(mat), fun = Rtsne::Rtsne,perplexity = i,max_iter=1000,dims=2)# uwot::tumap, uwot::umap, uwot::lvish, Rtsne::Rtsne, stats::prcomp
    autoplot(x, aes(fill=group), metadata = sample_info[, "group", drop=FALSE] ) +
      geom_point(size=5,alpha=0.99,shape=21,color="white")+
      theme_dr()+
      scale_fill_manual(name="Group",values = color_list[['group']])+ 
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
    ggsave(filename = paste0(outdir,"/","tSNE_sig_",TPM_cutoff,"_",TSS_cutoff,"_",Num,"_",method,"_perplexity",i,".pdf"), width =  18, height = 12)
    },
    error = function(e){
      message("seem err in tSNE")
    }
  )
}


## plot UMAP
#dim(t(mat))
umap.tmp <- umap::umap(t(mat))
#all(sample_info$sample==colnames(mat)) # all T!
sample_info[[paste0("Dim1")]] = umap.tmp$layout[,1]
sample_info[[paste0("Dim2")]] = umap.tmp$layout[,2]
sample_info[["Color"]] <- sample_info[["group"]]
#library(ggplot2)
plotDotSize <- 5
plotDotAlpha <- 0.95
figWid <- 18
figHigh <- 12
legend.position <- "none"
my_theme <- ggplot2::theme(aspect.ratio=1,
                    plot.title = ggplot2::element_text(size = 32,color="black",hjust = 0.5),
                    axis.title = ggplot2::element_text(size = 32,color ="black"),
                    axis.text = ggplot2::element_text(size= 30,color = "black"), #,face="bold
                    panel.grid.minor.y = ggplot2::element_blank(),
                    panel.grid.minor.x = ggplot2::element_blank(),
                    axis.text.x = ggplot2::element_text( hjust = 0.5 ), # angle = 45,
                    axis.text.y = ggplot2::element_text( hjust = 0 ), # angle = 45,
                    panel.grid=ggplot2::element_blank(),
                    legend.position = legend.position,#c(0.5,0.3),
                    legend.text = ggplot2::element_text(size= 20,color = "black"),
                    legend.title= ggplot2::element_text(size= 32,color = "black"))
p.tmp <- ggplot2::ggplot(data=sample_info,ggplot2::aes(x=Dim1,y=Dim2)) +
  ggplot2::geom_point(size=plotDotSize,alpha=plotDotAlpha, ggplot2::aes(color=Color)) + # ,fill=Fill,shape=Shape
  ggplot2::labs(title = paste0("UMAP.",nrow(mat))) +
  ggplot2::xlab(paste("Dim1")) +
  ggplot2::ylab(paste("Dim2")) +
  ggplot2::theme_minimal()+
  ggplot2::scale_color_manual(name="Color", values = color_list[["group"]])+
  # scale_fill_manual(name="Fill", values = plotGrpFill)+
  # scale_shape_discrete(name="Shape", values = plotGrpShape )+
  my_theme
ggsave(filename = paste0(outdir,"/","UMAP_sig_",TPM_cutoff,"_",TSS_cutoff,"_",Num,"_",method,".png"),plot = p.tmp,width = figWid,height = figHigh)
#data.table::fwrite(data.frame("sample"=rownames(sample.table),"UMAP_Dim1"=sample.table$Dim1,"UMAP_Dim2"=sample.table$Dim2),paste0(outFile,".UMAP.txt"),quote = F,sep = "\t",row.names = F,col.names = T)
