# test OSCC diff
# last 220322 by bpf
# b.p.f@qq.com

setwd("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/")
options(scipen = 4,stringsAsFactors = F,digits = 5)


# load all func.
source("./util.R")


setwd("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/")
options(scipen = 4,stringsAsFactors = F,digits = 5)
dst <- "lulab_oscc_plasma_diff"

## read meta
sample.table0 <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/lulab_oscc_plasma_diff/sample_table.txt",check.names = F,sep = "\t",header = T)

## read mat
count <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/lulab_oscc_plasma_diff/call_peak_dedup/count_matrix/cfpeakCNN_b5_d50_p1.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F)
#cpm0 <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/lulab_oscc_plasma_diff/call_peak_dedup/count_matrix/expeak_b5_d50_p1_CPM.txt",check.names = F,header = T)
rownames(count) <- count$feature
mat <- count[,2:ncol(count)]


### Meta vs. Local
colnames(sample.table0)
sample.table <- sample.table0[,c("Patient ID","Metastasis","Origin")]
colnames(sample.table) <- c("sample","group","origin")
rownames(sample.table) <- sample.table$sample
sample.table <- sample.table[sample.table$sample %in% colnames(count),]
sample.table <- sample.table[order(sample.table$group),]
positive_samples <- sample.table[sample.table$group=="Yes","sample"]
negative_samples <- sample.table[sample.table$group=="No","sample"]
samples <- c(positive_samples, negative_samples)
sample.table <- sample.table[samples,]
table(sample.table$group)

mat <- mat[,samples]

group <- c(rep("positive",length(positive_samples)),rep("negative",length(negative_samples)))
method <- "edger_glmlrt"
norm_method <- "TMM"

#res <- diff.v2(mat = mat, samples = samples, group = group, method = method, norm_method = norm_method, filterType = "NULL", featureType = "domain")
res <- diff.v2.dcb(mat = mat, samples = samples, group = group, method = method, norm_method = norm_method, filterType = "NULL", featureType = "domain")
#write.table(res,"./output/lulab/smallDomain_diff_EVvsCF_pair.txt",quote = F,sep = "\t",row.names = T,col.names = T)
dir.create(paste0("./output/",dst))
# write.table(res[["normMat"]],paste0("./output/",dst,"/allSmp_exPeak_smallDomain_diff_CRCvsNC.cpm"),quote = F,sep = "\t",row.names = T,col.names = T) # allSmp_
# write.table(res[["diffTable"]],paste0("./output/",dst,"/allSmp_exPeak_smallDomain_diff_CRCvsNC.diff"),quote = F,sep = "\t",row.names = T,col.names = T) # allSmp_
write.table(res[["normMat"]],paste0("./output/",dst,"/cfPeakCNN_smallDomain_diff_MvsL.cpm"),quote = F,sep = "\t",row.names = T,col.names = T) # allSmp_
write.table(res[["diffTable"]],paste0("./output/",dst,"/cfPeakCNN_smallDomain_diff_MvsL.diff"),quote = F,sep = "\t",row.names = T,col.names = T) # allSmp_


## select peak by top 
peak.list <- list()
seqs <- c(50,100,200,500,1000)
for (top in seqs){
  peak.list[[paste0("top",top)]] <- list()
}
disease <- "metastasis"
disease.label <- "M"
normal.label <- "L"

tmp <- read.table(paste0("./output/",dst,"/cfPeakCNN_smallDomain_diff_",disease.label,"vs",normal.label,".diff"),sep = "\t",header = T,check.names = F,stringsAsFactors = F)
tmp <- tmp[order(as.numeric(tmp$pvalue),decreasing = F),]
rownames(tmp) <- unlist(sapply(strsplit(rownames(tmp),"|",fixed=T),"[",4))

# op1: fixed cutoff
#tmp.sig <- tmp[tmp$pvalue<0.00001,] 

# op2: fixed top
for (i in seqs){
  peak.id <- rownames(tmp)[1:i]
  # print(length(peak.id))
  peak.list[[paste0("top",i)]][[disease]] <- peak.id
}
peak.ids <- list()
for (top in seqs){
  peak.ids[[paste0("top",top)]] <- unique(do.call("c",peak.list[[paste0("top",top)]])) 
  print(length( peak.ids[[paste0("top",top)]] ))
}



## PCA plot
dst <- "lulab_oscc_plasma_diff"
disease.label <- "M"
normal.label <- "L"
diff <- read.table(paste0("./output/",dst,"/cfPeakCNN_smallDomain_diff_",disease.label,"vs",normal.label,".diff"), header = T, sep="\t",check.names = F)
logcpm <- read.table(paste0("./output/",dst,"/cfPeakCNN_smallDomain_diff_",disease.label,"vs",normal.label,".cpm"), header = T, sep="\t",check.names = F)

logcpm <- logcpm[matrixStats::rowSds(as.matrix(logcpm))!=0,] # rm 0 std, or meet error for pca

# logcpm <- logcpm[unlist(sapply(strsplit(rownames(logcpm),"|",fixed=T),"[",4)) %in% peak.ids[["top100"]],] # rownames(diff)[diff$pvalue<=0.01]
#table(diff$padj<=0.01)
dim(logcpm)

# sample.table <- sample.table[sample.table$sample %in% colnames(logcpm),]
# table(sample.table$group)

pca <- stats::prcomp(t(logcpm), scale=TRUE)  ###prcomp
#choose top2 PC
pca.var <- pca$sdev^2  ## sdev
pca.var.per <- round(pca.var/sum(pca.var)*ncol(logcpm), 1)  #
barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")  ##
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2],
                       Z=pca$x[,3],
                       A=pca$x[,4])

library(ggplot2)
my_theme <-   theme(plot.title = element_text(size = 28,color="black",hjust = 0.5),
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
pca.plot1 <-
  ggplot(data=pca.data,aes(x=Y,y=X,fill=sample.table.pair$group))+
  geom_point(size=6,shape=21,color="white") +
  xlab(paste("")) +
  ylab(paste("PC1(",pca.var.per[1],"%","Var)",sep="")) +
  theme_bw()+
  # ggsci::scale_color_d3(name="Source")+
  ggsci::scale_fill_d3(name="Group")+
  #scale_color_discrete()+
  scale_shape_discrete(name="Disease")+
  my_theme
pca.plot2 <-
  ggplot(data=pca.data,aes(x=A,y=X,fill=sample.table.pair$group))+
  geom_point(size=6,shape=21,color="white") +
  xlab(paste("")) +
  ylab(paste("")) +
  theme_bw()+
  # ggsci::scale_color_d3(name="Source")+
  ggsci::scale_fill_d3(name="Group")+
  #scale_color_discrete()+
  scale_shape_discrete(name="Disease")+
  my_theme
pca.plot3 <-
  ggplot(data=pca.data,aes(x=Y,y=Z,fill=sample.table.pair$group))+
  geom_point(size=6,shape=21,color="white") +
  xlab(paste("PC2(",pca.var.per[2],"%","Var)",sep="")) +
  ylab(paste("PC3(",pca.var.per[3],"%","Var)",sep="")) +
  theme_bw()+
  # ggsci::scale_color_d3(name="Source")+
  ggsci::scale_fill_d3(name="Group")+
  #scale_color_discrete()+
  scale_shape_discrete(name="Disease")+
  my_theme
pca.plot4 <-
  ggplot(data=pca.data,aes(x=A,y=Z,fill=sample.table.pair$group))+
  geom_point(size=6,shape=21,color="white") +
  xlab(paste("PC4(",pca.var.per[4],"%","Var)",sep="")) +
  ylab(paste("")) +
  theme_bw()+
  # ggsci::scale_color_d3(name="Source")+
  ggsci::scale_fill_d3(name="Group")+
  #scale_color_discrete()+
  scale_shape_discrete(name="Disease")+
  my_theme
pca.plot <- ggpubr::ggarrange(align = "hv",legend = "right",common.legend = TRUE,ncol = 2,nrow = 2,plotlist = list(pca.plot1,pca.plot2,pca.plot3,pca.plot4))
# ggsave(filename = paste0("./output/",dst,"/allSmp_top100_cfPeakCNN_smallDomain_diff_",disease.label,"vs",normal.label,".cpm.PCA.pdf"),plot = pca.plot,width = 15,height = 12)
ggsave(filename = paste0("./output/",dst,"/allSmp_cfPeakCNN_smallDomain_diff_",disease.label,"vs",normal.label,".cpm.PCA.pdf"),plot = pca.plot,width = 15,height = 12)




## heatmap
ref <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt") #,header=T
ref <- as.data.frame(ref)
rownames(ref) <- ref$transcript_id
ref[1:3,1:3]


dst <- "lulab_oscc_plasma_diff"
disease.label <- "M"
normal.label <- "L"
diff <- read.table(paste0("./output/",dst,"/cfPeakCNN_smallDomain_diff_",disease.label,"vs",normal.label,".diff"), header = T, sep="\t",check.names = F)
logcpm0 <- read.table(paste0("./output/",dst,"/cfPeakCNN_smallDomain_diff_",disease.label,"vs",normal.label,".cpm"), header = T, sep="\t",check.names = F)
all(rownames(diff)==rownames(logcpm0)) # T !!!


seqs <- 1000
for(top in paste0("top",seqs)){
  print(top)
  # for(trend in c("+","-","")){
    # print(trend)
    tmp.id <- unlist(sapply(strsplit(rownames(logcpm0),"|",fixed=T),"[",4))
    # if(trend=="+"){
      # filter1 <- tmp.id %in% peak.ids[[top]] & diff$log2FoldChange>0
    # }else if(trend=="-"){
      # filter2 <- tmp.id %in% peak.ids[[top]] & diff$log2FoldChange<0
    # }else if(trend==""){
    #   filter <- tmp.id %in% peak.ids[[top]]
    # }
    logcpm <- logcpm0[tmp.id %in% peak.ids[[top]],] # rownames(diff)[diff$pvalue<=0.01]
#table(diff$pvalue<=0.01)

annotation_col <- sample.table[colnames(logcpm),]
annotation_col <- annotation_col[,c("group","origin"),drop=F]

txid <- unlist(sapply(strsplit(rownames(logcpm),"|",fixed=T),"[",3))
annotation_row <- data.frame( RNA = ref[txid,"transcript_type"], row.names = rownames(logcpm))
annotation_row$RNA <- gsub("_for|_rev","",annotation_row$RNA,perl = T)
tmp.id2 <- unlist(sapply(strsplit(rownames(logcpm),"|",fixed=T),"[",4))
annotation_row$trend <- ""
annotation_row$trend[diff[rownames(logcpm),]$log2FoldChange>0] <- "up"
annotation_row$trend[diff[rownames(logcpm),]$log2FoldChange<0] <- "down"
annotation_row$trend <- factor(annotation_row$trend,levels = c("up","down"))
table(annotation_row$trend)
# filter1 <- tmp.id %in% peak.ids[[top]] & diff$log2FoldChange>0
# filter2 <- tmp.id %in% peak.ids[[top]] & diff$log2FoldChange<0
tmp <- as.data.frame(table(annotation_row$RNA))
tmp$lab <- paste0(tmp$Var1," (n=",tmp$Freq,")")
# var1 <- tmp$Var1
annotation_row$RNA.lab <-tmp$lab[match(annotation_row$RNA,tmp$Var1)]
annotation_row$RNA <- factor(annotation_row$RNA, levels = c(rna,dna))
annotation_row <- annotation_row[order(annotation_row$RNA),]
table(annotation_row$RNA)


RNA_colors <- list()
for(i in 1:length(c(rna,dna))){
  j <- c(rna,dna)[i]
  RNA_colors[[j]] <- c(pal_nejm_adaptive()(15)[1:14],"#11838D")[i] # pal_d3_adaptive()(15)
}
RNA_colors <- do.call("c",RNA_colors)
names(RNA_colors) <- tmp$lab[match(names(RNA_colors),tmp$Var1)]
RNA_colors <- RNA_colors[!is.na(names(RNA_colors))]
ann_colors <- list(Metastasis=c("Yes"="firebrick2","No"="orange2"),Origin=c("cavity"="seagreen4","tongue"="orange4","lip"="chocolate1"),Regulation=c("up"="firebrick","down"="steelblue4"),Precursor=RNA_colors) # ,"salmon",seagreen,orange,chocolate

annotation_row <- annotation_row[order(annotation_row$trend,annotation_row$RNA,decreasing = F),]
annotation_row$RNA <- NULL
logcpm <- logcpm[rownames(annotation_row),] #positive_samples
colnames(annotation_row) <- c("Regulation","Precursor")
colnames(annotation_col) <- c("Metastasis","Origin")
annotation_col$Origin <- NULL
pheatmap::pheatmap(mat = logcpm,
                   annotation_col = annotation_col, #data frame that specifies the annotations shown on left side of the heatmap.
                   annotation_row = annotation_row,
                   annotation_colors = ann_colors, # RColorBrewer::brewer.pal(6, "Set3")
                   border_color="NA",
                   scale = "row", # row,none
                   #labels_col = 3, labels_row = 6,
                   cluster_cols = T,cluster_rows = F,
                   gaps_row = sum(annotation_row$trend=="up"), # invalid if cluster row/col
                   #cutree_cols = 2,cutree_rows = 3,
                   show_colnames=F, show_rownames=F,
                   fontsize = 12,
                   height = 12, width = 7,
                   color = rev( colorRampPalette( RColorBrewer::brewer.pal(n = 10, name = "RdBu"))(100) ),  # steelblue/dodgerblue4-firebrick2, viridis::viridis_pal()(100), #
                   #fontsize_row = 5,
                   filename = paste0("./output/",dst,"/allSmp_",top,"_cfPeakCNN_smallDomain_diff_",disease.label,"vs",normal.label,".cpm.heatmap3.pdf")
                   # filename = paste0("./output/",dst,"/allSmpP01_exPeak_smallDomain_diff_",disease.label,"vs",normal.label,".cpm.heatmap.pdf")
                   )
}


#prepare igv tbl 
# for (peak in rownames(annotation_row_merge)){
#   #peak <- "ENST00000385214_____1_30_49_+|pri_miRNA|ENST00000385214_____1|peak_6181|ENST00000385214_____1|30|49"  # "T295851_13316_13332_+|tucpRNA|T295851|peak_57431|T295851|13316|13332"
#   # peak <- "L2a__chr16___3509264____3509549_neg_84_99_+|repeats_rev|L2a__chr16___3509264____3509549_neg|peak_48841|L2a__chr16___3509264____3509549_neg|84|99"
#   for (dataset in unique(annotation_col_merge$dataset)){
#     #dataset <- "GSE110381 plasma"
#     tmp1.annotation_col_merge <- annotation_col_merge[annotation_col_merge$dataset==dataset,]
#     tmp1 <- as.data.frame(cpm.merge[peak,annotation_col_merge$dataset==dataset,drop=F])
#     for (group in unique(tmp1.annotation_col_merge$group)){
#       #group <- "CRC plasma"
#       # annotation_row_merge[peak,paste0(dataset,"_",group)] <- ""
#       
#       tmp2.annotation_col_merge <- tmp1.annotation_col_merge[tmp1.annotation_col_merge$group==group,]
#       tmp2 <- as.data.frame(t(tmp1[,tmp1.annotation_col_merge$group==group,drop=F]))
#       if(grepl("CRC",group)){
#         tmp2 <- tmp2[order(tmp2[,1],decreasing = T),,drop=F]
#       } else {
#         tmp2 <- tmp2[order(tmp2[,1],decreasing = F),,drop=F]
#       }
#       top3 <- rownames(tmp2)[1:3]
#       annotation_row_merge[peak,paste0(dataset,"_",group,"_",1:3)] <- top3
#     }
#   }
# }
# annotation_row_merge$peak <- unlist(sapply(strsplit(rownames(annotation_row_merge),"|",fixed=T),"[",4))
# annotation_row_merge$chr <- unlist(sapply(strsplit(rownames(annotation_row_merge),"|",fixed=T),"[",5))
# annotation_row_merge$start <- unlist(sapply(strsplit(rownames(annotation_row_merge),"|",fixed=T),"[",6))
# annotation_row_merge$end <- unlist(sapply(strsplit(rownames(annotation_row_merge),"|",fixed=T),"[",7))
# annotation_row_merge$width <- as.numeric(annotation_row_merge$end) - as.numeric(annotation_row_merge$start)
# annotation_row_merge$idx <- 1:nrow(annotation_row_merge)
# annotation_row_merge <- annotation_row_merge[,c((ncol(annotation_row_merge)-5):ncol(annotation_row_merge),1:(ncol(annotation_row_merge)-6))]
# #use annotation_row_merge for IGV eg.
# annotation_row_merge$enst <- ""
# enst.idx <- which(grepl("ENST",annotation_row_merge$chr))
# annotation_row_merge$enst[enst.idx] <- unlist(sapply(strsplit(annotation_row_merge$chr[enst.idx],"_"),"[",1))
# tmp <- as.data.frame(mygene::queryMany(annotation_row_merge$enst[enst.idx],scopes="ensembl.transcript",fields=c("ensembl.gene","symbol"),species="human"))
# dim(tmp)
# table(duplicated(tmp$ensembl.gene))
# table(duplicated(tmp$query))
# library(dplyr)
# library(biomaRt)
# mart <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")
# ann <- biomaRt::getBM(useCache = FALSE,c("hgnc_symbol","description","chromosome_name","band","strand","start_position","end_position","ensembl_gene_id"),"ensembl_gene_id", tmp$ensembl.gene, mart)
# dim(ann)
# table(duplicated(ann$hgnc_symbol))
# table(duplicated(ann$ensembl_gene_id))
# tmp$func <- ann$description[match(tmp$ensembl.gene,ann$ensembl_gene_id)]
# length(enst.idx)==nrow(tmp) # T !!!
# annotation_row_merge$func <- ""
# annotation_row_merge$func[enst.idx] <- tmp$func
# annotation_row_merge$symbol <- ""
# annotation_row_merge$symbol[enst.idx] <- tmp$symbol
# write.table(x = annotation_row_merge, file = "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/tmp/PRNNA-TCGA3_annotation_row_merge.txt",quote = T,sep = "\t",row.names = F,col.names = T)







## plot diff peak RNA type (pie plot)
df <- data.frame(peak.id = unlist(sapply(strsplit(rownames(mat),"|",fixed=T),"[",4)),
                 RNA = unlist(sapply(strsplit(rownames(mat),"|",fixed=T),"[",2)))
df$RNA <- gsub("_rev|_for","",df$RNA,perl=T)
rna <- c("rRNA","pri_miRNA","lncRNA", "mRNA", "piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA')
dna <- c("intron", "promoter", "enhancer", "repeats")
df$RNA <- factor(df$RNA,levels = c(rna,dna))
table(df$RNA)

df <- df[df$peak.id %in% peak.ids[['top1000']],] # 50
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
ggsave("lulab_oscc_plasma_diff_pie.pdf",width = 12,height = 7)





#diff volcano plot
dst <- "lulab_oscc_plasma_diff"
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
                    legend.text = element_text(size= 24,color = "black"),
                    legend.title= element_text(size= 26,color = "black"))
diff.plot <- list()
# tmp <- read.table(paste0("./output/",dst,"/exPeak_smallDomain_diff_",disease.label,"vs",normal.label,".diff"),sep = "\t",header = T,check.names = F)
tmp <- read.table(paste0("./output/",dst,"/cfPeakCNN_smallDomain_diff_",disease.label,"vs",normal.label,".diff"),sep = "\t",header = T,check.names = F,stringsAsFactors = F)

tmp$minus.log10.pvalue <- -log10(tmp$pvalue)
tmp$log10.baseMean <- log10(tmp$baseMean)

tmp$alpha <- maxmin.normalize.vec(tmp$minus.log10.pvalue)
# ggplot(tmp, aes(x=log2FoldChange,y=minus.log10.pvalue)) + 
#   geom_point() 
# ggplot(tmp, aes(x=baseMean.pos,y=baseMean.neg)) + 
#   geom_point() 
diff.plot[[disease.label]] <- ggplot(tmp, aes(x=log10.baseMean,y=log2FoldChange, fill=minus.log10.pvalue)) + 
  geom_point(aes(alpha=alpha),shape=21,color="white") + 
  labs(title = paste0(disease.label," vs. ",normal.label)) +
  ggraph::scale_fill_viridis(direction = -1) +
  geom_hline(yintercept = 0,linetype="dashed",color="salmon") + 
  theme_bw()+
  my_theme

# p <- ggpubr::ggarrange(plotlist = diff.plot, rel_widths = c(1,1,1), nrow = 1,  common.legend = T,align = "hv", axis = "b", legend = "right") #  labels = c("RBPs", "", "EV", "", "G4"),
ggsave(plot = diff.plot[[disease.label]], paste0("./lulab_oscc_plasma_diff_basemean.pdf"), width=9, height=9) # "_",sample








## 5. GO/KEGG ORA func. enrichment
options(stringsAsFactors = F)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(DOSE))
suppressPackageStartupMessages(library(enrichplot))
suppressPackageStartupMessages(library(org.Hs.eg.db))

### 5.1 prepare ref & func.
enst2entrezDF <- function(enst){
  library(mygene)
  tmp <- as.data.frame(mygene::queryMany(enst, scopes="ensembl.transcript",fields=c("entrezgene"),species="human")) #$entrezgene
  return(tmp)
}
enst2ensgDF <- function(enst){
  library(mygene)
  tmp <- as.data.frame(mygene::queryMany(enst, scopes="ensembl.transcript",fields=c("ensembl.gene"),species="human")) #$entrezgene
  return(tmp)
}

ora.kegg.offline <- function(gene.top.up,gene.top.down,gene.list,
                             m_t,
                             qValue=0.2, pValue=0.1, adjustPval="BH"){

  # 1.Over-Representation Analysis (ncbiid or otherId in TERM2GENE)
  # qValue <- 0.2
  # pValue <- 0.1
  # adjustPval <- "BH"
  # 1.1 ORA all gene
  # enrich.KEGG <- enricher(gene.top.ncbiid,
  #                         #universe = names(ccle.gene.spearman$g.l), 
  #                         TERM2GENE=m_t2e,
  #                         pvalueCutoff = pValue,
  #                         pAdjustMethod = adjustPval,
  #                         qvalueCutoff = qValue
  #                         #, maxGSSize = 350
  # )
  # 1.2 ORA up gene
  enrich.KEGG.up <- enricher(gene.top.up,
                             universe = gene.list, #names(ccle.gene.spearman$g.l),
                             TERM2GENE=m_t,
                             pvalueCutoff = pValue,
                             pAdjustMethod = adjustPval,
                             qvalueCutoff = qValue
                             #, maxGSSize = 350
  )
  enrich.KEGG.up@result$trend <- "up"
  enrich.KEGG.up@result$minus.log10.pvalue <- -log10(enrich.KEGG.up@result$pvalue)
  enrich.KEGG.up@result$minus.log10.padj <- -log10(enrich.KEGG.up@result$p.adjust)
  enrich.KEGG.up@result$minus.log10.qvalue <- -log10(enrich.KEGG.up@result$qvalue)
  
  # 1.3 ORA down gene
  enrich.KEGG.down <- enricher(gene.top.down,
                               universe = gene.list, #names(ccle.gene.spearman$g.l),
                               TERM2GENE=m_t,
                               pvalueCutoff = pValue,
                               pAdjustMethod = adjustPval,
                               qvalueCutoff = qValue
                               #, maxGSSize = 350
  )
  enrich.KEGG.down@result$trend <- "dw"
  enrich.KEGG.down@result$minus.log10.pvalue <- -log10(enrich.KEGG.down@result$pvalue)
  enrich.KEGG.down@result$minus.log10.padj <- -log10(enrich.KEGG.down@result$p.adjust)
  enrich.KEGG.down@result$minus.log10.qvalue <- -log10(enrich.KEGG.down@result$qvalue)
  
  return(as.data.frame(rbind(enrich.KEGG.up@result,enrich.KEGG.down@result)))
}

##offline
m_t <- read.table("/BioII/lulab_b/baopengfei/shared_reference/geneset/PATH_ID_NAME_KEGGplusHallmark.txt",sep = "\t",header = T, stringsAsFactors = F,row.names = 1)
m_t2g <- m_t[,c("short_DESCRPTION","ensembl_gene_id")]
m_t2e <- m_t[,c("short_DESCRPTION","ENTREZID")]

ref <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/transcript_table/all_newTxID.txt") #,header=T
ref <- as.data.frame(ref)
colnames(ref)
#rownames(ref) <- ref$transcript_id
ref[1:3,]
ref$gene_id2 <- unlist(sapply(strsplit(ref$gene_id,".",fixed=T),"[",1)) 
rna <- c("rRNA","pri_miRNA","lncRNA", "mRNA", "piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA')
dna <- c("intron","promoter", "enhancer", "repeats") 
ref$transcript_type <- gsub("_for|_rev","",ref$transcript_type,perl = T)
ref$transcript_type <- factor(ref$transcript_type,levels = c(rna,dna))
table(ref$transcript_type)
ref$transcript_id2 <- gsub("_____",".",ref$transcript_id)


### 5.2 read diff tbl
disease.label <- "M"
normal.label <- "L"
diff <- read.table(paste0("./output/",dst,"/cfPeakCNN_smallDomain_diff_",disease.label,"vs",normal.label,".diff"), header = T, sep="\t",check.names = F)
colnames(diff)
diff <- diff[,c("log2FoldChange","pvalue","baseMean")]
diff$peak.id <- unlist(sapply(strsplit(rownames(diff),"|",fixed=T),"[",4))
diff$RNA <- unlist(sapply(strsplit(rownames(diff),"|",fixed=T),"[",2))
diff$RNA <- gsub("_for|_rev","",diff$RNA,perl = T)
diff$RNA <- factor(diff$RNA,levels = c(rna,dna))
diff$enst <- unlist(sapply(strsplit(rownames(diff),"|",fixed=T),"[",3)) 
table(unique(diff$enst) %in% unique(ref$transcript_id)) # all TRUE
# diff$enst <- gsub("_____",".",diff$enst,fixed = T)

enst.idx <- (grepl("ENST",diff$enst))
table(enst.idx)
diff$enst[!enst.idx] <- ""
diff$enst2 <- substr(diff$enst,1,15)
diff$ensg[enst.idx] <- ref$gene_id[match(diff$enst[enst.idx],ref$transcript_id)]
diff$ensg2 <- ""
diff$ensg2[enst.idx] <- unlist(sapply(strsplit(diff$ensg[enst.idx],".",fixed=T),"[",1)) 

# tmp <- enst2entrezDF(unique(diff$enst2[enst.idx]))
# tmp <- tmp[,c("query","entrezgene")]
# tmp <- tmp[!duplicated(tmp$query),]
# diff$entrez <- ""
# diff$entrez[enst.idx] <- tmp$entrezgene[match(diff$enst2[enst.idx],tmp$query)]
# 
# tmp <- enst2ensgDF(unique(diff$enst2[enst.idx]))
# tmp <- tmp[is.na(tmp$notfound),]
# table(is.na(tmp$ensembl))
# tmp <- tmp[,c("query","ensembl")]
# tmp$ensembl2 <- ""
# for(i in 1:nrow(tmp)){
#   # i <- 58
#   tmp$ensembl2[i] <- as.character(tmp$ensembl[i][[1]])
# }
# tmp$ensembl2[grepl("c",tmp$ensembl2)] <- unlist(sapply(strsplit(tmp$ensembl2[grepl("c",tmp$ensembl2)],"\"",fixed=T),"[",2)) 
# tmp <- tmp[!duplicated(tmp$query),]
# diff$ensg <- ""
# diff$ensg[enst.idx] <- tmp$ensembl2[match(diff$enst2[enst.idx],tmp$query)]


### 5.3 get KEGG ORA tbl
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
ref2 <- list()
diff2 <- list()
univ <- list()


getORAtbl <- function(peak.ids,diff2,univ,enst.lab){
  #names(peak.ids)
  df.list <- list()
  for (seq in c(500,1000)){ # df.list
    # seq <- 1000
    print(seq)
    print(enst.lab)
    if(enst.lab=="DNA"){
      #table(is.na(diff2[[enst.lab]]$nearest.ensg))
      gene.top.up <- unique(diff2[[enst.lab]]$nearest.ensg[(diff2[[enst.lab]]$peak.id %in% peak.ids[[paste0("top",seq)]]) & (diff2[[enst.lab]]$log2FoldChange>0)])
      gene.top.down <- unique(diff2[[enst.lab]]$nearest.ensg[(diff2[[enst.lab]]$peak.id %in% peak.ids[[paste0("top",seq)]]) & (diff2[[enst.lab]]$log2FoldChange<0)])
    }else{
      gene.top.up <- unique(diff2[[enst.lab]]$ensg2[(diff2[[enst.lab]]$peak.id %in% peak.ids[[paste0("top",seq)]]) & (diff2[[enst.lab]]$log2FoldChange>0)])
      gene.top.down <- unique(diff2[[enst.lab]]$ensg2[(diff2[[enst.lab]]$peak.id %in% peak.ids[[paste0("top",seq)]]) & (diff2[[enst.lab]]$log2FoldChange<0)])
    }
    gene.list <- univ[[enst.lab]]
    gene.top.up <- gene.top.up[!is.na(gene.top.up)]
    gene.top.down <- gene.top.down[!is.na(gene.top.down)]
    gene.list <- gene.list[!is.na(gene.list)]
    
    # table(is.na(gene.list))
    tmp <- ora.kegg.offline(gene.top.up = gene.top.up, gene.top.down = gene.top.down, gene.list = gene.list, 
                            m_t = m_t2g )
    tmp$top <- seq
    df.list[[seq]] <- tmp
    # return(tmp)
  }
  
  df <- as.data.frame(do.call(rbind,df.list))
  df$lab <- enst.lab
  # table(df$trend)
  # table(df$top)
  return(df)
}



## 5.3.1. mRNA,lncRNA, pri_miRNA, #snoRNA,snRNA,srpRNA,Y_RNA
# lncRNA      mRNA pri_miRNA    snoRNA     snRNA    srpRNA     Y_RNA 
# 17030      4753      1231        47        87        16        82 
# for(j in c( "mRNA","lncRNA", "pri_miRNA","snoRNA","snRNA","srpRNA","Y_RNA")){
#   print(j)
#   diff2 <- diff[diff$RNA %in% c(j),]
#   print(table(diff2$ensg2 %in% m_t2g$ensembl_gene_id))
# }
# #seem only mRNA,pri_miRNA included in KEGG annotation !!

#gene universe list (better use all in annotation, instead of all peak enst)
enst.lab <- "mRNAmiRNA"
enst.types <- c("mRNA","pri_miRNA")
ref2[[enst.lab]] <- ref[ref$transcript_type %in% enst.types,]
diff2[[enst.lab]] <- diff[diff$RNA %in% enst.types,]
#table(ref.tmp$transcript_id %in% diff$enst)
univ[[enst.lab]] <- unique(ref2[[enst.lab]]$gene_id2)

df <- getORAtbl(peak.ids,diff2,univ,enst.lab)
write.table(df, paste0(pre,"/tmp/",dst,"-",enst.lab,"_ORA.txt"), quote = F, sep = "\t", row.names = F, col.names = T)




# 5.3.2. DNA
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
library(bedtoolsr)
options(bedtools.path = "/BioII/lulab_b/baopengfei/biosoft")
#For each feature in A, finds the closest feature (upstream or downstream) in B
#bedtools closest -t first -D ref -a -b
dna.sort <- bedtoolsr::bt.sort(i = "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/lulab_oscc_plasma_diff/call_peak_dedup/cfpeakCNN/b5_d50_p1_8DNA_gn.bed", 
                               g = "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/transcriptome_genome_sort_uniq_newTxID" )
nearest.mir <- bedtoolsr::bt.closest(t = "first", D = "ref",  # 
                                     a = dna.sort, 
                                     b = "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/bed_by_biotype/mRNAmiRNA.bed" )
nearest.mir <- nearest.mir[nearest.mir$V16!=".",]
# nearest.mir[1:3,]
# hist(nearest.mir$V25,breaks = 1000000)
# table(nearest.mir$V25==0) # half are overlapped !
# summary(nearest.mir$V25)
nearest.mir <- nearest.mir[abs(nearest.mir$V25)<=10000,] # keep record with nearest gene within 10kb
#table(unique(nearest.mir$V16) %in% unique(ref$transcript_id2))
nearest.mir$ensg <- ref$gene_id2[match(nearest.mir$V16,ref$transcript_id2)]

tmp <- data.table::fread( "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/bed_by_biotype/mRNAmiRNA.bed" )
#table(unique(tmp$V4) %in% unique(ref$transcript_id2)) # all T !!
#ref[10000:10003,]
tmp$ensg <- ref$gene_id2[match(tmp$V4,ref$transcript_id2)]

enst.lab <- "DNA"
enst.types <- dna
univ[[enst.lab]] <- unique(tmp$ensg) # unlist(sapply(strsplit(tmp$V4,".",fixed = T),"[",1)), currently use all, not expand 8DNA and find nearest genes !!!
ref2[[enst.lab]] <- ref[ref$transcript_type %in% enst.types,]
diff2[[enst.lab]] <- diff[diff$RNA %in% enst.types,]
diff2[[enst.lab]]$nearest.ensg <- nearest.mir$ensg[match(diff2[[enst.lab]]$peak.id,nearest.mir$V4)]

df <- getORAtbl(peak.ids,diff2,univ,enst.lab)
write.table(df, paste0(pre,"/tmp/",dst,"-",enst.lab,"_ORA.txt"), quote = F, sep = "\t", row.names = F, col.names = T)


# # 5.3.3. miR target (holding...)
# library(multiMiR)
# tmp <- multiMiR::get_multimir(mirna = s)@data
# tmp$cmb1 <- paste0(tmp$mature_mirna_id,":",tmp$target_symbol)
# tmp$cmb2 <- paste0(tmp$target_symbol,":",tmp$mature_mirna_id)
# grep("circ",tmp$target_symbol)
# #table(tmp$database) 
# my_cor_matrix$relation <- ""
# k <- paste0(my_cor_matrix$source,":",my_cor_matrix$target)
# my_cor_matrix$relation[(k %in% unique(tmp$cmb1)) | (k %in% unique(tmp$cmb2))] <- "miRNA-mRNA"
# table(my_cor_matrix$relation) # seems no records


# 5.4. dot plot
#str(df)
df <- list()
for(enst.lab in c("mRNAmiRNA","DNA")){
  df[[enst.lab]] <- read.table(paste0(pre,"/tmp/",dst,"-",enst.lab,"_ORA.txt"), sep = "\t", header = T)
}
df <- as.data.frame(do.call(rbind,df))
df$lab <- factor(df$lab, levels = c("mRNAmiRNA","DNA"))
df$trend <- factor(df$trend, levels = c("up","dw"))

df.sig <- df[df$pvalue<=0.05,]
df.sig$nchr <- nchar(as.character(df.sig$ID))
# df.sig$updw.level <- 0
# df.sig$updw.level[df.sig$trend == "up" & df.sig$lab == "mRNAmiRNA"] <- 1
# df.sig$updw.level[df.sig$trend == "up" & df.sig$lab == "DNA"] <- 2
# df.sig$updw.level[df.sig$trend == "dw" & df.sig$lab == "mRNAmiRNA"] <- 3
# df.sig$updw.level[df.sig$trend == "dw" & df.sig$lab == "DNA"] <- 4

# str(df.sig)
for(top in c(500,1000)){
  print(top)
  # table(df.sig$top == top)
  df.sig.tmp <- df.sig[df.sig$top == top,]
  df.sig.tmp <- df.sig.tmp[order(df.sig.tmp$trend,df.sig.tmp$lab,-df.sig.tmp$nchr,decreasing = T),] # ,df.sig.tmp$nchr
  df.sig.tmp$ID <- factor(df.sig.tmp$ID, levels = unique(df.sig.tmp$ID))
  
  ggplot(df.sig.tmp, aes(x=lab, y=ID) ) + # type2, type3
    geom_point(aes(color=minus.log10.pvalue, size=Count)) + # fill="grey",color="grey",  minus.log10.padj
    # scale_fill_manual(values = c("grey50","#ED782F")) + #"#ED782F","#B52B1A","#3E8C9D"
    # scale_color_manual(values = c("grey50","#ED782F")) +
    scale_color_continuous(type = "viridis") + # YlOrRd, viridis RColorBrewer::brewer.pal(n=8,name="Reds")
    scale_size_continuous(name = "Count") +
    # scale_fill_brewer(type = "seq", palette = RColorBrewer::brewer.pal(n=8,name="YlOrRd")) +
    ylab("") +
    scale_x_discrete(breaks = c("mRNAmiRNA","DNA"), labels = c("RNA","DNA")) +
    facet_grid(.~trend) + # lab
    theme_minimal() +  # base_size=12 void
    theme(#axis.ticks.x=element_blank(),
      #strip.text.y = element_blank(),
      # aspect.ratio = 10,
      strip.text = element_text(size=20),
      axis.title.x = element_text(size=20),
      axis.title.y = element_text(size=20),
      axis.text.x = element_text(size = 20, hjust = 1,vjust = 1, angle=45), #
      axis.text.y = element_text(size = 20, hjust = 1),
      plot.title = element_text(size=20),
      # strip.text = element_blank(),
      # legend.position = "none", #c(0.9,0.8),#,#
      legend.text = element_text(size= 16),
      legend.title= element_text(size= 16))
  ggsave(paste0("./",top,"_kegg_dot.pdf"), width = 10, height = 12)
}


