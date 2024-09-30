setwd("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/")
options(scipen = 4,stringsAsFactors = F,digits = 5)


# load all func.
source("./util.R")
#

library(bedtoolsr)
options(bedtools.path = "/BioII/lulab_b/baopengfei/biosoft")

dst <- "PRJNA540919_diff" # GSE110381_diff, PRJNA540919_diff
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"


diff <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",dst,"/cfPeakCNN_smallDomain_diff_CRCvsNC.diff"),sep = "\t",header = T,check.names = F,stringsAsFactors = F) # allSmp_
diff$feature <- unlist(sapply(strsplit(rownames(diff),"|",fixed = T),"[",4))
dim(diff)


for(p in c(100,200,500,1000)){ # 0.1,0.05,0.01,0.001
  # p <- 0.01
  #table(diff$log2FoldChange>0)
  #diff2 <- diff[diff$pvalue<=p & diff$log2FoldChange>0,]
  diff2 <- diff[diff$log2FoldChange>0,]
  diff2 <- diff2[order(diff2$pvalue,decreasing = F),]
  peakid <- unlist(sapply(strsplit(rownames(diff2)[1:p],"|",fixed = T),"[",4))
  peak2 <- peak[peak$V4 %in% peakid,]
  write.table(peak2,paste0(pre,"/tmp/",dst,"_CRChigh_p",p,".bed"),quote = F,sep = "\t",row.names = F,col.names = F) 
  
  #diff2 <- diff[diff$pvalue<=p & diff$log2FoldChange<0,]
  diff2 <- diff[diff$log2FoldChange<0,]
  diff2 <- diff2[order(diff2$pvalue,decreasing = F),]
  peakid <- unlist(sapply(strsplit(rownames(diff2)[1:p],"|",fixed = T),"[",4))
  peak2 <- peak[peak$V4 %in% peakid,]
  write.table(peak2,paste0(pre,"/tmp/",dst,"_NChigh_p",p,".bed"),quote = F,sep = "\t",row.names = F,col.names = F) 
}
#bg region bed
diff2 <- diff2[order(diff2$pvalue,decreasing = T),]
#diff2 <- diff2 # rm selected peakid
peakid <- unlist(sapply(strsplit(rownames(diff2)[sample(1:(as.integer(nrow(diff2)*0.8)),min(2000,as.integer(nrow(diff2)*0.5)),replace = F)],"|",fixed = T),"[",4))
peak2 <- peak[peak$V4 %in% peakid,]
write.table(peak2,paste0(pre,"/tmp/",dst,"_bg.bed"),quote = F,sep = "\t",row.names = F,col.names = F) 


{#bash (hub)
  pre="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
  outpre="/BioII/lulab_b/baopengfei/tmp"
    dst="PRJNA540919_diff" #"GSE110381_diff"
  for p in 100 200 500 1000 # 0.1 0.05 0.01 0.001
  do
  for g in CRC NC
  do
  i="$pre/tmp/${dst}_${g}high_p${p}.bed" 
  #i="${pre}/output/$dst/call_peak_all/${method}_11RNA.intersect.bed"
  echo $i
  
  mkdir -p $outpre/${dst}_${p}_${g}/
  cat $i | awk '$3-$2>=15' | grep -v "_chr" | cut -f1-6 > $outpre/${dst}_${p}_${g}/filterLen_11RNA.bed
  
  bedtools getfasta -name -s \
  -fi $pre/exSeek-dev/genome/hg38/fasta_newTxID/combine11.fa \
  -bed $outpre/${dst}_${p}_${g}/filterLen_11RNA.bed > \
  $outpre/${dst}_${p}_${g}/filterLen_11RNA.bed.fa
  
  bedtools getfasta -name -s \
  -fi $pre/exSeek-dev/genome/hg38/fasta_newTxID/combine11.fa \
  -bed $pre/tmp/${dst}_bg.bed > \
  $outpre/${dst}_${p}_${g}/bg.bed.fa
  
  #add bg sequence fa: --n
  sea_dir="$outpre/${dst}_${p}_${g}/sea" # sea, sea_selfBG
  mkdir -p $sea_dir
  /BioII/lulab_b/baopengfei/biosoft/meme_5-4-1/bin/sea --verbosity 1 \
  --p $outpre/${dst}_${p}_${g}/filterLen_11RNA.bed.fa \
  --n $outpre/${dst}_${p}_${g}/bg.bed.fa \
  --m /BioII/lulab_b/baopengfei/shared_reference/meme/motif_databases/RNA/POSTAR3_RBP.meme \
  -oc $sea_dir 
  done
  done
  #  --n $outpre/${dst}_${p}_${g}/bg.bed.fa \
}




#RPB gene universe list
postar <- readLines("/BioII/lulab_b/baopengfei/shared_reference/meme/motif_databases/RNA/POSTAR3_RBP.meme")
postar <- postar[grepl("^MOTIF",postar,perl = T)]
postar <- gsub("MOTIF| |human_eCLIP_MEME_|human_CLIPdb_MEME_|_1|_2|_3","",postar,perl = T)


crc.high <- read.table("/BioII/lulab_b/baopengfei/shared_reference/tmp/crc_up_rbp.txt",header = F)$V1 # annotation from ref


## get TCGA tpm and meta matrix
tumors <- c("LAML","LIHC","KIRC","PRAD","PAAD","COAD")
tcga.expr.list <- list()
for(tumor in tumors){
  # tumor <- "LAML"
  print(tumor)

tcga.expr <- data.table::fread(paste0("/BioII/lulab_b/baopengfei/projects/tcga/data/GDC-",tumor,"/RNA/TCGA-",tumor,".htseq_fpkm.tsv.gz")) # gz sometimes read truncated
print(dim(tcga.expr))
print((tcga.expr[1:3,1:3]))
#max(tcga.expr$`TCGA-DM-A288-01A`) # log2(fpkm+1)
tcga.mat <- 2^tcga.expr[,2:ncol(tcga.expr)]-1
# tcga.mat[1:3,1:3]
tcga.mat2 <- apply(tcga.mat,2,fpkmToTpm)
tcga.expr[,2:ncol(tcga.expr)] <- as.data.frame(tcga.mat2)

ensemb <- unlist(sapply(strsplit(tcga.expr$Ensembl_ID,".",fixed=T), "[", 1))
# tmp <- as.data.frame(mygene::queryMany(ensemb,scopes="ensembl.gene",fields=c("symbol"),species="human"))
# write.table(tmp,"/BioII/lulab_b/baopengfei/projects/tcga/data/GDC-COAD/RNA/ID.txt",quote = F,sep = "\t",col.names = T,row.names = T)
tmp <- read.table("/BioII/lulab_b/baopengfei/projects/tcga/data/GDC-COAD/RNA/ID.txt")
#tmp[1:3,]
tcga.expr$Ensembl_ID <- tmp$symbol[match(ensemb,tmp$query)]
table(duplicated(tcga.expr$Ensembl_ID))
tcga.expr <- tcga.expr[!duplicated(tcga.expr$Ensembl_ID) & !is.na(tcga.expr$Ensembl_ID),]
print(dim(tcga.expr))
tmp2 <- tcga.expr$Ensembl_ID
tcga.expr$Ensembl_ID <- NULL
#max(tcga.expr$`TCGA-DM-A288-01A`) 


# tcga.meta <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/TCGA_small/COAD_meta.txt",sep = "\t",header = T)
# tcga.meta <- readr::read_tsv("/BioII/lulab_b/baopengfei/projects/tcga/data/GDC-LAML/TCGA-LAML.GDC_phenotype.tsv.gz") # COAD, LAML
# table(tcga.meta$sample_type.samples,tcga.meta$leukemia_specimen_cell_source_type)
# table(tcga.meta$sample_type.samples,tcga.meta$site_of_resection_or_biopsy.diagnoses)
# table(tcga.meta$sample_type.samples,tcga.meta$primary_diagnosis.diagnoses)
# table(tcga.meta$sample_type_id.samples)

tcga.meta <- readr::read_tsv(paste0("/BioII/lulab_b/baopengfei/projects/tcga/data/GDC-",tumor,"/TCGA-",tumor,".GDC_phenotype.tsv.gz")) # COAD
if(tumor=="LAML"){
  tcga.meta <- tcga.meta[tcga.meta$sample_type.samples %in% c("Primary Blood Derived Cancer - Peripheral Blood"),] # Solid Tissue Normal,Primary Blood Derived Cancer - Peripheral Blood
}else{
  tcga.meta <- tcga.meta[tcga.meta$sample_type.samples %in% c("Primary Tumor"),]
}
#table(colnames(tcga.expr) %in% tcga.meta$submitter_id.samples)
tcga.expr <- as.data.frame(tcga.expr)

tcga.expr <- tcga.expr[,colnames(tcga.expr) %in% tcga.meta$submitter_id.samples]
tcga.expr <- log2(tcga.expr+1) # logTPM
tcga.expr <- as.data.frame(tcga.expr)
rownames(tcga.expr) <- tmp2
tcga.expr.list[[tumor]] <- tcga.expr
}
#tcga.expr <- NULL


## read GTEx tpm and meta matrix
meta <- data.table::fread("/BioII/lulab_b/baopengfei/projects/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",sep = "\t",header = T,check.names = F,stringsAsFactors = F)
table(meta$SMAFRZE)
meta <- meta[meta$SMAFRZE=="RNASEQ",] # no "Bone Marrow"
table(meta$SMTS)
meta <- meta[meta$SMTS %in% c("Blood","Colon","Small Intestine","Liver","Kidney","Lung","Spleen","Stomach","Prostate","Pancreas"),]
meta$SMTS <- factor(meta$SMTS,levels = c("Blood","Colon","Small Intestine","Liver","Kidney","Lung","Spleen","Stomach","Prostate","Pancreas"))
meta <- meta[order(meta$SMTS),]
expr <- data.table::fread("/BioII/lulab_b/baopengfei/projects/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct",skip = 2,sep = "\t",header = T,check.names = F,stringsAsFactors = F)
max(expr$`GTEX-1117F-0226-SM-5GZZ7`)




# loop for diff cutoff
outpre <- "/BioII/lulab_b/baopengfei/tmp"
run.sea <- function(x){
  # x <- "../WCHSU-FTC/output/i-pico/call_domain_bowtie2_noRepeats_MQ30_localmaxMin2_noQvalue/domains_localmax/domains_bothabundant_expand_sea/sea.tsv"
  print(x)
  sea <- read.table(x,sep = '\t',header = T,check.names = F,stringsAsFactors = F)
  sea$method <- "sea"
  sea$path <- x
  return(sea)
}

for(top in c(200, 500,1000)){ # 0.1,0.05,0.01,0.001;  100,200,500
  # top <- 500 # final version !
  p <- top
  
  rbp <- list()
  for(g in c("CRC","NC")){
    #g <- "NC"
    sea <- run.sea(paste0(outpre,"/",dst,"_",p,"_",g,"/sea/sea.tsv")) # sea, sea_selfBG: little records
    #filter Q-value?
    #filter top5 rank?
    sea$pvalue <- p
    sea$group <- g
    sea$RBP <- gsub("human_eCLIP_MEME_|human_CLIPdb_MEME_|_1|_2|_3","",sea$ID,perl = T)
    rbp[[paste(p,g)]] <- sea
    # table(duplicated(sea$RBP))
    print(length(unique(sea$RBP)))
  }
  df <- as.data.frame(do.call(rbind,rbp))
  length(unique(df$RBP))


#only cal for RBP
expr2 <- as.data.frame(expr[expr$Description %in% df$RBP,])
rownames(expr2) <- expr2$Description
dim(expr2)


#cal mean TPM in all each group of GTEx
smps <- intersect(meta$SAMPID,colnames(expr2))
meta <- meta[meta$SAMPID %in% smps,]
expr2 <- expr2[,meta$SAMPID]
expr2 <- log2(expr2+1) # logTPM
expr2$tx <- rownames(expr2)
expr2 <- reshape2::melt(expr2,id.vars = "tx")
colnames(expr2) <- c("tx","sample","logTPM")
expr2$group <- meta$SMTS[match(expr2$sample,meta$SAMPID)]
table(expr2$group)
expr3 <- as_tibble(expr2) %>% 
  dplyr::group_by(tx,group) %>% 
  dplyr::summarise(mean.logTPM=mean(logTPM,trim=0.05))
expr4 <- reshape2::acast(expr3,tx~group)


#select top
df <- df[order(df$group,df$RANK,decreasing = F),]
top.CRC <- unique(df$RBP[df$group=="CRC"]) #  & df$RANK %in% 1:500
# top.CRC <- top.CRC[(top.CRC %in% rownames(expr4)) & (top.CRC %in% rownames(tcga.expr.list[["COAD"]]))] # 33
top.NC <- unique(df$RBP[df$group=="NC"]) #  & df$RANK %in% 1:500
top.inter <- intersect(top.CRC,top.NC)
top.CRC <- top.CRC[!(top.CRC %in% top.inter) & top.CRC %in% rownames(expr4) & top.CRC %in% rownames(tcga.expr.list[["COAD"]])] # 33
# top.NC <- top.NC[!(top.NC %in% top.inter) & top.NC %in% rownames(expr4) & top.NC %in% rownames(tcga.expr[["COAD"]])] # 32
# # top.CRC <- top.CRC[1:10] # only show top10
# # top.NC <- top.NC[1:10] # only show top10

intersect(top.CRC,crc.high)
intersect(top.NC,crc.high)
# intersect(top.NC,crc.high)




## plot TCGA heatmap
tcga.expr2.df <- list()
for ( tumor in tumors){
  print(tumor)
  
  #only cal for RBP
  expr2 <- as.data.frame(tcga.expr.list[[tumor]][rownames(tcga.expr.list[[tumor]]) %in% df$RBP,])
  expr2$tx <- rownames(expr2)
  expr2 <- reshape2::melt(expr2,id.vars = "tx")
  colnames(expr2) <- c("tx","sample","logTPM")
  expr2$group <- tumor #meta$SMTS[match(expr2$sample,meta$SAMPID)]
  # table(expr2$group)
  expr3 <- as_tibble(expr2) %>% 
    dplyr::group_by(tx,group) %>% 
    dplyr::summarise(mean.logTPM=mean(logTPM,trim=0.05))
  tcga.expr2.df[[tumor]] <- reshape2::acast(expr3,tx~group)
}
tcga.expr2.df <- as.data.frame(do.call(cbind,tcga.expr2.df))


mat <- as.data.frame(tcga.expr2.df)[c(top.CRC),c("LIHC","KIRC","PRAD","PAAD","COAD")] # c("Blood","Colon") top.NC
# annotation_row <- data.frame(row.names = c(top.CRC,top.NC), group=c(rep("CRC plasma",length(top.CRC)),rep("NC plasma",length(top.NC))))
# ann_colors <- list(group=c("CRC plasma"="firebrick","NC plasma"="grey70")) # 
# annotation_row <- annotation_row[annotation_row$group=='CRC plasma',,drop=F] # only show CRC
rbp.order <- stats::hclust(dist(t(scale(t(mat))),method = "euclidean"),method = "ward.D2")  # scale by rows
rbp.order <- rbp.order$labels
# [1] "ZC3H11A"  "WDR3"     "PCBP2"    "DDX59"    "NONO"     "EIF3H"    "DDX51"    "IGF2BP3"  "APOBEC3C" "SFPQ"    
# [11] "PHF6"     "GEMIN5"   "AQR"      "AGGF1"    "DROSHA"   "QKI"
p <- pheatmap::pheatmap(mat = mat[rbp.order,], # [rownames(annotation_row),],
                   # annotation_col = annotation_col_merge, #data frame that specifies the annotations shown on left side of the heatmap.
                   # annotation_row = annotation_row,
                   # annotation_colors = ann_colors, # RColorBrewer::brewer.pal(6, "Set3")
                   border_color="grey20", # grey20
                   scale = "row", # row
                   #labels_col = 3, labels_row = 6,
                   cluster_cols = F,cluster_rows = F, treeheight_row = 20,
                   # gaps_col = 84,
                   gaps_row = length(top.CRC),
                   #cutree_cols = 2,cutree_rows = 3,
                   show_colnames=T, show_rownames=T,
                   fontsize = 12,
                   # height = 12,width =5,  
                   height = 4,width = 4,
                   color = viridis::viridis_pal()(100) #colorRampPalette(c("blue","white","red"))(100),  # steelblue/dodgerblue4-firebrick2
                   #fontsize_row = 5,
                   # filename =paste0("RBP_TCGA_heatmap_5all_",top,".pdf")
                   )
p$gtable$grobs[[3]]$gp=grid::gpar(col="salmon2", fontsize=12)# assuming that the xlabels are in the third grob
ggsave(plot = p, filename = paste0("RBP_TCGA_heatmap_CRC_5Cancer_",top,".pdf"),height = 4,width = 4 )
# #




## plot GTEx heatmap
mat <- as.data.frame(expr4)[c(top.CRC   ),c("Liver","Kidney","Prostate","Pancreas","Colon")] # c("Blood","Colon")
#top.NC
# "Blood","Liver","Spleen","Lung","Kidney","Prostate","Pancreas","Colon","Small Intestine")
# annotation_row <- data.frame(row.names = c(top.CRC,top.NC), group=c(rep("CRC plasma",length(top.CRC)),rep("NC plasma",length(top.NC))))
# ann_colors <- list(group=c("CRC plasma"="firebrick","NC plasma"="grey70")) # 
# annotation_row <- annotation_row[annotation_row$group=='CRC plasma',,drop=F] # only show CRC
p <- pheatmap::pheatmap(mat = mat[rbp.order,], # [rownames(annotation_row),]
                   # annotation_col = annotation_col_merge, #data frame that specifies the annotations shown on left side of the heatmap.
                   # annotation_row = annotation_row,
                   # annotation_colors = ann_colors, # RColorBrewer::brewer.pal(6, "Set3")
                   border_color="grey20", # grey20
                   scale = "row", # row
                   #labels_col = 3, labels_row = 6,
                   cluster_cols = F,cluster_rows = F,treeheight_row = 20,
                   # gaps_col = 84,
                   gaps_row = length(top.CRC),
                   #cutree_cols = 2,cutree_rows = 3,
                   show_colnames=T, show_rownames=T,
                   fontsize = 12,
                   # height = 12,width =5,
                   height = 4,width = 4,
                   color = viridis::viridis_pal()(100) #colorRampPalette(c("blue","white","red"))(100),  # steelblue/dodgerblue4-firebrick2
                   #fontsize_row = 5,
                   # filename = paste0("RBP_GTEx_heatmap_CRC_5Tissue_",top,".pdf")
                   )
p$gtable$grobs[[3]]$gp=grid::gpar(col="salmon2", fontsize=12)# assuming that the xlabels are in the third grob
ggsave(plot = p, filename = paste0("RBP_GTEx_heatmap_CRC_5Tissue_",top,".pdf"), height = 4,width = 4)
#









## ssGSEA
# df <- df[order(df$group,df$RANK,decreasing = F),]
# top.CRC <- unique(df$RBP[df$group=="CRC"]) #  & df$RANK %in% 1:50
# top.NC <- unique(df$RBP[df$group=="NC"]) #  & df$RANK %in% 1:50
# top.inter <- intersect(top.CRC,top.NC) # 32
# # top.CRC <- top.CRC[!(top.CRC %in% top.inter) & top.CRC %in% rownames(expr4) & top.CRC %in% rownames(tcga.expr2)] # 42
# top.CRC <- top.CRC[!(top.CRC %in% top.inter) & top.CRC %in% rownames(expr4) & top.CRC %in% rownames(tcga.expr2[["COAD"]])] # 16
# top.NC <- top.NC[!(top.NC %in% top.inter) & top.NC %in% rownames(expr4) & top.NC %in% rownames(tcga.expr2[["COAD"]])] # 79
# #top.CRC <- top.CRC[1:10]
# #top.NC <- top.NC[1:10]

#make geneset tab
top.list <- list()
top.list[["CRC"]] <- top.CRC
top.list[["NC"]] <- top.NC
method <- "gsva" # ssgsea,gsva,zscore,plage



#GTEx
for(g in c("CRC")){ # ,"NC"
  # g <- "CRC"
  print(g)
  geneset <- data.frame("geneset"="obs","gene"=top.list[[g]])
  #shuf 200 times
  for(seed in 1000:3000){ # 100:600
    set.seed(seed = seed)
    tmp <- data.frame("geneset"=paste0("shuf_",seed),
                      "gene"=postar[sample(1:length(postar),length(top.list[[g]]),replace = F)])
    geneset <- rbind(geneset,tmp)
  }
  dim(geneset) # 42*102
  #geneset <- geneset[!duplicated(geneset),]
  #geneset <- geneset[geneset$gene!="" & geneset$geneset!="",]
  geneset.l <-  base::split(geneset$gene,geneset$geneset)  
  
  table(expr$Description %in% postar)
  exp <- as.data.frame(expr[expr$Description %in% postar,])
  rownames(exp) <- exp$Description # no dup
  if(g=="CRC"){tissue <- "Colon"}else if(g=="NC"){tissue <- "Blood"}
  exp <- exp[,colnames(exp) %in% meta$SAMPID[meta$SMTS==tissue]]
  #table(meta$SMTS)
  if(max(as.matrix(exp))>50){
    exp <- log2(exp+1) # Gaussian need logTPM/CPM
    print("convert to log")
  }
  #max(as.matrix(exp))
  #exp[1:4,1:5]
  
  re <- GSVA::gsva(as.matrix(exp), 
                   geneset.l,
                   method=method,  # method=c("gsva", "ssgsea", "zscore", "plage"),
                   verbose=T,
                   kcdf="Gaussian", # Gaussian: microarray fluorescent units in logarithmic scale, RNA-seq log-CPMs, log-RPKMs or log-TPMs
                   parallel.sz=6)  
  re[1:3,1:3]
  data.table::fwrite(as.data.frame(re),paste0("tmp/GTEx_",method,"_",g,".txt"),quote = F,row.names = T,col.names = T,sep = "\t")
  re <- read.table(paste0("tmp/GTEx_",method,"_",g,".txt"),header = T,sep = "\t", check.names = F,row.names = 1)
  
  re2 <- as.data.frame(t(re))
  re2 <- reshape2::melt(re2)
  re2$group <- ""
  re2$group[grepl("obs",re2$variable)] <- "Observation"
  re2$group[grepl("shuf",re2$variable)] <- "Shuffle"
  re2$group <- factor(re2$group,levels=c("Shuffle","Observation"))
  table(re2$group)
  
  #plot ssGSEA box/violin
  if(g=="CRC"){col <- "#f87669"}else if(g=="NC"){col <- "#2fa1dd"}
  ggplot(re2,aes(x = group, y = value, fill = group)) +
    geom_violin()+
    geom_boxplot(width=0.1)+
    # facet_wrap( ~ sample) +
    # ylim(c(0,1))+
    labs(title = tissue ) + 
    ylab(paste0(method," score"))+
    scale_fill_manual(name="Type",values = c('grey50',col)) + #NC:"#2fa1dd"
    theme_bw() + 
    ggpubr::stat_compare_means(    
      label.x.npc = 0.3,
      label.y.npc = 0.8,
      hide.ns=F,
      size =4,#paired = TRUE,  
      #aes(group = Group,label = p.format), # p.signif,p.format
      label = "p.format",
      # symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05,1),symbols = c( "****","***", "**", "*", "ns")),
      method = "wilcox.test"
    )+
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
      legend.position = "none", #c(0.9,0.8),#,#
      legend.text = element_text(size= 16),
      legend.title= element_text(size= 16))
  ggsave(filename = paste0("GTEx-colon_",g,"_enrich_RBP_",method,"_box_",top,".pdf"), width = 5, height = 5)
}
#


#TCGA_COAD
g <- "CRC"
print(g)
geneset <- data.frame("geneset"="obs","gene"=top.list[[g]])
#shuf 200 times
for(seed in 1000:3000){ # 600
  set.seed(seed = seed)
  tmp <- data.frame("geneset"=paste0("shuf_",seed),
                    "gene"=postar[sample(1:length(postar),length(top.list[[g]]),replace = F)])
  geneset <- rbind(geneset,tmp)
}
dim(geneset) # 42*102
#geneset <- geneset[!duplicated(geneset),]
#geneset <- geneset[geneset$gene!="" & geneset$geneset!="",]
geneset.l <-  base::split(geneset$gene,geneset$geneset)  

tcga.expr.tmp <- tcga.expr.list[["COAD"]]
table(rownames(tcga.expr.tmp) %in% postar)
exp <- as.data.frame(tcga.expr.tmp[rownames(tcga.expr.tmp) %in% postar,])
tissue <- "CRC"
#table(meta$SMTS)
if(max(as.matrix(exp))>50){
  exp <- log2(exp+1) # Gaussian need logTPM/CPM
  print("convert to log")
}
re <- GSVA::gsva(as.matrix(exp), 
                 geneset.l,
                 method=method,  # method=c("gsva", "ssgsea", "zscore", "plage"), 
                 verbose=T,
                 kcdf="Gaussian", # Gaussian: microarray fluorescent units in logarithmic scale, RNA-seq log-CPMs, log-RPKMs or log-TPMs
                 parallel.sz=6)  
re[1:3,1:3]
data.table::fwrite(as.data.frame(re),paste0("tmp/TCGA_COAD_",method,"_",g,".txt"),quote = F,row.names = T,col.names = T,sep = "\t")
re <- read.table(paste0("tmp/TCGA_COAD_",method,"_",g,".txt"),header = T,sep = "\t", check.names = F,row.names = 1)

re2 <- as.data.frame(t(re))
re2 <- reshape2::melt(re2)
re2$group <- ""
re2$group[grepl("obs",re2$variable)] <- "Observation"
re2$group[grepl("shuf",re2$variable)] <- "Shuffle"
re2$group <- factor(re2$group,levels=c("Shuffle","Observation"))
table(re2$group)

#plot ssGSEA box/violin
if(g=="CRC"){col <- "#f87669"}else if(g=="NC"){col <- "#2fa1dd"}
ggplot(re2,aes(x = group, y = value, fill = group)) +
  geom_violin()+
  geom_boxplot(width=0.1)+
  # facet_wrap( ~ sample) +
  # ylim(c(0,1))+
  labs(title = tissue ) + 
  ylab(paste0(method," score"))+
  scale_fill_manual(name="Type",values = c('grey50',col)) + #NC:"#2fa1dd"
  theme_bw() + 
  ggpubr::stat_compare_means(    
    label.x.npc = 0.3,
    label.y.npc = 0.8,
    hide.ns=F,
    size =4,#paired = TRUE,  
    #aes(group = Group,label = p.format), # p.signif,p.format
    label = "p.format",
    method = "wilcox.test"
  )+
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
    legend.position = "none", #c(0.9,0.8),#,#
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16))
ggsave(filename = paste0("TCGA_COAD",g,"_enrich_RBP_",method,"_box_",top,".pdf"), width = 5, height = 5)


}

