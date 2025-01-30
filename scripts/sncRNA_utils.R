# sncRNA utils

runDcbPeakDiff <- function(count,sample.table,batch="NULL",pos.grp="Y",neg.grp="N",method="edger_glmlrt",norm_method="TMM",libSizePath="FALSE",outDir,outSurfix=""){
  # pos.grp <- "Y"
  # neg.grp <- "N"
  # libSizePath <- c("s1"=220000,"s2"=2010100) # named vec or tsv path
  # count <- read.table( count,header = T,sep = "\t",stringsAsFactors = F,check.names = F)

  # count = count
  # sample.table = sample.table
  # libSizePath = libSizePath
  # pos.grp = "Y"
  # neg.grp = "N"
  # method="edger_glmlrt"
  # norm_method="TMM"
  # outDir = outDir
  # outSurfix=""

  
  # if(file.exists(count)){
  if(!(nrow(count)>2 & ncol(count)>2)){
    message("read mat from path")
    count <- data.table::fread( count,header = T,sep = "\t",stringsAsFactors = F,check.names = F,data.table = F)
    #count <- read.table( paste0("/data2/lulab1/bpf/projects/WCHSU-FTC/output/",dst,"/call_peak_dedup/count_matrix/cfpeakCNN_b5_d50_p1_small.txt"),header = T,sep = "\t",stringsAsFactors = F,check.names = F)
    
    if("feature" %in% colnames(count)){
      rownames(count) <- count$feature
      count$feature <- NULL
    }else if("gene_id" %in% colnames(count)){
      rownames(count) <- count$gene_id
      count$gene_id <- NULL
    }

  }

  
  
  ### Meta vs. Local
  # colnames(sample.table0)
  # sample.table <- sample.table0
  sample.table <- sample.table[sample.table$sample %in% colnames(count),]
  sample.table <- sample.table[order(sample.table$group),]
  # table(sample.table$group)
  
  positive_samples <- sample.table[sample.table$group==pos.grp,"sample"]  # 
  negative_samples <- sample.table[sample.table$group==neg.grp,"sample"] #
  samples <- c(positive_samples, negative_samples)
  sample.table <- sample.table[samples,]
  mat <- count[,samples]
  
  ## read libSize
  if(libSizePath!="FALSE" & file.exists(libSizePath)){
    message("read libSize from path")
    tmp <- data.table::fread( libSizePath, header = F, sep = "\t", check.names = F, data.table = F )
    libSize <- tmp$V2
    names(libSize) <- tmp$V1
    libSize <- libSize[samples]
  }else if(libSizePath!="FALSE" & exists("libSizePath") & nrow(libSizePath)>1 & ncol(libSizePath)>1){
    message("use existing libSize var")
    libSize <- libSizePath[samples]
  }else{
    libSize <- "FALSE"
  }
  # message(libSize)
  
  group <- c(rep("positive",length(positive_samples)),rep("negative",length(negative_samples)))
  # method <- "edger_glmlrt"
  # norm_method <- "TMM"
  
  #res <- diff.v2(mat = mat, samples = samples, group = group, method = method, norm_method = norm_method, filterType = "NULL", featureType = "domain")
  if(batch=="NULL"){
    message("no batch")
    res <- diff.v2.dcb(mat = mat, samples = samples, group = group, 
                       method = method, norm_method = norm_method, filterType = "NULL", featureType = "domain", 
                       getCountGt1SmpFrac = TRUE, getPercExpByGrp = TRUE, useExternalLibSize = libSize)
  }else{
    message("has batch")
    res <- diff.v2.dcb(mat = mat, samples = samples, group = group, batch = as.factor(sample.table$batch), 
                       method = method, norm_method = norm_method, filterType = "NULL", featureType = "domain", 
                       getCountGt1SmpFrac = TRUE, getPercExpByGrp = TRUE, useExternalLibSize = libSize)
  }
  
  dir.create(outDir,recursive = T,showWarnings = F)
  # write.table(res[["normMat"]],paste0(outDir,"/cfPeakCNN_smallDomain_diff_",pos.grp,"vs",neg.grp,"",outSurfix,".cpm"),quote = F,sep = "\t",row.names = T,col.names = T) # 
  write.table(res[["diffTable"]],paste0(outDir,"/cfPeakCNN_smallDomain_diff_",pos.grp,"vs",neg.grp,"",outSurfix,".diff"),quote = F,sep = "\t",row.names = T,col.names = T) # 
}





# #TODO
# plotUMAP <- function(sample.table, logcpm ){
# }


# plot PCA
plotPCA <- function(sample.table, logcpm, topVarNum=1000, outFile, plotGrpColor, plotDotSize = 2, plotDotAlpha = 1,legend.position='none',figWid=15, figHigh=12, plotUMAP="NULL"){ # ,  plotGrpFill,plotGrpShape
  # topVarNum <- 1000
  # plotUMAPpcaNum <- 4
  message(nrow(logcpm))
  dir.create(dirname(outFile),recursive = T,showWarnings = F)
  # library(data.table)
  library(ggplot2)
  my_theme <-   theme(aspect.ratio=1,
                      plot.title = element_text(size = 28,color="black",hjust = 0.5),
                      axis.title = element_text(size = 28,color ="black"),
                      axis.text = element_text(size= 24,color = "black"), #,face="bold
                      panel.grid.minor.y = element_blank(),
                      panel.grid.minor.x = element_blank(),
                      axis.text.x = element_text( hjust = 0.5 ), # angle = 45,
                      axis.text.y = element_text( hjust = 0 ), # angle = 45,
                      panel.grid=element_blank(),
                      legend.position = legend.position,#c(0.5,0.3),
                      legend.text = element_text(size= 26,color = "black"),
                      legend.title= element_text(size= 26,color = "black"))
  
  logcpm <- logcpm[,sample.table$sample] 
  # sample.table <- sample.table[colnames(logcpm),]
  
  tmp.sd <- matrixStats::rowSds(as.matrix(logcpm))
  names(tmp.sd) <- rownames(logcpm)
  tmp.sd.sort <- sort(tmp.sd,decreasing = T)
  #head(tmp.sd.sort)
  
  # which(tmp.sd>0)[1:min(topVarNum,sum(tmp.sd>0))]
  logcpm <- logcpm[tmp.sd!=0,] # rm 0 std, or meet error for pca
  
  # table(duplicated(tmp.sd.sort))
  # 
  tmp.idx <- tmp.sd.sort[1:min(length(tmp.sd.sort),topVarNum)]
  tmp.idx <- tmp.idx[tmp.idx>0]
  logcpm <- logcpm[names(tmp.idx),] # use top 1k var features
  
  pca <- stats::prcomp(t(logcpm), scale=TRUE)  ###prcomp函数的横行必须是样本，所以倒置一下
  # dim(logcpm)
  pca <- summary(pca) # append sum info (importance et al.)
  # all(tmp$sdev==pca$sdev)
  # all(tmp$x==pca$x)
  # all(tmp$rotation==pca$rotation)
  #tmp <- tmp$importance
  #length(pca.var.per)
  #cor.test(pca.var.per,tmp[1,],method = "spearman") # 0.99
  #cor.test(pca$importance[1,],pca$importance[2,],method = "spearman") # 0.99
  
  # #choose top2 PC
  pca.var <- pca$sdev^2  ## sdev是标准偏差，十个样本，就有十个标准偏差，平方是避免负数的干扰
  pca.var.per <- round(pca.var/sum(pca.var)*ncol(logcpm), 1)  ##求每个样本的variation
  ElbowPoint <- PCAtools::findElbowPoint(pca.var.per)
  message("ElbowPoint: ",ElbowPoint)
  plotPCAnum <- max(7,min(50,ElbowPoint)) #limit to [7,50]
  exportPCAnum <- max(10,plotPCAnum) # min export PC1-10 table
  #barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")  ##用柱状图可视化
  #all(sample.table$sample==rownames(pca$x)) #T
  
  if(plotUMAP!="NULL"){
    #combine multiple PC to 2D plot
    library(umap) # or use uwot
    umap.lab <- sample.table$Color
    umap.df <- pca$x[,1:exportPCAnum]
    
    umap.tmp <- umap::umap(umap.df[,1:plotPCAnum])
    sample.table[[paste0("Dim1")]] = umap.tmp$layout[,1]
    sample.table[[paste0("Dim2")]] = umap.tmp$layout[,2]
    
    p.tmp <- ggplot(data=sample.table,aes(x=Dim1,y=Dim2)) +
      geom_point(size=plotDotSize,alpha=plotDotAlpha, aes(color=Color)) + # ,fill=Fill,shape=Shape
      labs(title = paste0("UMAP.topPC",plotPCAnum)) +
      xlab(paste("Dim1")) +
      ylab(paste("Dim2")) +
      theme_bw()+
      scale_color_manual(name="Color", values = plotGrpColor)+
      # scale_fill_manual(name="Fill", values = plotGrpFill)+
      # scale_shape_discrete(name="Shape", values = plotGrpShape )+
      my_theme
    ggsave(filename = paste0(outFile,".UMAP.png"),plot = p.tmp,width = figWid,height = figHigh)
    data.table::fwrite(data.frame("sample"=rownames(sample.table),"UMAP_Dim1"=sample.table$Dim1,"UMAP_Dim2"=sample.table$Dim2),paste0(outFile,".UMAP.txt"),quote = F,sep = "\t",row.names = F,col.names = T)
    #umap::plot.iris(umap.tmp, umap.lab)
  }
  
  pca.data <- data.frame(Sample=rownames(pca$x),
                         # Fill=sample.table$Fill,
                         Color=sample.table$Color,
                         # Shape=sample.table$Shape,
                         X=pca$x[,1],
                         Y=pca$x[,2],
                         Z=pca$x[,3],
                         A=pca$x[,4])
  
  pca.plot1 <-
    ggplot(data=pca.data,aes(x=Y,y=X))+
    geom_point(size=plotDotSize,alpha=plotDotAlpha, aes(color=Color)) + # ,fill=Fill,shape=Shape
    xlab(paste("")) +
    ylab(paste("PC1(",round(pca$importance[2,1]*100,digits = 2),"%","Var)",sep="")) +
    theme_bw()+
    scale_color_manual(name="Color", values = plotGrpColor)+
    # scale_fill_manual(name="Fill", values = plotGrpFill)+
    # scale_shape_discrete(name="Shape", values = plotGrpShape )+
    my_theme
  pca.plot2 <-
    ggplot(data=pca.data,aes(x=A,y=X))+
    geom_point(size=plotDotSize,alpha=plotDotAlpha, aes(color=Color)) + # ,fill=Fill,shape=Shape
    xlab(paste("")) +
    ylab(paste("")) +
    theme_bw()+
    scale_color_manual(name="Color", values = plotGrpColor)+
    # scale_fill_manual(name="Fill", values = plotGrpFill)+
    # scale_shape_discrete(name="Shape", values = plotGrpShape )+
    my_theme
  pca.plot3 <-
    ggplot(data=pca.data,aes(x=Y,y=Z))+
    geom_point(size=plotDotSize,alpha=plotDotAlpha, aes(color=Color)) + # ,fill=Fill,shape=Shape
    xlab(paste("PC2(",round(pca$importance[2,2]*100,digits = 2),"%","Var)",sep="")) +
    ylab(paste("PC3(",round(pca$importance[2,3]*100,digits = 2),"%","Var)",sep="")) +
    theme_bw()+
    scale_color_manual(name="Color", values = plotGrpColor)+
    # scale_fill_manual(name="Fill", values = plotGrpFill)+
    # scale_shape_discrete(name="Shape", values = plotGrpShape )+
    my_theme
  pca.plot4 <-
    ggplot(data=pca.data,aes(x=A,y=Z))+
    geom_point(size=plotDotSize,alpha=plotDotAlpha, aes(color=Color)) + # ,fill=Fill,shape=Shape
    xlab(paste("PC4(",round(pca$importance[2,4]*100,digits = 2),"%","Var)",sep="")) +
    ylab(paste("")) +
    theme_bw()+
    scale_color_manual(name="Color", values = plotGrpColor)+
    # scale_fill_manual(name="Fill", values = plotGrpFill)+
    # scale_shape_discrete(name="Shape", values = plotGrpShape )+
    my_theme
  pca.plot <- ggpubr::ggarrange(align = "hv",legend = "right",common.legend = TRUE,ncol = 2,nrow = 2,plotlist = list(pca.plot1,pca.plot2,pca.plot3,pca.plot4))
  # ggsave(filename = paste0("./output/",dst,"/allSmp_top100_cfPeakCNN_smallDomain_diff_",disease.label,"vs",normal.label,".cpm.PCA.pdf"),plot = pca.plot,width = 15,height = 12)
  ggsave(filename = outFile,plot = pca.plot,width = figWid,height = figHigh)
  
  # data.table::fwrite(as.data.frame(pca$rotation),paste0(outFile,".PCArotation.txt"),quote = F,sep = "\t",row.names = T,col.names = T)
  tmp <- as.data.frame(pca$importance[,1:exportPCAnum]) #round(digits=1) # 2nd row: importance
  data.table::fwrite(as.data.frame(cbind('sample'=names(tmp),t(tmp))),paste0(outFile,".PCAimportance.txt"),quote = F,sep = "\t",row.names = F,col.names = T)
  
  
  for (i in 1:exportPCAnum){
    pca.data[[paste0("PC",i)]] <- pca$x[,i]
  }
  colnames(pca.data)[colnames(pca.data)=="Sample"] <- "sample"
  data.table::fwrite(pca.data[,c("sample",paste0("PC",1:exportPCAnum))],paste0(outFile,".PCAdata.txt"),quote = F,sep = "\t",row.names = F,col.names = T)
}




# # plot heatmap
# p <- pheatmap::pheatmap(mat = mat, # [rownames(annotation_row),]
#                         # annotation_col = annotation_col_merge, #data frame that specifies the annotations shown on left side of the heatmap.
#                         # annotation_row = annotation_row,
#                         # annotation_colors = ann_colors, # RColorBrewer::brewer.pal(6, "Set3")
#                         border_color="grey20", # grey20
#                         scale = "row", # row
#                         #labels_col = 3, labels_row = 6,
#                         cluster_cols = F,cluster_rows = T,treeheight_row = 20,
#                         # gaps_col = 84,
#                         gaps_row = length(top.CRC),
#                         #cutree_cols = 2,cutree_rows = 3,
#                         show_colnames=T, show_rownames=T,
#                         fontsize = 12,
#                         # height = 12,width =5,
#                         height = 4,width = 4,
#                         color = viridis::viridis_pal()(100) #colorRampPalette(c("blue","white","red"))(100),  # steelblue/dodgerblue4-firebrick2
#                         #fontsize_row = 5,
#                         # filename = paste0("RBP_GTEx_heatmap_CRC_5Tissue_",top,".pdf")
# )
# 
# p$gtable$grobs[[3]]$gp=grid::gpar(col="salmon", fontsize=20)# assuming that the xlabels are in the third grob
# p$gtable$grobs[[4]]$gp=grid::gpar(col="salmon", fontsize=20)# assuming that the ylabels are in the fourth grob
# p$gtable$grobs[[1]]$gp=grid::gpar(col="salmon", lwd=2) # change the color of the dendrogram and set the linewidth to 2
# p$gtable$grobs[[5]]$gp=grid::gpar(col="salmon", fontsize="20", just="center") 
# p

runPCAtools <- function(mat,pheno.table,topN=1000, removeVar=0.1, col=rev( colorRampPalette( RColorBrewer::brewer.pal(n = 10, name = "RdBu"))(100) ), plotRsquared=F, scale=T, corFUN='spearman',main="log2CPM matrix",outFile="./pcatools-eigencorplot.pdf",width=7,height=7){
  # #test
  # mat = mat
  # pheno.table = pheno.table
  # topN=1000
  # removeVar=0.1
  # col <- rev( colorRampPalette( RColorBrewer::brewer.pal(n = 10, name = "RdBu"))(100) )
  # #col=c("#053061","#EDDED5","#67001F")
  # plotRsquared=F
  # scale=T
  # corFUN='spearman'
  # main="log2CPM matrix"
  # outFile="./pcatools-eigencorplot.pdf"
  # width=7
  # height=7

  set.seed(1234)
  library(tidyverse)
  library(PCAtools)
  # pheno.table <- sample.table0
  # colnames(pheno.table)
  # pheno.table <- pheno.table[,(colnames(pheno.table) %in% varColumns)]
  pheno.table[pheno.table==""] <- "-"
  pheno.table[pheno.table=="NO_INFO"] <- "-"
  pheno.table[is.na(pheno.table)] <- "-"
  #pheno.table[pheno.table$group=="NC",c("cfMeDIP_kit","WGS_Kit")] <- "QIAGEN"
  pheno.table <- pheno.table %>% mutate(across(where(is.character), factor))
  
  vst <- mat
  s <- apply(vst,1,sd)
  #table(s>0)
  names(s) <- rownames(vst)
  # s <- s[order(s,decreasing=T)]  # [1:200]
  #s <- t(scale(t(s)))
  s1 <- s[order(s,decreasing=T)]
  # if(length(s1))
  #message(length(s1))
  top <- names(s1)[1:max(length(s1),topN)]
  x <- vst[top,]
  x <- x[,rownames(pheno.table)]
  all(colnames(x) %in% rownames(pheno.table))
  all(colnames(x) == rownames(pheno.table))
  # dim(x)
  
  p <- pca(x, metadata = pheno.table, removeVar = removeVar)
  #screeplot(p, axisLabSize = 18, titleLabSize = 22)
  # p$metadata$lineage
  # biplot(p,shape = "group",legendPosition = "right", colby = "group") # , lab = NULL
  #pairsplot(p,colby = "seq_date")
  #plotloadings(p, labSize = 3)
  pdf(outFile,width = width,height = height) #
  p1 <- eigencorplot(p, plotRsquared = plotRsquared, col = col,
                     main = main,
                     signifSymbols = c('****', '***', '**', '*', ''),
                     scale = scale,
                     rotLabX = 0,
                     signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                     corFUN = corFUN, # spearman, pearson
                     metavars = colnames(pheno.table)  # colnames(pheno.table)[!(colnames(pheno.table) %in% c("cluster","group"))
  )
  print(p1)
  # ggsave(outFile,width = width,height = height)
  dev.off()
}


getSurvPval <- function(exprVec,sample.table,time.col="OS_time",event.col="OS"){ #
  library(survival)
  library(survminer)
  #name and value in same order
  #sample.table need has OS, OS_time et al.
  # mat2 <- mat[tmp.id == pID,]
  # exprVec <- mat2[1,]
  # time.col="OS_time"
  # event.col="OS"
  # getSurvPval(exprVec = as.numeric(mat[tmp.id == pID,]), sample.table = sample.table.pair, time.col = "OS_time", event.col = "OS")
  # exprVec <- 0
  # tryCatch(
  #   {
      if(is.numeric(exprVec[1])){
        exprVec.cut <- cut(as.numeric(exprVec), labels = c("low","high"),breaks =  as.numeric(quantile(exprVec, probs = c(0,0.5,1))), include.lowest = T) # need handle zero expr
      }else if(is.character(exprVec[1]) | is.factor(exprVec[1])){
        exprVec.cut <- as.factor(exprVec)
      }
    # if(exists("exprVec.cut") ){ # | is.nan(exprVec.cut[1])
    #   break(NaN)
    # }
    # summary(mat2)
    sample.table.tmp <- as_tibble(sample.table) %>% 
      dplyr::mutate("Group"=exprVec.cut, "event.col"=sample.table[[event.col]], "time.col"=sample.table[[time.col]]) %>% 
      # clean_names() %>%
      # mutate(churn = ifelse(churn == "Yes", 1, 0)) %>%
      dplyr::mutate_if(is.character, as_factor) %>% 
      as.data.frame # need convert to df, tbl will meet error
    # time.col <- time.col
    # event.col <- event.col
    # message(time.col)
    sfit_tmp <- survminer::surv_fit(Surv(time.col, event.col) ~ Group, data = sample.table.tmp) # why time.col only take from envir ? https://github.com/kassambara/survminer/issues/342
    # ?survfit
    ## summary(sfit_tmp)
    ## eval(parse(text = "x"))
    ## get(sample.table.tmp[["OS_time"]])
    res.suvr <- (survminer::surv_pvalue(sfit_tmp)$pval)
    
  #   },
  #   error = function(cond) {
  #     message("too many zero expr? skip")
  #     res.suvr <- NaN
  #     # exprVec.cut=NaN
  #   }
  #   
  # )
  return(res.suvr)
  
}


plotSurv <- function(exprVec,sample.table,time.col="OS_time",event.col="OS",palette="aaas",facet="NULL",outFile="surv.pdf", width = 7, height = 7){
  library(survival)
  library(survminer)
  #x <- "peak_21548"
  #x <- ENST00000364469.1_0_35_+|Y_RNA|ENST00000364469.1|peak_449|ENST00000364469.1|0|35
  # outFile <- paste0("surv",x,".pdf")
  # exprVec <- mat[tmp.id == x,]
  # time.col="OS_time"
  # event.col="OS"
  # palette="aaas"
  # facet="NULL"
  # outFile="surv.pdf"
  # 
  # print(x)
  # exprVec <- 1:4
  # exprVec <- c("a","b")
  if(is.numeric(exprVec[1])){
    exprVec.cut <- cut(as.numeric(exprVec), labels = c("low","high"),breaks =  as.numeric(quantile(exprVec, probs = c(0,0.5,1))), include.lowest = T) # need handle zero expr
  }else if(is.character(exprVec[1]) | is.factor(exprVec[1])){
    exprVec.cut <- as.factor(exprVec)
  }
  # if(exists("exprVec.cut") ){ # | is.nan(exprVec.cut[1])
  #   break(NaN)
  # }
  # summary(mat2)
  sample.table.tmp <- as_tibble(sample.table) %>% 
    dplyr::mutate("Group"=exprVec.cut, "event.col"=sample.table[[event.col]], "time.col"=sample.table[[time.col]]) %>% 
    # clean_names() %>%
    # mutate(churn = ifelse(churn == "Yes", 1, 0)) %>%
    dplyr::mutate_if(is.character, as_factor) %>% 
    as.data.frame # need convert to df, tbl will meet error
  # time.col <- time.col
  # event.col <- event.col
  # message(time.col)
  sfit_tmp <- survminer::surv_fit(Surv(time.col, event.col) ~ Group, data = sample.table.tmp) # why time.col only take from envir ? https://github.com/kassambara/survminer/issues/342
  # ?survfit
  # summary(sfit_2)
  # surv_pvalue(sfit_2)
  # survminer::surv_summary(sfit_2)
  tmp.theme <- theme_classic() +
    theme(plot.title = element_text(size = 16,color="black",hjust = 0.5,face="bold"),
          axis.title = element_text(size = 24,color ="black"),
          axis.text = element_text(size= 20,color = "black"), # ,face="bold"
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5),
          axis.text.y = element_text(angle = 0, hjust = 0.9, vjust = 0.5),
          panel.grid=element_blank(),
          legend.position = "right",
          axis.line = element_line(colour = "black"),
          legend.text = element_text(size= 24),
          legend.title= element_text(size= 24),
          strip.text = element_text(size = 20),
          strip.background = element_blank(),
          panel.spacing = unit(1.5, "lines"),
          panel.border = element_rect(colour = "black", fill=NA, size=2))
  if(facet=="NULL"){
    p <- ggsurvplot(
      sfit_tmp,
      # color = c("#246BAE","#B51F2E"), #col2(10),  "#67001F" "#B51F2E" "#DC6F58" "#F7B698" "#FDEBDF" "#E5F0F6" "#A7CFE4" "#549EC9" "#246BAE" "#053061"
      palette=palette,
      conf.int = F,pval = T,
      data = sample.table.tmp
    ) + xlab("Survival time (mon.)")+
      ylab(event.col)+
      labs(title = "") #+ tmp.theme
    # print(p)
  }else if(facet %in% colnames(sample.table.tmp)){
    p <- ggsurvplot_facet( # ggsurvplot_facet
      sfit_tmp,
      palette=palette,
      conf.int = F, #pval.coord = c(0.1,0.2), pval.method=T,
      data  = sample.table.tmp,
      panel.labs.font = list(face = NULL, color = NULL, size = 28, angle = NULL),
      facet.by = c(facet), # race, tumor_stage
      palette = "jco", pval = TRUE #, nrow = 1
    ) + xlab("Survival time (mon.)")+
      ylab(event.col)+
      labs(title = "") # + tmp.theme
    #   theme_tq() +
    #   # scale_fill_tq() +
    #   # scale_color_tq() +
  }else(
    stop("facet var not found in smp tbl")
  )
  # ggsave(plot = p, filename = outFile)
  pdf(file = outFile, width = width, height = width)
  print(p)
  dev.off()
}



