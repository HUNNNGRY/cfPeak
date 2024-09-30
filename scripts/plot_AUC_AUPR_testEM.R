# test for gnStarEM
# last 2406 by b.p.f@qq.com

library(tidyverse)
library(dplyr)
library(ggpubr)

source("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/scripts/util.R")
#
library(pROC)  # For AUROC
library(PRROC) # For AUPR

calculate_auroc_aupr <- function(TP, TN, FP, FN) {
  # Calculate TPR, FPR, Precision, and Recall
  TPR <- TP / (TP + FN)          # True Positive Rate (Recall)
  FPR <- FP / (FP + TN)          # False Positive Rate
  Precision <- TP / (TP + FP)    # Precision
  
  # Create the input for ROC and PR curves
  labels <- c(rep(1, length(TP)), rep(0, length(FP))) # Positive and negative labels
  scores <- c(TP / (TP + FN), FP / (FP + TN))         # Scores based on TPR and FPR
  
  # AUROC calculation
  roc_curve <- roc(labels, scores)
  auroc <- auc(roc_curve)
  
  # AUPR calculation
  pr_curve <- pr.curve(scores.class0 = scores[labels == 1], scores.class1 = scores[labels == 0], curve = TRUE)
  aupr <- pr_curve$auc.integral
  
  return(list(AUROC = auroc, AUPR = aupr))
}

# # Example usage
# # TP, TN, FP, FN vectors under different cutoffs
# TP <- c(100, 90, 80)
# TN <- c(90, 80, 70)
# FP <- c(10, 20, 30)
# FN <- c(20, 30, 40)
# 
# results <- calculate_auroc_aupr(TP, TN, FP, FN)
# print(results)

#

# calculate confusion matrix (deprecated, test-only) -----------------------------------------------------------------------
#setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/") 
#
#mypal <- pal_d3_adaptive("category10", alpha = 0.9)(9)
black.col <- alpha("grey90",0.3)
methCols <- c("Piranha"="#1F77B4E5","CLIPper"="#FF7F0EE5","CLAM"="#2CA02CE5","cfPeak"="red","cfPeak.CNN"="firebrick","blockbuster"="#9467BDE5","blockbuster.k100"="#9C564BE5","blockbuster.EM"="#E377C2E5")


MAP="tx" # gn tx
cfpeakMAP="tx" # gn tx
TYPE="all4" # all  all2 all3 all4  called
message(MAP," mode")
dir.create(paste0("./figure/",MAP,"_",cfpeakMAP,"_",TYPE))

#dsts <- c("TCGA-COAD_small_NC_NCpool") # "GSE71008","GSE94533", "GSE110381", "WSQ_SMARTer_NEB","GSE94582" ,"GSE50676"
dsts <- c("CNP0003091_urine_NCpool","GSE129255","GSE112343_NCpool","GSE56686", "snc_pandora_hsa_HeLa_NCpool", "encode_small_colon_NCpool")
# dsts <- c("GSE71008_NCpool","GSE94533_NCpool","GSE123972_NCpool", "GSE94582_NCpool","WSQ_SMARTer_NEB_NCpool","GSE148861_GSE148862_NCpool", "Phospho-RNA-seq_NCpool", "GSE110381_NCpool", "AGO2_IP_NCpool")
#sample_ids <- c("NCpool", "TCGA-4T-AA8H-01A-11H-A41D-13_mirna_gdc_realn") # "NCpool", # TCGA-4T-AA8H-01A-11H-A41D-13_mirna_gdc_realn
#"GSE50676_NCpool"

if(MAP == "tx"){
  methods <- c("Piranha","CLIPper","CLAM","bowtie2_cfpeak") #StarEM_win50,Star_k100,bowtie2_cfpeak,,Piranha,CLIPper,CLAM  ,"bowtie2_cfpeakCNN" "StarEM_win50"
  methods.labs <- c("Piranha","CLIPper","CLAM","cfPeak") #blockbuster.k100 cfPeak.CNN
  #StarEM_win2000
  # exclud.meth <- c("blockbuster.k100","blockbuster.EM") # "blockbuster.k100" ,"cfPeak.CNN"
  # methods.labs <- methods.labs[!methods.labs %in% exclud.meth]
}else if(MAP == "gn"){
  #sample_ids <- c("NCpool") #c("SRR2105340","SRR2105341","SRR2105342","SRR2105346","SRR2105335")  #c("SRR5230634","SRR5230661","SRR5230662","SRR5230667") # c("SRR2105336","SRR2105128","SRR2105129","SRR2105335","SRR2105336"), c("NC_PKU-2392860_smart_1","NC_PKU-2392860_smart_PNK_1","NC_PKU-2392860_1"),c(CleanTag_Lab5,NEBNext_Lab1,NEBNext_Lab2,TruSeq_Lab1,TruSeq_Lab2,N4_A_Lab5,N4_B_Lab2),c("GM12878","GM12891")
  # methods <- c("piranha_by_sample/b5_p01","clipper_by_sample/b5_p05","clam_by_sample/b5_p005","expeakCNN_by_sample/b5_d50_p1") # ,"domains_localmax_significant/b5_d05_p05"
  #methods <- c("StarEM_win50","Star_k100","Star_default","bowtie2_cfpeak") # ,"domains_localmax_significant/b5_d05_p05"
  methods <- c("Star_default","Piranha","CLIPper","CLAM","bowtie2_cfpeak") #StarEM_win50,Star_k100,bowtie2_cfpeak,,Piranha,CLIPper,CLAM  ,"bowtie2_cfpeakCNN" "StarEM_win50",
  methods.labs <- c("blockbuster","Piranha","CLIPper","CLAM","cfPeak") #blockbuster.k100
  #StarEM_win2000
  # exclud.meth <- c("blockbuster.k100","blockbuster.EM","cfPeak.CNN") # "blockbuster.k100"
  # methods.labs <- methods.labs[!methods.labs %in% exclud.meth]
}

for (dst in dsts){
  # dst <- "test"
  # print(dst)
  # dst <- "AGO2_IP_NCpool"
  # tmp <- tmp_df[tmp_df$dst==dst,]
  # sample_ids <- read.table(paste0("data/",dst,"/sample_ids_NC_test.txt"))$V1
  sample_ids <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/",dst,"/sample_ids_raw.txt"))$V1
  #sample_ids <- "TCGA-A6-2685-11A-01R-1757-13_mirna_gdc_realn"
  
  for (sample_id in sample_ids){
    # sample_id <- "CRC"
    # print(sample_id)
    
    ## get matrix in R
    res <- list()
    # for (frac in c("")){ # "0.1","0.5"
      frac <- ""
    #   for (top in c("")){ # "500","1000","2000","5000","10564"; '9426', '7077', '4893','3419'
        top <- ""

        for (method in methods){
          # method <- "Star_default"
          if(grepl("cfpeak",tolower(method))){
            message("cfpeak")
            MAP2 <- cfpeakMAP
          }else{
            message("other")
            MAP2 <- MAP
          }
          
          inputFile <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/test/",dst,"_",method,"_",sample_id,"_",MAP2,"_",TYPE,".txt") # OPT1: .txt6   OPT2:.txt7
          method.name <- method
          
          tmp <- read.table(inputFile,header = T,sep = "")
          tmp$specificity <- 1-tmp$fpr
          tmp$replicate <- ""
          colnames(tmp)[colnames(tmp)=="replicate"] <- "sample"
          colnames(tmp)[colnames(tmp)=="counts"] <- "cutoff"
          tmp$method <- method.name
          tmp$dst <- dst
          tmp$top <- top
          tmp$frac <- frac
          tmp$smp <- sample_id
          
          res[[paste0(dst,sample_id,method.name,frac,top)]] <- tmp
        }


    res.df <- as.data.frame(do.call(rbind,res))
    table(res.df$dst, res.df$method)
    res.df$method[res.df$method=="Star_default"] <- "blockbuster"
    res.df$method[res.df$method=="Star_k100"] <- "blockbuster.k100"
    res.df$method[res.df$method=="StarEM_win50"] <- "blockbuster.EM"
    res.df$method[res.df$method=="bowtie2_cfpeak"] <- "cfPeak"
    res.df$method[res.df$method=="bowtie2_cfpeakCNN"] <- "cfPeak.CNN"
    table(res.df$method)
    #res.df[1:3]
    # summary(res.df$fpr)
    # omit any record that not in [0,1]
    res.df <- res.df[res.df$precision>=0 & res.df$precision<=1 & res.df$recall>=0 & res.df$recall<=1 & res.df$fpr>=0 & res.df$fpr<=1 & res.df$f1>=0 & res.df$f1<=1 & res.df$specificity>=0 & res.df$specificity<=1 ,]
    res.df <- res.df[res.df$method %in% methods.labs,]
    
    ## plot standard-overlap/recall/precison pie/bar plot 
    res.df2 <- res.df
    res.df2 <- dplyr::as_tibble(res.df2) %>% 
      dplyr::group_by(dst,method,sample,smp) %>%
      # dplyr::summarise() %>% 
      dplyr::arrange(cutoff) %>% # only select min cutoff (all peaks)
      dplyr::distinct(method, sample, dst, .keep_all = TRUE) %>% 
      tidyr::pivot_longer(cols = precision:specificity, names_to = "group",values_to = "number") %>%
      dplyr::mutate(total.peak.num=TP+FP+TN+FN, no.number=1-number) %>% # TN+FN=0
      tidyr::pivot_longer(cols = c("number","no.number"), names_to = "group.fill",values_to = "number") %>%
      dplyr::mutate(method2=paste0(method,".",group.fill)) %>% 
      dplyr::group_by(dst,method,method2,group,group.fill,cutoff,sample) %>%
      dplyr::summarise(number_mean = mean(number, na.rm=T), #,trim = 0.05, na.rm=T
                       number_sd = sd(number),
                       number_upper = Rmisc::CI(number,ci = 0.95)[1],
                       number_lower = Rmisc::CI(number,ci = 0.95)[3]
      ) %>% 
      dplyr::mutate(number=number_mean)
    
    table(res.df2$group.fill)
    res.df2$method2 <- factor(res.df2$method2)
    res.df2$method <- factor(res.df2$method,levels = methods.labs)
    res.df2$label <- base::format(res.df2$number,digits = 3)
    res.df2$label[res.df2$group.fill=="no.number"] <- ""
    #res.df2$number <- res.df2$number
    write.table(res.df2,paste0("./figure/",MAP,"_",cfpeakMAP,"_",TYPE,"/",dst,"_benchmark.txt"),quote = F,row.names = F,col.names = T,sep = '\t')
    
    
    ### plot test bar
    this.dst <- dst #"TCGA-COAD_small"  # GSE71008_NCpool
    for (g in unique(res.df2$group)){
      # g <- "precision"
      print(g)
      tmp <- res.df2[res.df2$group==g & res.df2$dst==this.dst,]
      tmp <- tmp[tmp$group.fill=="number",]
      ggpubr::ggbarplot(data = tmp, x = "method", y = "number", fill = "method", position = position_stack(), width = 0.5) + # , add="mean_se", , facet.by=c("group")
      scale_fill_manual(values = methCols) +
      geom_text(aes(x=method,y=0.15*max(tmp$number),label=label),size=6,color="black",angle=90)+
        labs(title="",x="", y = tocapital(g))+
        theme_minimal() +  # base_size=12
        theme(#axis.ticks.x=element_blank(),
          #strip.text.y = element_blank(),
          aspect.ratio = 1,
          strip.text = element_text(size=20),
          axis.title.x = element_text(size=20),
          axis.title.y = element_text(size=24),
          axis.text.x = element_text(size = 20,hjust = 1,vjust = 0.5,angle = 90), #
          axis.text.y = element_text(size = 20),
          plot.title = element_text(size=20),
          # strip.text = element_blank(),
          legend.position = "none", #c(0.9,0.8),#,#
          legend.text = element_text(size= 16),
          legend.title= element_text(size= 16))
      ggsave(paste0("./figure/",MAP,"_",cfpeakMAP,"_",TYPE,"/",this.dst,"_",sample_id,"_",g,".png"),width = 8,height = 6)
    }
    #
    
    
    
    ## plot ROC,PR (optional)
    {
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
    tmp_df <- res.df %>% 
      # dplyr::filter(smp=="TCGA-A6-2685-11A-01R-1757-13_mirna_gdc_realn") %>% 
      # dplyr::group_by(dst,method,cutoff,sample) %>%
      # dplyr::summarise(fpr_mean = mean(fpr, na.rm=T), #,trim = 0.05, na.rm=T
      #                  fpr_sd = sd(fpr),
      #                  fpr_upper = Rmisc::CI(fpr,ci = 0.95)[1],
      #                  fpr_lower = Rmisc::CI(fpr,ci = 0.95)[3],
      #                  
      #                  recall_mean = mean(recall, na.rm=T), #,trim = 0.05, na.rm=T
      #                  recall_sd = sd(recall),
      #                  recall_upper = Rmisc::CI(recall,ci = 0.95)[1],
      #                  recall_lower = Rmisc::CI(recall,ci = 0.95)[3],
      #                  
      #                  precision_mean = mean(precision, na.rm=T), #,trim = 0.05, na.rm=T
      #                  precision_sd = sd(precision),
      #                  precision_upper = Rmisc::CI(precision,ci = 0.95)[1],
      #                  precision_lower = Rmisc::CI(precision,ci = 0.95)[3]
      # ) %>% 
      # dplyr::mutate(fpr=fpr_mean,recall=recall_mean,precision=precision_mean) %>% 
      as.data.frame() #%>%
    tmp_df <- tmp_df[order(tmp_df$sample,tmp_df$method,tmp_df$cutoff),]
    #head(tmp_df)
    tmp_df$method <- factor(tmp_df$method,levels = methods.labs)
    write.table(tmp_df,paste0("./figure/",MAP,"_",cfpeakMAP,"_",TYPE,"/",dst,"_benchmark2.txt"),quote = F,row.names = F,col.names = T,sep = '\t')
    # tmp1 <- calculate_auroc_aupr(TP = tmp_df$TP, TN = tmp_df$TN, FP = tmp_df$FP, FN = tmp_df$FN)
    # message(tmp1)
    # plot ROC and PR
    my_theme <- theme_minimal() +  # base_size=12
      theme(#axis.ticks.x=element_blank(),
        #strip.text.y = element_blank(),
        aspect.ratio = 1,
        strip.text = element_text(size=20),
        axis.title.x = element_text(size=24),
        axis.title.y = element_text(size=24),
        axis.text.x = element_text(size = 20), # angle = 90,vjust = 0.5,hjust = 1
        axis.text.y = element_text(size = 20),
        plot.title = element_text(size=20),
        # strip.text = element_blank(),
        legend.position = "none", #c(0.9,0.8),#,#
        legend.text = element_text(size= 16),
        legend.title= element_text(size= 16))
    }
    
    {
    #tmp_df$method <- factor(tmp_df$method,levels = c("Piranha","CLIPper","CLAM","exPeak")) # ,"LocalMax"
    #str(tmp_df)
    tmp_roc <- ggplot(tmp_df, aes(x = fpr, y = recall, color = method)) +
      # geom_point() +
      geom_path(size = 1) +
      xlim(c(0,1)) +
      # ylim(c(0,1)) +
      xlab("False Positive Rate") +
      ylab("Recall") +
      scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1.00), limits = c(0,1),labels = c("0.00","0.25","0.50","0.75","1.00")) +
      geom_abline(slope = 1, linetype = 2) +
      # facet_wrap(sample ~ ., nrow = 1,) +
      # ggsci::scale_color_d3()+
        scale_color_manual(values = methCols) +
      # scale_color_brewer(palette = "RdYlBu", direction = -1) +
      # theme_minimal() +
      my_theme
      # theme(legend.position = "")
    # theme(plot.background = element_rect(fill = "white"))
    # tmp_roc
    tmp_pr <- ggplot(tmp_df, aes(x = recall, y = precision, color = method)) +
      # geom_point()+
      geom_path(size = 1) +
      xlim(c(0,1)) +
      xlab("Recall") +
      ylab("Precision") +
      scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1.00), limits = c(0,1),labels = c("0.00","0.25","0.50","0.75","1.00")) +
      # facet_wrap(sample ~ ., nrow = 1, scales = "free_y") +
      # ggsci::scale_color_d3()+
      scale_color_manual(values = methCols) +
      # scale_color_brewer(palette = "RdYlBu", direction = -1) +
      # theme_minimal() +
      my_theme
    # theme(plot.background = element_rect(fill = "white"))
    # tmp_pr
    ggpubr::ggarrange(plotlist = list(tmp_roc,tmp_pr), heights=c(1,1), ncol = 1, nrow = 2, align = "hv",common.legend = T,legend = "none")
    }
    ggsave(paste0("./figure/",MAP,"_",cfpeakMAP,"_",TYPE,"/",this.dst,"_",sample_id,"_roc_pr.png"), width = 5, height = 10, dpi = 600)
    # dir.create("data/figures/roc", recursive = TRUE)
    # dir.create("data/figures/pr", recursive = TRUE)
    # ggsave("./roc.pdf", tmp_roc, width = 10, height = 5, dpi = 600)
    # ggsave("./pr.pdf", tmp_pr, width = 10, height = 5, dpi = 600)
    # }
    #
    
    
    
    
    
    # # plot mult-dataset AUROC/AUPR barplot (fig3, calculate area under curve manually, deprecated)
    # # devtools::install_github('smin95/smplot2')
    # library(smplot2)
    # res2.df <- res.df 
    # res2.list <- list()
    # for (frac in c("")){ # "0.1","0.5"
    #   for (top in c("")){ # "500","1000","2000","5000","10564"; '9426', '7077', '4893','3419'
    #     for (dst in unique(res2.df$dst)){
    #       # dst <- "WSQ_SMARTer_NEB_NCpool"
    #       print(dst)
    #       tmp <- res2.df[res2.df$dst==dst & res2.df$frac==frac & res2.df$top==top,]
    #       for (smp in unique(tmp$sample)){
    #         # smp <- "NCpool_smart"
    #         print(smp)
    #         tmp2 <- tmp[tmp$sample==smp,]
    #         for (method in unique(tmp2$method)){
    #           method.name <- method
    #           print(method.name)
    #           tmp3 <- tmp2[tmp2$method==method,]
    #           res2.list[[paste0(frac,top,dst,smp,method.name)]] <- data.frame("dataset"=dst,"method"=method.name,"sample"=smp,"frac"=frac,"top"=top,
    #                                                                           "AUROC"=smplot2::sm_auc(tmp3$fpr, tmp3$recall),
    #                                                                           "AUPR"=smplot2::sm_auc(tmp3$recall, tmp3$precision))
    #         }
    #       }
    #     }
    #   }
    # }
    # res2.list.df <- do.call("rbind",res2.list)
    # res2.list.df$method <- factor(res2.list.df$method, levels = methods.labs)
    # ggpubr::ggbarplot(data = res2.list.df, x = "dataset", y = "AUROC", fill = "method", position = position_dodge(width = 0.8)) + #, add="mean_se", facet.by=c("dataset")
    #   # scale_fill_d3()+
    #   scale_color_manual(values = methCols) +
    #   scale_fill_manual(values = methCols) +
    #   labs(title="",x="", y = "AUROC")+
    #   # xlim(c(0,1))+   
    #   ylim(c(0,1))+
    #   theme_minimal() +  # base_size=12
    #   theme(#axis.ticks.x=element_blank(),
    #     #strip.text.y = element_blank(),
    #     strip.text = element_text(size=20),
    #     axis.title.x = element_text(size=20),
    #     axis.title.y = element_text(size=20),
    #     axis.text.x = element_text(size = 20,hjust = 1,vjust = 0.5,angle = 90), #
    #     axis.text.y = element_text(size = 20),
    #     plot.title = element_text(size=20),
    #     # strip.text = element_blank(),
    #     legend.position = "right", #c(0.9,0.8),#,#
    #     legend.text = element_text(size= 16),
    #     legend.title= element_text(size= 16)) + 
    #   facet_grid(frac~top)
    # #ggsave("multi_AUROC.pdf",width = 12,height = 8)
    # ggsave("./figure/",MAP,"_",cfpeakMAP,"_",TYPE,"/multi_AUROC.pdf",width = 24,height = 10)
    # #table(res2.list.df$X)
    # ggpubr::ggbarplot(data = res2.list.df, x = "dataset", y = "AUPR", fill = "method", position = position_dodge(width = 0.8)) + #, add="mean_se", facet.by=c("dataset")
    #   # scale_fill_d3()+
    #   scale_color_manual(values = methCols) +
    #   scale_fill_manual(values = methCols) +
    #   labs(title="",x="", y = "AUPR")+
    #   # xlim(c(0,1))+   
    #   ylim(c(0,1))+
    #   theme_minimal() +  # base_size=12
    #   theme(#axis.ticks.x=element_blank(),
    #     #strip.text.y = element_blank(),
    #     strip.text = element_text(size=20),
    #     axis.title.x = element_text(size=20),
    #     axis.title.y = element_text(size=20),
    #     axis.text.x = element_text(size = 20,hjust = 1,vjust = 0.5,angle = 90), #
    #     axis.text.y = element_text(size = 20),
    #     plot.title = element_text(size=20),
    #     # strip.text = element_blank(),
    #     legend.position = "right", #c(0.9,0.8),#,#
    #     legend.text = element_text(size= 16),
    #     legend.title= element_text(size= 16)) +
    #   facet_grid(frac~top)
    # #ggsave("multi_AUPR.pdf",width = 12,height = 8)
    # ggsave(paste0("./figure/",MAP,"_",cfpeakMAP,"_",TYPE,"/multi_AUPR.pdf"),width = 24,height = 10)
    
  }
}
#




# calculate confusion matrix (multi-dst) -----------------------------------------------------------------------
black.col <- alpha("grey90",0.3)
methCols <- c("Piranha"="#1F77B4E5","CLIPper"="#FF7F0EE5","CLAM"="#2CA02CE5","cfPeak"="red","cfPeak.CNN"="firebrick","blockbuster"="#9467BDE5","blockbuster.k100"="#9C564BE5","blockbuster.EM"="#E377C2E5")

blood.dsts <- c("GSE71008_NCpool2","GSE94533_NCpool2","GSE123972_NCpool2", "GSE110381_NCpool2", "GSE94582_NCpool2", "Phospho-RNA-seq_NCpool2", "AGO2_IP_NCpool2") # TGIRT-seq
blood.dsts.lab <- c("GSE71008","GSE94533","GSE123972", "GSE110381", "GSE94582", "GSE126051", "AGO2_IP") # TGIRT-seq
names(blood.dsts.lab) <- blood.dsts
nonblood.dsts <- c("CNP0003091_urine_NCpool","GSE129255","GSE112343_NCpool","GSE56686")
nonblood.dsts.lab <- c("CNP0003091_urine","GSE129255_urine_EV","GSE112343_bile","GSE56686_semen")
names(nonblood.dsts.lab) <- nonblood.dsts
tissue.dsts <- c("TCGA-COAD_small_NC_NCpool","snc_pandora_hsa_HeLa_NCpool", "encode_small_colon_NCpool") # 
tissue.dsts.lab <- c("TCGA_COAD","PANDORA_HeLa", "ENCODE_colon") # TCGA-COAD_small_NC_NCpool
names(tissue.dsts.lab) <- tissue.dsts
dsts <- c(tissue.dsts,nonblood.dsts,blood.dsts)
#dsts <- c("TCGA-COAD_small_NC_NCpool") # "GSE71008","GSE94533", "GSE110381", "WSQ_SMARTer_NEB","GSE94582" ,"GSE50676"
# dsts <- c("CNP0003091_urine_NCpool","GSE129255","GSE112343_NCpool","GSE56686", "snc_pandora_hsa_HeLa_NCpool", "encode_small_colon_NCpool")
# dsts <- c("GSE71008_NCpool","GSE94533_NCpool","GSE123972_NCpool", "GSE94582_NCpool","WSQ_SMARTer_NEB_NCpool","GSE148861_GSE148862_NCpool", "Phospho-RNA-seq_NCpool", "GSE110381_NCpool", "AGO2_IP_NCpool")
#sample_ids <- c("NCpool", "TCGA-4T-AA8H-01A-11H-A41D-13_mirna_gdc_realn") # "NCpool", # TCGA-4T-AA8H-01A-11H-A41D-13_mirna_gdc_realn
#"GSE50676_NCpool"

for(MAP in c("tx")){ # gn
  for(cfpeakMAP in c("tx")){ # "gn"
    for(TYPE in c("all4")){ # called,  all,  all2, all3, all4

      
      #MAP="gn" # gn tx
# cfpeakMAP="tx" # gn (depre) tx
# TYPE="called" # "all"  called
message(MAP," mode")
#dir.create(paste0("./figure/",MAP,"_",cfpeakMAP,"_",TYPE))

if(MAP == "tx"){
  methods <- c("Piranha","CLIPper","CLAM","bowtie2_cfpeak") #StarEM_win50,Star_k100,bowtie2_cfpeak,,Piranha,CLIPper,CLAM  ,"bowtie2_cfpeakCNN" "StarEM_win50"
  methods.labs <- c("Piranha","CLIPper","CLAM","cfPeak") #blockbuster.k100 cfPeak.CNN
  #StarEM_win2000
  # exclud.meth <- c("blockbuster.k100","blockbuster.EM") # "blockbuster.k100" ,"cfPeak.CNN"
  # methods.labs <- methods.labs[!methods.labs %in% exclud.meth]
}else if(MAP == "gn"){
  #sample_ids <- c("NCpool") #c("SRR2105340","SRR2105341","SRR2105342","SRR2105346","SRR2105335")  #c("SRR5230634","SRR5230661","SRR5230662","SRR5230667") # c("SRR2105336","SRR2105128","SRR2105129","SRR2105335","SRR2105336"), c("NC_PKU-2392860_smart_1","NC_PKU-2392860_smart_PNK_1","NC_PKU-2392860_1"),c(CleanTag_Lab5,NEBNext_Lab1,NEBNext_Lab2,TruSeq_Lab1,TruSeq_Lab2,N4_A_Lab5,N4_B_Lab2),c("GM12878","GM12891")
  # methods <- c("piranha_by_sample/b5_p01","clipper_by_sample/b5_p05","clam_by_sample/b5_p005","expeakCNN_by_sample/b5_d50_p1") # ,"domains_localmax_significant/b5_d05_p05"
  #methods <- c("StarEM_win50","Star_k100","Star_default","bowtie2_cfpeak") # ,"domains_localmax_significant/b5_d05_p05"
  methods <- c("Piranha","CLIPper","CLAM","bowtie2_cfpeak") #StarEM_win50,Star_k100,bowtie2_cfpeak,,Piranha,CLIPper,CLAM  ,"bowtie2_cfpeakCNN" "StarEM_win50", Star_default
  methods.labs <- c("Piranha","CLIPper","CLAM","cfPeak") #blockbuster.k100  blockbuster
  #StarEM_win2000
  # exclud.meth <- c("blockbuster.k100","blockbuster.EM","cfPeak.CNN") # "blockbuster.k100"
  # methods.labs <- methods.labs[!methods.labs %in% exclud.meth]
}


## get matrix in R
res <- list()
for (dst in dsts){
  # dst <- "test"
  # print(dst)
  # dst <- "AGO2_IP_NCpool"
  # tmp <- tmp_df[tmp_df$dst==dst,]
  # sample_ids <- read.table(paste0("data/",dst,"/sample_ids_NC_test.txt"))$V1
  sample_ids <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/",dst,"/sample_ids_raw.txt"))$V1
  #sample_ids <- "TCGA-A6-2685-11A-01R-1757-13_mirna_gdc_realn"
  message(dst)
  for (sample_id in sample_ids){
    # sample_id <- "CRC"
    # print(sample_id)
    
    ## get matrix in R

    # for (frac in c("")){ # "0.1","0.5"
    frac <- ""
    #   for (top in c("")){ # "500","1000","2000","5000","10564"; '9426', '7077', '4893','3419'
    top <- ""
    
    for (method in methods){
      # method <- "Star_default"
      if(grepl("cfpeak",tolower(method))){
        message("cfpeak")
        MAP2 <- cfpeakMAP
      }else{
        message("other")
        MAP2 <- MAP
      }
      
      inputFile <- paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/test/",dst,"_",method,"_",sample_id,"_",MAP2,"_",TYPE,".txt") # OPT1: .txt6   OPT2:.txt7
      method.name <- method
      
      tmp <- read.table(inputFile,header = T,sep = "")
      tmp$specificity <- 1-tmp$fpr
      tmp$replicate <- ""
      colnames(tmp)[colnames(tmp)=="replicate"] <- "sample"
      colnames(tmp)[colnames(tmp)=="counts"] <- "cutoff"
      tmp$method <- method.name
      tmp$dst <- dst
      tmp$top <- top
      tmp$frac <- frac
      tmp$smp <- sample_id
      
      res[[paste0(dst,sample_id,method.name,frac,top)]] <- tmp
    }
  }
}

res.df <- as.data.frame(do.call(rbind,res))
table(res.df$dst, res.df$method)
res.df$method[res.df$method=="Star_default"] <- "blockbuster"
res.df$method[res.df$method=="Star_k100"] <- "blockbuster.k100"
res.df$method[res.df$method=="StarEM_win50"] <- "blockbuster.EM"
res.df$method[res.df$method=="bowtie2_cfpeak"] <- "cfPeak"
res.df$method[res.df$method=="bowtie2_cfpeakCNN"] <- "cfPeak.CNN"
table(res.df$method)
unique(res.df$dst)
#res.df[1:3]
# summary(res.df$fpr)
# omit any record that not in [0,1]
res.df <- res.df[res.df$precision>=0 & res.df$precision<=1 & res.df$recall>=0 & res.df$recall<=1 & res.df$fpr>=0 & res.df$fpr<=1 & res.df$f1>=0 & res.df$f1<=1 & res.df$specificity>=0 & res.df$specificity<=1 ,]
res.df <- res.df[res.df$method %in% methods.labs,]


## plot standard-overlap/recall/precison pie/bar plot 
# pval_seq2 <- c(0,10^seq(-15,-2,5),0.05,0.1)
# res.df2 <- res.df[res.df$cutoff %in% pval_seq2,]
res.df2 <- res.df #[res.df$dst=="CNP0003091_urine_NCpool" & res.df$smp=="NCpool" & res.df$method=="blockbuster",]
# calculate_auroc_aupr(res.df2$TP, res.df2$TN, res.df2$FP, res.df2$FN)
res.df2$sample <- res.df2$smp
#res.df2$recurrence <- res.df2$TP/(res.df2$TP+res.df2$FP+res.df2$TN+res.df2$FN)
res.df2 <- dplyr::as_tibble(res.df2) %>% 
  dplyr::group_by(dst,method,sample) %>%
  # dplyr::mutate(AUROC=calculate_auroc_aupr(TP, TN, FP, FN)$AUROC,AUPR=calculate_auroc_aupr(TP, TN, FP, FN)$AUPR) %>%
  dplyr::arrange(cutoff) %>% # only select min cutoff (all peaks)
  dplyr::distinct(method, sample, dst, .keep_all = TRUE) %>% 
  tidyr::pivot_longer(cols = precision:specificity, names_to = "group",values_to = "number") %>%
  dplyr::mutate(total.peak.num=TP+FP+TN+FN, no.number=1-number) %>% # TN+FN=0
  tidyr::pivot_longer(cols = c("number","no.number"), names_to = "group.fill",values_to = "number") %>%
  dplyr::mutate(method2=paste0(method,".",group.fill))
# table(res.df2$group.fill)


res.df2$method2 <- factor(res.df2$method2, levels = as.vector(outer(methods.labs, c(".no.number",".number"), paste0)))
res.df2$method <- factor(res.df2$method, levels = methods.labs)
# table(res.df2$method)
#table(res.df2$method2 )
#str(res.df2)
#unique(res.df2$method2)
# mypal <- pal_d3_adaptive("category10", alpha = 0.9)(4)
# black.col <- alpha("grey90",0.3)
#scales::show_col(mypal)
#table(res.df2$group)
#tmp <- res.df2[res.df2$group=="TP",]
# res.df2 <- res.df2[res.df2$group!="fpr",]
res.df2$label <- base::format(res.df2$number,digits = 1)
res.df2$label[res.df2$group.fill=="no.number"] <- ""
table(res.df2$label)
table(res.df2$group)
table(res.df2$dst)
table(res.df2$smp)


# ### plot single dst (GSE71008_NCpool) (fig3)
# this.dst <- "TCGA-COAD_small_NC_NCpool"  # GSE71008_NCpool
# for (g in unique(res.df2$group)){
#   # g <- "precision"
#   print(g)
#   #table(res.df2$dst==this.dst)
#   #unique(res.df2$dst)
#   tmp <- res.df2[res.df2$group==g & res.df2$dst==this.dst,]
#   tmp <- tmp[tmp$group.fill=="number",]
#   ggpubr::ggbarplot(data = tmp, x = "method", y = "number", fill = "method", position = position_stack(), width = 0.5) + # , add="mean_se", , facet.by=c("group")
#   scale_fill_manual(values = methCols) +
#     geom_text(aes(x=method,y=0.15*max(tmp$number),label=label),size=6,color="black",angle=90)+
#     labs(title="",x="", y = g)+
#     # xlim(c(0,1))+
#     # ylim(c(0,1))+
#     theme_minimal() +  # base_size=12
#     theme(#axis.ticks.x=element_blank(),
#       #strip.text.y = element_blank(),
#       aspect.ratio = 1,
#       strip.text = element_text(size=20),
#       axis.title.x = element_text(size=20),
#       axis.title.y = element_text(size=20),
#       axis.text.x = element_text(size = 20,hjust = 1,vjust = 0.5,angle = 90), #
#       axis.text.y = element_text(size = 20),
#       plot.title = element_text(size=20),
#       # strip.text = element_blank(),
#       legend.position = "none", #c(0.9,0.8),#,#
#       legend.text = element_text(size= 16),
#       legend.title= element_text(size= 16))
#   ggsave(paste0("./figure/",this.dst,"_",g,".pdf"),width = 5,height = 5)
# }





### plot multi dsts (GSE71008_NCpool...) (suppl fig3 ?)
#dsts <- c("GSE71008_NCpool", "GSE94533_NCpool")  # , "GSE110381_NCpool", "GSE94582_NCpool"
#"GSE123972_NCpool", "Phospho-RNA-seq_NCpool", "GSE110381_NCpool","GSE94582_NCpool","WSQ_SMARTer_NEB_NCpool","GSE148861_GSE148862_NCpool","AGO2_IP_NCpool","GSE50676_NCpool"
# res.df2$X <- res.df2$dst
# res.df2$X <- gsub("_NCpool|_NCpool2","",res.df2$X,perl = T)
# res.df2$X <- gsub("_small|hsa_|snc_","",res.df2$X,perl = T)
# res.df2$X <- gsub("_NEBNext","",res.df2$X )
# res.df2$X <- gsub("WSQ_SMARTer_NEB","in-house",res.df2$X )
# res.df2$X <- gsub("Phospho-RNA-seq","GSE126051",res.df2$X )
# res.df2$X <- gsub("GSE148861_GSE148862","GSE148861",res.df2$X )
# # res.df2$X <- factor(res.df2$X, levels = c("GSE71008", "GSE94533","GSE123972", "GSE110381","GSE94582", "GSE126051", "in-house","AGO2_IP","GSE148861"))
# table(res.df2$X)
# res.df2$X2 <- paste0(res.df2$X,".",res.df2$sample)
# table(res.df2$X2)
# res.df2$X2 <- gsub("NCpool_|NCpool2_","",res.df2$X2,perl = T)
# res.df2$X2 <- gsub(".NCpool","",res.df2$X2,fixed = T)
all.dst <- c(tissue.dsts.lab,nonblood.dsts.lab,blood.dsts.lab)
res.df2$X <- all.dst[res.df2$dst]
res.df2$X <- factor(res.df2$X,levels = all.dst)
# table(res.df2$X2)
res.df2$sample <- gsub("NCpool_","",res.df2$sample,fixed = T)
res.df2$sample <- gsub("NCpool","",res.df2$sample,fixed = T)
res.df2 <- res.df2[order(res.df2$X,res.df2$sample,decreasing = F),]
res.df2 <- res.df2[!res.df2$sample %in% c("TCGA-A6-2683-11A-01R-1757-13_mirna_gdc_realn",
                                          "TCGA-A6-2684-11A-01R-1757-13_mirna_gdc_realn"),]
res.df2$sample <- gsub("TCGA-4T-AA8H-01A-11H-A41D-13_mirna_gdc_realn","01A",res.df2$sample,fixed = T)
res.df2$sample <- gsub("TCGA-A6-2685-11A-01R-1757-13_mirna_gdc_realn","11A",res.df2$sample,fixed = T)
res.df2$X2 <- paste0(res.df2$X," ",res.df2$sample)
res.df2$X2 <- factor(res.df2$X2,levels = unique(res.df2$X2))
#res.df2$Specificity <- 1-res.df2$fpr

#str(res.df2)
for (g in unique(res.df2$group)){
  # g <- "precision"
  # if(g == "fpr"){
  #   res.df2
  # }
  print(g)
  tmp <- res.df2[res.df2$group==g,]
  tmp <- tmp[tmp$group.fill=="number",]
  maxY <- max(tmp$number)
  if(g=="fpr"){
    g="FPR"
  }else{
    g <- tocapital(g)
  }
  ggpubr::ggbarplot(data = tmp, x = "X2", y = "number", fill = "method", position = position_dodge(width = 0.8)) + #, add="mean_se", facet.by=c("dataset")
    # scale_fill_d3()+
    scale_fill_manual(values = methCols) +
    labs(title="",x="", y = g)+ # stringr::str_to_title(g)
    # xlim(c(0,1))+   
    # ylim(c(0,1))+
    scale_y_continuous(breaks = seq(0,maxY,round(maxY/3,digits = 1))) +
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
      legend.position = "none", #c(0.9,0.8),#,#
      legend.text = element_text(size= 16),
      legend.title= element_text(size= 16)) #+
  # facet_grid(frac~top)
  #ggsave("multi_AUPR.pdf",width = 12,height = 8)
  ggsave(paste0("./figure/multi_",MAP,"_",cfpeakMAP,"_",TYPE,"_",g,".png"),width = 12,height = 5)
}
    }
  }
}
#



### plot single dst (GSE71008_NCpool2) (fig3)
this.dst <- "GSE71008_NCpool2"  # 
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
  ggsave(paste0("./figure/",this.dst,"_",g,".png"),width = 5,height = 5)
}









