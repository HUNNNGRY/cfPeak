#coverage meta plot
#test EnrichedHeatmap + deeptools + R
# last 20221207 by b.p.f@qq.com

# def. func. ------------------------------------------------
source("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/scripts/util.R")
#











# comapre Fig1 coverage meta     plot (enrichedHeatmap + R, old, 221212) ------------------------------------------------
#plot tx meta-plot and Heatmap (long reads & peaks, on CDS, UTR)

## get matrix by EnrichedHeatmap in R
#library(EnrichedHeatmap)
# setwd("/Users/baopengfei/Desktop/lulab/tmp/")
setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC")
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggsci)
library(Rmisc)
# library(data.table)
conflict_prefer("select", "dplyr")


#tx
up=50
down=50
bin=5 # 10 for test
ratio=0.3
#tbed_11RNA_primary for clipper; tbed_RNA_EM for expeak !!!!

# use NCpool
plasma.dst <- "GSE71008_NCpool"  # "GSE110381_NCpool"
cmb <- data.frame(signal=c(paste0("output/",plasma.dst,"/call_peak_all/tbed_11RNA_primary/","NCpool",".bed.gz"), # tbed_RNA_EM
                           # paste0("output/","TCGA_small_NCpool","/call_peak_all/tbed_RNA_EM/","NCpool",".bed.gz"),
                           paste0("output/","GSE148861_GSE148862_NCpool","/call_peak_all/tbed_11RNA_primary/","NCpool",".bed.gz"),
                           paste0("output/","GSE50676_NCpool","/call_peak_all/tbed_11RNA_primary/","NCpool",".bed.gz")
                           # paste0("output/","AGO2_IP_NCpool","/call_peak_dedup/tbed_RNA_EM/","NCpool",".bed.gz")
),
region=c(paste0("output/",plasma.dst,"/call_peak_all/clipper_by_sample/b5_p05/","NCpool",".bed"),
         paste0("output/","GSE148861_GSE148862_NCpool","/call_peak_all/clipper_by_sample/b5_p05/","NCpool",".bed"),
         #paste0("output/","TCGA_small_NCpool","/call_peak_all/clipper_by_sample/b5_p05/","NCpool",".bed"),
         paste0("output/","GSE50676_NCpool","/call_peak_all/clipper_by_sample/b5_p05/","NCpool",".bed")
         # paste0("output/","AGO2_IP_NCpool","/call_peak_dedup/domains_clipper_by_sample/b5_p05/","NCpool",".bed")
),
signal.label=c("extracellular small","cellular small","cellular CLIP"), # ,"ex RIP"
region.label=c("CLIPper","CLIPper","CLIPper") # ,"CLIPper","exPeak"
)

## use all smps (in coordance with other fig)
#need use read bed as signal, thus only NC pool used 


#one way using sapply
res.list <- lapply(1:nrow(cmb), function(j) {get.mat(signal=as.character(cmb$signal[j]), region=as.character(cmb$region[j]),
                                                     signal.label=as.character(cmb$signal.label[j]),region.label=as.character(cmb$region.label[j]),
                                                     up = up, down = down, bin = bin, ratio = ratio)} )
res.list.df <- do.call(rbind,res.list)


# plot in R
res.list.df$signal <- factor(res.list.df$signal,levels = c("cellular CLIP", "cellular small","extracellular small"))
# res.list.df$region <- CLIPper
# table(res.list.df$signal)

#(500+500)/20+10=60,再乘以1个样本,再乘以1个区域,为60列
clean_refp <- dplyr::as_tibble(res.list.df)
#clean_refp <- refp #%>% select(c(-1:-3,-5,-6))
colnames(clean_refp)[c(1,3)] <- c('V4','sample')

#max-min scale
clean_refp[,4:ncol(clean_refp)] <- t(maxmin.normalize(t(clean_refp[,4:ncol(clean_refp)])))

long_refp <- reshape2::melt(clean_refp,id.vars = c('V4','sample','region'),value.name = 'signal') %>%
  select(-variable)
#add sample name
#long_refp$sample <- rep(c("FTA-10","csFTA-10"),
#                        each = nrow(refp)*22) # (100 + 100 + 20)/10

# add x position
long_refp$pos <- rep(c(1:(ncol(res.list.df)-3)),each = nrow(res.list.df)) # ,times = nrow(cmb)
# dim(long_refp)

# calculate means with CI, plot line
filnal_scaler <- long_refp %>%
  # remove na
  tidyr::drop_na() %>%
  dplyr::group_by(pos,sample,region) %>% # sample,
  # mean and 95% interval confidence
  dplyr::summarise(mean_signal = mean(signal), #,trim = 0.05, na.rm=T
                   sd = sd(signal),
                   upper = CI(signal,ci = 0.95)[1],
                   lower = CI(signal,ci = 0.95)[3]
  )

#no 2nd max-min norm here (not needed for this fig)

# plot
# p <- ggplot(filnal_scaler,aes(x = pos,y = mean_signal)) +
#   geom_line(size = 1) + #aes(color = sample),
#   theme_classic(base_size = 14) +
#   scale_color_d3(name = '') +
#   # x label
#   scale_x_continuous(breaks = c(0,25,35,ncol(res.list.df)),
#                      labels = c('-500b','Start','End','+500b')) +
#   xlab('') + ylab('Normalized signal') +
#   theme(aspect.ratio = 0.8)
#p
# add CI
up_bins <- as.integer(up/bin) # ceiling
down_bins <- as.integer(down/bin) # ceiling
target_bins <- as.integer((up_bins + down_bins)*(ratio/(1-ratio)))
ggplot(filnal_scaler,aes(x = pos,y = mean_signal)) +
  # add 0.95 interval
  geom_ribbon(aes(ymin = lower,
                  ymax = upper,
                  fill = sample), # sample
              alpha = 0.5) +
  geom_line(aes(color = sample), size = 1) + # 
  theme_classic(base_size = 16) +
  # ylim(c(0,1.2)) +
  # scale_color_nejm(name = 'Data type') +
  scale_fill_manual(name = '', values = c("#FD6905","#843692","#2C68A9")) +
  scale_color_manual(name = 'Data type', values = c("#FD6905","#843692","#2C68A9")) +
  # geom_vline(xintercept = c(up_bins,up_bins+target_bins),color="grey50",linetype="dashed")+
  geom_vline(xintercept = c((1+up_bins+target_bins+down_bins)*0.5),color="grey50",linetype="dashed")+ # need +1 for 1-based coord
  # x label
  scale_x_continuous(breaks = c(1,up_bins+1,up_bins+target_bins,up_bins+target_bins+down_bins), # need +1 for 1-based coord
                     labels = c(paste0("-",up,"nt"),'Start','End',paste0("+",down,"nt")) ) +
  scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1.00), limits = c(0,1),labels = c("0","0.25","0.50","0.75","1.00")) +
  xlab('') + ylab('Min-max scaled depth') +
  theme(aspect.ratio = 0.4,
        strip.background = element_rect(color = NA,fill = 'grey'),
        strip.text = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20,angle = 0,hjust = 0.5),
        axis.text.y = element_text(size=20),
        legend.position =  "right",#c(0.9,0.8),
        legend.text = element_text(size= 16),
        legend.title= element_text(size= 16)) 
#ggsave("cov_meta.pdf",width = 9,height = 6)



filnal_scaler$sample <- factor(filnal_scaler$sample, levels = c("extracellular small", "cellular small","cellular CLIP"))
filnal_scaler$y <- as.numeric(filnal_scaler$sample) * 0.5
library(ggridges)
ggplot(filnal_scaler, aes(x = pos, y=y, height=mean_signal, fill=sample))+ #
  geom_ridgeline_gradient() + # stat = "summary"
  # theme_ridges() + 
  # theme(legend.position = "none")
  scale_fill_manual(name = '', values = alpha(colour = c("#2C68A9","#843692","#FD6905"),alpha = 0.8)) +
  # scale_color_manual(name = 'Data type', values = c("#FD6905","#843692","#2C68A9")) + 
  geom_vline(xintercept = c((1+up_bins+target_bins+down_bins)*0.5),color="black",linetype="dashed") + # need +1 for 1-based coord
  # x label
  scale_x_continuous(breaks = c(1,up_bins+1,up_bins+target_bins,up_bins+target_bins+down_bins), # need +1 for 1-based coord
                     labels = c(paste0("-",up,"nt"),'Start','End',paste0("+",down,"nt")) ) +
  # scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1.00), limits = c(0,1),labels = c("0","0.25","0.50","0.75","1.00")) +
  xlab('') + ylab('Min-max scaled depth') +
  # theme_bw() + 
  # ylim(c(1.7,4)) + # need change if change dst !!!
  theme_classic() +
  theme(aspect.ratio = 0.4,
        axis.ticks.y=element_blank(), 
        strip.background = element_rect(color = NA,fill = 'grey'),
        strip.text = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20,angle = 0,hjust = 0.5),
        axis.text.y = element_blank(), #element_text(size=20),
        legend.position =  "right",#c(0.9,0.8),
        legend.text = element_text(size= 16),
        legend.title= element_text(size= 16)) 
ggsave("cov_meta_ridge.pdf",width = 9,height = 6)




# heatmap (same region order for 3 datatype, not like fig4 enrich re-order each type)
# plot (slowly)
# long_refp$log10.signal <- log10(long_refp$signal+0.00001)
long_refp_sort <- long_refp %>%
  # dplyr::filter(pos %in% 25:35) %>% # only rank/sort by central regions
  dplyr::group_by(V4,sample) %>%  # sample,
  #mean of each row/region
  dplyr::summarise(mean_signal_region = mean(signal,na.rm=T)
  )
long_refp_sort <- long_refp_sort[order(long_refp_sort$mean_signal_region,decreasing = F),]
long_refp2 <- long_refp
long_refp2$V4 <- factor(long_refp2$V4,levels = unique(long_refp_sort$V4))

clr <- c("#2C68A9","#843692","#FD6905")
datatypes <- as.character(unique(long_refp2$sample)) # "extracellular small" "cellular small"      "cellular CLIP"  
for (i in 1:length(datatypes)){
  datatype <- datatypes[i]
  print(datatype)
  long_refp2_tmp <- long_refp2[long_refp2$sample==datatype,]
  tmp.p <- ggplot(long_refp2_tmp, # %>% filter(sample == 'H3K27ac')
                  aes(x = pos,y = V4)) +
    geom_tile(aes(fill = signal)) +
    theme_bw() +
    # coord_cartesian(expand = 0) +
    # scale_x_continuous(breaks = c(0,10,12,22),
    #                    labels = c('-100b','Start','End','+100b')) +
    scale_x_continuous(breaks = c(1,up_bins+1,up_bins+target_bins,up_bins+target_bins+down_bins), # need +1 for 1-based coord
                       labels = c(paste0("-",up,"nt"),'Start','End',paste0("+",down,"nt")) ) +
    scale_fill_gradient(low = 'white',high = clr[i]) + # 'salmon'
    # scale_fill_viridis_c() +
    # scale_fill_continuous("blues") + 
    ylab('') + xlab('') +
    theme(aspect.ratio = 0.15,
          strip.background = element_blank(), #element_rect(color = NA,fill = 'grey'),
          strip.text = element_text(size=20),
          axis.title.x = element_text(size=20),
          axis.title.y = element_text(size=20),
          axis.text.x = element_text(size=20,angle = 0,hjust = 0.5),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position =  "right",#c(0.9,0.8),
          legend.text = element_text(size= 16),
          legend.title= element_text(size= 16)) #+
  # theme(aspect.ratio = 1.5,
  
  #       panel.grid = element_blank(),
  #       axis.text.x = element_text(size=20,angle = 0,hjust = 0.5),
  #       # axis.text.y = element_text(size=20),
  #       legend.position = "right",
  #       legend.text = element_text(size= 16),
  #       legend.title= element_text(size= 16))+
  # facet_grid(sample ~ . , scales = "free" )
  ggsave(plot = tmp.p,filename = paste0(datatype,"_cov_heat.pdf"),width = 12,height = 3)
}




# comapre Fig4 overlap  meta     plot (enrichedHeatmap + R, deprecated, 221212) ------------------------------------------------
setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC")
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggsci)
library(Rmisc)
#library(data.table)
conflict_prefer("select", "dplyr")


#tx
up=50
down=50
bin=5 # 10 for test
ratio=0.3

dst <- "GSE71008" # GSE71008,GSE110381,GSE123972

## NCpool (change to 11RNA next time ?)
# cmb <- expand.grid(signal=c("/BioII/lulab_b/baopengfei/shared_reference/RBP/splitByRBP/all_tx_score1.bed",
#                             "/BioII/lulab_b/baopengfei/shared_reference/structure/quadratlas_g4grinder_tx_score1.bed",
#                             "/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/lulab/11RNAsmallDomain_EVenrich_score1.bed"),
# region=c(paste0("output/","GSE71008_NCpool","/call_peak_all/expeak_by_sample/b5_d50_p1/","NCpool",".bed"),
#          paste0("output/","GSE71008_NCpool","/call_peak_all/clipper_by_sample/b5_p05/","NCpool",".bed"),
#          paste0("output/","GSE71008_NCpool","/call_peak_all/clam_by_sample/b5_p005/","NCpool",".bed"),
#          paste0("output/","GSE71008_NCpool","/call_peak_all/piranha_by_sample/b5_p01/","NCpool",".bed")
#          )
# )
## sample-wise consensus
cmb <- expand.grid(signal=c("/BioII/lulab_b/baopengfei/shared_reference/RBP/splitByRBP/all_tx_score1.bed",
                            "/BioII/lulab_b/baopengfei/shared_reference/structure/quadratlas_g4grinder_tx_score1.bed",
                            "/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/lulab/11RNAsmallDomain_EVenrich_score1.bed"),
                   region=c(paste0("output/",dst,"/call_peak_all/expeak/b5_d50_p1_11RNA.bed"),
                            paste0("output/",dst,"/call_peak_all/clipper/b5_p05_11RNA.bed"),
                            paste0("output/",dst,"/call_peak_all/clam/b5_p005_11RNA.bed"),
                            paste0("output/",dst,"/call_peak_all/piranha/b5_p01_11RNA.bed")
                   )
)


cmb$signal.label <- ""
cmb$signal.label[grepl("all_tx",cmb$signal)] <- c("RBPs") #,"small","CLIP","ex RIP"
cmb$signal.label[grepl("g4grinder",cmb$signal)] <- c("G4") #,"small","CLIP","ex RIP"
cmb$signal.label[grepl("EVenrich",cmb$signal)] <- c("EV") #,"small","CLIP","ex RIP"
# table(cmb$signal.label)
cmb$region.label <- "Piranha"
cmb$region.label[grepl("clam",cmb$region)] <- "CLAM"
cmb$region.label[grepl("clipper",cmb$region)] <- "CLIPper"
cmb$region.label[grepl("expeak",cmb$region)] <- "exPeak"


#one way using sapply
res.list <- lapply(1:nrow(cmb), function(j) {get.mat(signal=as.character(cmb$signal[j]), region=as.character(cmb$region[j]),
                                                     signal.label=as.character(cmb$signal.label[j]),region.label=as.character(cmb$region.label[j]),
                                                     up = up, down = down, bin = bin, ratio = ratio, )} )
res.list.df <- do.call(rbind,res.list)
write.table(res.list.df,"tmp/comapre_Fig5_overlap_meta.txt",row.names = T,col.names = T,sep = "\t",quote = F)
res.list.df <- read.table("tmp/comapre_Fig5_overlap_meta.txt",sep = "\t",header = T)
# res.list.df[1:3,]
# table(res.list.df$signal)


# plot in R
res.list.df$signal <- gsub("G4iM","G4",res.list.df$signal)
res.list.df$signal <- factor(res.list.df$signal,levels = c("RBPs","EV","G4"))
res.list.df$region <- factor(res.list.df$region,levels = c("Piranha","CLIPper","CLAM","exPeak"))
# table(res.list.df$region)


#(500+500)/20+10=60,再乘以1个样本,再乘以1个区域,为60列
clean_refp <- dplyr::as_tibble(res.list.df)
#clean_refp <- refp #%>% select(c(-1:-3,-5,-6))
colnames(clean_refp)[c(1,3)] <- c('V4','sample')
#plot(x=1:(ncol(clean_refp)-3), y=as.numeric(clean_refp[1,4:ncol(clean_refp)]))
table(clean_refp$region,clean_refp$sample)

clean_refp <- clean_refp[which( apply(clean_refp[,4:ncol(clean_refp)] ,1, sd)>=0),] # filter low sd
#~ half num of peaks have no ref overlap: sd==0

#max-min scale
clean_refp[,4:ncol(clean_refp)] <- t(maxmin.normalize(t(clean_refp[,4:ncol(clean_refp)])))
# clean_refp.tmp <- clean_refp[clean_refp$sample=="G4iM",]
# summary(as.numeric(clean_refp.tmp[3,4:ncol(clean_refp.tmp)]))

# clean_refp[1:6,1:6]
long_refp <- reshape2::melt(clean_refp,id.vars = c('V4','sample','region'),value.name = 'signal') %>%
  select(-variable)

# add x position
long_refp$pos <- rep(c(1:(ncol(clean_refp)-3)),each = nrow(clean_refp)) # ,times = nrow(cmb)
# dim(long_refp)


#lowess model 
correct.loess <- function(long_refp){
  lo = stats::loess(signal ~ pos, data = long_refp)
  long_refp$signal = mean(long_refp$signal, na.rm = TRUE) * long_refp$bc/stats::predict(lo, newdata = long_refp[,c("signal","pos")])
  long_refp$signal = round(long_refp$signal, digits = 2)
  if (any(long_refp$signal < 0, na.rm = TRUE))
    long_refp$signal[long_refp$signal < 0] = 0
  return(long_refp)
}
get.loess <- function(long_refp){
  lo = stats::loess(signal ~ pos, data = long_refp, span = 0.5) # default: span=0.75, higher is smoother
  long_refp$signal = stats::predict(lo, newdata = long_refp_tmp$pos)
  
  long_refp$signal = round(long_refp$signal, digits = 2)
  if (any(long_refp$signal < 0, na.rm = TRUE))
    long_refp$signal[long_refp$signal < 0] = 0
  return(long_refp)
}

l <- list()
for(smp in as.character(unique(long_refp$sample))){
  #smp <- "RBPs"
  print(smp)
  for (region in as.character(unique(long_refp$region))){
    #region <- "exPeak"
    print(region)
    long_refp_tmp <- long_refp[long_refp$sample==smp & long_refp$region==region,]
    #plot(long_refp_tmp$pos, long_refp_tmp$signal)
    long_refp_tmp2 <- get.loess(long_refp_tmp[,c("pos","signal")])
    long_refp_tmp$signal <- long_refp_tmp2$signal
    colnames(long_refp_tmp)[colnames(long_refp_tmp)=="signal"] <- "mean_signal"
    l[[paste0(smp,region)]] <- long_refp_tmp
    #table(long_refp_tmp$pos==long_refp_tmp2$pos)
  }
}
filnal_scaler <- do.call(rbind,l)

# # calculate means with CI, plot line
# filnal_scaler <- long_refp %>%
#   # remove na
#   tidyr::drop_na() %>%
#   dplyr::group_by(pos,sample,region) %>% # sample,
#   # mean and 95% interval confidence
#   dplyr::summarise(mean_signal = mean(signal, na.rm=T), #,trim = 0.05, na.rm=T
#                    sd = sd(signal),
#                    upper = CI(signal,ci = 0.95)[1],
#                    lower = CI(signal,ci = 0.95)[3],
#                    max = max(signal),
#                    min = min(signal)
#   )

up_bins <- as.integer(up/bin) # ceiling
down_bins <- as.integer(down/bin) # ceiling
target_bins <- round((up_bins + down_bins)*(ratio/(1-ratio)),digits = 0) #as.integer((up_bins + down_bins)*(ratio/(1-ratio)))
table(filnal_scaler$sample)
table(filnal_scaler$region)
filnal_scaler$region <- factor(filnal_scaler$region,levels = c("Piranha","CLIPper","CLAM","exPeak"))

#2nd max-min scale (only display mean line, no ribbon)
filnal_scaler <- filnal_scaler %>% 
  dplyr::group_by(sample,region) %>%
  dplyr::mutate(mean_signal_norm= (mean_signal - min(mean_signal))/(max(mean_signal)-min(mean_signal))^as.logical(sd(mean_signal)) )

ggplot(filnal_scaler,aes(x = pos,y = mean_signal_norm)) + # 
  # add 0.95 interval
  # geom_ribbon(aes(ymin = lower,
  #                 ymax = upper,
  #                 fill = region), # sample
  #             alpha = 0.5) +
  geom_line(aes(color = region), size = 2) + # 
  theme_classic(base_size = 16) +
  scale_color_d3(name = 'Data type') +
  scale_fill_d3(name = '') +
  geom_vline(xintercept = c(up_bins+1,up_bins+target_bins),color="grey50",linetype="dashed")+
  # geom_vline(xintercept = c((1+up_bins+target_bins+down_bins)*0.5),color="grey50",linetype="dashed")+ # need +1 for 1-based coord
  # x label
  scale_x_continuous(breaks = c(1,up_bins+1,up_bins+target_bins,up_bins+target_bins+down_bins), # need +1 for 1-based coord
                     labels = c(paste0("-",up,"nt"),'Start','End',paste0("+",down,"nt")) ) +
  scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1.0), limits = c(0,1), labels = c("0","0.25","0.50","0.75","1.00")) +
  xlab('') + ylab('Scaled depth') +
  theme(aspect.ratio = 0.6,
        strip.background = element_rect(color = NA,fill = 'grey'),
        strip.text = element_text(size=30),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=30),
        axis.text.x = element_blank(), #  element_text(size=20,angle = 0,hjust = 0.5)
        axis.text.y = element_text(size=30),
        legend.position =  "right",#c(0.9,0.8),
        legend.text = element_text(size= 24),
        legend.title= element_text(size= 30)) + 
  facet_grid(#.~sample, 
    cols = vars(sample), 
    scales = "free_y")
ggsave("cov_meta.pdf",width = 20,height = 8)


# heatmap (need further optimization)
# # plot (slowly)
# # long_refp$log10.signal <- log10(long_refp$signal+0.00001)
# long_refp_sort <- long_refp %>%
#   # dplyr::filter(pos %in% 25:35) %>% # only rank/sort by central regions
#   dplyr::group_by(V4,sample) %>%  # sample,
#   #mean of each row/region
#   dplyr::summarise(mean_signal_region = mean(signal,na.rm=T)
#   )
# long_refp_sort <- long_refp_sort[order(long_refp_sort$mean_signal_region,decreasing = F),]
long_refp2 <- long_refp
# long_refp2$V4 <- factor(long_refp2$V4,levels = unique(long_refp_sort$V4))

# table(long_refp2$region)
clr <- c("#C9121E","#FC6910","#269321","#1B61A5")
methods <- as.character(unique(long_refp2$region)) # "extracellular small" "cellular small"      "cellular CLIP"  
#methods <- "exPeak"
for (i in 1:length(methods)){
  #i <- 1
  method <- methods[i]
  print(method)
  long_refp2_tmp <- long_refp2[long_refp2$region==method,]
  
  plot.list <- list()
  for (sample in c("RBPs","EV","G4")){ # as.character(unique(long_refp2_tmp$sample))
    # sample <- "RBPs"
    print(sample)
    long_refp2_tmp2 <- long_refp2_tmp[long_refp2_tmp$sample==sample,]
    long_refp_sort <- long_refp2_tmp2 %>%
      # dplyr::filter(pos %in% 25:35) %>% # only rank/sort by central regions
      dplyr::group_by(V4) %>%  # sample,
      #mean of each row/region
      dplyr::summarise(mean_signal_region = mean(signal,na.rm=T)
      )
    long_refp_sort <- long_refp_sort[order(long_refp_sort$mean_signal_region,decreasing = F),]
    long_refp2_tmp3 <- long_refp2_tmp2
    long_refp2_tmp3$V4 <- factor(long_refp2_tmp3$V4,levels = unique(long_refp_sort$V4))
    table(long_refp2_tmp3$sample)
    #table(is.na(long_refp2_tmp3$signal))
    table(long_refp2_tmp3$pos)
    tmp.p <- ggplot(long_refp2_tmp3, aes(x = pos,y = V4)) +
      geom_tile(aes(fill = signal)) + # 
      theme_minimal() +
      # coord_cartesian(expand = 0) +
      scale_x_continuous(breaks = c(1,up_bins+0.5,up_bins+target_bins+0.5,up_bins+target_bins+down_bins), # need +1 for 1-based coord
                         limits = c(1,up_bins+target_bins+down_bins),
                         labels = c(paste0("-",up,"nt"),'Start','End',paste0("+",down,"nt"))) +
      scale_fill_gradient(low = 'white',high = clr[i]) + # , na.value = "white")
      ylab('') + xlab('') +
      theme(
        aspect.ratio = 0.4,
        plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
        # plot.background = element_rect(fill = "white"),
        # panel.background = element_rect(fill = "white", colour = "white"),
        # strip.background = element_rect(color = NA,fill = 'white'),
        strip.text = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20,angle = 0,hjust = 0.5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position =  "none",#c(0.9,0.8),
        legend.text = element_text(size= 16),
        legend.title= element_text(size= 16)) #+
    # theme(aspect.ratio = 1.5,
    # plot.margin = unit(c(5.5, 5.5, 0, 5.5), "pt")
    #       panel.grid = element_blank(),
    #       axis.text.x = element_text(size=20,angle = 0,hjust = 0.5),
    #       # axis.text.y = element_text(size=20),
    #       legend.position = "right",
    #       legend.text = element_text(size= 16),
    #       legend.title= element_text(size= 16))+
    # facet_grid( . ~ sample ,  scales = "free_y" ) # region, space = "free_y",
    plot.list[[paste0(method,"_",sample)]] <- tmp.p
  }
  # p <- cowplot::plot_grid(plotlist = plot.list, rel_widths = c(1,1,1), nrow = 1, align = "hv", axis = "b")
  p <- ggarrange(plotlist = plot.list, rel_widths = c(1,1,1), nrow = 1,  common.legend = T,align = "hv", axis = "b", legend = "right") #  labels = c("RBPs", "", "EV", "", "G4"),
  ggsave(plot = p, paste0(method,"_cov_heat.pdf"), width=15, height=3) # "_",sample
}
# ggsave(paste0("cov_heat.pdf"), width=15, height=3) # "_",sample






# comapre Fig4 overlap  meta + WPS/COV plot sample-wise (enrichedHeatmap + R, newest, 230424) ------------------------------------------------
setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC")
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggsci)
library(Rmisc)
library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("rename", "dplyr")

## read ref
ref <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)
table(ref$transcript_type)


#tx
up=50
down=50
bin=5 # 10 for test
ratio=0.3 # 20/(50*2+20)

dst <- "GSE71008" # GSE71008,GSE110381,GSE123972
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
smps <- read.table(paste0(pre,"/exSeek-dev/data/",dst,"/sample_ids_NC_test15.txt"))$V1


## pre-filter boundary peak (sample-wise)
#smps <- smps[sample(1:length(smps),10,replace = F)]
for(i in 1:length(smps)){
  #smp <- "SAMN03863458"
  smp <- smps[i]
  print(smp)
  for(m in c("piranha_by_sample/b5_p01","clipper_by_sample/b5_p05","clam_by_sample/b5_p005","expeakCNN_by_sample/b5_d50_p1")){
    print(m)
    df <- read.table(paste0("output/",dst,"/call_peak_all/",m,"/intersect/",smp,".bed6"))
    
    df$txLen <- ref$tx.length[match(df$V1,ref$transcript_id)]
    #peakLen <- 20
    flankLen <- 50
    df <- df[(df$V2-flankLen)>=0 & (df$V3+flankLen)<=df$txLen,] 
    write.table(df,paste0("output/",dst,"/call_peak_all/",m,"/intersect/",smp,"_rmBoundary.bed"),quote = F,sep = "\t",row.names = F,col.names = F) # _noBoundary
  }
}



# cal meta-cov table
#use _noBoundary.bed6 seem not good
#smps <- smps[6:15]
for(i in 1:length(smps)){
  smp <- smps[i]
  print(smp)
  cmb <- expand.grid(signal=c("/BioII/lulab_b/baopengfei/shared_reference/RBP/splitByRBP/all_tx_score1.bed",
                              "/BioII/lulab_b/baopengfei/shared_reference/structure/quadratlas_g4grinder_tx_score1.bed",
                              "/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/lulab/11RNAsmallDomain_EVenrich_score1.bed"),
                     region=c(paste0("output/",dst,"/call_peak_all/piranha_by_sample/b5_p01/intersect/",smp,"_rmBoundary.bed"), # _noBoundary
                              paste0("output/",dst,"/call_peak_all/clipper_by_sample/b5_p05/intersect/",smp,"_rmBoundary.bed"), # _noBoundary
                              paste0("output/",dst,"/call_peak_all/clam_by_sample/b5_p005/intersect/",smp,"_rmBoundary.bed"), # _noBoundary
                              paste0("output/",dst,"/call_peak_all/expeakCNN_by_sample/b5_d50_p1/intersect/",smp,"_rmBoundary.bed")) # _noBoundary
  )
  
  
  cmb$signal.label <- ""
  cmb$signal.label[grepl("all_tx",cmb$signal)] <- c("RBPs") #,"small","CLIP","ex RIP"
  cmb$signal.label[grepl("g4grinder",cmb$signal)] <- c("G4") #,"small","CLIP","ex RIP"
  cmb$signal.label[grepl("EVenrich",cmb$signal)] <- c("EV") #,"small","CLIP","ex RIP"
  # table(cmb$signal.label)
  cmb$region.label <- "Piranha"
  cmb$region.label[grepl("clam",tolower(cmb$region))] <- "CLAM"
  cmb$region.label[grepl("clipper",tolower(cmb$region))] <- "CLIPper"
  cmb$region.label[grepl("expeak",tolower(cmb$region))] <- "exPeak"
  
  # need prefilter boundary peaks (EnrichedHeatmap::normalizeToMatrix seem not take chr.size info when extend region)
  
  res.list <- lapply(1:nrow(cmb), function(j) {get.mat(signal=as.character(cmb$signal[j]), region=as.character(cmb$region[j]),
                                                       signal.label=as.character(cmb$signal.label[j]),region.label=as.character(cmb$region.label[j]),
                                                       up = up, down = down, bin = bin, ratio = ratio )} )
  res.list.df <- do.call(rbind,res.list)
  # write.table(res.list.df,paste0("tmp/comapre_Fig5_overlap_meta_",smp,".txt2"),row.names = T,col.names = T,sep = "\t",quote = F)
  write.table(res.list.df,paste0("tmp/comapre_Fig5_overlap_meta_",smp,".txt3"),row.names = T,col.names = T,sep = "\t",quote = F)
  
  # res.list.df[1:3,]
}



## read meta-cov table
l <- list()
# smps <- Sys.glob("tmp/comapre_Fig5_overlap_meta_SAMN*.txt")
# smps <- basename(smps)
# smps <- gsub("comapre_Fig5_overlap_meta","",smps)
# smps <- gsub(".txt","",smps,fixed = T)
for  (smp in smps){
  print(smp)
  # tmp <- read.table(paste0("tmp/comapre_Fig5_overlap_meta_",smp,".txt2"),sep = "\t",header = T)
  tmp <- read.table(paste0("tmp/comapre_Fig5_overlap_meta_",smp,".txt3"),sep = "\t",header = T)
  tmp$smp <- smp
  print(dim(tmp))
  l[[smp]] <- tmp
}
res.list.df <- do.call(rbind,l)
# table(res.list.df$signal)
# for(i in 1:length(l)){
#   print(i);
#   print(dim(l[[i]]))
# }


# plot in R
res.list.df$signal <- gsub("G4iM","G4",res.list.df$signal)
res.list.df$signal <- factor(res.list.df$signal,levels = c("RBPs","EV","G4"))
res.list.df$region <- factor(res.list.df$region,levels = c("Piranha","CLIPper","CLAM","exPeak"))
# table(res.list.df$region)


#(500+500)/20+10=60,再乘以1个样本,再乘以1个区域,为60列
clean_refp <- dplyr::as_tibble(res.list.df)
#clean_refp$smp <- NULL
clean_refp <- clean_refp[,c(1:3,ncol(clean_refp),4:(ncol(clean_refp)-1))]
#clean_refp <- refp #%>% select(c(-1:-3,-5,-6))
colnames(clean_refp)[c(1,3)] <- c('V4','sample')
#plot(x=1:(ncol(clean_refp)-3), y=as.numeric(clean_refp[1,4:ncol(clean_refp)]))
table(clean_refp$region,clean_refp$sample)

clean_refp <- clean_refp[which( apply(clean_refp[,5:ncol(clean_refp)] ,1, sd)>0),] # filter low sd
#~ half num of peaks have no ref overlap: sd==0

#1st max-min scale
clean_refp[,5:ncol(clean_refp)] <- t(maxmin.normalize(t(clean_refp[,5:ncol(clean_refp)])))
#clean_refp.tmp <- clean_refp[clean_refp$sample=="G4iM",]
#summary(as.numeric(clean_refp.tmp[3,4:ncol(clean_refp.tmp)]))

# clean_refp[1:6,1:6]
long_refp <- reshape2::melt(clean_refp,id.vars = c('V4','sample','region','smp'),value.name = 'signal') %>%
  select(-variable)
long_refp[1:3,]
#summary(long_refp$signal)

#add x position
long_refp$pos <- rep(c(1:(ncol(clean_refp)-4)),each = nrow(clean_refp)) # ,times = nrow(cmb)
# dim(long_refp)


# run loess model 
# job::job({
l <- list()
for(sample in as.character(unique(long_refp$sample))){
  #sample <- "RBPs"
  print(sample)
  for (region in as.character(unique(long_refp$region))){
    #region <- "exPeak"
    print(region)
    for (smp in as.character(unique(long_refp$smp))){
      # smp <- "SAMN03863400"
      print(smp)
      # print(Sys.time())
      long_refp_tmp <- long_refp[long_refp$sample==sample & long_refp$region==region & long_refp$smp==smp,] # not merge smps
      #plot(long_refp_tmp$pos, long_refp_tmp$signal)
      # long_refp_tmp <- get.loess(long_refp_tmp) #[,c("pos","signal")] # get.lowess, get.loess
      long_refp_tmp <- as.data.frame(get.loess.smooth(long_refp_tmp)) #[,c("pos","signal")] # get.lowess, get.loess
      # long_refp_tmp$signal[long_refp_tmp$signal < 0] <- 0
      # long_refp_tmp$signal[long_refp_tmp$signal > 1] <- 1
      long_refp_tmp$smp <- smp
      long_refp_tmp$sample <- sample
      long_refp_tmp$region <- region
      colnames(long_refp_tmp)[colnames(long_refp_tmp)=="signal"] <- "mean_signal"
      l[[paste0(sample,region,smp)]] <- long_refp_tmp
      #table(long_refp_tmp$pos==long_refp_tmp2$pos)
      # write.table(long_refp_tmp,paste0("tmp/comapre_Fig5_overlap_meta_final_scaler_",sample,"_",region,".txt2"),row.names = F,col.names = T,sep = "\t",quote = F)
      write.table(long_refp_tmp,paste0("tmp/comapre_Fig5_overlap_meta_final_scaler_",sample,"_",region,"_",smp,".txt3"),row.names = F,col.names = T,sep = "\t",quote = F)
    }
  }
}
# })


## read in lowess
l <- list()
for(sample in as.character(unique(long_refp$sample))){
  #sample <- "RBPs"
  print(sample)
  for (region in as.character(unique(long_refp$region))){
    #region <- "exPeak"
    print(region)
    
    for (smp in as.character(unique(long_refp$smp))){
      # smp <- "SAMN03863400"
      print(smp)
      # tmp <- read.table( paste0("tmp/comapre_Fig5_overlap_meta_final_scaler_",sample,"_",region,".txt2"), sep = "\t",header = T)
      tmp <- read.table( paste0("tmp/comapre_Fig5_overlap_meta_final_scaler_",sample,"_",region,"_",smp,".txt3"), sep = "\t",header = T)
      # tmp$smp <- smp
      l[[paste0(sample,region,smp)]] <- tmp
    }
  }
}
filnal_scaler <- do.call(rbind,l)
table(filnal_scaler$sample)
table(filnal_scaler$region)
table(filnal_scaler$smp)
summary(filnal_scaler$mean_signal)


up_bins <- as.integer(up/bin) # ceiling
down_bins <- as.integer(down/bin) # ceiling
target_bins <- round((up_bins + down_bins)*(ratio/(1-ratio)),digits = 0) #as.integer((up_bins + down_bins)*(ratio/(1-ratio)))
table(filnal_scaler$sample)
table(filnal_scaler$region)
filnal_scaler$region <- factor(filnal_scaler$region,levels = c("Piranha","CLIPper","CLAM","exPeak"))

#2nd max-min scale (only display mean line, no CI/ribbon)
filnal_scaler <- filnal_scaler %>% 
  dplyr::group_by(sample,region) %>% # add sample ?
  dplyr::mutate(mean_signal_norm= (mean_signal - min(mean_signal))/(max(mean_signal)-min(mean_signal))^as.logical(sd(mean_signal)) )

#2nd loess (if no CI/ribbon implemented)
#add CI/ribbon
filnal_scaler <- filnal_scaler %>%
  tidyr::drop_na() %>% # remove na
  dplyr::group_by(pos,sample,region) %>% # sample,
  # mean and 95% interval confidence
  dplyr::summarise(mean_signal_norm2 = mean(mean_signal_norm, na.rm=T), #,trim = 0.05, na.rm=T
                   sd = sd(mean_signal_norm),
                   upper = CI(mean_signal_norm,ci = 0.95)[1],
                   lower = CI(mean_signal_norm,ci = 0.95)[3],
                   max = max(mean_signal_norm),
                   min = min(mean_signal_norm)
  )


filnal_scaler2 <- filnal_scaler[filnal_scaler$region %in% "exPeak",]
filnal_scaler2$sample <- factor(filnal_scaler2$sample, levels = c("RBPs","EV","G4"))
ggplot(filnal_scaler2,aes(x = pos,y = mean_signal_norm2)) +
  # add 0.95 interval
  geom_ribbon(aes(ymin = lower,
                  ymax = upper,
                  fill = "salmon"), # sample
              alpha = 0.5) +
  geom_line(aes(color = region), size = 2) + # 
  theme_classic(base_size = 16) +
  scale_color_manual(values = pal_d3()(4)[4]) + 
  # scale_color_d3(name = 'Data type') +
  # scale_fill_d3(name = '') +
  geom_vline(xintercept = c(up_bins+1,up_bins+target_bins),color="grey50",linetype="dashed")+
  # geom_vline(xintercept = c((1+up_bins+target_bins+down_bins)*0.5),color="grey50",linetype="dashed")+ # need +1 for 1-based coord
  # x label
  scale_x_continuous(breaks = c(1,up_bins+1,up_bins+target_bins,up_bins+target_bins+down_bins), # need +1 for 1-based coord
                     labels = c(paste0("-",up,"nt"),'Start','End',paste0("+",down,"nt")) ) +
  scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1.0), limits = c(0,1), labels = c("0","0.25","0.50","0.75","1.00")) +
  xlab('') + ylab('Scaled depth') +
  theme(aspect.ratio = 0.6,
        strip.background = element_rect(color = NA,fill = 'white'),
        strip.text = element_text(size=30),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=30),
        axis.text.x = element_blank(), #  element_text(size=20,angle = 0,hjust = 0.5)
        axis.text.y = element_text(size=30),
        legend.position =  "right",#c(0.9,0.8),
        legend.text = element_text(size= 24),
        legend.title= element_text(size= 30) ) + 
  facet_grid(#.~sample, 
    cols = vars(sample), 
    scales = "free_y")
ggsave("cov_meta1_v2.pdf",width = 20,height = 8)

filnal_scaler3 <- filnal_scaler[!(filnal_scaler$region %in% "exPeak"),]
table(filnal_scaler3$sample,filnal_scaler3$region)
filnal_scaler3$sample <- factor(filnal_scaler3$sample, levels = c("RBPs","EV","G4"))
ggplot(filnal_scaler3,aes(x = pos,y = mean_signal_norm2)) +
  # add 0.95 interval
  geom_ribbon(aes(ymin = lower,
                  ymax = upper,
                  fill = region), # sample
              alpha = 0.5) +
  geom_line(aes(color = region), size = 2) + # 
  theme_classic(base_size = 16) +
  scale_color_d3(name = 'Data type') +
  scale_fill_d3(name = '') +
  geom_vline(xintercept = c(up_bins+1,up_bins+target_bins),color="grey50",linetype="dashed")+
  # geom_vline(xintercept = c((1+up_bins+target_bins+down_bins)*0.5),color="grey50",linetype="dashed")+ # need +1 for 1-based coord
  # x label
  scale_x_continuous(breaks = c(1,up_bins+1,up_bins+target_bins,up_bins+target_bins+down_bins), # need +1 for 1-based coord
                     labels = c(paste0("-",up,"nt"),'Start','End',paste0("+",down,"nt")) ) +
  scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1.0), limits = c(0,1), labels = c("0","0.25","0.50","0.75","1.00")) +
  xlab('') + ylab('Scaled depth') +
  theme(aspect.ratio = 0.6,
        strip.background = element_rect(color = NA,fill = 'white'),
        strip.text = element_text(size=30),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=30),
        axis.text.x = element_blank(), #  element_text(size=20,angle = 0,hjust = 0.5)
        axis.text.y = element_text(size=30),
        legend.position =  "right",#c(0.9,0.8),
        legend.text = element_text(size= 24),
        legend.title= element_text(size= 30)) + 
  facet_grid(#.~sample, 
    cols = vars(sample), 
    scales = "free_y")
ggsave("cov_meta2_v2.pdf",width = 20,height = 8)



# heatmap (need further optimization?)
# # plot (slowly)
# # long_refp$log10.signal <- log10(long_refp$signal+0.00001)
# long_refp_sort <- long_refp %>%
#   # dplyr::filter(pos %in% 25:35) %>% # only rank/sort by central regions
#   dplyr::group_by(V4,sample) %>%  # sample,
#   #mean of each row/region
#   dplyr::summarise(mean_signal_region = mean(signal,na.rm=T)
#   )
# long_refp_sort <- long_refp_sort[order(long_refp_sort$mean_signal_region,decreasing = F),]
long_refp2 <- long_refp[long_refp$region=="exPeak",]
# long_refp2$V4 <- factor(long_refp2$V4,levels = unique(long_refp_sort$V4))

# table(long_refp2$region)
clr <- c("#C9121E","#FC6910","#269321","#1B61A5")
methods <- as.character(unique(long_refp2$region)) # 
#methods <- "exPeak"
for (i in 1:length(methods)){
  #i <- 1
  method <- methods[i]
  print(method)
  long_refp2_tmp <- long_refp2[long_refp2$region==method,]
  
  plot.list <- list()
  for (sample in c("RBPs","EV","G4")){ # as.character(unique(long_refp2_tmp$sample))
    # sample <- "RBPs"
    print(sample)
    long_refp2_tmp2 <- long_refp2_tmp[long_refp2_tmp$sample==sample,]
    long_refp_sort <- long_refp2_tmp2 %>%
      # dplyr::filter(pos %in% 25:35) %>% # only rank/sort by central regions
      dplyr::group_by(V4) %>%  # sample,
      #mean of each row/region
      dplyr::summarise(mean_signal_region = mean(signal,na.rm=T)
      )
    long_refp_sort <- long_refp_sort[order(long_refp_sort$mean_signal_region,decreasing = F),]
    long_refp2_tmp3 <- long_refp2_tmp2
    long_refp2_tmp3$V4 <- factor(long_refp2_tmp3$V4,levels = unique(long_refp_sort$V4))
    table(long_refp2_tmp3$sample)
    #table(is.na(long_refp2_tmp3$signal))
    table(long_refp2_tmp3$pos)
    tmp.p <- ggplot(long_refp2_tmp3, aes(x = pos,y = V4)) +
      geom_tile(aes(fill = signal)) + # 
      theme_minimal() +
      # coord_cartesian(expand = 0) +
      scale_x_continuous(breaks = c(1,up_bins+0.5,up_bins+target_bins+0.5,up_bins+target_bins+down_bins), # need +1 for 1-based coord
                         limits = c(1,up_bins+target_bins+down_bins),
                         labels = c(paste0("-",up,"nt"),'Start','End',paste0("+",down,"nt"))) +
      scale_fill_gradient(low = 'white',high = clr[i]) + # , na.value = "white")
      ylab('') + xlab('') +
      theme(
        aspect.ratio = 0.4,
        plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
        # plot.background = element_rect(fill = "white"),
        # panel.background = element_rect(fill = "white", colour = "white"),
        # strip.background = element_rect(color = NA,fill = 'white'),
        strip.text = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20,angle = 0,hjust = 0.5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position =  "none",#c(0.9,0.8),
        legend.text = element_text(size= 16),
        legend.title= element_text(size= 16)) #+
    # theme(aspect.ratio = 1.5,
    # plot.margin = unit(c(5.5, 5.5, 0, 5.5), "pt")
    #       panel.grid = element_blank(),
    #       axis.text.x = element_text(size=20,angle = 0,hjust = 0.5),
    #       # axis.text.y = element_text(size=20),
    #       legend.position = "right",
    #       legend.text = element_text(size= 16),
    #       legend.title= element_text(size= 16))+
    # facet_grid( . ~ sample ,  scales = "free_y" ) # region, space = "free_y",
    plot.list[[paste0(method,"_",sample)]] <- tmp.p
  }
  # p <- cowplot::plot_grid(plotlist = plot.list, rel_widths = c(1,1,1), nrow = 1, align = "hv", axis = "b")
  p <- ggpubr::ggarrange(plotlist = plot.list, rel_widths = c(1,1,1), nrow = 1,  common.legend = T,align = "hv", axis = "b", legend = "right") #  labels = c("RBPs", "", "EV", "", "G4"),
  ggsave(plot = p, paste0(method,"_cov_heat_v2.pdf"), width=15, height=3) # "_",sample
}
# ggsave(paste0("cov_heat.pdf"), width=15, height=3) # "_",sample

long_refp2 <- long_refp[!(long_refp$region=="exPeak"),]
# long_refp2$V4 <- factor(long_refp2$V4,levels = unique(long_refp_sort$V4))

# table(long_refp2$region)
#clr <- c("#FC6910","#269321","#1B61A5") # 
clr <- c("#1B61A5","#FC6910","#269321")

#plot(x=1:3,y=1:3,col=clr)
methods <- as.character(unique(long_refp2$region)) # 
#methods <- "exPeak"
for (i in 1:length(methods)){
  #i <- 1
  method <- methods[i]
  print(method)
  long_refp2_tmp <- long_refp2[long_refp2$region==method,]
  
  plot.list <- list()
  for (sample in c("RBPs","EV","G4")){ # as.character(unique(long_refp2_tmp$sample))
    # sample <- "RBPs"
    print(sample)
    long_refp2_tmp2 <- long_refp2_tmp[long_refp2_tmp$sample==sample,]
    long_refp_sort <- long_refp2_tmp2 %>%
      # dplyr::filter(pos %in% 25:35) %>% # only rank/sort by central regions
      dplyr::group_by(V4) %>%  # sample,
      #mean of each row/region
      dplyr::summarise(mean_signal_region = mean(signal,na.rm=T)
      )
    long_refp_sort <- long_refp_sort[order(long_refp_sort$mean_signal_region,decreasing = F),]
    long_refp2_tmp3 <- long_refp2_tmp2
    long_refp2_tmp3$V4 <- factor(long_refp2_tmp3$V4,levels = unique(long_refp_sort$V4))
    table(long_refp2_tmp3$sample)
    #table(is.na(long_refp2_tmp3$signal))
    table(long_refp2_tmp3$pos)
    tmp.p <- ggplot(long_refp2_tmp3, aes(x = pos,y = V4)) +
      geom_tile(aes(fill = signal)) + # 
      theme_minimal() +
      # coord_cartesian(expand = 0) +
      scale_x_continuous(breaks = c(1,up_bins+0.5,up_bins+target_bins+0.5,up_bins+target_bins+down_bins), # need +1 for 1-based coord
                         limits = c(1,up_bins+target_bins+down_bins),
                         labels = c(paste0("-",up,"nt"),'Start','End',paste0("+",down,"nt"))) +
      scale_fill_gradient(low = 'white',high = clr[i]) + # , na.value = "white")
      ylab('') + xlab('') +
      theme(
        aspect.ratio = 0.4,
        plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
        # plot.background = element_rect(fill = "white"),
        # panel.background = element_rect(fill = "white", colour = "white"),
        # strip.background = element_rect(color = NA,fill = 'white'),
        strip.text = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20,angle = 0,hjust = 0.5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position =  "none",#c(0.9,0.8),
        legend.text = element_text(size= 16),
        legend.title= element_text(size= 16)) #+
    # theme(aspect.ratio = 1.5,
    # plot.margin = unit(c(5.5, 5.5, 0, 5.5), "pt")
    #       panel.grid = element_blank(),
    #       axis.text.x = element_text(size=20,angle = 0,hjust = 0.5),
    #       # axis.text.y = element_text(size=20),
    #       legend.position = "right",
    #       legend.text = element_text(size= 16),
    #       legend.title= element_text(size= 16))+
    # facet_grid( . ~ sample ,  scales = "free_y" ) # region, space = "free_y",
    plot.list[[paste0(method,"_",sample)]] <- tmp.p
  }
  # p <- cowplot::plot_grid(plotlist = plot.list, rel_widths = c(1,1,1), nrow = 1, align = "hv", axis = "b")
  p <- ggarrange(plotlist = plot.list, rel_widths = c(1,1,1), nrow = 1,  common.legend = T,align = "hv", axis = "b", legend = "right") #  labels = c("RBPs", "", "EV", "", "G4"),
  ggsave(plot = p, paste0(method,"_cov_heat_v2.pdf"), width=15, height=3) # "_",sample
}







# WPS
setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC")
dst.expeak <- "GSE71008_NCpool"
dst <- "GSE71008_NCpool_min30_prot15_full20"  #GSE71008_min25_prot15_full20, "GSE71008_NCpool_min15", "GSE71008_NCpool_min30_prot15_full20"
#smps <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/",dst.expeak,"/sample_ids_NC_test5.txt"))$V1
smps <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/",dst.expeak,"/sample_ids.txt"))$V1 # NCpool
pre <- paste0("/BioII/lulab_b/baopengfei/gitsoft/cfDNA-2016/results/intermediate/",dst)

#lowess model
# for(smp in smps){
#
# #smp <- "SAMN03863458"
#   print(smp)
#   print(Sys.time())
for (type in c("COV", "WPS", "STARTS")){
  #type <- "COV"
  print(type)
  # print(Sys.time())
  for (method in c("exPeak","CLIPper","CLAM","Piranha")){
    #method <- "exPeak"
    print(method)
    # print(Sys.time())
    
    l0 <- list()
    for(smp in smps){
      # print(smp)
      wps <- data.table::fread(paste0(pre,"/table/GRCh38/target/",method,"--",smp,"_",type,".csv.gz")) # CLAM sample-wise raw peak contain ","
      wps$smp <- smp
      wps$region <- method
      wps$sample <- type
      l0[[type]] <- wps
    }
    wps <- do.call(rbind,l0)
    
    clean_refp <- dplyr::as_tibble(wps)
    clean_refp <- clean_refp[,c(1,ncol(clean_refp),ncol(clean_refp)-1,ncol(clean_refp)-2,2:(ncol(clean_refp)-3))]
    #clean_refp <- refp #%>% select(c(-1:-3,-5,-6))
    colnames(clean_refp)[c(1:4)] <- c('V4','sample','region','smp') # sample:RBP,EV; region: expeak,clipper; smp: NCpool
    colnames(clean_refp)[c(5:ncol(clean_refp))] <- paste0("pos_",1:(ncol(clean_refp)-4))
    #plot(x=1:(ncol(clean_refp)-3), y=as.numeric(clean_refp[1,4:ncol(clean_refp)]))
    # table(clean_refp$region,clean_refp$sample)
    # clean_refp[1:3,1:7]
    clean_refp <- clean_refp[which( apply(clean_refp[,5:ncol(clean_refp)] ,1, sd)>0),] # filter low sd
    #~ half num of peaks have no ref overlap: sd==0
    
    #max-min scale
    clean_refp[,5:ncol(clean_refp)] <- t(maxmin.normalize(t(clean_refp[,5:ncol(clean_refp)])))
    #summary(wps[,2:(ncol(wps)-4)] )
    
    long_refp <- reshape2::melt(clean_refp,id.vars = c('V4','sample','region','smp'),value.name = 'signal') %>%
      dplyr::select(-variable)
    
    # add x position
    long_refp$pos <- rep(c(1:(ncol(clean_refp)-4)),each = nrow(clean_refp))
    # not smooth to keep diff peaks/records as CI/ribbon
    # long_refp$signal <- get.loess(long_refp)$signal # [,c("pos","signal")] # lowess seem poor curve
    # long_refp$signal[long_refp$signal < 0] <- 0
    # long_refp$signal[long_refp$signal > 1] <- 1
    colnames(long_refp)[colnames(long_refp)=="signal"] <- "mean_signal"
    # write.table(long_refp,paste0(pre,"/",dst.expeak,"_",method,"_",type,".txt"),quote = F,sep = "\t",row.names = F,col.names = T)
    write.table(long_refp,paste0(pre,"/",dst.expeak,"_",method,"_",type,".txt2"),quote = F,sep = "\t",row.names = F,col.names = T)
  }
}


# read in
tmpList <- list()
for(type in c("COV", "WPS", "STARTS")){
  #type <- "WPS"
  print(type)
  for (method in c("exPeak","CLIPper","CLAM","Piranha")){
    #method <- "exPeak"
    print(method)
    # tmp <- data.table::fread(paste0(pre,"/",dst.expeak,"_",method,"_",type,".txt"),check.names = F,sep = "\t",header = T)
    tmp <- data.table::fread(paste0(pre,"/",dst.expeak,"_",method,"_",type,".txt2"),check.names = F,sep = "\t",header = T)
    
    # #rm dup (peak not needed, though diff)
    # tmp$V4 <- NULL
    # tmp <- tmp[!duplicated(paste(tmp$sample,tmp$region,tmp$smp,tmp$mean_signal,tmp$pos)),]
    
    tmpList[[paste0(type,"_",method)]] <- tmp
  }
}
filnal_scaler <- do.call(rbind,tmpList)
filnal_scaler <- filnal_scaler[!(filnal_scaler$sample %in% c("STARTS","WPS_v2")),] #
filnal_scaler[1:2,]


up <- 50
down <- 50
up_bins <- 50 # ceiling
down_bins <- 50 # ceiling
target_bins <- 20 #as.integer((up_bins + down_bins)*(ratio/(1-ratio)))
table(filnal_scaler$sample)
table(filnal_scaler$region,filnal_scaler$sample)
#filnal_scaler$region <- factor(filnal_scaler$region,levels = c("Piranha","CLIPper","CLAM","exPeak"))
#table((filnal_scaler$V4))
#summary(filnal_scaler$mean_signal)


#2nd max-min scale (only display mean line, no CI/ribbon)
table(filnal_scaler$smp)
filnal_scaler <- filnal_scaler %>%
  dplyr::group_by(V4,sample,region,smp) %>% # for NCpool, 
  dplyr::mutate(mean_signal_norm= (mean_signal - min(mean_signal))/(max(mean_signal)-min(mean_signal))^as.logical(sd(mean_signal)) )

#add CI/ribbon
filnal_scaler <- filnal_scaler %>%
  tidyr::drop_na() %>% # remove na
  dplyr::group_by(pos,sample,smp,region) %>% # V4/peaks/records shown as CI/ribbon
  # mean and 95% interval confidence
  dplyr::summarise(mean_signal_norm2 = mean(mean_signal_norm, na.rm=T), #,trim = 0.05, na.rm=T
                   sd = sd(mean_signal_norm),
                   upper = CI(mean_signal_norm,ci = 0.99)[1],
                   lower = CI(mean_signal_norm,ci = 0.99)[3],
                   max = max(mean_signal_norm),
                   min = min(mean_signal_norm)
  )
# #2nd loess (loess after max-min-scale to aviod ridges-line)
# filnal_scaler <- filnal_scaler %>%
#   dplyr::group_by(sample,region,smp) %>%
#   dplyr::mutate(mean_signal_norm2 = get.loess.vec(pos,mean_signal_norm))
# filnal_scaler$mean_signal_norm2[filnal_scaler$mean_signal_norm2>1] <- 1
# filnal_scaler$mean_signal_norm2[filnal_scaler$mean_signal_norm2<0] <- 0

#expeak
filnal_scaler2 <- filnal_scaler[filnal_scaler$region == "exPeak",]
filnal_scaler2$sample <- factor(filnal_scaler2$sample,levels = c("COV","WPS")) # WPS_v2
ggplot(filnal_scaler2,aes(x = pos,y = mean_signal_norm2)) + # mean_signal, mean_signal_norm
  # add 0.95 interval
  geom_ribbon(aes(ymin = lower,
                  ymax = upper,
                  fill = sample), # sample
              alpha = 0.5) +
  geom_line(aes(color = sample), size = 2) + #
  scale_color_manual(values = c("salmon","chocolate")) +
  scale_fill_manual(values = alpha(c("salmon","chocolate"),0.5)) +
  # scale_color_d3(name = 'Data type') +
  # scale_fill_d3(name = '') +
  geom_vline(xintercept = c(up_bins+target_bins*0.5+1),color="grey50",linetype="dashed")+
  geom_hline(yintercept = 0,color="grey30",linetype="solid")+
  # geom_vline(xintercept = c((1+up_bins+target_bins+down_bins)*0.5),color="grey50",linetype="dashed")+ # need +1 for 1-based coord
  # x label
  scale_x_continuous(breaks = c(1,up_bins+target_bins*0.5+1,up_bins+target_bins+down_bins+1), # need +1 for 1-based coord
                     labels = c(paste0("-",up,""),'Center',paste0("+",down,"")) ) +
  scale_y_continuous(breaks = c(0,1.0), limits = c(0,1), labels = c("0","1")) +
  xlab('') + ylab('Scaled value') +
  # facet_grid(sample~.,scales = "free_y") +  #,rows = 4
  theme_classic(base_size = 16) +
  theme(aspect.ratio = 0.6,
        strip.background = element_rect(color = NA,fill = 'white'),
        strip.text = element_text(size=30),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=30),
        axis.text.x = element_text(size=30), #  element_text(size=20,angle = 0,hjust = 0.5)
        axis.text.y = element_text(size=30),
        legend.position =  "right",#c(0.9,0.8),
        legend.text = element_text(size= 24),
        legend.title= element_text(size= 30))
ggsave("wps_meta_expeak_v2.pdf",width = 9,height = 9)


#other 3 methods
filnal_scaler3 <- filnal_scaler[filnal_scaler$region != "exPeak",]
filnal_scaler3$sample <- factor(filnal_scaler3$sample,levels = c("COV","WPS")) 
ggplot(filnal_scaler3,aes(x = pos,y = mean_signal_norm2)) + # mean_signal, mean_signal_norm
  # add 0.95 interval
  geom_ribbon(aes(ymin = lower,
                  ymax = upper,
                  fill = region), # sample
              alpha = 0.5) +
  geom_line(aes(color = region), size = 2) + #
  scale_color_manual(values = pal_d3()(4)[c(3,2,1)]) +
  scale_fill_manual(values = alpha(pal_d3()(4)[c(3,2,1)],0.5)) +
  # scale_color_d3(name = 'Data type') +
  # scale_fill_d3(name = '') +
  geom_vline(xintercept = c(up_bins+target_bins*0.5+1),color="grey50",linetype="dashed")+
  geom_hline(yintercept = 0,color="grey30",linetype="solid")+
  # geom_vline(xintercept = c((1+up_bins+target_bins+down_bins)*0.5),color="grey50",linetype="dashed")+ # need +1 for 1-based coord
  # x label
  scale_x_continuous(breaks = c(1,up_bins+target_bins*0.5+1,up_bins+target_bins+down_bins+1), # need +1 for 1-based coord
                     labels = c(paste0("-",up,""),'Center',paste0("+",down,"")) ) +
  scale_y_continuous(breaks = c(0,1.0), limits = c(-0.1,1.1), labels = c("0","1")) +
  xlab('') + ylab('Scaled value') +
  facet_grid(sample~.,scales = "free_y") +  #,rows = 4
  theme_classic(base_size = 16) +
  theme(aspect.ratio = 0.6,
        strip.background = element_rect(color = NA,fill = 'white'),
        strip.text = element_text(size=30),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=30),
        axis.text.x = element_text(size=30), #  element_text(size=20,angle = 0,hjust = 0.5)
        axis.text.y = element_text(size=30),
        legend.position =  "right",#c(0.9,0.8),
        legend.text = element_text(size= 24),
        legend.title= element_text(size= 30))
ggsave("wps_meta_3methods_v2.pdf",width = 12,height = 12)


#






# comapre Fig4 overlap  meta plot sample-wise (i-pico, 230624) ------------------------------------------------
setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC")
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggsci)
library(Rmisc)
library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("rename", "dplyr")

## read ref
ref <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)
table(ref$transcript_type)


#tx
up=50
down=50
bin=5 # 10 for test
ratio=0.3 # 20/(50*2+20)

dst <- "i-pico" # GSE71008,GSE110381,GSE123972
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
smps <- read.table(paste0(pre,"/exSeek-dev/data/",dst,"/sample_ids.txt"))$V1

#/usr/bin/Rscript /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/scripts/gbed2tbed.R -i /BioII/lulab_b/baopengfei/shared_reference/RBP/splitByRBP/APOBEC3C.bed -o /BioII/lulab_b/baopengfei/shared_reference/RBP/splitByRBP/APOBEC3C_tx.bed
#for i in DDX3X SRSF3 ABCF1; do bedtools sort -i /BioII/lulab_b/baopengfei/shared_reference/RBP/splitByRBP/${i}_tx.bed  | bedtools merge -s -o distinct,mean,distinct -c 4,5,6 -i stdin > /BioII/lulab_b/baopengfei/shared_reference/RBP/splitByRBP/${i}_tx_merge.bed; done
#for i in `cat data/i-pico/sample_ids.txt`; do zcat ../output/i-pico/call_peak_dedup/tbed_11RNA_primary/${i}.bed.gz  | awk '{{if(($3-$2)<65) {{print $0}}}}' > ../output/i-pico/call_peak_dedup/tbed_11RNA_primary/${i}_LenLT65.bed;done
for(rbp in c("DDX3X", "SRSF3", "ABCF1")){
  #rbp <- "DDX3X"
  print(rbp)
  df <- read.table(paste0("/BioII/lulab_b/baopengfei/shared_reference/RBP/splitByRBP/",rbp,"_tx_merge.bed"))
  df$txLen <- ref$tx.length[match(df$V1,ref$transcript_id)]
  #df[1:3,]
  df <- df[(df$V3-df$V2)>10 & (df$V3-df$V2)<100,]
  #peakLen <- 20
  flankLen <- 50
  df <- df[(df$V2-flankLen)>=0 & (df$V3+flankLen)<=df$txLen,]
  write.table(df,paste0("/BioII/lulab_b/baopengfei/shared_reference/RBP/splitByRBP/",rbp,"_tx_merge_rmBoundary.bed"),quote = F,sep = "\t",row.names = F,col.names = F) # _noBoundary
}
#for i in DDX3X SRSF3 ABCF1 ; do /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/scripts/shufPeakBed.R -i /BioII/lulab_b/baopengfei/shared_reference/RBP/splitByRBP/${i}_tx_merge_rmBoundary.bed; done


# cal meta-cov table
#use _noBoundary.bed6 seem not good
#smps <- smps[6:15]
# rbp <- "DDX3X" # DDX3X, SRSF3, ABCF1
for(rbp in c("DDX3X",  "ABCF1")){ # "SRSF3",
for(i in 1:length(smps)){
  smp <- smps[i]
  print(smp)
  cmb <- expand.grid(signal=c(paste0("output/",dst,"/call_peak_dedup/tbed_11RNA_primary/",smp,"_LenLT45.bed"),
                              paste0("output/",dst,"/call_peak_dedup/tbed_11RNA_primary/",smp,"_LenLT65.bed"),
                              paste0("output/",dst,"/call_peak_dedup/tbed_11RNA_primary/",smp,"_LenLT85.bed")), # .bed.gz  _LenLT45.bed
                     region=c(paste0("/BioII/lulab_b/baopengfei/shared_reference/RBP/splitByRBP/",rbp,"_tx_merge_rmBoundary.bed"), 
                              paste0("/BioII/lulab_b/baopengfei/shared_reference/RBP/splitByRBP/",rbp,"_tx_merge_rmBoundary.bed.shuffle")
                              # "/BioII/lulab_b/baopengfei/shared_reference/RBP/splitByRBP/SRSF3_tx.bed",
                              # "/BioII/lulab_b/baopengfei/shared_reference/RBP/splitByRBP/ABCF1_tx.bed"
                              ) 
  )
  cmb$signal.label <- "bed"
  cmb$signal.label[grepl("45",cmb$signal)] <- c("45") #,"small","CLIP","ex RIP"
  cmb$signal.label[grepl("65",cmb$signal)] <- c("65") #,"small","CLIP","ex RIP"
  cmb$signal.label[grepl("85",cmb$signal)] <- c("85") #,"small","CLIP","ex RIP"
  # table(cmb$signal.label)
  cmb$region.label <- ""
  cmb$region.label[grepl(rbp,(cmb$region))] <- "RBP"
  cmb$region.label[grepl("shuffle",(cmb$region))] <- "shuffle"
  # cmb$region.label[grepl("ABCF1",(cmb$region))] <- "ABCF1"
  
  # need prefilter boundary peaks (EnrichedHeatmap::normalizeToMatrix seem not take chr.size info when extend region)
  
  res.list <- lapply(1:nrow(cmb), function(j) {get.mat(signal=as.character(cmb$signal[j]), region=as.character(cmb$region[j]),
                                                       signal.label=as.character(cmb$signal.label[j]),region.label=as.character(cmb$region.label[j]),
                                                       up = up, down = down, bin = bin, ratio = ratio )} )
  res.list.df <- do.call(rbind,res.list)
  # write.table(res.list.df,paste0("tmp/comapre_Fig5_overlap_meta_",smp,".txt2"),row.names = T,col.names = T,sep = "\t",quote = F)
  write.table(res.list.df,paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/i-pico/diffgene-CFEV_ExonUTR_from_longestTx/i-pico_comapre_Fig5_overlap_meta_",smp,"_",rbp,".txt4"),row.names = T,col.names = T,sep = "\t",quote = F)
}
}


rbp <- "DDX3X"
## read meta-cov table
l <- list()
for  (smp in smps){
  print(smp)
  # tmp <- read.table(paste0("tmp/comapre_Fig5_overlap_meta_",smp,".txt2"),sep = "\t",header = T)
  tmp <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/i-pico/diffgene-CFEV_ExonUTR_from_longestTx/i-pico_comapre_Fig5_overlap_meta_",smp,"_",rbp,".txt4"),sep = "\t",header = T)
  tmp$smp <- smp
  print(dim(tmp))
  l[[smp]] <- tmp
}
res.list.df <- do.call(rbind,l)


# plot in R
#res.list.df$signal <- gsub("G4iM","G4",res.list.df$signal)
table(res.list.df$signal)
#res.list.df$signal <- factor(res.list.df$signal,levels = c("RBPs","EV","G4"))
table(res.list.df$region)
#res.list.df$region <- factor(res.list.df$region,levels = c("Piranha","CLIPper","CLAM","exPeak"))
# table(res.list.df$region)


#(500+500)/20+10=60,再乘以1个样本,再乘以1个区域,为60列
clean_refp <- dplyr::as_tibble(res.list.df)
#clean_refp$smp <- NULL
#clean_refp[1:3,]
clean_refp <- clean_refp[,c(1:3,ncol(clean_refp),4:(ncol(clean_refp)-1))]
#clean_refp <- refp #%>% select(c(-1:-3,-5,-6))
colnames(clean_refp)[c(1,3)] <- c('V4','sample')
#plot(x=1:(ncol(clean_refp)-3), y=as.numeric(clean_refp[1,4:ncol(clean_refp)]))
table(clean_refp$region,clean_refp$sample)

clean_refp <- clean_refp[which( apply(clean_refp[,5:ncol(clean_refp)] ,1, sd)>0),] # filter low sd （noise）, seem not approriate when compare with bg regions
#~ half num of peaks have no ref overlap: sd==0

#1st max-min scale
clean_refp[,5:ncol(clean_refp)] <- t(maxmin.normalize(t(clean_refp[,5:ncol(clean_refp)])))
#clean_refp.tmp <- clean_refp[clean_refp$sample=="G4iM",]
#summary(as.numeric(clean_refp.tmp[3,4:ncol(clean_refp.tmp)]))

# clean_refp[1:6,1:6]
long_refp <- reshape2::melt(clean_refp,id.vars = c('V4','sample','region','smp'),value.name = 'signal') %>%
  select(-variable)
long_refp[1:3,]
#summary(long_refp$signal)

#add x position
long_refp$pos <- rep(c(1:(ncol(clean_refp)-4)),each = nrow(clean_refp)) # ,times = nrow(cmb)
# dim(long_refp)


table(long_refp$sample)
# run loess model 
# job::job({
l <- list()
for(sample in as.character(unique(long_refp$sample))){
  #sample <- 45
  print(sample)
  for (region in as.character(unique(long_refp$region))){
    #region <- "shuffle"
    print(region)
    for (smp in as.character(unique(long_refp$smp))){
      # smp <- "L39-5-BC3"
      print(smp)
      # print(Sys.time())
      long_refp_tmp <- long_refp[long_refp$sample==sample & long_refp$region==region & long_refp$smp==smp,] # not merge smps
      #plot(long_refp_tmp$pos, long_refp_tmp$signal)
      # long_refp_tmp <- get.loess(long_refp_tmp) #[,c("pos","signal")] # get.lowess, get.loess
      long_refp_tmp <- as.data.frame(get.loess.smooth(long_refp_tmp)) #[,c("pos","signal")] # get.lowess, get.loess
      # long_refp_tmp$signal[long_refp_tmp$signal < 0] <- 0
      # long_refp_tmp$signal[long_refp_tmp$signal > 1] <- 1
      long_refp_tmp$smp <- smp
      long_refp_tmp$sample <- sample
      long_refp_tmp$region <- region
      colnames(long_refp_tmp)[colnames(long_refp_tmp)=="signal"] <- "mean_signal"
      l[[paste0(sample,region,smp)]] <- long_refp_tmp
      #table(long_refp_tmp$pos==long_refp_tmp2$pos)
      # write.table(long_refp_tmp,paste0("tmp/comapre_Fig5_overlap_meta_final_scaler_",sample,"_",region,"_",smp,".txt3"),row.names = F,col.names = T,sep = "\t",quote = F)
      write.table(long_refp_tmp,paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/i-pico/diffgene-CFEV_ExonUTR_from_longestTx/i-pico_comapre_Fig5_overlap_meta_final_scaler_",sample,"_",region,"_",smp,"_",rbp,".txt4"),row.names = F,col.names = T,sep = "\t",quote = F)
    }
  }
}
# })


## read in lowess
l <- list()
for(sample in as.character(unique(long_refp$sample))){
  #sample <- "RBPs"
  print(sample)
  for (region in as.character(unique(long_refp$region))){
    #region <- "exPeak"
    print(region)
    
    for (smp in as.character(unique(long_refp$smp))){
      # smp <- "SAMN03863400"
      print(smp)
      # tmp <- read.table( paste0("tmp/comapre_Fig5_overlap_meta_final_scaler_",sample,"_",region,"_",smp,".txt3"), sep = "\t",header = T)
      tmp <- read.table( paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/i-pico/diffgene-CFEV_ExonUTR_from_longestTx/i-pico_comapre_Fig5_overlap_meta_final_scaler_",sample,"_",region,"_",smp,"_",rbp,".txt4"), sep = "\t",header = T)
      # tmp$smp <- smp
      l[[paste0(sample,region,smp)]] <- tmp
    }
  }
}
filnal_scaler <- do.call(rbind,l)
table(filnal_scaler$sample)
table(filnal_scaler$region)
table(filnal_scaler$smp)
summary(filnal_scaler$mean_signal)


up_bins <- as.integer(up/bin) # ceiling
down_bins <- as.integer(down/bin) # ceiling
target_bins <- round((up_bins + down_bins)*(ratio/(1-ratio)),digits = 0) #as.integer((up_bins + down_bins)*(ratio/(1-ratio)))
table(filnal_scaler$sample)
table(filnal_scaler$region)
#filnal_scaler$region <- factor(filnal_scaler$region,levels = c("Piranha","CLIPper","CLAM","exPeak"))

ggplot(filnal_scaler,aes(x = pos,y = mean_signal)) +
  # add 0.95 interval
  # geom_ribbon(aes(ymin = lower,
  #                 ymax = upper, fill=region), # sample
  #             alpha = 0.5) +
  geom_line(aes(color = region), size = 2) + # 
  theme_classic(base_size = 16) +
  scale_color_manual(values = c("steelblue","grey50") ) +  # [4]
  scale_fill_manual(values = c(alpha("steelblue",0.5),"grey90") ) +  # [4]
  # scale_color_d3(name = 'Data type') +
  # scale_fill_d3(name = '') +
  geom_vline(xintercept = c(up_bins+1,up_bins+target_bins),color="grey50",linetype="dashed")+
  # geom_vline(xintercept = c((1+up_bins+target_bins+down_bins)*0.5),color="grey50",linetype="dashed")+ # need +1 for 1-based coord
  # x label
  scale_x_continuous(breaks = c(1,up_bins+1,up_bins+target_bins,up_bins+target_bins+down_bins), # need +1 for 1-based coord
                     labels = c(paste0("-",up,"nt"),'Start','End',paste0("+",down,"nt")) ) +
  scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1.0), limits = c(0,1), labels = c("0","0.25","0.50","0.75","1.00")) +
  xlab('') + ylab('Scaled depth') +
  theme(aspect.ratio = 0.6,
        strip.background = element_rect(color = NA,fill = 'white'),
        strip.text = element_text(size=30),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=30),
        axis.text.x = element_blank(), #  element_text(size=20,angle = 0,hjust = 0.5)
        axis.text.y = element_text(size=30),
        legend.position =  "right",#c(0.9,0.8),
        legend.text = element_text(size= 24),
        legend.title= element_text(size= 30) ) + 
  facet_grid(#.~sample, 
    cols = vars(sample), 
    scales = "free_y")

#2nd max-min scale (only display mean line, no CI/ribbon)
filnal_scaler <- filnal_scaler %>% 
  dplyr::group_by(sample,region) %>% # add sample ?
  dplyr::mutate(mean_signal_norm= (mean_signal - min(mean_signal))/(max(mean_signal)-min(mean_signal))^as.logical(sd(mean_signal)) )

ggplot(filnal_scaler,aes(x = pos,y = mean_signal_norm)) +
  # add 0.95 interval
  # geom_ribbon(aes(ymin = lower,
  #                 ymax = upper, fill=region), # sample
  #             alpha = 0.5) +
  geom_line(aes(color = region), size = 2) + # 
  theme_classic(base_size = 16) +
  scale_color_manual(values = c("steelblue","grey50") ) +  # [4]
  scale_fill_manual(values = c(alpha("steelblue",0.5),"grey90") ) +  # [4]
  # scale_color_d3(name = 'Data type') +
  # scale_fill_d3(name = '') +
  geom_vline(xintercept = c(up_bins+1,up_bins+target_bins),color="grey50",linetype="dashed")+
  # geom_vline(xintercept = c((1+up_bins+target_bins+down_bins)*0.5),color="grey50",linetype="dashed")+ # need +1 for 1-based coord
  # x label
  scale_x_continuous(breaks = c(1,up_bins+1,up_bins+target_bins,up_bins+target_bins+down_bins), # need +1 for 1-based coord
                     labels = c(paste0("-",up,"nt"),'Start','End',paste0("+",down,"nt")) ) +
  scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1.0), limits = c(0,1), labels = c("0","0.25","0.50","0.75","1.00")) +
  xlab('') + ylab('Scaled depth') +
  theme(aspect.ratio = 0.6,
        strip.background = element_rect(color = NA,fill = 'white'),
        strip.text = element_text(size=30),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=30),
        axis.text.x = element_blank(), #  element_text(size=20,angle = 0,hjust = 0.5)
        axis.text.y = element_text(size=30),
        legend.position =  "right",#c(0.9,0.8),
        legend.text = element_text(size= 24),
        legend.title= element_text(size= 30) ) + 
  facet_grid(#.~sample, 
    cols = vars(sample), 
    scales = "free_y")

#2nd loess (if no CI/ribbon implemented)
#add CI/ribbon
filnal_scaler <- filnal_scaler %>%
  tidyr::drop_na() %>% # remove na
  dplyr::group_by(pos,sample,region) %>% # sample,
  # mean and 95% interval confidence
  dplyr::summarise(mean_signal_norm2 = mean(mean_signal_norm, na.rm=T), #,trim = 0.05, na.rm=T
                   sd = sd(mean_signal_norm),
                   upper = CI(mean_signal_norm,ci = 0.95)[1],
                   lower = CI(mean_signal_norm,ci = 0.95)[3],
                   max = max(mean_signal_norm),
                   min = min(mean_signal_norm)
  )

table(filnal_scaler$region)
filnal_scaler2 <- filnal_scaler[(filnal_scaler$region %in% "RBP") & (filnal_scaler2$sample==85),] # RBP, shuffle
#filnal_scaler2$sample <- factor(filnal_scaler2$sample, levels = c("45","65","85"))
ggplot(filnal_scaler2,aes(x = pos,y = mean_signal_norm2)) +
  # add 0.95 interval
  geom_ribbon(aes(ymin = lower,
                  ymax = upper, fill=region), # sample
              alpha = 0.5) +
  geom_line(aes(color = region), size = 2) + # 
  theme_classic(base_size = 16) +
  scale_color_manual(values = c("steelblue","grey50") ) +  # [4]
  scale_fill_manual(values = c(alpha("steelblue",0.5),"grey90") ) +  # [4]
  # scale_color_d3(name = 'Data type') +
  # scale_fill_d3(name = '') +
  geom_vline(xintercept = c(up_bins+1,up_bins+target_bins),color="grey50",linetype="dashed")+
  # geom_vline(xintercept = c((1+up_bins+target_bins+down_bins)*0.5),color="grey50",linetype="dashed")+ # need +1 for 1-based coord
  # x label
  scale_x_continuous(breaks = c(1,up_bins+1,up_bins+target_bins,up_bins+target_bins+down_bins), # need +1 for 1-based coord
                     labels = c(paste0("-",up,"nt"),'Start','End',paste0("+",down,"nt")) ) +
  scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1.0), limits = c(0,1.0), labels = c("0","0.25","0.50","0.75","1.00")) +
  xlab('') + ylab('Scaled depth') +
  theme(aspect.ratio = 0.6,
        strip.background = element_rect(color = NA,fill = 'white'),
        strip.text = element_text(size=30),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=30),
        axis.text.x = element_blank(), #  element_text(size=20,angle = 0,hjust = 0.5)
        axis.text.y = element_text(size=30),
        legend.position =  "right",#c(0.9,0.8),
        legend.text = element_text(size= 24),
        legend.title= element_text(size= 30) ) + 
  facet_grid(#.~sample, 
    cols = vars(sample), 
    scales = "free_y")
# ggsave("cov_meta1_v2.pdf",width = 20,height = 8)



# heatmap (need further optimization?)
long_refp2 <- long_refp # [!(long_refp$region=="exPeak"),]
clr <- c("#1B61A5","#FC6910","#269321")

#plot(x=1:3,y=1:3,col=clr)
methods <- as.character(unique(long_refp2$region)) # 
#methods <- "exPeak"
for (i in 1:length(methods)){
  #i <- 1
  method <- methods[i]
  print(method)
  long_refp2_tmp <- long_refp2[long_refp2$region==method,]
  
  plot.list <- list()
  for (sample in c("85")){ # as.character(unique(long_refp2_tmp$sample))
    # sample <- "RBPs"
    print(sample)
    long_refp2_tmp2 <- long_refp2_tmp[long_refp2_tmp$sample==sample,]
    long_refp_sort <- long_refp2_tmp2 %>%
      # dplyr::filter(pos %in% 25:35) %>% # only rank/sort by central regions
      dplyr::group_by(V4) %>%  # sample,
      #mean of each row/region
      dplyr::summarise(mean_signal_region = mean(signal,na.rm=T)
      )
    long_refp_sort <- long_refp_sort[order(long_refp_sort$mean_signal_region,decreasing = F),]
    long_refp2_tmp3 <- long_refp2_tmp2
    long_refp2_tmp3$V4 <- factor(long_refp2_tmp3$V4,levels = unique(long_refp_sort$V4))
    table(long_refp2_tmp3$sample)
    #table(is.na(long_refp2_tmp3$signal))
    table(long_refp2_tmp3$pos)
    tmp.p <- ggplot(long_refp2_tmp3, aes(x = pos,y = V4)) +
      geom_tile(aes(fill = signal)) + # 
      theme_minimal() +
      # coord_cartesian(expand = 0) +
      scale_x_continuous(breaks = c(1,up_bins+0.5,up_bins+target_bins+0.5,up_bins+target_bins+down_bins), # need +1 for 1-based coord
                         limits = c(1,up_bins+target_bins+down_bins),
                         labels = c(paste0("-",up,"nt"),'Start','End',paste0("+",down,"nt"))) +
      scale_fill_gradient(low = 'white',high = clr[i]) + # , na.value = "white")
      ylab('') + xlab('') +
      theme(
        aspect.ratio = 0.4,
        plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
        # plot.background = element_rect(fill = "white"),
        # panel.background = element_rect(fill = "white", colour = "white"),
        # strip.background = element_rect(color = NA,fill = 'white'),
        strip.text = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20,angle = 0,hjust = 0.5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position =  "none",#c(0.9,0.8),
        legend.text = element_text(size= 16),
        legend.title= element_text(size= 16)) #+
    # theme(aspect.ratio = 1.5,
    # plot.margin = unit(c(5.5, 5.5, 0, 5.5), "pt")
    #       panel.grid = element_blank(),
    #       axis.text.x = element_text(size=20,angle = 0,hjust = 0.5),
    #       # axis.text.y = element_text(size=20),
    #       legend.position = "right",
    #       legend.text = element_text(size= 16),
    #       legend.title= element_text(size= 16))+
    # facet_grid( . ~ sample ,  scales = "free_y" ) # region, space = "free_y",
    plot.list[[paste0(method,"_",sample)]] <- tmp.p
  }
  # p <- cowplot::plot_grid(plotlist = plot.list, rel_widths = c(1,1,1), nrow = 1, align = "hv", axis = "b")
  p <- ggarrange(plotlist = plot.list, rel_widths = c(1,1,1), nrow = 1,  common.legend = T,align = "hv", axis = "b", legend = "right") #  labels = c("RBPs", "", "EV", "", "G4"),
  ggsave(plot = p, paste0(method,"_cov_heat_v2.pdf"), width=15, height=3) # "_",sample
}





#



# compare Fig4 RSS meta (icSHaPE, why unexpected) ------------
## read ref
ref <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)
table(ref$transcript_type)
ref <- ref[ref$transcript_type=="mRNA",]
ref$transcript_id2 <- unlist(sapply(strsplit(ref$transcript_id,"_____"),"[",1))
ref[10000:10003,]


#todo: expand region of selectedBaseReactivities

## read icShape data
file_path <- "/BioII/lulab_b/baopengfei/projects/RNA-structure/data/cell-2016-icShape/GSE74353_HS_293T_icSHAPE_InVivo_BaseReactivities_filterLen.txt"
icshape <- readLines(file_path)
icshape[1:2]
length(icshape)
icshape.l <- sapply(strsplit(icshape,"\t"),"[")
length(icshape.l[[4]])
icshape.df <- do.call("rbind",icshape.l)
icshape.df <- as.data.frame(icshape.df)
# dim(icshape.df)
# table(duplicated(icshape.df$V1))
g <- icshape.df$V1
gl <- as.numeric(icshape.df$V2)
icshape.df <- icshape.df[,4:ncol(icshape.df)]
icshape.df[1:4,1:10]
# table(icshape.df=="NULL")
# FALSE     TRUE 
# 39936730 59468710
icshape.df[icshape.df=="NULL" | icshape.df=="NA"] <- NaN
r <- nrow(icshape.df)
icshape.df <- matrix(as.numeric(as.matrix(icshape.df)),nrow=r)
#table(icshape.df>1)
#hist(icshape.df,breaks = 100000,xlim = c(0,2))
icshape.df[icshape.df>1] <- 1
rownames(icshape.df) <- g

#only keep the same mRNA tx version/length
table(g %in% ref$transcript_id2)
# FALSE  TRUE 
# 1359  8801 
table(paste(g,gl) %in% paste(ref$transcript_id2,ref$tx.length))
# FALSE  TRUE 
# 9124  1036
icshape.df <- icshape.df[ paste(g,gl) %in% paste(ref$transcript_id2,ref$tx.length) , ]
g2 <- g[paste(g,gl) %in% paste(ref$transcript_id2,ref$tx.length)]
gl2 <- gl[paste(g,gl) %in% paste(ref$transcript_id2,ref$tx.length)]
dim( icshape.df )

icshape.df <- as.data.frame(icshape.df)
icshape.df$tx <- g2
icshape.df.long <- reshape2::melt(icshape.df, id.var="tx", value.name = "value")
head(icshape.df.long)
#table(is.nan(icshape.df.long$value))
icshape.df.long <- icshape.df.long[!is.nan(icshape.df.long$value),]
icshape.df.long$start <- as.numeric(gsub("V","",icshape.df.long$variable))-1
table(icshape.df.long$start>=0)
icshape.df.long$end <- icshape.df.long$start+1
icshape.df.long$chr <- ref$transcript_id[match(icshape.df.long$tx,ref$transcript_id2)]
icshape.df.long$name <- "X"
icshape.df.long$score <- icshape.df.long$value
icshape.df.long$strand <- "+"
icshape.df.long <- icshape.df.long[,c("chr","start","end","name","score","strand")]
icshape.df.long <- icshape.df.long[!is.na(icshape.df.long$score),]
#table(is.na(icshape.df.long$score))
dim(icshape.df.long)
icshape.df.long <- icshape.df.long[order(icshape.df.long$chr,as.numeric(icshape.df.long$start),decreasing=F),]
head(icshape.df.long,n=10)
data.table::fwrite(icshape.df.long,"/BioII/lulab_b/baopengfei/projects/RNA-structure/data/cell-2016-icShape/GSE74353_HS_293T_icSHAPE_InVivo_BaseReactivities_filterLen_exPeak_mRNA.bed",quote = F,sep = "\t",row.names = F,col.names = F)
#1036*9784

#cat /BioII/lulab_b/baopengfei/projects/RNA-structure/data/cell-2016-icShape/GSE74353_HS_293T_icSHAPE_InVivo_BaseReactivities_filterLen_exPeak_mRNA.bed|bedtools sort > /BioII/lulab_b/baopengfei/projects/RNA-structure/data/cell-2016-icShape/GSE74353_HS_293T_icSHAPE_InVivo_BaseReactivities_filterLen_exPeak_mRNA_sort.bed
#genomeCoverageBed -g genome/hg38/chrom_sizes/transcriptome_genome_sort_uniq_newTxID -i /BioII/lulab_b/baopengfei/projects/RNA-structure/data/cell-2016-icShape/GSE74353_HS_293T_icSHAPE_InVivo_BaseReactivities_filterLen_exPeak_mRNA_sort.bed -bg -split > /BioII/lulab_b/baopengfei/projects/RNA-structure/data/cell-2016-icShape/GSE74353_HS_293T_icSHAPE_InVivo_BaseReactivities_filterLen_exPeak_mRNA_sort.bg
# bed <- rtracklayer::import.bed("/BioII/lulab_b/baopengfei/projects/RNA-structure/data/cell-2016-icShape/GSE74353_HS_293T_icSHAPE_InVivo_BaseReactivities_filterLen_exPeak_mRNA_sort.bed")
# rtracklayer::export.bedGraph(object=bed, con="/BioII/lulab_b/baopengfei/projects/RNA-structure/data/cell-2016-icShape/GSE74353_HS_293T_icSHAPE_InVivo_BaseReactivities_filterLen_exPeak_mRNA_sort.bg")

# icshape.df.long2 <- icshape.df.long[icshape.df.long$chr=="ENST00000053468_____3",]
# tmp <- rle(icshape.df.long2$score)
# length(tmp$lengths)
##seem not decrease much record number

# bedtools unionbedg -i /BioII/lulab_b/baopengfei/projects/RNA-structure/data/cell-2016-icShape/GSE74353_HS_293T_icSHAPE_InVivo_BaseReactivities_filterLen_exPeak_mRNA_sort.bg \
#   test.bg > /BioII/lulab_b/baopengfei/projects/RNA-structure/data/cell-2016-icShape/GSE74353_HS_293T_icSHAPE_InVivo_BaseReactivities_filterLen_exPeak_mRNA_sort.bg.unionbedg 
# 






setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC")
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggsci)
library(Rmisc)
library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("which", "base")

#tx
up=50
down=50
bin=5 # 10 for test, 5
ratio=0.3 # 20/(50*2+20)

dst <- "GSE71008" # GSE71008,GSE110381,GSE123972
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
smps <- read.table(paste0(pre,"/exSeek-dev/data/",dst,"/sample_ids_NC_test15.txt"))$V1


## pre-filter boundary peak (sample-wise)
#smps <- smps[sample(1:length(smps),10,replace = F)]
for(i in 1:length(smps)){
  #smp <- "SAMN03863458"
  smp <- smps[i]
  print(smp)
  for(m in c("piranha_by_sample/b5_p01","clipper_by_sample/b5_p05","clam_by_sample/b5_p005","expeakCNN_by_sample/b5_d50_p1")){
    print(m)
    df <- read.table(paste0("output/",dst,"/call_peak_all/",m,"/intersect/",smp,".bed6"))
    
    df$txLen <- ref$tx.length[match(df$V1,ref$transcript_id)]
    #peakLen <- 20
    flankLen <- 50
    df <- df[(df$V2-flankLen)>=0 & (df$V3+flankLen)<=df$txLen,] 
    write.table(df,paste0("output/",dst,"/call_peak_all/",m,"/intersect/",smp,"_rmBoundary.bed"),quote = F,sep = "\t",row.names = F,col.names = F) # _noBoundary
  }
}


# cal meta-cov table
#use _noBoundary.bed6 seem not good
#smps <- smps[6:15]
for(i in 1:length(smps)){
  smp <- smps[i]
  #smp <- "SAMN03863400"
  print(smp)
  cmb <- expand.grid(signal=c("/BioII/lulab_b/baopengfei/projects/RNA-structure/data/cell-2016-icShape/GSE74353_HS_293T_icSHAPE_InVivo_BaseReactivities_filterLen_exPeak_mRNA.bed"),
                     region=c(#paste0("output/",dst,"/call_peak_all/piranha_by_sample/b5_p01/intersect/",smp,"_rmBoundary.bed"), # _noBoundary
                              paste0("output/",dst,"/call_peak_all/clipper_by_sample/b5_p05/intersect/",smp,"_rmBoundary.bed"), # _noBoundary
                              #paste0("output/",dst,"/call_peak_all/clam_by_sample/b5_p005/intersect/",smp,"_rmBoundary.bed"), # _noBoundary
                              paste0("output/",dst,"/call_peak_all/expeakCNN_by_sample/b5_d50_p1/intersect/",smp,"_rmBoundary.bed")) # _noBoundary
  )
  
  
  cmb$signal.label <- "icshape"
  # cmb$signal.label[grepl("all_tx",cmb$signal)] <- c("RBPs") #,"small","CLIP","ex RIP"
  # cmb$signal.label[grepl("g4grinder",cmb$signal)] <- c("G4") #,"small","CLIP","ex RIP"
  # cmb$signal.label[grepl("EVenrich",cmb$signal)] <- c("EV") #,"small","CLIP","ex RIP"
  # table(cmb$signal.label)
  cmb$region.label <- "Piranha"
  cmb$region.label[grepl("clam",tolower(cmb$region))] <- "CLAM"
  cmb$region.label[grepl("clipper",tolower(cmb$region))] <- "CLIPper"
  cmb$region.label[grepl("expeak",tolower(cmb$region))] <- "exPeak"
  
  # need prefilter boundary peaks (EnrichedHeatmap::normalizeToMatrix seem not take chr.size info when extend region)
  
  res.list <- lapply(1:nrow(cmb), function(j) {get.mat(signal=as.character(cmb$signal[j]), region=as.character(cmb$region[j]),
                                                       signal.label=as.character(cmb$signal.label[j]),region.label=as.character(cmb$region.label[j]),
                                                       up = up, down = down, bin = bin, ratio = ratio, 
                                                       mean_mode = "coverage" )} ) # mean_mode doesnot change much
  res.list.df <- do.call(rbind,res.list)
  write.table(res.list.df,paste0("tmp/comapre_Fig5_overlap_meta_icshape_",smp,".txt3"),row.names = T,col.names = T,sep = "\t",quote = F)
  
  # res.list.df[1:3,]
}



## read meta-cov table
l <- list()
# smps <- Sys.glob("tmp/comapre_Fig5_overlap_meta_SAMN*.txt")
# smps <- basename(smps)
# smps <- gsub("comapre_Fig5_overlap_meta","",smps)
# smps <- gsub(".txt","",smps,fixed = T)
for  (smp in smps){
  print(smp)
  # tmp <- read.table(paste0("tmp/comapre_Fig5_overlap_meta_",smp,".txt2"),sep = "\t",header = T)
  tmp <- read.table(paste0("tmp/comapre_Fig5_overlap_meta_icshape_",smp,".txt3"),sep = "\t",header = T)
  tmp$smp <- smp
  print(dim(tmp))
  l[[smp]] <- tmp
}
res.list.df <- do.call(rbind,l)
# table(res.list.df$signal)
# for(i in 1:length(l)){
#   print(i);
#   print(dim(l[[i]]))
# }


# plot in R
res.list.df$region <- factor(res.list.df$region,levels = c("CLIPper","exPeak")) # c("Piranha","CLIPper","CLAM","exPeak")
# table(res.list.df$region)


#(500+500)/20+10=60,再乘以1个样本,再乘以1个区域,为60列
clean_refp <- dplyr::as_tibble(res.list.df)
#clean_refp$smp <- NULL
clean_refp <- clean_refp[,c(1:3,ncol(clean_refp),4:(ncol(clean_refp)-1))]
#clean_refp <- refp #%>% select(c(-1:-3,-5,-6))
colnames(clean_refp)[c(1,3)] <- c('V4','sample')
#plot(x=1:(ncol(clean_refp)-3), y=as.numeric(clean_refp[1,4:ncol(clean_refp)]))
table(clean_refp$region,clean_refp$sample)

clean_refp <- clean_refp[which( apply(clean_refp[,5:ncol(clean_refp)] ,1, sd)>0),] # filter low sd
#~ half num of peaks have no ref overlap: sd==0

#1st max-min scale
clean_refp[,5:ncol(clean_refp)] <- t(maxmin.normalize(t(clean_refp[,5:ncol(clean_refp)])))
#clean_refp.tmp <- clean_refp[clean_refp$sample=="G4iM",]
#summary(as.numeric(clean_refp.tmp[3,4:ncol(clean_refp.tmp)]))

# clean_refp[1:6,1:6]
long_refp <- reshape2::melt(clean_refp,id.vars = c('V4','sample','region','smp'),value.name = 'signal') %>%
  select(-variable)
long_refp[1:3,]
#summary(long_refp$signal)

#add x position
long_refp$pos <- rep(c(1:(ncol(clean_refp)-4)),each = nrow(clean_refp)) # ,times = nrow(cmb)
# dim(long_refp)


# run loess model 
# job::job({
l <- list()
for(sample in as.character(unique(long_refp$sample))){
  #sample <- "RBPs"
  print(sample)
  for (region in as.character(unique(long_refp$region))){
    #region <- "exPeak"
    print(region)
    for (smp in as.character(unique(long_refp$smp))){
      # smp <- "SAMN03863400"
      print(smp)
      # print(Sys.time())
      long_refp_tmp <- long_refp[long_refp$sample==sample & long_refp$region==region & long_refp$smp==smp,] # not merge smps
      #plot(long_refp_tmp$pos, long_refp_tmp$signal)
      # long_refp_tmp <- get.loess(long_refp_tmp) #[,c("pos","signal")] # get.lowess, get.loess
      long_refp_tmp <- as.data.frame(get.loess.smooth(long_refp_tmp)) #[,c("pos","signal")] # get.lowess, get.loess
      # long_refp_tmp$signal[long_refp_tmp$signal < 0] <- 0
      # long_refp_tmp$signal[long_refp_tmp$signal > 1] <- 1
      long_refp_tmp$smp <- smp
      long_refp_tmp$sample <- sample
      long_refp_tmp$region <- region
      colnames(long_refp_tmp)[colnames(long_refp_tmp)=="signal"] <- "mean_signal"
      l[[paste0(sample,region,smp)]] <- long_refp_tmp
      #table(long_refp_tmp$pos==long_refp_tmp2$pos)
      # write.table(long_refp_tmp,paste0("tmp/comapre_Fig5_overlap_meta_final_scaler_",sample,"_",region,".txt2"),row.names = F,col.names = T,sep = "\t",quote = F)
      write.table(long_refp_tmp,paste0("tmp/comapre_Fig5_overlap_meta_icshape_final_scaler_",sample,"_",region,"_",smp,".txt3"),row.names = F,col.names = T,sep = "\t",quote = F)
    }
  }
}
# })


## read in lowess
l <- list()
for(sample in as.character(unique(long_refp$sample))){
  #sample <- "RBPs"
  print(sample)
  for (region in as.character(unique(long_refp$region))){
    #region <- "exPeak"
    print(region)
    
    for (smp in as.character(unique(long_refp$smp))){
      # smp <- "SAMN03863400"
      print(smp)
      # tmp <- read.table( paste0("tmp/comapre_Fig5_overlap_meta_final_scaler_",sample,"_",region,".txt2"), sep = "\t",header = T)
      tmp <- read.table( paste0("tmp/comapre_Fig5_overlap_meta_icshape_final_scaler_",sample,"_",region,"_",smp,".txt3"), sep = "\t",header = T)
      # tmp$smp <- smp
      l[[paste0(sample,region,smp)]] <- tmp
    }
  }
}
filnal_scaler <- do.call(rbind,l)
table(filnal_scaler$sample)
table(filnal_scaler$region)
table(filnal_scaler$smp)
summary(filnal_scaler$mean_signal)


up_bins <- as.integer(up/bin) # ceiling
down_bins <- as.integer(down/bin) # ceiling
target_bins <- round((up_bins + down_bins)*(ratio/(1-ratio)),digits = 0) #as.integer((up_bins + down_bins)*(ratio/(1-ratio)))
table(filnal_scaler$sample)
table(filnal_scaler$region)
filnal_scaler$region <- factor(filnal_scaler$region,levels = c("CLIPper","exPeak")) # c("Piranha","CLIPper","CLAM","exPeak")

#2nd max-min scale (only display mean line, no CI/ribbon)
filnal_scaler <- filnal_scaler %>% 
  dplyr::group_by(sample,region) %>% # add sample ?
  dplyr::mutate(mean_signal_norm= (mean_signal - min(mean_signal))/(max(mean_signal)-min(mean_signal))^as.logical(sd(mean_signal)) )

#2nd loess (if no CI/ribbon implemented)
#add CI/ribbon
filnal_scaler <- filnal_scaler %>%
  tidyr::drop_na() %>% # remove na
  dplyr::group_by(pos,sample,region) %>% # sample,
  # mean and 95% interval confidence
  dplyr::summarise(mean_signal_norm2 = mean(mean_signal_norm, na.rm=T), #,trim = 0.05, na.rm=T
                   sd = sd(mean_signal_norm),
                   upper = CI(mean_signal_norm,ci = 0.95)[1],
                   lower = CI(mean_signal_norm,ci = 0.95)[3],
                   max = max(mean_signal_norm),
                   min = min(mean_signal_norm)
  )


filnal_scaler2 <- filnal_scaler[filnal_scaler$region %in% "exPeak",]
# filnal_scaler2$sample <- factor(filnal_scaler2$sample, levels = c("RBPs","EV","G4"))
ggplot(filnal_scaler2,aes(x = pos,y = mean_signal_norm2)) +
  # add 0.95 interval
  geom_ribbon(aes(ymin = lower,
                  ymax = upper,
                  fill = "salmon"), # sample
              alpha = 0.5) +
  geom_line(aes(color = region), size = 2) + # 
  theme_classic(base_size = 16) +
  scale_color_manual(values = pal_d3()(4)[4]) + 
  # scale_color_d3(name = 'Data type') +
  # scale_fill_d3(name = '') +
  geom_vline(xintercept = c(up_bins+1,up_bins+target_bins),color="grey50",linetype="dashed")+
  # geom_vline(xintercept = c((1+up_bins+target_bins+down_bins)*0.5),color="grey50",linetype="dashed")+ # need +1 for 1-based coord
  # x label
  scale_x_continuous(breaks = c(1,up_bins+1,up_bins+target_bins,up_bins+target_bins+down_bins), # need +1 for 1-based coord
                     labels = c(paste0("-",up,"nt"),'Start','End',paste0("+",down,"nt")) ) +
  scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1.0), limits = c(0,1), labels = c("0","0.25","0.50","0.75","1.00")) +
  xlab('') + ylab('Scaled depth') +
  theme(aspect.ratio = 0.6,
        strip.background = element_rect(color = NA,fill = 'white'),
        strip.text = element_text(size=30),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=30),
        axis.text.x = element_blank(), #  element_text(size=20,angle = 0,hjust = 0.5)
        axis.text.y = element_text(size=30),
        legend.position =  "right",#c(0.9,0.8),
        legend.text = element_text(size= 24),
        legend.title= element_text(size= 30) ) + 
  facet_grid(#.~sample, 
    cols = vars(sample), 
    scales = "free_y")
#ggsave("cov_meta1_v2.pdf",width = 20,height = 8)

filnal_scaler3 <- filnal_scaler[!(filnal_scaler$region %in% "exPeak"),]
table(filnal_scaler3$sample,filnal_scaler3$region)
filnal_scaler3$sample <- factor(filnal_scaler3$sample, levels = c("RBPs","EV","G4"))
ggplot(filnal_scaler3,aes(x = pos,y = mean_signal_norm2)) +
  # add 0.95 interval
  geom_ribbon(aes(ymin = lower,
                  ymax = upper,
                  fill = region), # sample
              alpha = 0.5) +
  geom_line(aes(color = region), size = 2) + # 
  theme_classic(base_size = 16) +
  scale_color_d3(name = 'Data type') +
  scale_fill_d3(name = '') +
  geom_vline(xintercept = c(up_bins+1,up_bins+target_bins),color="grey50",linetype="dashed")+
  # geom_vline(xintercept = c((1+up_bins+target_bins+down_bins)*0.5),color="grey50",linetype="dashed")+ # need +1 for 1-based coord
  # x label
  scale_x_continuous(breaks = c(1,up_bins+1,up_bins+target_bins,up_bins+target_bins+down_bins), # need +1 for 1-based coord
                     labels = c(paste0("-",up,"nt"),'Start','End',paste0("+",down,"nt")) ) +
  scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1.0), limits = c(0,1), labels = c("0","0.25","0.50","0.75","1.00")) +
  xlab('') + ylab('Scaled depth') +
  theme(aspect.ratio = 0.6,
        strip.background = element_rect(color = NA,fill = 'white'),
        strip.text = element_text(size=30),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=30),
        axis.text.x = element_blank(), #  element_text(size=20,angle = 0,hjust = 0.5)
        axis.text.y = element_text(size=30),
        legend.position =  "right",#c(0.9,0.8),
        legend.text = element_text(size= 24),
        legend.title= element_text(size= 30)) + 
  facet_grid(#.~sample, 
    cols = vars(sample), 
    scales = "free_y")
#ggsave("cov_meta2_v2.pdf",width = 20,height = 8)



# heatmap (need further optimization?)
# # plot (slowly)
# # long_refp$log10.signal <- log10(long_refp$signal+0.00001)
# long_refp_sort <- long_refp %>%
#   # dplyr::filter(pos %in% 25:35) %>% # only rank/sort by central regions
#   dplyr::group_by(V4,sample) %>%  # sample,
#   #mean of each row/region
#   dplyr::summarise(mean_signal_region = mean(signal,na.rm=T)
#   )
# long_refp_sort <- long_refp_sort[order(long_refp_sort$mean_signal_region,decreasing = F),]
long_refp2 <- long_refp[long_refp$region=="exPeak",]
# long_refp2$V4 <- factor(long_refp2$V4,levels = unique(long_refp_sort$V4))

# table(long_refp2$region)
clr <- c("#C9121E","#FC6910","#269321","#1B61A5")
methods <- as.character(unique(long_refp2$region)) # 
#methods <- "exPeak"
for (i in 1:length(methods)){
  #i <- 1
  method <- methods[i]
  print(method)
  long_refp2_tmp <- long_refp2[long_refp2$region==method,]
  
  plot.list <- list()
  for (sample in c("icshape")){ # as.character(unique(long_refp2_tmp$sample))
    # sample <- "RBPs"
    print(sample)
    long_refp2_tmp2 <- long_refp2_tmp[long_refp2_tmp$sample==sample,]
    long_refp_sort <- long_refp2_tmp2 %>%
      # dplyr::filter(pos %in% 25:35) %>% # only rank/sort by central regions
      dplyr::group_by(V4) %>%  # sample,
      #mean of each row/region
      dplyr::summarise(mean_signal_region = mean(signal,na.rm=T)
      )
    long_refp_sort <- long_refp_sort[order(long_refp_sort$mean_signal_region,decreasing = F),]
    long_refp2_tmp3 <- long_refp2_tmp2
    long_refp2_tmp3$V4 <- factor(long_refp2_tmp3$V4,levels = unique(long_refp_sort$V4))
    table(long_refp2_tmp3$sample)
    #table(is.na(long_refp2_tmp3$signal))
    table(long_refp2_tmp3$pos)
    tmp.p <- ggplot(long_refp2_tmp3, aes(x = pos,y = V4)) +
      geom_tile(aes(fill = signal)) + # 
      theme_minimal() +
      # coord_cartesian(expand = 0) +
      scale_x_continuous(breaks = c(1,up_bins+0.5,up_bins+target_bins+0.5,up_bins+target_bins+down_bins), # need +1 for 1-based coord
                         limits = c(1,up_bins+target_bins+down_bins),
                         labels = c(paste0("-",up,"nt"),'Start','End',paste0("+",down,"nt"))) +
      scale_fill_gradient(low = 'white',high = clr[i]) + # , na.value = "white")
      ylab('') + xlab('') +
      theme(
        aspect.ratio = 0.4,
        plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
        # plot.background = element_rect(fill = "white"),
        # panel.background = element_rect(fill = "white", colour = "white"),
        # strip.background = element_rect(color = NA,fill = 'white'),
        strip.text = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20,angle = 0,hjust = 0.5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position =  "none",#c(0.9,0.8),
        legend.text = element_text(size= 16),
        legend.title= element_text(size= 16)) #+
    # theme(aspect.ratio = 1.5,
    # plot.margin = unit(c(5.5, 5.5, 0, 5.5), "pt")
    #       panel.grid = element_blank(),
    #       axis.text.x = element_text(size=20,angle = 0,hjust = 0.5),
    #       # axis.text.y = element_text(size=20),
    #       legend.position = "right",
    #       legend.text = element_text(size= 16),
    #       legend.title= element_text(size= 16))+
    # facet_grid( . ~ sample ,  scales = "free_y" ) # region, space = "free_y",
    plot.list[[paste0(method,"_",sample)]] <- tmp.p
  }
  # p <- cowplot::plot_grid(plotlist = plot.list, rel_widths = c(1,1,1), nrow = 1, align = "hv", axis = "b")
  p <- ggpubr::ggarrange(plotlist = plot.list, rel_widths = c(1,1,1), nrow = 1,  common.legend = T,align = "hv", axis = "b", legend = "right") #  labels = c("RBPs", "", "EV", "", "G4"),
  ggsave(plot = p, paste0(method,"_icshape_heat_v2.pdf"), width=15, height=3) # "_",sample
}
# ggsave(paste0("cov_heat.pdf"), width=15, height=3) # "_",sample

long_refp2 <- long_refp[!(long_refp$region=="exPeak"),]
# long_refp2$V4 <- factor(long_refp2$V4,levels = unique(long_refp_sort$V4))

# table(long_refp2$region)
#clr <- c("#FC6910","#269321","#1B61A5") # 
clr <- c("#1B61A5","#FC6910","#269321")

#plot(x=1:3,y=1:3,col=clr)
methods <- as.character(unique(long_refp2$region)) # 
#methods <- "exPeak"
for (i in 1:length(methods)){
  #i <- 1
  method <- methods[i]
  print(method)
  long_refp2_tmp <- long_refp2[long_refp2$region==method,]
  
  plot.list <- list()
  for (sample in c("icshape")){ # as.character(unique(long_refp2_tmp$sample))
    # sample <- "RBPs"
    print(sample)
    long_refp2_tmp2 <- long_refp2_tmp[long_refp2_tmp$sample==sample,]
    long_refp_sort <- long_refp2_tmp2 %>%
      # dplyr::filter(pos %in% 25:35) %>% # only rank/sort by central regions
      dplyr::group_by(V4) %>%  # sample,
      #mean of each row/region
      dplyr::summarise(mean_signal_region = mean(signal,na.rm=T)
      )
    long_refp_sort <- long_refp_sort[order(long_refp_sort$mean_signal_region,decreasing = F),]
    long_refp2_tmp3 <- long_refp2_tmp2
    long_refp2_tmp3$V4 <- factor(long_refp2_tmp3$V4,levels = unique(long_refp_sort$V4))
    table(long_refp2_tmp3$sample)
    #table(is.na(long_refp2_tmp3$signal))
    table(long_refp2_tmp3$pos)
    tmp.p <- ggplot(long_refp2_tmp3, aes(x = pos,y = V4)) +
      geom_tile(aes(fill = signal)) + # 
      theme_minimal() +
      # coord_cartesian(expand = 0) +
      scale_x_continuous(breaks = c(1,up_bins+0.5,up_bins+target_bins+0.5,up_bins+target_bins+down_bins), # need +1 for 1-based coord
                         limits = c(1,up_bins+target_bins+down_bins),
                         labels = c(paste0("-",up,"nt"),'Start','End',paste0("+",down,"nt"))) +
      scale_fill_gradient(low = 'white',high = clr[i]) + # , na.value = "white")
      ylab('') + xlab('') +
      theme(
        aspect.ratio = 0.4,
        plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
        # plot.background = element_rect(fill = "white"),
        # panel.background = element_rect(fill = "white", colour = "white"),
        # strip.background = element_rect(color = NA,fill = 'white'),
        strip.text = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20,angle = 0,hjust = 0.5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position =  "none",#c(0.9,0.8),
        legend.text = element_text(size= 16),
        legend.title= element_text(size= 16)) #+
    # theme(aspect.ratio = 1.5,
    # plot.margin = unit(c(5.5, 5.5, 0, 5.5), "pt")
    #       panel.grid = element_blank(),
    #       axis.text.x = element_text(size=20,angle = 0,hjust = 0.5),
    #       # axis.text.y = element_text(size=20),
    #       legend.position = "right",
    #       legend.text = element_text(size= 16),
    #       legend.title= element_text(size= 16))+
    # facet_grid( . ~ sample ,  scales = "free_y" ) # region, space = "free_y",
    plot.list[[paste0(method,"_",sample)]] <- tmp.p
  }
  # p <- cowplot::plot_grid(plotlist = plot.list, rel_widths = c(1,1,1), nrow = 1, align = "hv", axis = "b")
  p <- ggarrange(plotlist = plot.list, rel_widths = c(1,1,1), nrow = 1,  common.legend = T,align = "hv", axis = "b", legend = "right") #  labels = c("RBPs", "", "EV", "", "G4"),
  ggsave(plot = p, paste0(method,"_icshape_heat_v2.pdf"), width=15, height=3) # "_",sample
}

#


# comapre Fig1 length, cpm density plot (ggplot, new, 221212) ------------------------------------------------
setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/")
options(stringsAsFactors = F)
options(bedtools.path = "/BioII/lulab_b/baopengfei/biosoft") # /BioII/lulab_b/baopengfei/anaconda3/envs/exvariance/bin
suppressPackageStartupMessages(library(bedtoolsr))
library(tidyverse)
library(dplyr)

plasma.dst <- "GSE71008_NCpool" # "GSE110381_NCpool" 
region <- c(paste0("output/",plasma.dst,"/call_peak_all/clipper_by_sample/b5_p05/","NCpool",".bed"),  # GSE71008
            # paste0("output/","TCGA_small_NCpool","/call_peak_all/clipper_by_sample/b5_p05/","NCpool",".bed"),
            paste0("output/","GSE148861_GSE148862_NCpool","/call_peak_all/clipper_by_sample/b5_p05/","NCpool",".bed"),
            paste0("output/","GSE50676_NCpool","/call_peak_all/clipper_by_sample/b5_p05/","NCpool",".bed")
            # paste0("output/","AGO2_IP_NCpool","/call_peak_dedup/domains_clipper_by_sample/b5_p05/","NCpool",".bed")
)

# pre <- "GSE71008/call_peak_all" # "AGO2_IP/call_peak_dedup_bk", "GSE71008/call_peak_all", "GSE133684/call_peak_dedupByPos"  
# in1 <- Sys.glob(paste0("output/",pre,"/domains/b5_p*_8DNA.bed"))

bed.list <- list()
sum.len <- function(x){
  # x <- "output/AGO2_IP/call_domain/domains/5/50_gn_intersectAGO2.bed"
  # x <- "output/GSE71008_NCpool/call_peak_all/domains_clipper_by_sample/b5_p05/NCpool.bed"
  print(x)
  
  # #x <- "output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax_by_sample/interectAGO2/FTA-12_1.bed"
  rbp <- data.table::fread(paste0(x),data.table = F, header = F,sep = "\t",stringsAsFactors = F)
  
  bed <- as_tibble(rbp)
  bed <- bed[,1:6]
  colnames(bed) <- c("chr","start","end","peak","score","strand")
  summary(bed$end-bed$start)
  # bed <- bed %>% 
  #   filter(end-start<=200 , end-start>=10)
  bed$sample <- basename(x)
  bed$path <- x
  bed.list[[x]] <- bed
  # return(bed)
}
# res2.df2 <- do.call(rbind,mclapply(c(in1,in2,in3,  in5), sum.len, mc.cores = 5) )  # 11RNA
res2.df2 <- do.call(rbind,mclapply(region, sum.len, mc.cores = 2) )  # 8DNA
res2.df2[1:3,]

res2.df2$dataset <- "Null"
res2.df2$dataset[grepl(plasma.dst,res2.df2$path)] <- "extracellular small"
res2.df2$dataset[grepl("TCGA_small|GSE148861",res2.df2$path,perl = T)] <- "cellular small"
# res2.df2$dataset[grepl("AGO2_IP",res2.df2$path)] <- "ex RIP"
res2.df2$dataset[grepl("GSE50676",res2.df2$path)] <- "cellular CLIP"
table(res2.df2$dataset)
res2.df2$dataset <- factor(res2.df2$dataset,levels = c("cellular CLIP","cellular small","extracellular small"))

res2.df2$width <- res2.df2$end-res2.df2$start
# hist(res2.df2$width)
p2 <- ggplot(res2.df2, aes(x=width, color=dataset, fill=dataset)) +
  geom_density(alpha=0.2) +
  # scale_fill_nejm()+
  # scale_color_nejm()+
  scale_fill_manual(name = '', values = c("#FD6905","#843692","#2C68A9")) +
  scale_color_manual(name = 'Data type', values = c("#FD6905","#843692","#2C68A9")) +
  # scale_x_discrete(label=c("Domain Ratio","MIR Ratio"))+
  labs(title="",x="", y = "Peak length density")+
  # facet_grid(method~.)+
  # ylim(c(0,50))+
  # xlim(c(0,1))+
  theme_classic(base_size=12) + 
  theme(#axis.ticks.x=element_blank(),  
    aspect.ratio = 0.5,
    #strip.text.y = element_blank(),
    #strip.text.x = element_text(face="bold",family="arial",size=20),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size = 20), # ,hjust = 1,vjust = 0.5
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size=20),
    strip.text = element_text(size = 20),
    legend.position = "right", #c(0.9,0.8),#,# 
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16))
p2
#ggsave("length_meta.pdf",width = 9,height = 6)


res2.df2[1:3,]
res2.df2$dataset <- factor(res2.df2$dataset, levels = c("extracellular small", "cellular small","cellular CLIP"))
library(ggridges)
ggplot(res2.df2, aes(x=width, y=dataset, fill=dataset))+ #
  geom_density_ridges_gradient(scale = 1.5, rel_min_height = 0.01, gradient_lwd = 3) + 
  # theme_ridges() + 
  # theme(legend.position = "none")
  scale_fill_manual(name = '', values = alpha(colour = c("#2C68A9","#843692","#FD6905"),alpha = 0.8)) +
  xlab('') + ylab('Peak length density') +
  # theme_bw() + 
  xlim(c(10,80))+
  # ylim(c(1,4))+
  theme_classic() +
  theme(aspect.ratio = 0.4,
        strip.background = element_rect(color = NA,fill = 'grey'),
        strip.text = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20,angle = 0,hjust = 0.5),
        axis.text.y = element_blank(), #element_text(size=20),
        legend.position =  "right",#c(0.9,0.8),
        legend.text = element_text(size= 16),
        legend.title= element_text(size= 16)) 
ggsave("cov_meta_ridge.pdf",width = 9,height = 6)



# # plot peak num
# df0 <- as.data.frame(table(res2.df2$dataset))
# ggplot(df0, aes(x=Freq, y=Var1, fill=Var1))+ #
#   geom_bar(stat = "identity") + 
#   scale_fill_manual(name = '', values = alpha(colour = c("#2C68A9","#843692","#FD6905"),alpha = 0.8)) +
#   xlab('') + ylab('Peak number') +
#   theme_classic()
# 
# # plot peak num per tx
# df <- as.data.frame(table(res2.df2$chr,res2.df2$dataset))
# df <- df[df$Freq!=0,]
# df$Var2 <- factor(df$Var2, levels = c("extracellular small", "cellular small","cellular CLIP"))
# #library(ggridges)
# ggplot(df, aes(x=Freq, y=Var2, fill=Var2))+ #
#   geom_boxplot()+
#   # geom_density_ridges_gradient(scale = 1.5, rel_min_height = 0.01, gradient_lwd = 3) + 
#   # theme_ridges() + 
#   # theme(legend.position = "none")
#   scale_fill_manual(name = '', values = alpha(colour = c("#2C68A9","#843692","#FD6905"),alpha = 0.8)) +
#   xlab('') + ylab('Peak number per tx') +
#   # theme_bw() + 
#   xlim(c(0,20))+
#   # ylim(c(1,4))+
#   theme_classic() +
#   theme(aspect.ratio = 0.4,
#         strip.background = element_rect(color = NA,fill = 'grey'),
#         strip.text = element_text(size=20),
#         axis.title.x = element_text(size=20),
#         axis.title.y = element_text(size=20),
#         axis.text.x = element_text(size=20,angle = 0,hjust = 0.5),
#         axis.text.y = element_blank(), #element_text(size=20),
#         legend.position =  "right",#c(0.9,0.8),
#         legend.text = element_text(size= 16),
#         legend.title= element_text(size= 16)) 
# 
# 
# # uniq tx num with peak
# df2 <- as.data.frame(table(df$Var2))
# 
# 
# plot peak full-length ratio
#...


# plot cpm
type <- "CPM" # CPM,RPKM
method <- "clipper_b5_p05" # clipper_b5_p05, expeak_b5_d50_p1
rpkms <- c(paste0("output/",plasma.dst,"/call_peak_all/count_matrix/",method,"_",type,".txt"),  # GSE71008, _CPMtmm.txt
           # paste0("output/","TCGA_small_NCpool","/call_peak_all/count_matrix/",method,"_",type,".txt"),
           paste0("output/","GSE148861_GSE148862_NCpool","/call_peak_all/count_matrix/",method,"_",type,".txt"),
           paste0("output/","GSE50676_NCpool","/call_peak_all/count_matrix/",method,"_",type,".txt")
)
rpkm.list <- list()
sum.rpkm <- function(x){
  # x <- "output/GSE71008_NCpool/call_peak_all/count_matrix/clipper_b5_p05_RPKM.txt"
  # print(x)
  # y <- "output/GSE71008_NCpool/call_peak_all/count_matrix/clipper_b5_p05.txt"
  # #x <- "output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax_by_sample/interectAGO2/FTA-12_1.bed"
  rpkm <- data.table::fread(paste0(x),data.table = F, header = T,check.names = F, sep = "\t",stringsAsFactors = F)
  # y <- data.table::fread(paste0(x,".txt"),data.table = F, header = T,check.names = F, sep = "\t",stringsAsFactors = F)
  # cpm$gene_id <- y$feature
  
  rpkm$path <- x
  rpkm.list[[x]] <- rpkm
  # return(bed)
}
# res2.df2 <- do.call(rbind,mclapply(c(in1,in2,in3,  in5), sum.len, mc.cores = 5) )  # 11RNA
res2.df3 <- do.call(rbind,mclapply(rpkms, sum.rpkm, mc.cores = 1) )  # 8DNA
res2.df3[1:3,]

table(res2.df3$path)
res2.df3$dataset <- "extracellular small"
res2.df3$dataset[grepl("GSE50676",res2.df3$path)] <- "cellular CLIP"
res2.df3$dataset[grepl("TCGA_small|GSE148861",res2.df3$path,perl = T)] <- "cellular small"
res2.df3$dataset <- factor(res2.df3$dataset, levels = c("extracellular small", "cellular small","cellular CLIP"))
table(res2.df3$dataset)
library(ggridges)
#table(is.na(res2.df3$Sample1))
res2.df3$NCpool.log <- log10(res2.df3$NCpool+1)
ggplot(res2.df3, aes(x=NCpool.log, y=dataset, fill=dataset))+ #
  geom_density_ridges_gradient(scale = 1.5,na.rm = T, rel_min_height = 0.01, gradient_lwd = 3) + 
  # theme_ridges() + 
  # theme(legend.position = "none")
  scale_fill_manual(name = '', values = alpha(colour = c("#2C68A9","#843692","#FD6905"),alpha = 0.8)) +
  xlab('') + ylab('Peak CPM density') +
  scale_x_continuous(breaks = c(0,2,4),labels = c("0","100","10000"), limits = c(0,4)) +
  # theme_bw() + 
  # xlim(c(0,4))+
  # ylim(c(1,4))+
  theme_classic() +
  theme(aspect.ratio = 0.5,
        strip.background = element_rect(color = NA,fill = 'grey'),
        strip.text = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20,angle = 0,hjust = 0.5),
        axis.text.y = element_blank(), #element_text(size=20),
        legend.position =  "right",#c(0.9,0.8),
        legend.text = element_text(size= 16),
        legend.title= element_text(size= 16)) 
ggsave("cpm_mat_ridge.pdf",width = 9,height = 6)








# comapre Fig1 full-len ratio density plot and peak RNA type ratio barplot (ggplot, new, 221212) ------------------
options(bedtools.path = "/BioII/lulab_b/baopengfei/biosoft") # /BioII/lulab_b/baopengfei/anaconda3/envs/exvariance/bin
suppressPackageStartupMessages(library(bedtoolsr))
suppressPackageStartupMessages(library(Rsamtools))
library(tidyverse)
library(dplyr)

## read ref
ref <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt",header = T,sep = "\t", stringsAsFactors = F)
ref <- ref[!is.na(ref$tx.length),]

## read peak
#sample <- 'L38-4-BC3'
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
region <- c(paste0(pre,"/output/",plasma.dst,"/call_peak_all/clipper_by_sample/b5_p05/","NCpool",".bed"),  # GSE71008
            # paste0(pre,"/output/","TCGA_small_NCpool","/call_peak_all/clipper_by_sample/b5_p05/","NCpool",".bed"),
            paste0(pre,"/output/","GSE148861_GSE148862_NCpool","/call_peak_all/clipper_by_sample/b5_p05/","NCpool",".bed"),
            paste0(pre,"/output/","GSE50676_NCpool","/call_peak_all/clipper_by_sample/b5_p05/","NCpool",".bed")
            # paste0("output/","AGO2_IP_NCpool","/call_peak_dedup/domains_clipper_by_sample/b5_p05/","NCpool",".bed")
)
samples <- region
readDomain <- function(sample){
  #sample <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008_NCpool/call_peak_all/clipper_by_sample/b5_p05/NCpool.bed"
  print(sample)
  # peak <- rtracklayer::import.bed(as.character(paste0("../WCHSU-FTC/output/WSQ_SMARTer_NEB/call_domain_withRepeats_all/domains_localmax_by_sample_EM/b5_d05_p01/",sample,".bed")))
  peak <- data.table::fread(sample,data.table = F,sep = '\t',header = F,stringsAsFactors = F)
  peak <- as.data.frame(peak)
  peak <- peak[,1:6]
  colnames(peak) <- c("chr","start","end","name","score","strand")
  #peak[1:3,]
  #x <- as.data.frame(cbind(id=as.vector(GenomeInfoDb::seqnames(peak)),value=as.numeric(IRanges::width(peak))))
  peak$sample <- sample
  return(peak)
}

tmp <- lapply(samples,readDomain)
x <- do.call("rbind",tmp)

x <- x[!grepl("_chr",x$chr),]# rm 8DNA (necessary for expeak, not for clipper)
x$tx.length <- ref$tx.length[match(x$chr,ref$transcript_id)]
x$RNA <- ref$transcript_type[match(x$chr,ref$transcript_id)]
x$RNA <- gsub("\\.for|\\.rev|_for|_rev","",x$RNA,perl = T)
rna <- c("rRNA","pri_miRNA","lncRNA", "mRNA","piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA')
dna <- c("intron", "promoter", "enhancer","repeats") # , 
x$RNA <- factor(x$RNA,levels = c(rna)) # dna
#table(x$RNA)

x <- as.data.frame(x)
x$length <- as.numeric(x$end-x$start) #/x$tx.length # /x$length
x$length.ratio <- x$length/x$tx.length

## plot peak length(ratio) ridge density
#x[1:3,]
# df2 <- dplyr::as_tibble(x) %>% 
#   dplyr::select(sample,chr,name,group,lib,group.lib,RNA,length,length.ratio) #%>% 
# # dplyr::group_by(group.lib,RNA) %>% 
# # dplyr::mutate(mean.length.ratio=mean(length.ratio,na.rm=T,trim=0.1) )
# df2[1:3,]
x$dataset <- "extracellular small"
x$dataset[grepl("GSE50676",x$sample)] <- "cellular CLIP"
x$dataset[grepl("TCGA_small|GSE148861",x$sample,perl = T)] <- "cellular small"
x$dataset <- factor(x$dataset, levels = c("extracellular small", "cellular small","cellular CLIP"))
table(x$dataset)

library(ggridges)
ggplot(x, aes(x=length.ratio, y=dataset, fill=dataset))+ #
  geom_density_ridges_gradient(scale = 1.5,na.rm = T, rel_min_height = 0.01, gradient_lwd = 3) + 
  # theme_ridges() + 
  # theme(legend.position = "none")
  scale_fill_manual(name = '', values = alpha(colour = c("#2C68A9","#843692","#FD6905"),alpha = 0.8)) +
  xlab('') + ylab('Peak full-tx ratio') +
  # scale_x_continuous(breaks = c(0,2,4),labels = c("0","100","10000"), limits = c(0,4)) +
  # theme_bw() + 
  xlim(c(0,1))+
  # ylim(c(1,4))+
  # facet_grid(RNA~.) +
  theme_classic() +
  theme(aspect.ratio = 0.5,
        strip.background = element_rect(color = NA,fill = 'grey'),
        strip.text = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20,angle = 0,hjust = 0.5),
        axis.text.y = element_blank(), #element_text(size=20),
        legend.position =  "right",#c(0.9,0.8),
        legend.text = element_text(size= 16),
        legend.title= element_text(size= 16)) 
ggsave("full-len-ratio_ridge.pdf",width = 9,height = 6)


## plot RNA type ratio
x.sum <- as.data.frame(table(x$dataset,x$RNA))
ggplot(x.sum,aes(x=Freq,y=Var1,fill=Var2))+ 
  geom_bar(position = "fill",width = 0.8,linewidth=0.3, stat = "identity") + #  , color="white" 
  #geom_errorbar()+
  xlab("")+
  ylab(paste0(""))+
  geom_vline(xintercept = 0,linetype="dashed",color="grey10")+
  #ggsci::scale_fill_d3()+
  scale_fill_nejm_adaptive(name="Species")+
  # scale_fill_manual(values = c("grey30","steelblue","firebrick","grey90","steelblue2","firebrick2","steelblue1","salmon"))+
  theme_classic() +
  theme(aspect.ratio = 0.5,
        strip.background = element_rect(color = NA,fill = 'grey'),
        strip.text = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20,angle = 0,hjust = 0.5),
        axis.text.y = element_blank(), #element_text(size=20),
        legend.position =  "right",#c(0.9,0.8),
        legend.text = element_text(size= 16),
        legend.title= element_text(size= 16)) 
# theme_void() + 
# theme(plot.title = element_text(size = 24,color="black",hjust = 0.5,face="bold"),
#       axis.title = element_text(size = 24,color ="black"), 
#       axis.text = element_text(size= 20,color = "black"),
#       axis.text.x = element_text(size= 20,color = "black",vjust = 0.5), # angle = 90,,hjust = 1
#       #panel.grid=element_blank(),
#       panel.grid.major.x=element_blank(),
#       panel.grid.minor.x = element_blank(),
#       panel.grid.major.y = element_line(color = "grey50",linetype = "dashed"), #size= 1,
#       #panel.grid.minor.y = element_blank(),
#       panel.border = element_blank(),
#       #axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, color = c(rep("grey50",length(other)),rep("steelblue",length(rna)),rep("firebrick",length(dna)))), #  
#       legend.position = "right",#c(.25,.6),
#       legend.text = element_text(size= 16),
#       legend.title= element_text(size= 16))
ggsave("./3data_peak_RNA_ratio.pdf",width = 9,height = 6)
#RColorBrewer::display.brewer.all()



## tx coverage ratio(tx density)
x <- do.call("rbind",tmp)

x <- x[!grepl("_chr",x$chr),]# rm 8DNA (necessary for expeak, not for clipper)
x$tx.length <- ref$tx.length[match(x$chr,ref$transcript_id)]
x$RNA <- ref$transcript_type[match(x$chr,ref$transcript_id)]
x$RNA <- gsub("\\.for|\\.rev|_for|_rev","",x$RNA,perl = T)
rna <- c("rRNA","pri_miRNA","lncRNA", "mRNA","piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA')
dna <- c("intron", "promoter", "enhancer","repeats") # , 
x$RNA <- factor(x$RNA,levels = c(rna)) # dna
#x[1:3,]
#str(x)
#table(x$RNA)
x <- as_tibble(x) %>% 
  group_by(chr,tx.length,sample) %>% 
  summarise(sumLen=sum(end-start))
summary(x$sumLen)

x <- as.data.frame(x)
#x$length <- as.numeric(x$end-x$start) #/x$tx.length # /x$length
x$span.ratio <- x$sumLen/x$tx.length

## plot peak span(ratio) ridge density
table(x$sample)
x$dataset <- "extracellular small"
x$dataset[grepl("GSE50676",x$sample)] <- "cellular CLIP"
x$dataset[grepl("TCGA_small|GSE148861",x$sample,perl = T)] <- "cellular small"
x$dataset <- factor(x$dataset, levels = c("extracellular small", "cellular small","cellular CLIP"))
table(x$dataset)

library(ggridges)
ggplot(x, aes(x=span.ratio, y=dataset, fill=dataset))+ #span.ratio, sumLen
  geom_density_ridges_gradient(scale = 1.5,na.rm = T, rel_min_height = 0.01, gradient_lwd = 3) + 
  # theme_ridges() + 
  # theme(legend.position = "none")
  scale_fill_manual(name = '', values = alpha(colour = c("#2C68A9","#843692","#FD6905"),alpha = 0.8)) +
  xlab('') + ylab('Peak full-tx span ratio') +
  # scale_x_continuous(breaks = c(0,2,4),labels = c("0","100","10000"), limits = c(0,4)) +
  # theme_bw() + 
  xlim(c(0,1))+
  # xlim(c(0,200))+
  # ylim(c(1,4))+
  # facet_grid(RNA~.) +
  theme_classic() +
  theme(aspect.ratio = 0.5,
        strip.background = element_rect(color = NA,fill = 'grey'),
        strip.text = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20,angle = 0,hjust = 0.5),
        axis.text.y = element_blank(), #element_text(size=20),
        legend.position =  "right",#c(0.9,0.8),
        legend.text = element_text(size= 16),
        legend.title= element_text(size= 16)) 
ggsave("full-len-span_ridge.pdf",width = 9,height = 6)

ggplot(x, aes(x=sumLen, y=dataset, fill=dataset))+ #span.ratio, sumLen
  geom_density_ridges_gradient(scale = 1.5,na.rm = T, rel_min_height = 0.01, gradient_lwd = 3) + 
  # theme_ridges() + 
  # theme(legend.position = "none")
  scale_fill_manual(name = '', values = alpha(colour = c("#2C68A9","#843692","#FD6905"),alpha = 0.8)) +
  xlab('') + ylab('Peak full-tx span ratio') +
  # scale_x_continuous(breaks = c(0,2,4),labels = c("0","100","10000"), limits = c(0,4)) +
  # theme_bw() + 
  # xlim(c(0,1))+
  xlim(c(0,200))+
  # ylim(c(1,4))+
  # facet_grid(RNA~.) +
  theme_classic() +
  theme(aspect.ratio = 0.5,
        strip.background = element_rect(color = NA,fill = 'grey'),
        strip.text = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20,angle = 0,hjust = 0.5),
        axis.text.y = element_blank(), #element_text(size=20),
        legend.position =  "right",#c(0.9,0.8),
        legend.text = element_text(size= 16),
        legend.title= element_text(size= 16)) 
ggsave("full-len-sumlen_ridge.pdf",width = 9,height = 6)




# suppl fig 4: tRNA ridge/complexheatmap: (230103) ------------------------------
setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC")
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggsci)
library(Rmisc)
library(GenomicRanges)
library(IRanges)
# library(data.table)
conflicted::conflict_prefer("select", "dplyr")

## get tRNA gtf
#chrM: 26266,Leu_tRNA
gtf <- rtracklayer::readGFF("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/gtf_by_biotype/tRNA.gtf")
gtf <- as.data.frame(gtf)
gtf <- gtf[,c(-2,-3)]
gtf$start <- gtf$start-1  # convert 1-based to 0-based
gtf$transcript_type <- gsub("_tRNA","",gtf$transcript_type)




## 1.peak cov plot
#tx
up=0
down=0
bin=3
ratio=1

cmb <- data.frame(signal=c(paste0("output/","GSE71008_NCpool","/call_peak_all/expeak/b5_d50_p1.bed")),
                  region=c(paste0("exSeek-dev/genome/hg38/tbed/tRNA.bed")),
                  signal.label=c("ex sm"),
                  region.label=c("tRNA") 
)


#one way using sapply
res.list <- lapply(1:nrow(cmb), function(j) {get.mat(signal=as.character(cmb$signal[j]), region=as.character(cmb$region[j]),
                                                     signal.label=as.character(cmb$signal.label[j]),region.label=as.character(cmb$region.label[j]),
                                                     up = up, down = down, bin = bin, ratio = ratio)} )
res.list.df <- do.call(rbind,res.list)
rownames(res.list.df) <- res.list.df$bed


#plot in R
#res.list.df$signal <- factor(res.list.df$signal,levels = c("cellular CLIP", "cellular small","extracellular small"))
#(500+500)/20+10=60,再乘以1个样本,再乘以1个区域,为60列
clean_refp <- dplyr::as_tibble(res.list.df)
#clean_refp <- refp #%>% select(c(-1:-3,-5,-6))
colnames(clean_refp)[c(1,3)] <- c('V4','sample')

#merge same tRNA
# tx.id <- unlist(sapply(strsplit(clean_refp$V4,"_"),"[",2))
# tx.type <- gtf$transcript_type[match(tx.id,gtf$transcript_id)]
# table(tx.type)
# clean_refp <- clean_refp[tx.type=="Pro",]

#max-min scale
clean_refp[,4:ncol(clean_refp)] <- t(maxmin.normalize(t(clean_refp[,4:ncol(clean_refp)])))
long_refp <- reshape2::melt(clean_refp,id.vars = c('V4','sample','region'),value.name = 'signal') %>%
  select(-variable)
#add sample name
#long_refp$sample <- rep(c("FTA-10","csFTA-10"),
#                        each = nrow(refp)*22) # (100 + 100 + 20)/10
# add x position
long_refp$pos <- rep(c(1:(ncol(clean_refp)-3)),each = nrow(clean_refp)) # ,times = nrow(cmb)
long_refp[1:3,]
long_refp$tx.id <- unlist(sapply(strsplit(long_refp$V4,"_"),"[",2))
long_refp$tx.type <- gtf$transcript_type[match(long_refp$tx.id,gtf$transcript_id)]

#aver each amino acid. type
filnal_scaler <- long_refp %>%
  # remove na
  tidyr::drop_na() %>%
  dplyr::group_by(pos,sample,region,tx.type) %>% # sample,
  # mean and 95% interval confidence
  dplyr::summarise(mean_signal = mean(signal), #,trim = 0.05, na.rm=T
                   sd = sd(signal),
                   upper = CI(signal,ci = 0.95)[1],
                   lower = CI(signal,ci = 0.95)[3]
  )

#2nd max-min scale (only display mean line, no ribbon)
filnal_scaler <- filnal_scaler %>% 
  dplyr::group_by(sample,region,tx.type) %>%
  dplyr::mutate(mean_signal_norm= (mean_signal - min(mean_signal))/(max(mean_signal)-min(mean_signal))^as.logical(sd(mean_signal)) )
table(filnal_scaler$tx.type)
filnal_scaler[1:3,]

#heatmap (need further optimization)
up_bins <- as.integer(up/bin) # ceiling
down_bins <- as.integer(down/bin) # ceiling
target_bins <- 20 #as.integer((up_bins + down_bins)*(ratio/(1-ratio)))

#filnal_scaler_wide[1:3,1:9]
filnal_scaler_wide <- filnal_scaler %>% 
  tidyr::pivot_wider(id_cols = c("sample","region","tx.type"),names_from = "pos",values_from = "mean_signal_norm")

mat <- as.data.frame(filnal_scaler_wide[,4:ncol(filnal_scaler_wide)])
rownames(mat) <- filnal_scaler_wide$tx.type
hc2 <- hclust(dist(mat))
#peak.order <- hc$order
table(rownames(mat)==hc2$labels)
pdf("tRNA_hclust.pdf",width = 36,height = 6)
plot(hc2, hang = -1)
dev.off()
#factor: left->right

#long_refp2 <- long_refp
#table(long_refp2$V4 %in% hc2$labels[hc2$order])
filnal_scaler$tx.type <- factor(filnal_scaler$tx.type,levels = hc2$labels[hc2$order])
#length(hc2$labels)
# clr <- c("#2C68A9","#843692","#FD6905")
# datatypes <- as.character(unique(long_refp2$sample)) # "extracellular small" "cellular small"      "cellular CLIP"  
# for (i in 1:length(datatypes)){
#   datatype <- datatypes[i]
#   print(datatype)
clr <- "firebrick" # purple,salmon,firebrick
tmp.p2 <- ggplot(filnal_scaler, # %>% filter(sample == 'H3K27ac')
                 aes(x = pos,y = tx.type)) +
  geom_tile(aes(fill = mean_signal_norm),color="white") +
  scale_x_continuous(breaks = c(up_bins+1,up_bins+target_bins), # need +1 for 1-based coord
                     labels = c("Start 5'","End 3'") ) +
  scale_fill_gradient(low = 'white',high = clr) + # 'salmon'
  # scale_fill_viridis_c() +
  # scale_fill_continuous("blues") +
  ylab('') + xlab('') +
  theme_void() + 
  theme(aspect.ratio = 1,
        strip.background = element_blank(), #element_rect(color = NA,fill = 'grey'),
        strip.text = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20,angle = 0,hjust = 0.5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position =  "right",#c(0.9,0.8),
        legend.text = element_text(size= 16),
        legend.title= element_text(size= 16)) #+
#factor: down->up
tmp.p2
ggsave(plot = tmp.p2,filename = paste0("tRNA_heat_peak.pdf"),width = 8,height = 8)

#tx.id <- unlist(sapply(strsplit(hc2$labels[hc2$order],"_"),"[",2))
tx.type <- hc2$labels[hc2$order] #gtf$transcript_type[match(tx.id,gtf$transcript_id)]
lab.df <- data.frame("tx.type"=tx.type,"pos"=1:length(tx.type),"x"=1)
lab.df$tx.type <- factor(lab.df$tx.type,levels = unique(tx.type))
tmp.p1 <- ggplot(lab.df, aes(x = x,y = pos)) +
  geom_tile(aes(fill = tx.type)) +
  scale_fill_nejm_adaptive()+
  theme_void() +
  ylab('') + xlab('') +
  theme(aspect.ratio = 15,
        strip.background = element_blank(), #element_rect(color = NA,fill = 'grey'),
        strip.text = element_blank(),
        axis.title = element_blank(),
        # axis.title.y = element_text(size=20),
        # axis.text.x = element_text(size=20,angle = 0,hjust = 0.5),
        axis.text = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position =  "right",#c(0.9,0.8),
        legend.text = element_text(size= 16),
        legend.title= element_text(size= 16))
tmp.p1
#factor: down->up
ggsave(plot = tmp.p1,filename = paste0("tRNA_heatLab_peak.pdf"),width = 6,height = 10)
#cowplot::plot_grid(plotlist = list(tmp.p1,tmp.p2),rel_widths = c(0.1,1), nrow = 1, align = "hv")




## 2.read cov plot
#tx
up=0
down=0
bin=3
ratio=1

cmb <- data.frame(signal=c(paste0("output/","GSE71008_NCpool","/call_peak_all/tbed_RNA_EM/","NCpool",".bed.gz")),
                  region=c(paste0("exSeek-dev/genome/hg38/tbed/tRNA.bed")),
                  signal.label=c("ex sm"),
                  region.label=c("tRNA") 
)


#one way using sapply
res.list <- lapply(1:nrow(cmb), function(j) {get.mat(signal=as.character(cmb$signal[j]), region=as.character(cmb$region[j]),
                                                     signal.label=as.character(cmb$signal.label[j]),region.label=as.character(cmb$region.label[j]),
                                                     up = up, down = down, bin = bin, ratio = ratio)} )
res.list.df <- do.call(rbind,res.list)
rownames(res.list.df) <- res.list.df$bed


#plot in R
#res.list.df$signal <- factor(res.list.df$signal,levels = c("cellular CLIP", "cellular small","extracellular small"))
#(500+500)/20+10=60,再乘以1个样本,再乘以1个区域,为60列
clean_refp <- dplyr::as_tibble(res.list.df)
#clean_refp <- refp #%>% select(c(-1:-3,-5,-6))
colnames(clean_refp)[c(1,3)] <- c('V4','sample')

#merge same tRNA
# tx.id <- unlist(sapply(strsplit(clean_refp$V4,"_"),"[",2))
# tx.type <- gtf$transcript_type[match(tx.id,gtf$transcript_id)]
# table(tx.type)
# clean_refp <- clean_refp[tx.type=="Pro",]

#max-min scale
clean_refp[,4:ncol(clean_refp)] <- t(maxmin.normalize(t(clean_refp[,4:ncol(clean_refp)])))
long_refp <- reshape2::melt(clean_refp,id.vars = c('V4','sample','region'),value.name = 'signal') %>%
  select(-variable)
#add sample name
#long_refp$sample <- rep(c("FTA-10","csFTA-10"),
#                        each = nrow(refp)*22) # (100 + 100 + 20)/10
# add x position
long_refp$pos <- rep(c(1:(ncol(clean_refp)-3)),each = nrow(clean_refp)) # ,times = nrow(cmb)
long_refp[1:3,]
long_refp$tx.id <- unlist(sapply(strsplit(long_refp$V4,"_"),"[",2))
long_refp$tx.type <- gtf$transcript_type[match(long_refp$tx.id,gtf$transcript_id)]

#aver each amino acid. type
filnal_scaler <- long_refp %>%
  # remove na
  tidyr::drop_na() %>%
  dplyr::group_by(pos,sample,region,tx.type) %>% # sample,
  # mean and 95% interval confidence
  dplyr::summarise(mean_signal = mean(signal), #,trim = 0.05, na.rm=T
                   sd = sd(signal),
                   upper = CI(signal,ci = 0.95)[1],
                   lower = CI(signal,ci = 0.95)[3]
  )

#2nd max-min scale (only display mean line, no ribbon)
filnal_scaler <- filnal_scaler %>% 
  dplyr::group_by(sample,region,tx.type) %>%
  dplyr::mutate(mean_signal_norm= (mean_signal - min(mean_signal))/(max(mean_signal)-min(mean_signal))^as.logical(sd(mean_signal)) )
table(filnal_scaler$tx.type)
filnal_scaler[1:3,]

#heatmap (need further optimization)
up_bins <- as.integer(up/bin) # ceiling
down_bins <- as.integer(down/bin) # ceiling
target_bins <- 20 #as.integer((up_bins + down_bins)*(ratio/(1-ratio)))

#filnal_scaler_wide[1:3,1:9]
filnal_scaler_wide <- filnal_scaler %>% 
  tidyr::pivot_wider(id_cols = c("sample","region","tx.type"),names_from = "pos",values_from = "mean_signal_norm")

mat <- as.data.frame(filnal_scaler_wide[,4:ncol(filnal_scaler_wide)])
rownames(mat) <- filnal_scaler_wide$tx.type
hc <- hclust(dist(mat))
#peak.order <- hc$order
table(rownames(mat)==hc$labels)
pdf("tRNA_hclust.pdf",width = 36,height = 6)
plot(hc, hang = -1)
dev.off()
#factor: left->right
#rev( hc2$labels[hc2$order])
#long_refp2 <- long_refp
#table(long_refp2$V4 %in% hc2$labels[hc2$order])
filnal_scaler$tx.type <- factor(filnal_scaler$tx.type,levels = hc2$labels[hc2$order])
# # clr <- c("#2C68A9","#843692","#FD6905")
# # datatypes <- as.character(unique(long_refp2$sample)) # "extracellular small" "cellular small"      "cellular CLIP"  
# # for (i in 1:length(datatypes)){
# #   datatype <- datatypes[i]
# #   print(datatype)
# clr <- "salmon" # purple,salmon,firebrick
# table(long_refp2$sample)
# tmp.p2 <- ggplot(long_refp2, # %>% filter(sample == 'H3K27ac')
#                 aes(x = pos,y = V4)) +
#   geom_tile(aes(fill = signal)) +
#   theme_void() +
#   scale_x_continuous(breaks = c(up_bins+1,up_bins+target_bins), # need +1 for 1-based coord
#                      labels = c('Start 5','End 3') ) +
#   scale_fill_gradient(low = 'white',high = clr) + # 'salmon'
#   # scale_fill_viridis_c() +
#   # scale_fill_continuous("blues") +
#   ylab('') + xlab('') +
#   theme(aspect.ratio = 1,
#         strip.background = element_blank(), #element_rect(color = NA,fill = 'grey'),
#         strip.text = element_text(size=20),
#         axis.title.x = element_text(size=20),
#         axis.title.y = element_text(size=20),
#         axis.text.x = element_text(size=20,angle = 0,hjust = 0.5),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         legend.position =  "right",#c(0.9,0.8),
#         legend.text = element_text(size= 16),
#         legend.title= element_text(size= 16)) #+
# #factor: down->up
# tmp.p2
library(ggridges)
#long_refp2[1:3,]
#table(filnal_scaler$pos)
tmp.p2 <- ggplot(filnal_scaler, aes(x = pos, y=tx.type, height=mean_signal_norm, fill=..x..))+ #
  # geom_ridgeline_gradient(width=0.01) + # stat = "summary"
  geom_density_ridges_gradient(stat = "identity",scale = 1.3, na.rm = T, rel_min_height = 0.001, gradient_lwd = 0, color="white", alpha=0.6) +  # 
  # theme_ridges() +
  # theme(legend.position = "none")
  
  ggraph::scale_fill_viridis(option = "C", direction = 1, alpha = 0.8, begin = 0, end = 0.6) + # name = "value",
  # scale_fill_manual(name = '', values = alpha(colour = c("#2C68A9","#843692","#FD6905"),alpha = 0.8)) +
  # scale_color_manual(name = 'Data type', values = c("#FD6905","#843692","#2C68A9")) +
  geom_vline(xintercept = c((1+up_bins+target_bins+down_bins)*0.5),color="black",linetype="dashed") + # need +1 for 1-based coord
  # x label
  # scale_x_continuous(breaks = c(1,up_bins+1,up_bins+target_bins,up_bins+target_bins+down_bins), # need +1 for 1-based coord
  #                    labels = c(paste0("-",up,"nt"),'Start','End',paste0("+",down,"nt")) ) +
  scale_x_continuous(breaks = c(up_bins+1,up_bins+target_bins), # need +1 for 1-based coord
                     labels = c("Start 5'","End 3'") ) +
  # scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1.00), limits = c(0,1),labels = c("0","0.25","0.50","0.75","1.00")) +
  xlab('') + ylab('Min-max scaled depth') +
  # theme_bw() +
  # ylim(c(1.7,4)) + # need change if change dst !!!
  theme_classic() +
  theme(aspect.ratio = 1,
        axis.ticks.y=element_blank(),
        strip.background = element_rect(color = NA,fill = 'grey'),
        strip.text = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20,angle = 0,hjust = 0.5),
        axis.text.y = element_blank(), #element_text(size=20),
        legend.position =  "right",#c(0.9,0.8),
        legend.text = element_text(size= 16),
        legend.title= element_text(size= 16))
tmp.p2
ggsave(plot = tmp.p2,filename = paste0("tRNA_heat.pdf"),width = 8,height = 8)

#tx.id <- unlist(sapply(strsplit(hc$labels[hc$order],"_"),"[",2))
tx.type <- hc2$labels[hc2$order] #  gtf$transcript_type[match(tx.id,gtf$transcript_id)]
lab.df <- data.frame("tx.type"=tx.type,"pos"=1:length(tx.type),"x"=1)
lab.df$tx.type <- factor(lab.df$tx.type,levels = unique(tx.type))
lab.df[1:3,]
tmp.p1 <- ggplot(lab.df, aes(x = x,y = pos)) +
  geom_tile(aes(fill = tx.type)) +
  scale_fill_nejm_adaptive()+
  theme_void() +
  ylab('') + xlab('') +
  theme(aspect.ratio = 15,
        strip.background = element_blank(), #element_rect(color = NA,fill = 'grey'),
        strip.text = element_blank(),
        axis.title = element_blank(),
        # axis.title.y = element_text(size=20),
        # axis.text.x = element_text(size=20,angle = 0,hjust = 0.5),
        axis.text = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position =  "right",#c(0.9,0.8),
        legend.text = element_text(size= 16),
        legend.title= element_text(size= 16))
tmp.p1
#factor: down->up
ggsave(plot = tmp.p1,filename = paste0("tRNA_heatLab.pdf"),width = 6,height = 10)
#cowplot::plot_grid(plotlist = list(tmp.p1,tmp.p2),rel_widths = c(0.1,1), nrow = 1, align = "hv")










## EM reads length distribution (ridge plot)
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
region.gr <- data.table::fread(as.character(paste0(pre,"/output/","GSE71008_NCpool","/call_peak_all/tbed_RNA_EM/","NCpool",".bed.gz")))
region.gr <- GRanges(region.gr$V1, IRanges(region.gr$V2, region.gr$V3), name=region.gr$V4, score=region.gr$V5, strand=region.gr$V6)
df <- as.data.frame(region.gr)
#0-based
df$score <- 1

df$tx.type <- gtf$transcript_type[match(df$seqnames,gtf$transcript_id)]
#table(is.na(df$tx.type) )
df <- df[!is.na(df$tx.type),]
df$width <- df$end-df$start
#
tmp.mat <- as_tibble(df) %>% 
  dplyr::group_by(width,tx.type) %>% 
  dplyr::summarise(value=sum(score)) %>% 
  pivot_wider(names_from = "width",values_from = "value")
tmp.mat[is.na(tmp.mat)] <- 0
tmp.mat <- as.data.frame(tmp.mat)
rownames(tmp.mat) <- tmp.mat[,"tx.type"]
tmp.mat <- tmp.mat[,-1]

hc2 <- hclust(dist(tmp.mat))
pdf("readLength_hclust.pdf",width = 12,height = 6)
plot(hc2, hang = -1)
dev.off()

library(ggridges)
library(ggplot2)
library(fontawesome)
table(df$tx.type)
str(df)
df$tx.type <- factor(df$tx.type,levels = hc2$labels[hc2$order])
p1 <- ggplot(df, aes(x = width,y=tx.type, fill = ..x..))+ #
  geom_density_ridges_gradient(scale = 1.5, rel_min_height = 0.01, gradient_lwd = 0)+ #stat="binline",  bins = 400,
  #geom_ridgeline_gradient()
  xlab("")+
  ylab(paste0("Reads length Density "))+
  #geom_vline(xintercept = 0,linetype="dashed",color="grey10")+
  #xlim(c(0,500))+
  # xlim(c(0,500))+
  #geom_hline(yintercept = c(0))+
  ggraph::scale_fill_viridis(option = "C",direction = -1) + # name = "value",
  scale_color_manual(values = "white") +
  #?geom_density_ridges_gradient
  theme_bw() + 
  theme(plot.title = element_text(size = 24,color="black",hjust = 0.5,face="bold"),
        axis.title = element_text(size = 16,color ="black",face="bold"), 
        axis.text = element_text(size= 20,color = "black"),
        #panel.grid=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "black"), #size= 1,
        #panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(), # angle = 45, hjust = 1 
        legend.position = "none",
        legend.text = element_text(size= 16),
        legend.title= element_text(size= 16))
p1
ggplot2::ggsave(paste0("readLength_ridge",".pdf"),plot = p1, width = 5, height = 8)





## unitas tRFs reads meta distribution
tmp <- read.table("/BioII/lulab_b/baopengfei/biosoft/unitas/UNITAS_03-01-2023_tRNA.sam_#1/unitas.tRF-table.txt",skip = 3, check.names = F, sep = "\t",header = T,row.names = 1)
tmp <- tmp[,!grepl("absolute",colnames(tmp))]
#colnames(tmp) <- gsub(" (fractionated)","",colnames(tmp),perl = F)
colnames(tmp) <- unlist(sapply(strsplit(colnames(tmp)," "),"[",1))
unique(rownames(tmp))
tmp$id <- rownames(tmp)
tmp$tx.type <- unlist(sapply(strsplit(tmp$id,"-"),"[",2))
table(tmp$tx.type)

colnames(tmp)
df2 <- as_tibble(tmp) %>% 
  tidyr::pivot_longer(cols = "5p-tR-halves":"misc-tRFs",names_to = "frag.type",values_to = "count") %>% 
  dplyr::group_by(tx.type,frag.type) %>% 
  dplyr::summarize(count=sum(count))
df22 <- df2[df2$tx.type %in% c("Gly","Glu"),]
sum(df22$count)
tmp.mat2 <- df2 %>% 
  tidyr::pivot_wider(names_from = "frag.type", values_from = "count")
tmp.mat2 <- as.data.frame(tmp.mat2)
rownames(tmp.mat2) <- tmp.mat2$tx.type
tmp.mat2 <- tmp.mat2[,-1]

hc3 <- hclust(dist(tmp.mat2))
pdf("tRNAfrag_hclust.pdf",width = 8,height = 4)
plot(hc3, hang = -1)
dev.off()

df2$tx.type <- factor(df2$tx.type,levels = hc3$labels[hc3$order])
df2$frag.type <- factor(df2$frag.type,levels = c("misc-tRFs","5p-tR-halves","3p-tR-halves","3p-CCA-tRFs","5p-tRFs","3p-tRFs","tRNA-leader","tRF-1"))
#unique(df2$frag.type)

ggplot(df2,aes(x=count,y=tx.type,fill=frag.type))+ 
  geom_bar(position = "fill",width = 0.8,linewidth=0.3, stat = "identity", color="black") +
  #geom_errorbar()+
  xlab("")+
  ylab(paste0(""))+
  geom_vline(xintercept = 0,linetype="dashed",color="grey10")+
  #ggsci::scale_fill_d3()+
  # scale_fill_nejm_adaptive(name="Species")+
  scale_fill_manual(values = c("grey30","steelblue","firebrick","grey90","steelblue2","firebrick2","steelblue1","salmon"))+
  theme_void() + 
  theme(plot.title = element_text(size = 24,color="black",hjust = 0.5,face="bold"),
        axis.title = element_text(size = 24,color ="black"), 
        axis.text = element_text(size= 20,color = "black"),
        axis.text.x = element_text(size= 20,color = "black",vjust = 0.5), # angle = 90,,hjust = 1
        #panel.grid=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey50",linetype = "dashed"), #size= 1,
        #panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        #axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, color = c(rep("grey50",length(other)),rep("steelblue",length(rna)),rep("firebrick",length(dna)))), #  
        legend.position = "right",#c(.25,.6),
        legend.text = element_text(size= 16),
        legend.title= element_text(size= 16))
ggsave("./unitas_frag_reads_ratio.pdf",width = 8,height = 8)
ggplot(df2,aes(x=count,y=tx.type))+ 
  geom_bar(width = 0.8,linewidth=0.3, stat = "identity", fill="grey50", color="grey50") + # position = "fill",
  #geom_errorbar()+
  xlab("")+
  ylab(paste0(""))+
  geom_vline(xintercept = 0,linetype="dashed",color="grey10")+
  #ggsci::scale_fill_d3()+
  # scale_fill_nejm_adaptive(name="Species")+
  scale_fill_manual(values = c("grey30","steelblue","firebrick","grey90","steelblue2","firebrick2","steelblue1","salmon"))+
  theme_void() + 
  theme(plot.title = element_text(size = 24,color="black",hjust = 0.5,face="bold"),
        axis.title = element_text(size = 24,color ="black"), 
        axis.text = element_text(size= 20,color = "black"),
        axis.text.x = element_text(size= 20,color = "black",vjust = 0.5), # angle = 90,,hjust = 1
        #panel.grid=element_blank(),
        panel.grid.major=element_blank(),
        # panel.grid.minor.x = element_blank(),
        # panel.grid.major.y = element_line(color = "grey50",linetype = "dashed"), #size= 1,
        #panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        #axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, color = c(rep("grey50",length(other)),rep("steelblue",length(rna)),rep("firebrick",length(dna)))), #  
        legend.position = "right",#c(.25,.6),
        legend.text = element_text(size= 16),
        legend.title= element_text(size= 16))
ggsave("./unitas_frag_reads_num.pdf",width = 5,height = 8)




## find most abundant fragment (2017 Tosor, NAR, dimer)
cat /BioII/lulab_b/baopengfei/biosoft/unitas/UNITAS_03-01-2023_tRNA.sam_\#1/unitas.full_annotation_matrix.txt | grep "Gly-GCC" | grep "5p-tR-half" \
| cut -f 1,2 | sort -k2 -n | tail
>Gly-GCC (5p-tR-half)
GCATTGGTGGTTCAGTGGTAGAATTCTCGCC # 24348 reads

cat /BioII/lulab_b/baopengfei/biosoft/unitas/UNITAS_03-01-2023_tRNA.sam_\#1/unitas.full_annotation_matrix.txt | grep "Glu-CTC" | grep "5p-tR-half" \
| cut -f 1,2 | sort -k2 -n | tail

>Glu-CUC (5p-tR-half)
TCCCTGGTGGTCTAGTGGTTAGGATTCGGCGCT # 43683 reads



#homo dimer in forna webserver
>Gly-GCC(5p-tR-half)_1
GCATTGGTGGTTCAGTGGTAGAATTCTCGCC
>Gly-GCC(5p-tR-half)_2
GCATTGGTGGTTCAGTGGTAGAATTCTCGCC
#(color)
>Gly-GCC(5p-tR-half)_1
2-17:steelblue 18-27:salmon 29-30:salmon 1,28,31:grey
>Gly-GCC(5p-tR-half)_2
2-17:steelblue 18-27:salmon 29-30:salmon 1,28,31:grey


#heter dimer in forna webserver
>Gly-GCC(5p-tR-half)
GCATTGGTGGTTCAGTGGTAGAATTCTCGCC
>Glu-CUC(5p-tR-half)
TCCCTGGTGGTCTAGTGGTTAGGATTCGGCGCT
#(color)
>Gly-GCC(5p-tR-half)
2-17:steelblue 18-19:salmon 21-27:salmon 29-31:salmon 1,20,28,31:grey
>Glu-CUC(5p-tR-half)
2-17:orange 18-19:salmon 21-27:salmon 29-31:salmon 1,20,28,31-33:grey


#






# suppl fig 3: mir seed cov in exPeak peaks ------------------------------
dst <- "GSE50676_NCpool"


## read ref
ref <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)
table(ref$transcript_type)

#op1: 
# ## get top expressed miR in blood cells
# top100 <- data.table::fread("/lulabdata/baopengfei/shared_reference/tissueAtlas2/hsa_snc_expression.csv",header = T, check.names=F, sep = ",",stringsAsFactors = F)
# # head(top100)
# # table(top100$organ)
# # table(top100$type)
# top100 <- top100[top100$organ=="lymph_node" & top100$type=="mirna",]
# top100 <- top100[order(top100$expression,decreasing = T),]
# top100 <- top100$acc[1:1000] # nrow(top20) # has great impact in results !!!
# 
# ## get seed seqs
# seed <- data.table::fread("/BioII/lulab_b/baopengfei/shared_reference/targetScan/miR_Family_Info.txt",header = T,sep = "\t",stringsAsFactors = F)
# seed <- as_tibble(seed)
# seed$`Seed+m8` <- gsub("U","T",seed$`Seed+m8`) # 9994 total records
# #seed[100:106,]
# sequences <- (unique(seed[seed$`MiRBase ID` %in% top100,"Seed+m8"])) 
# sequences <- sequences$`Seed+m8`
# sequences.rev <- sapply(sequences,function(x) seq_compl(seq_rev(x)) )  # add reverse complement, also reverse direction 
# sequences <- sequences.rev # unique(c(sequences,sequences.rev)) # sequences.rev
# seq.num <- length(sequences)
# sequences <- paste((sequences),collapse = "|")

#op2:
## use Lymphoblastoid CLIP ref (2014 Nature and  2012 plos pathogens)
#https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1002484#s5
utr <- readxl::read_xls(path = "/lulabdata/baopengfei/shared_reference/2012_plospathogen_lymphoblastoid_AGO_CLIP/ppat.1002484.s010.xls", skip = 1,)
cds <- readxl::read_xls(path = "/lulabdata/baopengfei/shared_reference/2012_plospathogen_lymphoblastoid_AGO_CLIP/ppat.1002484.s011.xls", skip = 1,)
seed <- as.data.frame(rbind(utr,cds))
table(seed$`Human miRs`=="NA")
table(seed$`Viral miRs`=="NA")
table(seed$`Site is specific to EBV infected cell`)
table(seed$`Transcript Location`)
seed <- seed[seed$`Human miRs`!="NA" & seed$`Site is specific to EBV infected cell`=="NO",]

expr <- data.table::fread("/lulabdata/baopengfei/shared_reference/2012_plospathogen_lymphoblastoid_AGO_CLIP/ppat.1002484.s015.csv",header = T,sep = ",",stringsAsFactors = F, skip = 1)
table(seed$`Sequence of miRNA-interaction site` %in% expr$ClusterSequence)

for (i in unique(expr$TranscriptLocation)){
  print(paste0(i,": ",length(unique(expr$SeedSequence[expr$TranscriptLocation==i]))))
}
# [1] "5UTR: 275"
# [1] "CDS: 1500"   !!!!!!!
# [1] "Intergenic: 2545"
# [1] "3UTR: 2124"
# [1] "intron: 1877"
# [1] "miRNA: 186"
# [1] "ncTranscript: 97"
# [1] ">250: 257
sequences <- unique(expr$SeedSequence[expr$TranscriptLocation=="3UTR"]) # todo: use union of 3UTR and CDS, filter topN seeds by count
sequences.rev <- sapply(sequences,function(x) seq_compl(seq_rev(x)) )
sequences <- sequences.rev # sequences
seq.num <- length(sequences)
sequences <- paste((sequences),collapse = "|")

## get seed into bed
sum.match.freq.intersect <- function(x){
  # x <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE50676_NCpool/call_peak_all/expeakCNN_by_sample/b5_d50_p1/intersect/NCpool.bed"
  print(x)
  bed <- data.table::fread(paste0(x),data.table = F, header = F,sep = "\t",stringsAsFactors = F)
  bed <- bed[,1:6]
  colnames(bed) <- c("chr","start","end","peak","score","strand")
  #summary(bed$end-bed$start)
  bed <- bed %>% 
    dplyr::filter(end-start<=200 , end-start>=10)
  
  #only keep mRNA
  bed$RNA <- ref$transcript_type[match(bed$chr,ref$transcript_id)]
  bed <- bed[bed$RNA=="mRNA",1:6]
  
  #expand
  bed.exp <- bedtoolsr::bt.slop(b = 150,s = T,
                                i=bed,
                                g = "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/transcriptome_genome_sort_uniq_newTxID") # 50 + seed length
  #rm boundary regions
  bed <- bed.exp[(bed.exp$V3-bed.exp$V2)==(bed$end-bed$start+300),1:6] #
  colnames(bed) <- c("chr","start","end","peak","score","strand")
  bed$score <- 1 # score unify to 1
  
  
  bedtoolsr::bt.getfasta(fo = paste0(x,".tmp.fa"),nameOnly = T, bed = bed, s = T, 
                         fi = "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/fasta_newTxID/combine19.fa")
  fasta <- rtracklayer::import(paste0(x,".tmp.fa"),format = "fasta")
  fasta <- as.data.frame(fasta)
  fasta$id <- rownames(fasta)
  fasta$id <- gsub("\\(-\\)|\\(\\+\\)","",fasta$id)
  #fasta[1:3,]
  
 
  #tmp.list <- list()
  ReadMotifFreqByRow <- function(i){
    # i <- 1
    tmp <- as.data.frame(stringr::str_locate_all(pattern = sequences, fasta$x[i]))
    if(nrow(tmp)!=0){
      tmp$chr <- fasta$id[i]
      # tmp.list[[as.character(i)]] <- tmp
      tmp$start <- bed$start[i] + tmp$start  #note: need convert to full-length tx coordinate !!!
      tmp$end <- bed$start[i] + tmp$end
      return(tmp)
    }
    # return(tmp)
  }
  tmp.df <- do.call(rbind,parallel::mclapply(1:nrow(fasta), ReadMotifFreqByRow, mc.cores = 1) )
  # for(i in 1:nrow(fasta)){
  #   if(i %% 1000 == 0){
  #     print(i)
  #   }
  #   # i <- 1
  #   tmp <- as.data.frame(stringr::str_locate_all(pattern = sequences, fasta$x[i]))
  #   if(nrow(tmp)!=0){
  #     tmp$chr <- fasta$id[i]
  #     tmp.list[[as.character(i)]] <- tmp
  #   }
  #   # return(tmp)
  # }
  tmp.df$seqnames <- bed$chr[match(tmp.df$chr,bed$peak)]
  tmp.df$score <- 1
  tmp.df$strand <- "+"
  tmp.df <- tmp.df[,c(4,1:3,5:6)]
  file.remove(paste0(x,".tmp.fa"))
  return(tmp.df)
}
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
seed.bed <- sum.match.freq.intersect(paste0(pre,"/output/GSE50676_NCpool/call_peak_all/expeakCNN_by_sample/b5_d50_p1/intersect/NCpool.bed"))
#data.table::fwrite(file = paste0(pre, "/output/",dst,"/call_peak_all/expeakCNN_by_sample/b5_d50_p1/intersect/NCpool.bed.miRseedFor.allPeak.top20"),x=seed.bed,quote = F,sep = "\t",row.names = F,col.names = F)
#data.table::fwrite(file = paste0(pre, "/output/",dst,"/call_peak_all/expeakCNN_by_sample/b5_d50_p1/intersect/NCpool.bed.miRseedFor.allPeak.top200"),x=seed.bed,quote = F,sep = "\t",row.names = F,col.names = F)
#data.table::fwrite(file = paste0(pre, "/output/",dst,"/call_peak_all/expeakCNN_by_sample/b5_d50_p1/intersect/NCpool.bed.miRseedRev.allPeak.top1000"),x=seed.bed,quote = F,sep = "\t",row.names = F,col.names = F)
#data.table::fwrite(file = paste0(pre, "/output/",dst,"/call_peak_all/expeakCNN_by_sample/b5_d50_p1/intersect/NCpool.bed.miRseedFor.allPeak.CDS1500"),x=seed.bed,quote = F,sep = "\t",row.names = F,col.names = F)
#data.table::fwrite(file = paste0(pre, "/output/",dst,"/call_peak_all/expeakCNN_by_sample/b5_d50_p1/intersect/NCpool.bed.miRseedRev.allPeak.CDS1500"),x=seed.bed,quote = F,sep = "\t",row.names = F,col.names = F)
data.table::fwrite(file = paste0(pre, "/output/",dst,"/call_peak_all/expeakCNN_by_sample/b5_d50_p1/intersect/NCpool.exp.bed.miRseedRev.mRNAPeak.3UTR2124"),x=seed.bed,quote = F,sep = "\t",row.names = F,col.names = F)



## enrichedHeatmap
setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC")
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggsci)
library(Rmisc)
library(conflicted)
conflict_prefer("select", "dplyr")

#tx
up=100
down=100
bin=5 # 10 for test
ratio=0.3

#use NCpool
expeak.all <- read.table(paste0(pre, "/output/",dst,"/call_peak_all/expeakCNN_by_sample/b5_d50_p1/intersect/NCpool.bed6"),header = F)
expeak.only <- read.table(paste0(pre, "/output/",dst,"/call_peak_all/expeakCNN_by_sample/b5_d50_p1/intersect/NCpool_expeakCNNOnly.bed6"),header = F)
table((expeak.all$V4 %in% expeak.only$V4))
expeak.other <- expeak.all[!(expeak.all$V4 %in% expeak.only$V4),]
expeak.other[1:3,]
expeak.other$V5 <- 1
expeak.only$V5 <- 1
data.table::fwrite(file = paste0(pre, "/output/",dst,"/call_peak_all/expeakCNN_by_sample/b5_d50_p1/intersect/NCpool_expeakCNNother_score1.bed6"),x=expeak.other,quote = F,sep = "\t",row.names = F,col.names = F)
data.table::fwrite(file = paste0(pre, "/output/",dst,"/call_peak_all/expeakCNN_by_sample/b5_d50_p1/intersect/NCpool_expeakCNNonly_score1.bed6"),x=expeak.only,quote = F,sep = "\t",row.names = F,col.names = F)

cmb <- data.frame(
  signal=c( paste0(pre,"/output/",dst,"/call_peak_all/expeakCNN_by_sample/b5_d50_p1/intersect/NCpool.exp.bed.miRseedRev.mRNAPeak.3UTR2124") ) ,
  region=c(paste0(pre, "/output/",dst,"/call_peak_all/expeakCNN_by_sample/b5_d50_p1/intersect/NCpool_expeakCNNother_score1.bed6"),
           paste0(pre, "/output/",dst,"/call_peak_all/expeakCNN_by_sample/b5_d50_p1/intersect/NCpool_expeakCNNonly_score1.bed6") ) ,
  signal.label=c("wbc"), # ,"ex RIP"
  region.label=c("exPeakOther","exPeakOnly") # ,"CLIPper","exPeak"
)

#one way using sapply
res.list <- lapply(1:nrow(cmb), function(j) {get.mat(signal=as.character(cmb$signal[j]), region=as.character(cmb$region[j]),
                                                     signal.label=as.character(cmb$signal.label[j]),region.label=as.character(cmb$region.label[j]),
                                                     up = up, down = down, bin = bin, ratio = ratio)} )
res.list.df <- do.call(rbind,res.list)
# write.table(res.list.df,"tmp/comapre_Fig3_supp_seedFor_allPeak_top20_overlap_meta.txt",row.names = T,col.names = T,sep = "\t",quote = F)
# res.list.df <- read.table("tmp/comapre_Fig3_supp_seedFor_allPeak_top20_overlap_meta.txt",sep = "\t",header = T)
#write.table(res.list.df,"tmp/comapre_Fig3_supp_seedFor_allPeak_top200_overlap_meta.txt",row.names = T,col.names = T,sep = "\t",quote = F)
#res.list.df <- read.table("tmp/comapre_Fig3_supp_seedFor_allPeak_top200_overlap_meta.txt",sep = "\t",header = T)
# write.table(res.list.df,"tmp/comapre_Fig3_supp_seedFor_allPeak_top1000_overlap_meta.txt",row.names = T,col.names = T,sep = "\t",quote = F)
# res.list.df <- read.table("tmp/comapre_Fig3_supp_seedFor_allPeak_top1000_overlap_meta.txt",sep = "\t",header = T)
# write.table(res.list.df,"tmp/comapre_Fig3_supp_seedRev_allPeak_top1000_overlap_meta.txt",row.names = T,col.names = T,sep = "\t",quote = F)
# res.list.df <- read.table("tmp/comapre_Fig3_supp_seedRev_allPeak_top1000_overlap_meta.txt",sep = "\t",header = T)
# write.table(res.list.df,"tmp/comapre_Fig3_supp_seedFor_allPeak_CDS1500_overlap_meta.txt",row.names = T,col.names = T,sep = "\t",quote = F)
# res.list.df <- read.table("tmp/comapre_Fig3_supp_seedFor_allPeak_CDS1500_overlap_meta.txt",sep = "\t",header = T)
# write.table(res.list.df,"tmp/comapre_Fig3_supp_seedRev_allPeak_CDS1500_overlap_meta.txt",row.names = T,col.names = T,sep = "\t",quote = F)
# res.list.df <- read.table("tmp/comapre_Fig3_supp_seedRev_allPeak_CDS1500_overlap_meta.txt",sep = "\t",header = T)
write.table(res.list.df,"tmp/comapre_Fig3_supp_seedRev_mRNAPeakExp_3UTR2124_overlap_meta.txt",row.names = T,col.names = T,sep = "\t",quote = F)
res.list.df <- read.table("tmp/comapre_Fig3_supp_seedRev_mRNAPeakExp_3UTR2124_overlap_meta.txt",sep = "\t",header = T)
res.list.df[1:3,]


## plot in R
table(res.list.df$region)
res.list.df[1:3,]
res.list.df$region <- factor(res.list.df$region,levels = c("exPeakOther","exPeakOnly"))
clean_refp <- dplyr::as_tibble(res.list.df)
colnames(clean_refp)[c(1,3)] <- c('V4','sample')
table(clean_refp$region,clean_refp$sample)

clean_refp <- clean_refp[base::which( apply(clean_refp[,4:ncol(clean_refp)] ,1, sd)>0),] # filter low sd

#1st max-min scale
clean_refp[,4:ncol(clean_refp)] <- t(maxmin.normalize(t(clean_refp[,4:ncol(clean_refp)])))

long_refp <- reshape2::melt(clean_refp,id.vars = c('V4','sample','region'),value.name = 'signal') %>%
  select(-variable)

#add x position
long_refp$pos <- rep(c(1:(ncol(clean_refp)-3)),each = nrow(clean_refp)) # ,times = nrow(cmb)
#dim(long_refp)

#1st loess model 
l <- list()
for(smp in as.character(unique(long_refp$sample))){
  #smp <- "wbc"
  print(smp)
  for (region in as.character(unique(long_refp$region))){
    #region <- "exPeak"
    print(region)
    long_refp_tmp <- long_refp[long_refp$sample==smp & long_refp$region==region,]
    # long_refp_tmp2 <- get.loess(long_refp_tmp[,c("pos","signal")])
    long_refp_tmp2 <- get.loess.smooth(long_refp_tmp[,c("pos","signal")])
    # long_refp_tmp$signal <- long_refp_tmp2$signal
    colnames(long_refp_tmp2)[colnames(long_refp_tmp2)=="signal"] <- "mean_signal"
    long_refp_tmp2$region <- region
    long_refp_tmp2$sample <- smp
    l[[paste0(smp,region)]] <- long_refp_tmp2
    #table(long_refp_tmp$pos==long_refp_tmp2$pos)
  }
}
filnal_scaler <- do.call(rbind,l)

up_bins <- as.integer(up/bin) # ceiling
down_bins <- as.integer(down/bin) # ceiling
target_bins <- round((up_bins + down_bins)*(ratio/(1-ratio)),digits = 0) #as.integer((up_bins + down_bins)*(ratio/(1-ratio)))
filnal_scaler$region <- gsub("exPeakOther","exPeak other",filnal_scaler$region)
filnal_scaler$region <- gsub("exPeakOnly","exPeak only",filnal_scaler$region)
filnal_scaler$region <- factor(filnal_scaler$region,levels = c("exPeak other","exPeak only"))
table(filnal_scaler$region)

#2nd max-min scale (only display mean line, no ribbon)
filnal_scaler <- filnal_scaler %>% 
  dplyr::group_by(sample,region) %>% # no need to add smp for NCpool
  dplyr::mutate(mean_signal_norm= (mean_signal - min(mean_signal))/(max(mean_signal)-min(mean_signal))^as.logical(sd(mean_signal)) )

#2nd loess (loess after max-min-scale to aviod ridges-line)
filnal_scaler <- filnal_scaler %>%
  dplyr::group_by(sample,region) %>%
  dplyr::mutate(mean_signal_norm2 = get.loess.vec(pos,mean_signal_norm)) # fast, no need for loess.smooth
hist(filnal_scaler$mean_signal_norm2)
filnal_scaler$mean_signal_norm2[filnal_scaler$mean_signal_norm2>1] <- 1
filnal_scaler$mean_signal_norm2[filnal_scaler$mean_signal_norm2<0] <- 0

#
table(filnal_scaler$sample)
filnal_scaler$sample <- factor(filnal_scaler$sample, levels = c("wbc"))
#filnal_scaler <- filnal_scaler[filnal_scaler$region %in% "exPeak",]
ggplot(filnal_scaler,aes(x = pos,y = mean_signal_norm2)) +
  # add 0.95 interval
  # geom_ribbon(aes(ymin = lower,
  #                 ymax = upper,
  #                 fill = region), # sample
  #             alpha = 0.5) +
  geom_line(aes(color = region), size = 2) + # 
  theme_classic(base_size = 16) +
  scale_color_manual(values = c("grey50", pal_d3()(4)[4])) + 
  # scale_color_d3(name = 'Data type') +
  scale_fill_d3(name = '') +
  geom_vline(xintercept = c(up_bins+1,up_bins+target_bins),color="grey50",linetype="dashed")+
  # geom_vline(xintercept = c((1+up_bins+target_bins+down_bins)*0.5),color="grey50",linetype="dashed")+ # need +1 for 1-based coord
  # x label
  scale_x_continuous(breaks = c(1,up_bins+1,up_bins+target_bins,up_bins+target_bins+down_bins), # need +1 for 1-based coord
                     labels = c(paste0("-",up,"nt"),'Start','End',paste0("+",down,"nt")) ) +
  scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1.0), limits = c(0,1), labels = c("0","0.25","0.50","0.75","1.00")) +
  xlab('') + ylab('Scaled depth') +
  theme(aspect.ratio = 0.6,
        strip.background = element_rect(color = NA,fill = 'white'),
        strip.text = element_text(size=30),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=30),
        axis.text.x = element_blank(), #  element_text(size=20,angle = 0,hjust = 0.5)
        axis.text.y = element_text(size=30),
        legend.position =  "right",#c(0.9,0.8),
        legend.text = element_text(size= 24),
        legend.title= element_text(size= 30) ) + 
  facet_grid(#.~sample, 
    cols = vars(sample), 
    scales = "free_y")
#ggsave("seedFor_allPeak_top20_cov_meta.pdf",width = 20,height = 8)
#ggsave("seedFor_allPeak_top200_cov_meta.pdf",width = 20,height = 8)
#ggsave("seedAll_allPeak_top1000_cov_meta.pdf",width = 20,height = 8)
#ggsave("seedFor_allPeak_top1000_cov_meta.pdf",width = 20,height = 8)
#ggsave("seedFor_allPeak_CDStop1500_cov_meta.pdf",width = 20,height = 8)
#ggsave("seedRev_allPeak_CDStop1500_cov_meta.pdf",width = 20,height = 8)
ggsave("exp_seedRev_mRNAPeak_3UTR2124_cov_meta.pdf",width = 20,height = 8)





## heatmap
long_refp2 <- long_refp 
hist(long_refp2$signal)

table(long_refp2$signal==0)
clr <- c("grey50", pal_d3()(4)[4])  #c("#1B61A5","#FC6910","#269321")
methods <- as.character(unique(long_refp2$region)) # 
long_refp2[1:3,]
table(long_refp2$pos)
#methods <- "exPeak"
for (i in 1:length(methods)){
  #i <- 1
  method <- methods[i]
  print(method)
  long_refp2_tmp <- long_refp2[long_refp2$region==method,]
  
  plot.list <- list()
  for (sample in as.character(unique(long_refp2$sample))){ # as.character(unique(long_refp2_tmp$sample))
    # sample <- "wbc"
    print(sample)
    long_refp2_tmp2 <- long_refp2_tmp[long_refp2_tmp$sample==sample,]
    long_refp_sort <- long_refp2_tmp2 %>%
      # dplyr::filter(pos %in% 20:37) %>% # only rank/sort by central regions
      dplyr::group_by(V4) %>%  # sample,
      #mean of each row/region
      dplyr::summarise(mean_signal_region = mean(signal,na.rm=T)
      )
    long_refp_sort <- long_refp_sort[order(long_refp_sort$mean_signal_region,decreasing = F),]
    long_refp2_tmp3 <- long_refp2_tmp2
    long_refp2_tmp3$V4 <- factor(long_refp2_tmp3$V4,levels = unique(long_refp_sort$V4))
    table(long_refp2_tmp3$sample)
    #table(is.na(long_refp2_tmp3$signal))
    table(long_refp2_tmp3$pos)
    tmp.p <- ggplot(long_refp2_tmp3, aes(x = pos,y = V4)) +
      geom_tile(aes(fill = signal)) + # 
      theme_minimal() +
      # coord_cartesian(expand = 0) +
      scale_x_continuous(breaks = c(1,up_bins+0.5,up_bins+target_bins+0.5,up_bins+target_bins+down_bins), # need +1 for 1-based coord
                         limits = c(1,up_bins+target_bins+down_bins),
                         labels = c(paste0("-",up,"nt"),'Start','End',paste0("+",down,"nt"))) +
      scale_fill_gradient(low = 'white',high = clr[i]) + # , na.value = "white")
      ylab('') + xlab('') +
      theme(
        aspect.ratio = 2.5,
        plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
        # plot.background = element_rect(fill = "white"),
        # panel.background = element_rect(fill = "white", colour = "white"),
        # strip.background = element_rect(color = NA,fill = 'white'),
        strip.text = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20,angle = 0,hjust = 0.5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position =  "none",#c(0.9,0.8),
        legend.text = element_text(size= 16),
        legend.title= element_text(size= 16)) #+
    # theme(aspect.ratio = 1.5,
    # plot.margin = unit(c(5.5, 5.5, 0, 5.5), "pt")
    #       panel.grid = element_blank(),
    #       axis.text.x = element_text(size=20,angle = 0,hjust = 0.5),
    #       # axis.text.y = element_text(size=20),
    #       legend.position = "right",
    #       legend.text = element_text(size= 16),
    #       legend.title= element_text(size= 16))+
    # facet_grid( . ~ sample ,  scales = "free_y" ) # region, space = "free_y",
    plot.list[[paste0(method,"_",sample)]] <- tmp.p
  }
  # p <- cowplot::plot_grid(plotlist = plot.list, rel_widths = c(1,1,1), nrow = 1, align = "hv", axis = "b")
  p <- ggpubr::ggarrange(plotlist = plot.list, rel_widths = c(1,1,1), nrow = 1,  common.legend = T,align = "hv", axis = "b", legend = "right") #  labels = c("RBPs", "", "EV", "", "G4"),
  ggsave(plot = p, paste0(method,"_exp_mRNA_lymMirSeedRev_cov_heat_v2.pdf"), width=12, height=15) # "_",sample
}
#




# suppl fig 3: mir seed cov in exPeak peaks (test, deprecated?) ------------------------------
dst <- "GSE50676_NCpool"


## read ref
ref <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)
table(ref$transcript_type)

#op1: 
# ## get top expressed miR in blood cells
# top100 <- data.table::fread("/lulabdata/baopengfei/shared_reference/tissueAtlas2/hsa_snc_expression.csv",header = T, check.names=F, sep = ",",stringsAsFactors = F)
# # head(top100)
# # table(top100$organ)
# # table(top100$type)
# top100 <- top100[top100$organ=="lymph_node" & top100$type=="mirna",]
# top100 <- top100[order(top100$expression,decreasing = T),]
# top100 <- top100$acc[1:1000] # nrow(top20) # has great impact in results !!!
# 
# ## get seed seqs
# seed <- data.table::fread("/BioII/lulab_b/baopengfei/shared_reference/targetScan/miR_Family_Info.txt",header = T,sep = "\t",stringsAsFactors = F)
# seed <- as_tibble(seed)
# seed$`Seed+m8` <- gsub("U","T",seed$`Seed+m8`) # 9994 total records
# #seed[100:106,]
# sequences <- (unique(seed[seed$`MiRBase ID` %in% top100,"Seed+m8"])) 
# sequences <- sequences$`Seed+m8`
# sequences.rev <- sapply(sequences,function(x) seq_compl(seq_rev(x)) )  # add reverse complement, also reverse direction 
# sequences <- sequences.rev # unique(c(sequences,sequences.rev)) # sequences.rev
# seq.num <- length(sequences)
# sequences <- paste((sequences),collapse = "|")

#op2:
## use Lymphoblastoid CLIP ref (2014 Nature and  2012 plos pathogens)
#https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1002484#s5
utr <- readxl::read_xls(path = "/lulabdata/baopengfei/shared_reference/2012_plospathogen_lymphoblastoid_AGO_CLIP/ppat.1002484.s010.xls", skip = 1,)
cds <- readxl::read_xls(path = "/lulabdata/baopengfei/shared_reference/2012_plospathogen_lymphoblastoid_AGO_CLIP/ppat.1002484.s011.xls", skip = 1,)
seed <- as.data.frame(rbind(utr,cds))
table(seed$`Human miRs`=="NA")
table(seed$`Viral miRs`=="NA")
table(seed$`Site is specific to EBV infected cell`)
table(seed$`Transcript Location`)
seed <- seed[seed$`Human miRs`!="NA" & seed$`Site is specific to EBV infected cell`=="NO",]

expr <- data.table::fread("/lulabdata/baopengfei/shared_reference/2012_plospathogen_lymphoblastoid_AGO_CLIP/ppat.1002484.s015.csv",header = T,sep = ",",stringsAsFactors = F, skip = 1)
#table(seed$`Sequence of miRNA-interaction site` %in% expr$ClusterSequence)
#expr$ClusterSequence seem to be cluster of target sequences
#expr$SeedSequence seem to be seed of expr$miRNAs
#expr interval seem to be hg19 of miR seed binding sites

for (i in unique(expr$TranscriptLocation)){
  print(paste0(i,": ",length(unique(expr$SeedSequence[expr$TranscriptLocation==i]))))
}
# [1] "5UTR: 275"
# [1] "CDS: 1500"   !!!!!!!
# [1] "Intergenic: 2545"
# [1] "3UTR: 2124"
# [1] "intron: 1877"
# [1] "miRNA: 186"
# [1] "ncTranscript: 97"
# [1] ">250: 257

# expr2 <- expr[expr$TranscriptLocation %in% c("3UTR","CDS"),]
# expr2 <- expr2[order(expr2$ReadCount,decreasing = T),]
# sequences <- head(unique(expr2$SeedSequence),1500) # use union of 3UTR and CDS, filter top1500 seeds by count, seem not as expected
sequences <- unique(expr$SeedSequence[expr$TranscriptLocation=="3UTR"])
sequences.rev <- sapply(sequences,function(x) seq_compl(seq_rev(x)) )
sequences <- sequences.rev # whether or not to use complementary reverse 
seq.num <- length(sequences)
sequences <- paste((sequences),collapse = "|")

## get seed into bed
sum.match.freq.intersect <- function(x){
  # x <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE50676_NCpool/call_peak_all/expeakCNN_by_sample/b5_d50_p1/intersect/NCpool.bed"
  print(x)
  bed <- data.table::fread(paste0(x),data.table = F, header = F,sep = "\t",stringsAsFactors = F)
  bed <- bed[,1:6]
  colnames(bed) <- c("chr","start","end","peak","score","strand")
  #summary(bed$end-bed$start)
  bed <- bed %>% 
    dplyr::filter(end-start<=200 , end-start>=10)

  #only keep mRNA
  bed$RNA <- ref$transcript_type[match(bed$chr,ref$transcript_id)]
  bed <- bed[bed$RNA=="mRNA",1:6] # miRNA lengt not allowed exp 50 nt both
  
  #expand
  bed.exp <- bedtoolsr::bt.slop(b = 58,s = T,
                                i=bed,
                                g = "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/transcriptome_genome_sort_uniq_newTxID") # 50 + seed length
  #rm boundary regions
  bed <- bed.exp[(bed.exp$V3-bed.exp$V2)==(bed$end-bed$start+116),1:6] #
  colnames(bed) <- c("chr","start","end","peak","score","strand")
  bed$score <- 1 # score unify to 1
  
  bedtoolsr::bt.getfasta(fo = paste0(x,".tmp.fa"),nameOnly = T, bed = bed, s = T, 
                         fi = "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/fasta_newTxID/combine19.fa")
  fasta <- rtracklayer::import(paste0(x,".tmp.fa"),format = "fasta")
  fasta <- as.data.frame(fasta)
  fasta$id <- rownames(fasta)
  fasta$id <- gsub("\\(-\\)|\\(\\+\\)","",fasta$id)
  
  
  #tmp.list <- list()
  ReadMotifFreqByRow <- function(i){
    # i <- 1
    tmp <- as.data.frame(stringr::str_locate_all(pattern = sequences, fasta$x[i]))
    if(nrow(tmp)!=0){
      tmp$chr <- fasta$id[i]
      # tmp.list[[as.character(i)]] <- tmp
      tmp$start <- bed$start[i] + tmp$start  #note: need convert to full-length tx coordinate !!!
      tmp$end <- bed$start[i] + tmp$end
      return(tmp)
    }
    # return(tmp)
  }
  tmp.df <- do.call(rbind,parallel::mclapply(1:nrow(fasta), ReadMotifFreqByRow, mc.cores = 1) )
  # for(i in 1:nrow(fasta)){
  #   if(i %% 1000 == 0){
  #     print(i)
  #   }
  #   # i <- 1
  #   tmp <- as.data.frame(stringr::str_locate_all(pattern = sequences, fasta$x[i]))
  #   if(nrow(tmp)!=0){
  #     tmp$chr <- fasta$id[i]
  #     tmp.list[[as.character(i)]] <- tmp
  #   }
  #   # return(tmp)
  # }
  tmp.df$seqnames <- bed$chr[match(tmp.df$chr,bed$peak)]
  tmp.df$score <- 1
  tmp.df$strand <- "+"
  tmp.df <- tmp.df[,c(4,1:3,5:6)]
  file.remove(paste0(x,".tmp.fa"))
  return(tmp.df)
}
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
seed.bed <- sum.match.freq.intersect(paste0(pre,"/output/GSE50676_NCpool/call_peak_all/expeakCNN_by_sample/b5_d50_p1/intersect/NCpool.bed"))
#data.table::fwrite(file = paste0(pre, "/output/",dst,"/call_peak_all/expeakCNN_by_sample/b5_d50_p1/intersect/NCpool.bed.miRseedFor.allPeak.top20"),x=seed.bed,quote = F,sep = "\t",row.names = F,col.names = F)
#data.table::fwrite(file = paste0(pre, "/output/",dst,"/call_peak_all/expeakCNN_by_sample/b5_d50_p1/intersect/NCpool.bed.miRseedFor.allPeak.top200"),x=seed.bed,quote = F,sep = "\t",row.names = F,col.names = F)
#data.table::fwrite(file = paste0(pre, "/output/",dst,"/call_peak_all/expeakCNN_by_sample/b5_d50_p1/intersect/NCpool.bed.miRseedRev.allPeak.top1000"),x=seed.bed,quote = F,sep = "\t",row.names = F,col.names = F)
#data.table::fwrite(file = paste0(pre, "/output/",dst,"/call_peak_all/expeakCNN_by_sample/b5_d50_p1/intersect/NCpool.bed.miRseedFor.allPeak.CDS1500"),x=seed.bed,quote = F,sep = "\t",row.names = F,col.names = F)
#data.table::fwrite(file = paste0(pre, "/output/",dst,"/call_peak_all/expeakCNN_by_sample/b5_d50_p1/intersect/NCpool.bed.miRseedRev.allPeak.CDS1500"),x=seed.bed,quote = F,sep = "\t",row.names = F,col.names = F)
#data.table::fwrite(file = paste0(pre, "/output/",dst,"/call_peak_all/expeakCNN_by_sample/b5_d50_p1/intersect/NCpool.bed.miRseedRev.allPeak.3UTR2124"),x=seed.bed,quote = F,sep = "\t",row.names = F,col.names = F)
#data.table::fwrite(file = paste0(pre, "/output/",dst,"/call_peak_all/expeakCNN_by_sample/b5_d50_p1/intersect/NCpool.bed.miRseedRev.allPeak.3UTRCDS1500"),x=seed.bed,quote = F,sep = "\t",row.names = F,col.names = F)
data.table::fwrite(file = paste0(pre, "/output/",dst,"/call_peak_all/expeakCNN_by_sample/b5_d50_p1/intersect/NCpool.exp.bed.miRseedRev.mRNAPeak.3UTR2124"),x=seed.bed,quote = F,sep = "\t",row.names = F,col.names = F)
#data.table::fwrite(file = paste0(pre, "/output/",dst,"/call_peak_all/expeakCNN_by_sample/b5_d50_p1/intersect/NCpool.bed.miRseedRev.mRNAPeak.3UTR2124"),x=seed.bed,quote = F,sep = "\t",row.names = F,col.names = F)

seed.bed.center <- seed.bed
table(seed.bed.center$start-seed.bed.center$end)
seed.bed.center$start <- as.integer((seed.bed.center$start+seed.bed.center$end) * 0.5)
seed.bed.center$end <- seed.bed.center$start + 1
# all seeds have length <12
seed.bed.center[1:3,]
seed.bed.center.exp <- bedtoolsr::bt.slop(b = 5,s = T, i=seed.bed.center,
                                      g = "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/transcriptome_genome_sort_uniq_newTxID")
seed.bed.center.exp[1:3,]
table(seed.bed.center.exp$V2>=0)
table(seed.bed.center.exp$V5)
#rm boundary regions
table((seed.bed.center.exp$V3-seed.bed.center.exp$V2)==(seed.bed.center$end-seed.bed.center$start+10))
seed.bed.center <- seed.bed.center.exp[(seed.bed.center.exp$V3-seed.bed.center.exp$V2)==(seed.bed.center$end-seed.bed.center$start+10),1:6] #
data.table::fwrite(file = paste0(pre, "/output/",dst,"/call_peak_all/expeakCNN_by_sample/b5_d50_p1/intersect/NCpool.exp.bed.miRseedRev.mRNAPeak.3UTR2124.center"),
                   x=seed.bed.center,quote = F,sep = "\t",row.names = F,col.names = F)


## enrichedHeatmap
setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC")
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggsci)
library(Rmisc)
# library(data.table)
library(conflicted)
conflict_prefer("select", "dplyr")

#tx
up=50
down=50
bin=5 # 10 for test
ratio=0.1 #11/(50*2+11)

#use NCpool
expeak.all <- read.table(paste0(pre, "/output/",dst,"/call_peak_all/expeakCNN_by_sample/b5_d50_p1/intersect/NCpool.bed6"),header = F)
expeak.only <- read.table(paste0(pre, "/output/",dst,"/call_peak_all/expeakCNN_by_sample/b5_d50_p1/intersect/NCpool_expeakCNNOnly.bed6"),header = F)
table((expeak.all$V4 %in% expeak.only$V4))
expeak.other <- expeak.all[!(expeak.all$V4 %in% expeak.only$V4),]
expeak.other[1:3,]
expeak.other$V5 <- 1
expeak.only$V5 <- 1
data.table::fwrite(file = paste0(pre, "/output/",dst,"/call_peak_all/expeakCNN_by_sample/b5_d50_p1/intersect/NCpool_expeakCNNother_score1.bed6"),x=expeak.other,quote = F,sep = "\t",row.names = F,col.names = F)
data.table::fwrite(file = paste0(pre, "/output/",dst,"/call_peak_all/expeakCNN_by_sample/b5_d50_p1/intersect/NCpool_expeakCNNonly_score1.bed6"),x=expeak.only,quote = F,sep = "\t",row.names = F,col.names = F)


cmb <- data.frame(
  signal=c(paste0(pre, "/output/",dst,"/call_peak_all/expeakCNN_by_sample/b5_d50_p1/intersect/NCpool_expeakCNNother_score1.bed6"),
           paste0(pre, "/output/",dst,"/call_peak_all/expeakCNN_by_sample/b5_d50_p1/intersect/NCpool_expeakCNNonly_score1.bed6") ) ,
  region=c( paste0(pre,"/output/",dst,"/call_peak_all/expeakCNN_by_sample/b5_d50_p1/intersect/NCpool.exp.bed.miRseedRev.mRNAPeak.3UTR2124.center") ) ,
  
  signal.label=c("exPeak","exPeakOnly"), # ,"CLIPper","exPeak"
  region.label=c("wbc") # ,"ex RIP"
)

#one way using sapply
res.list <- lapply(1:nrow(cmb), function(j) {get.mat(signal=as.character(cmb$signal[j]), region=as.character(cmb$region[j]),
                                                     signal.label=as.character(cmb$signal.label[j]),region.label=as.character(cmb$region.label[j]),
                                                     up = up, down = down, bin = bin, ratio = ratio)} )
res.list.df <- do.call(rbind,res.list)
# write.table(res.list.df,"tmp/comapre_Fig3_supp_seedFor_allPeak_top20_overlap_meta.txt",row.names = T,col.names = T,sep = "\t",quote = F)
# res.list.df <- read.table("tmp/comapre_Fig3_supp_seedFor_allPeak_top20_overlap_meta.txt",sep = "\t",header = T)
#write.table(res.list.df,"tmp/comapre_Fig3_supp_seedFor_allPeak_top200_overlap_meta.txt",row.names = T,col.names = T,sep = "\t",quote = F)
#res.list.df <- read.table("tmp/comapre_Fig3_supp_seedFor_allPeak_top200_overlap_meta.txt",sep = "\t",header = T)
# write.table(res.list.df,"tmp/comapre_Fig3_supp_seedFor_allPeak_top1000_overlap_meta.txt",row.names = T,col.names = T,sep = "\t",quote = F)
# res.list.df <- read.table("tmp/comapre_Fig3_supp_seedFor_allPeak_top1000_overlap_meta.txt",sep = "\t",header = T)
# write.table(res.list.df,"tmp/comapre_Fig3_supp_seedRev_allPeak_top1000_overlap_meta.txt",row.names = T,col.names = T,sep = "\t",quote = F)
# res.list.df <- read.table("tmp/comapre_Fig3_supp_seedRev_allPeak_top1000_overlap_meta.txt",sep = "\t",header = T)
# write.table(res.list.df,"tmp/comapre_Fig3_supp_seedFor_allPeak_CDS1500_overlap_meta.txt",row.names = T,col.names = T,sep = "\t",quote = F)
# res.list.df <- read.table("tmp/comapre_Fig3_supp_seedFor_allPeak_CDS1500_overlap_meta.txt",sep = "\t",header = T)
# write.table(res.list.df,"tmp/comapre_Fig3_supp_seedRev_allPeak_CDS1500_overlap_meta.txt",row.names = T,col.names = T,sep = "\t",quote = F)
# res.list.df <- read.table("tmp/comapre_Fig3_supp_seedRev_allPeak_CDS1500_overlap_meta.txt",sep = "\t",header = T)
# write.table(res.list.df,"tmp/comapre_Fig3_supp_seedRev_allPeak_3UTR2124_overlap_meta.txt",row.names = T,col.names = T,sep = "\t",quote = F)
# res.list.df <- read.table("tmp/comapre_Fig3_supp_seedRev_allPeak_3UTR2124_overlap_meta.txt",sep = "\t",header = T)
# write.table(res.list.df,"tmp/comapre_Fig3_supp_seedRev_allPeak_3UTRCDS1500_overlap_meta.txt",row.names = T,col.names = T,sep = "\t",quote = F)
# res.list.df <- read.table("tmp/comapre_Fig3_supp_seedRev_allPeak_3UTRCDS1500_overlap_meta.txt",sep = "\t",header = T)
# write.table(res.list.df,"tmp/comapre_Fig3_supp_seedRev_mRNAPeakExp_3UTRCDS1500_overlap_meta.txt",row.names = T,col.names = T,sep = "\t",quote = F)
# res.list.df <- read.table("tmp/comapre_Fig3_supp_seedRev_mRNAPeakExp_3UTRCDS1500_overlap_meta.txt",sep = "\t",header = T)
# write.table(res.list.df,"tmp/comapre_Fig3_supp_seedRev_mRNAPeakExp_3UTR2124_overlap_meta.txt",row.names = T,col.names = T,sep = "\t",quote = F)
# res.list.df <- read.table("tmp/comapre_Fig3_supp_seedRev_mRNAPeakExp_3UTR2124_overlap_meta.txt",sep = "\t",header = T)
# write.table(res.list.df,"tmp/comapre_Fig3_supp_seedRevCenter_mRNAPeak_3UTR2124_overlap_meta.txt",row.names = T,col.names = T,sep = "\t",quote = F)
# res.list.df <- read.table("tmp/comapre_Fig3_supp_seedRevCenter_mRNAPeak_3UTR2124_overlap_meta.txt",sep = "\t",header = T)
write.table(res.list.df,"tmp/comapre_Fig3_supp_seedRevCenter_mRNAPeak_3UTR2124_overlap_meta.txt",row.names = T,col.names = T,sep = "\t",quote = F)
res.list.df <- read.table("tmp/comapre_Fig3_supp_seedRevCenter_mRNAPeak_3UTR2124_overlap_meta.txt",sep = "\t",header = T)
# res.list.df[1:3,]


## plot in R
table(res.list.df$region)
table(res.list.df$signal)
res.list.df$signal[res.list.df$signal=='exPeak'] <- "exPeakOther"
res.list.df$signal <- factor(res.list.df$signal,levels = c("exPeakOther","exPeakOnly"))
#res.list.df$region[res.list.df$region=='exPeak'] <- "exPeakOther"
#res.list.df$region <- factor(res.list.df$region,levels = c("exPeakOther","exPeakOnly"))
clean_refp <- dplyr::as_tibble(res.list.df)
clean_refp[1:3,1:5]
colnames(clean_refp)[c(1,3)] <- c('V4','sample')
table(clean_refp$region,clean_refp$sample)

clean_refp <- clean_refp[base::which( apply(clean_refp[,4:ncol(clean_refp)] ,1, sd)>0),] # filter low sd

#1st max-min scale
clean_refp[,4:ncol(clean_refp)] <- t(maxmin.normalize(t(clean_refp[,4:ncol(clean_refp)])))

long_refp <- reshape2::melt(clean_refp,id.vars = c('V4','sample','region'),value.name = 'signal') %>%
  select(-variable)

#add x position
long_refp$pos <- rep(c(1:(ncol(clean_refp)-3)),each = nrow(clean_refp)) # ,times = nrow(cmb)
# dim(long_refp)


#1st loess model 
l <- list()
for(smp in as.character(unique(long_refp$sample))){
  #smp <- "exPeak"
  print(smp)
  for (region in as.character(unique(long_refp$region))){
    #region <- "wbc"
    print(region)
    long_refp_tmp <- long_refp[long_refp$sample==smp & long_refp$region==region,]
    # long_refp_tmp2 <- get.loess(long_refp_tmp[,c("pos","signal")])
    long_refp_tmp2 <- get.loess.smooth(long_refp_tmp[,c("pos","signal")])
    # long_refp_tmp$signal <- long_refp_tmp2$signal
    colnames(long_refp_tmp2)[colnames(long_refp_tmp2)=="signal"] <- "mean_signal"
    long_refp_tmp2$region <- region
    long_refp_tmp2$sample <- smp
    l[[paste0(smp,region)]] <- long_refp_tmp2
    #table(long_refp_tmp$pos==long_refp_tmp2$pos)
  }
}
filnal_scaler <- do.call(rbind,l)


up_bins <- as.integer(up/bin) # ceiling
down_bins <- as.integer(down/bin) # ceiling
target_bins <- round((up_bins + down_bins)*(ratio/(1-ratio)),digits = 0) #as.integer((up_bins + down_bins)*(ratio/(1-ratio)))
filnal_scaler$sample <- gsub("exPeakOther","exPeak other",filnal_scaler$sample)
filnal_scaler$sample <- gsub("exPeakOnly","exPeak only",filnal_scaler$sample)
filnal_scaler$sample <- factor(filnal_scaler$sample,levels = c("exPeak other","exPeak only"))
# filnal_scaler$region <- gsub("exPeakOther","exPeak other",filnal_scaler$region)
# filnal_scaler$region <- gsub("exPeakOnly","exPeak only",filnal_scaler$region)
# filnal_scaler$region <- factor(filnal_scaler$region,levels = c("exPeak other","exPeak only"))
table(filnal_scaler$region,filnal_scaler$sample)

#2nd max-min scale (only display mean line, no ribbon)
filnal_scaler <- filnal_scaler %>% 
  # dplyr::group_by(sample,region) %>% # no need to add smp for NCpool
  dplyr::group_by(sample,region) %>% 
  dplyr::mutate(mean_signal_norm= (mean_signal - min(mean_signal))/(max(mean_signal)-min(mean_signal))^as.logical(sd(mean_signal)) )

#2nd loess (loess after max-min-scale to aviod ridges-line)
filnal_scaler <- filnal_scaler %>%
  dplyr::group_by(sample,region) %>%
  dplyr::mutate(mean_signal_norm2 = get.loess.vec(pos,mean_signal_norm)) # fast, no need for loess.smooth
hist(filnal_scaler$mean_signal_norm2)
table(filnal_scaler$mean_signal_norm2<=0)
filnal_scaler$mean_signal_norm2[filnal_scaler$mean_signal_norm2>1] <- 1
filnal_scaler$mean_signal_norm2[filnal_scaler$mean_signal_norm2<0] <- 0

#
# filnal_scaler$sample <- factor(filnal_scaler$sample, levels = c("wbc"))
filnal_scaler$region <- factor(filnal_scaler$region, levels = c("wbc"))
table(filnal_scaler$region )
#table(filnal_scaler$pos)
#filnal_scaler <- filnal_scaler[filnal_scaler$region %in% "exPeak",]
ggplot(filnal_scaler,aes(x = pos,y = mean_signal_norm2)) +
  # add 0.95 interval
  # geom_ribbon(aes(ymin = lower,
  #                 ymax = upper,
  #                 fill = region), # sample
  #             alpha = 0.5) +
  # geom_line(aes(color = region), size = 2) + # 
  geom_line(aes(color = sample), size = 2) + # 
  theme_classic(base_size = 16) +
  scale_color_manual(values = c("grey50", pal_d3()(4)[4])) + 
  # scale_color_d3(name = 'Data type') +
  scale_fill_d3(name = '') +
  geom_vline(xintercept = c(up_bins+1,up_bins+target_bins),color="grey50",linetype="dashed")+
  # geom_vline(xintercept = c((1+up_bins+target_bins+down_bins)*0.5),color="grey50",linetype="dashed")+ # need +1 for 1-based coord
  # x label
  scale_x_continuous(breaks = c(1,up_bins+1,up_bins+target_bins,up_bins+target_bins+down_bins), # need +1 for 1-based coord
                     labels = c(paste0("-",up,"nt"),'Start','End',paste0("+",down,"nt")) ) +
  scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1.0), limits = c(0,1), labels = c("0","0.25","0.50","0.75","1.00")) +
  xlab('') + ylab('Scaled depth') +
  theme(aspect.ratio = 0.6,
        strip.background = element_rect(color = NA,fill = 'white'),
        strip.text = element_text(size=30),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=30),
        axis.text.x = element_blank(), #  element_text(size=20,angle = 0,hjust = 0.5)
        axis.text.y = element_text(size=30),
        legend.position =  "right",#c(0.9,0.8),
        legend.text = element_text(size= 24),
        legend.title= element_text(size= 30) ) + 
  facet_grid(#.~sample, 
    cols = vars(region), 
    scales = "free_y")

#why as a valley ????

ggsave("seedRev_mRNAPeakExp_3UTR2124_cov_meta.pdf",width = 20,height = 8)





# ## heatmap
# long_refp2 <- long_refp 
# clr <- c("grey50", pal_d3()(4)[4])  #c("#1B61A5","#FC6910","#269321")
# methods <- as.character(unique(long_refp2$region)) # 
# #methods <- "exPeak"
# for (i in 1:length(methods)){
#   #i <- 1
#   method <- methods[i]
#   print(method)
#   long_refp2_tmp <- long_refp2[long_refp2$region==method,]
#   
#   plot.list <- list()
#   for (sample in as.character(unique(long_refp2$sample))){ # as.character(unique(long_refp2_tmp$sample))
#     # sample <- "wbc"
#     print(sample)
#     long_refp2_tmp2 <- long_refp2_tmp[long_refp2_tmp$sample==sample,]
#     long_refp_sort <- long_refp2_tmp2 %>%
#       # dplyr::filter(pos %in% 25:35) %>% # only rank/sort by central regions
#       dplyr::group_by(V4) %>%  # sample,
#       #mean of each row/region
#       dplyr::summarise(mean_signal_region = mean(signal,na.rm=T)
#       )
#     long_refp_sort <- long_refp_sort[order(long_refp_sort$mean_signal_region,decreasing = F),]
#     long_refp2_tmp3 <- long_refp2_tmp2
#     long_refp2_tmp3$V4 <- factor(long_refp2_tmp3$V4,levels = unique(long_refp_sort$V4))
#     table(long_refp2_tmp3$sample)
#     #table(is.na(long_refp2_tmp3$signal))
#     table(long_refp2_tmp3$pos)
#     tmp.p <- ggplot(long_refp2_tmp3, aes(x = pos,y = V4)) +
#       geom_tile(aes(fill = signal)) + # 
#       theme_minimal() +
#       # coord_cartesian(expand = 0) +
#       scale_x_continuous(breaks = c(1,up_bins+0.5,up_bins+target_bins+0.5,up_bins+target_bins+down_bins), # need +1 for 1-based coord
#                          limits = c(1,up_bins+target_bins+down_bins),
#                          labels = c(paste0("-",up,"nt"),'Start','End',paste0("+",down,"nt"))) +
#       scale_fill_gradient(low = 'white',high = clr[i]) + # , na.value = "white")
#       ylab('') + xlab('') +
#       theme(
#         aspect.ratio = 2.5,
#         plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
#         # plot.background = element_rect(fill = "white"),
#         # panel.background = element_rect(fill = "white", colour = "white"),
#         # strip.background = element_rect(color = NA,fill = 'white'),
#         strip.text = element_text(size=20),
#         axis.title.x = element_text(size=20),
#         axis.title.y = element_text(size=20),
#         axis.text.x = element_text(size=20,angle = 0,hjust = 0.5),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         legend.position =  "none",#c(0.9,0.8),
#         legend.text = element_text(size= 16),
#         legend.title= element_text(size= 16)) #+
#     # theme(aspect.ratio = 1.5,
#     # plot.margin = unit(c(5.5, 5.5, 0, 5.5), "pt")
#     #       panel.grid = element_blank(),
#     #       axis.text.x = element_text(size=20,angle = 0,hjust = 0.5),
#     #       # axis.text.y = element_text(size=20),
#     #       legend.position = "right",
#     #       legend.text = element_text(size= 16),
#     #       legend.title= element_text(size= 16))+
#     # facet_grid( . ~ sample ,  scales = "free_y" ) # region, space = "free_y",
#     plot.list[[paste0(method,"_",sample)]] <- tmp.p
#   }
#   # p <- cowplot::plot_grid(plotlist = plot.list, rel_widths = c(1,1,1), nrow = 1, align = "hv", axis = "b")
#   p <- ggpubr::ggarrange(plotlist = plot.list, rel_widths = c(1,1,1), nrow = 1,  common.legend = T,align = "hv", axis = "b", legend = "right") #  labels = c("RBPs", "", "EV", "", "G4"),
#   ggsave(plot = p, paste0(method,"_lymMirSeed_cov_heat_v2.pdf"), width=12, height=15) # "_",sample
# }
# #



# test ------------------------------------------------------------
library(EnrichedHeatmap)
#BiocManager::install("EnrichedHeatmap")

set.seed(123)
load(system.file("extdata", "chr21_test_data.RData", package = "EnrichedHeatmap"))
#ls()
#rpkm
genes@elementType
tss = IRanges::promoters(genes, upstream = 0, downstream = 1)
tss[1:5]
#?promoters

mat1 = normalizeToMatrix(H3K4me3, tss, value_column = "coverage", 
                         extend = c(1000, 2000), mean_mode = "w0", w = 50,
                         background = 0, smooth = TRUE)
mat1


#EnrichedHeatmap(mat1, col = c("white", "red"), name = "H3K4me3")
#quantile(H3K4me3$coverage, c(0, 0.25, 0.5, 0.75, 0.99, 1))
#quantile(mat1, c(0, 0.25, 0.5, 0.75, 0.99, 1))
# mat1_trim = normalizeToMatrix(H3K4me3, tss, value_column = "coverage", 
#                               extend = 5000, mean_mode = "w0", w = 50, keep = c(0, 0.99))
# EnrichedHeatmap(mat1_trim, col = c("white", "red"), name = "H3K4me3")
# H3K4me3

library(circlize)
col_fun = colorRamp2(quantile(mat1, c(0.01, 0.99)), c("white", "red"))


set.seed(123)
EnrichedHeatmap(mat1, col = col_fun, name = "H3K4me3", 
                use_raster = TRUE, # reduce plot size
                width = 3,
                row_km = 3, 
                #cluster_rows = T,
                top_annotation = HeatmapAnnotation(enriched = anno_enriched(ylim = c(0, 10),
                                                                            axis_param = list(
                                                                              at = c(0, 5, 10),
                                                                              labels = c("zero", "five", "ten"),
                                                                              side = "left",
                                                                              facing = "outside"
                                                                            ),
                                                                            gp = gpar(col = 2:4, lty = 1:3))),
                column_title = "Enrichment of H3K4me3", 
                #column_title_gp,
                row_title_rot = 0
)


## meth
meth[1:3]
mat2 = normalizeToMatrix(meth, tss, value_column = "meth", mean_mode = "absolute",
                         extend = 5000, w = 50, background = NA)
meth_col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
EnrichedHeatmap(mat2, col = meth_col_fun, name = "methylation", column_title = "methylation near TSS")


## target ares regions
mat3 = normalizeToMatrix(meth, cgi, value_column = "meth", mean_mode = "absolute",
                         extend = c(0,5000), 
                         k=20, w = 50, target_ratio = 0.3,
                         background = NA, smooth = TRUE)
EnrichedHeatmap(mat3, col = meth_col_fun, name = "methylation", axis_name_rot = 90,
                column_title = "methylation near CGI")


## multiple heatmap
EnrichedHeatmap(mat1, col = col_fun, name = "H3K4me3",
                top_annotation = HeatmapAnnotation(enrich = anno_enriched(axis_param = list(side = "left")))) + 
  EnrichedHeatmap(mat2, col = meth_col_fun, name = "methylation") +
  ComplexHeatmap::Heatmap(log2(rpkm+1), col = c("white", "orange"), name = "log2(rpkm+1)", 
                          show_row_names = FALSE, width = unit(5, "mm"))
#?kmeans
#tmp <- as.data.frame(mat1)
color3 <- c("steelblue","grey50","firebrick")
col_fun = colorRamp2(quantile(mat1, c(0.01,0.99)), c("white","firebrick"))
meth_col_fun = colorRamp2(c(0, 0.5, 1), c("white","grey50","firebrick"))

partition = paste0("cluster", kmeans(mat1, centers = 3)$cluster)
lgd = Legend(at = c("cluster1", "cluster2", "cluster3"), title = "Clusters", 
             type = "lines", legend_gp = gpar(col = color3))
ht_list = Heatmap(partition, col = structure(color3, names = paste0("cluster", 1:3)), name = "partition",
                  show_row_names = FALSE, width = unit(2, "mm")) +
  EnrichedHeatmap(mat1, col = col_fun, name = "H3K4me3",
                  top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = color3))), 
                  column_title = "H3K4me3") + 
  EnrichedHeatmap(mat2, col = meth_col_fun, name = "methylation",
                  top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = color3))), 
                  column_title = "Methylation") +
  Heatmap(log2(rpkm+1), col = c("white", "orange"), name = "log2(rpkm+1)", 
          show_row_names = FALSE, width = unit(15, "mm"),
          top_annotation = HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = color3),
                                                                    outline = FALSE, axis_param = list(side = "right"))))
draw(ht_list, split = partition, annotation_legend_list = list(lgd), 
     ht_gap = unit(c(2, 8, 8), "mm"))
#?Heatmap


# ht_list = draw(ht_list)
# row_order(ht_list)


## Summarize from a list of matrices
#assume you have a list of histone modification signals for different samples 
#and you want to visualize the mean pattern across samples.
#you can first normalize histone mark signals for each sample and then calculate means values across all samples
mat_list = NULL
for(i in seq_along(hm_gr_list)) {
  mat_list[[i]] = normalizeToMatrix(hm_gr_list[[i]], tss, value_column = ...)
}
#mat_list:
#the first dimension corresponds to genes, 
#the second dimension corresponds to windows 
#the third dimension corresponds to samples

#the mean signal across all samples can be calculated on the third dimension: getSignalsFromList
mat_mean = getSignalsFromList(mat_list)
EnrichedHeatmap(mat_mean)



## pre your own matrix for EnrichedHeatmap
attributes(mat2) = NULL
dim(mat2) = dim(mat1)
as.normalizedMatrix(mat2, 
                    k_upstream = 100, 
                    k_downstream = 100, 
                    k_target = 0,
                    extend = c(5000, 5000), 
                    signal_name = "H3K4me3", 
                    target_name = "TSS"
)


## globally set heatmap params
ht_opt()
EnrichedHeatmap(...)
ht_opt(RESET = TRUE)





# test final EnrichedHeatmap --------------------------------------------------------------------
#BiocManager::install("EnrichedHeatmap")
library(EnrichedHeatmap)

#get data
load(system.file("extdata", "chr21_test_data.RData", package = "EnrichedHeatmap"))

tss = promoters(genes, upstream = 0, downstream = 1)
tss[1:5]
H3K4me3[1:5]

#pre mat
mat1 = normalizeToMatrix(H3K4me3, tss, 
                         value_column = "coverage", 
                         extend = c(1000, 2000), 
                         mean_mode = "w0", 
                         w = 50, # bins of extended windows
                         #k =20, # bins of target regions 
                         background = 0, 
                         smooth = TRUE
)
mat1
mat1.df <- as.data.frame(mat1)

#plot heatmap
library(circlize)
col_fun = colorRamp2(quantile(mat1, c(0.01, 0.99)), c("white", "red"))  # project color to value (0.01,0.99 can prevent outlier)

set.seed(123)
EnrichedHeatmap(mat1, col = col_fun, name = "H3K4me3", 
                #use_raster = TRUE, # reduce plot size
                width = 3, 
                row_km = 3, # order rows option1
                #cluster_rows = T, # order rows option2
                #row_split = sample(c("A", "B"), length(genes), replace = TRUE), # split rows by given partition vector/variables
                
                top_annotation = HeatmapAnnotation(enriched = anno_enriched(ylim = c(0, 10),
                                                                            axis_param = list(
                                                                              at = c(0, 5, 10),
                                                                              labels = c("zero", "five", "ten"),
                                                                              side = "left",
                                                                              facing = "outside"
                                                                            ),
                                                                            gp = gpar(col = 2:4, lty = 1:3))),
                column_title = "Enrichment of H3K4me3", 
                #column_title_gp,
                row_title_rot = 0
)


#plot multi-heatmap
color3 <- c("steelblue","grey50","firebrick")
col_fun = colorRamp2(quantile(mat1, c(0.01,0.99)), c("white","firebrick"))
meth_col_fun = colorRamp2(c(0, 0.5, 1), c("white","grey50","firebrick"))

partition = paste0("cluster", kmeans(mat1, centers = 3)$cluster)
lgd = Legend(at = c("cluster1", "cluster2", "cluster3"), title = "Clusters", 
             type = "lines", legend_gp = gpar(col = color3))
ht_list = Heatmap(partition, col = structure(color3, names = paste0("cluster", 1:3)), name = "partition",
                  show_row_names = FALSE, width = unit(2, "mm")) +
  EnrichedHeatmap(mat1, col = col_fun, name = "H3K4me3",
                  top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = color3))), 
                  column_title = "H3K4me3") + 
  EnrichedHeatmap(mat2, col = meth_col_fun, name = "methylation",
                  top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = color3))), 
                  column_title = "Methylation") +
  Heatmap(log2(rpkm+1), col = c("white", "orange"), name = "log2(rpkm+1)", 
          show_row_names = FALSE, width = unit(15, "mm"),
          top_annotation = HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = color3),
                                                                    outline = FALSE, axis_param = list(side = "right"))))
draw(ht_list, split = partition, annotation_legend_list = list(lgd), 
     ht_gap = unit(c(2, 8, 8), "mm"))



##Summarize from a list of matrices
#assume you have a list of histone modification signals for different samples 
#and you want to visualize the mean pattern across samples.
#you can first normalize histone mark signals for each sample and then calculate means values across all samples
mat_list = NULL
for(i in seq_along(hm_gr_list)) {
  mat_list[[i]] = normalizeToMatrix(hm_gr_list[[i]], tss, value_column = ...)
}
#mat_list:
#the first dimension corresponds to genes, 
#the second dimension corresponds to windows 
#the third dimension corresponds to samples

#the mean signal across all samples can be calculated on the third dimension: getSignalsFromList
mat_mean = getSignalsFromList(mat_list)
EnrichedHeatmap(mat_mean)



##pre your own matrix for EnrichedHeatmap
attributes(mat2) = NULL
dim(mat2) = dim(mat1)
as.normalizedMatrix(mat2, 
                    k_upstream = 100, 
                    k_downstream = 100, 
                    k_target = 0,
                    extend = c(5000, 5000), 
                    signal_name = "H3K4me3", 
                    target_name = "TSS"
)














# test deeptools --------------------------------------------------------------------
# setwd("/Users/baopengfei/Desktop/lulab/tmp/")
setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC")
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggsci)
library(Rmisc)
library(data.table)
#install.packages("Rmisc")
conflict_prefer("select", "dplyr")

full.len <- 150 # nt
smps <- c("FTA-10") # samples for signals

refp <- data.table::fread('./ftc_small_all_50k.mat.gz',header = F)
#refp <- fread('./ftc_small_AGO2_gn_point.mat',header = F)
head(refp[1:3,1:8])
dim(refp)
# [1] 181   606
#上下游各100b，加target 20b，加和除以10 binsize,再乘以2个样本,再乘以1个区域,加上前6列基因信息为44列
# (100 + 100 + 20)/10*2*1
#sample1.pos1,sample1.posN...; sampleN.posN,...,sampleN.posN

clean_refp <- refp %>% select(c(-1:-3,-5,-6))
#head(clean_refp[1:3,1:8])
long_refp <- melt(clean_refp,id.vars = 'V4',value.name = 'signal') %>%
  select(-variable)
#add sample name
# long_refp$sample <- rep(c("FTA-10"), # ,"csFTA-10"
#                         each = nrow(refp)*22) # (100 + 100 + 20)/10
long_refp$sample <- rep(smps, each = nrow(refp)*full.len*0.1) # (100 + 100 + 20)/10

# add x position
long_refp$pos <- rep(c(1:(full.len*0.1)),each = nrow(refp),times = length(smps))

# calculate means with CI, plot line
filnal_scaler <- long_refp %>%
  # remove na
  tidyr::drop_na() %>%
  dplyr::group_by(sample,pos) %>%
  # mean and 95% interval confidence
  dplyr::summarise(mean_signal = mean(signal),
                   sd = sd(signal),
                   upper = CI(signal,ci = 0.95)[1],
                   lower = CI(signal,ci = 0.95)[3]
  )


# plot
p <- ggplot(filnal_scaler,aes(x = pos,y = mean_signal)) +
  geom_line(aes(color = sample),size = 1) +
  theme_classic(base_size = 14) +
  scale_color_d3(name = '') +
  # x label
  scale_x_continuous(breaks = c(0,10,12,22),
                     labels = c('-100b','Start','End','+100b')) +
  xlab('') + ylab('Normalized signal') +
  theme(aspect.ratio = 0.8)

p
# facet
p + facet_wrap(~sample,scales = 'free_x',ncol = 2) +
  theme(strip.background = element_rect(color = NA,fill = 'grey'),
        axis.text.x = element_text(angle = 45,hjust = 1))

# add CI
ggplot(filnal_refp,aes(x = pos,y = mean_signal)) +
  # add 0.95 interval
  geom_ribbon(aes(ymin = lower,
                  ymax = upper,
                  fill = sample),
              alpha = 0.5) +
  geom_line(aes(color = sample),size = 1) +
  theme_classic(base_size = 16) +
  scale_color_d3(name = '') +
  scale_fill_d3(name = '') +
  # x label
  scale_x_continuous(breaks = c(0,10,12,22),
                     labels = c('-100b','Start','End','+100b')) +
  xlab('') + ylab('Normalized signal') +
  theme(aspect.ratio = 0.8)

# heatmap
# plot (slowly)
long_refp$log10.signal <- log10(long_refp$signal+0.00001)
long_refp_sort <- long_refp %>% 
  # dplyr::filter(pos %in% 10:12) %>% # only rank/sort by central regions
  dplyr::group_by(sample,V4) %>%
  #mean of each row/region
  dplyr::summarise(mean_signal_region = mean(signal,na.rm=T)
  )
long_refp_sort <- long_refp_sort[order(long_refp_sort$mean_signal_region,decreasing = F),]
long_refp$V4 <- factor(long_refp$V4,levels = unique(long_refp_sort$V4))
ggplot(long_refp, # %>% filter(sample == 'H3K27ac')
       aes(x = pos,y = V4)) + 
  geom_tile(aes(fill = log10.signal)) +
  theme_bw() +
  # coord_cartesian(expand = 0) +
  scale_x_continuous(breaks = c(0,10,12,22),
                     labels = c('-100b','Start','End','+100b')) +
  scale_fill_gradient(low = 'white',high = 'red') +
  ylab('') + xlab('') +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        aspect.ratio = 1.5)



# test enrichedHeatmap + R ------------------------------------------------
#plot gn Heatmap (long reads & peaks, on CDS, UTR)

## pre-run in bash to get peak and reads bed
#only consider long RNA-seq, small RNA not reasonable
#get reads bed
samtools sort -@4 -n  /Users/baopengfei/Desktop/lulab/tmp/elFTA-10_L4.gn.bam >  /Users/baopengfei/Desktop/lulab/tmp/elFTA-10_L4.gn.sort.bam; 
samtools view -bf 0x2 -q 0  /Users/baopengfei/Desktop/lulab/tmp/elFTA-10_L4.gn.sort.bam \
| bedtools bamtobed -i stdin -bedpe -mate1 \
| awk -v s="reverse" 'BEGIN{{FS=OFS="\t"}}{{if($9=="+"){{min=$2;max=$6}} else if($9=="-"){{min=$5;max=$3}} \
            {{if(s=="reverse") print($1,min,max,$7,$8,$10); else if(s=="forward") print($1,min,max,$7,$8,$9)}} }}' \
| pigz -c -p 4 >  /Users/baopengfei/Desktop/lulab/tmp/elFTA-10_L4.gn.bed.gz

samtools sort -@4 -n  /Users/baopengfei/Desktop/lulab/tmp/clFTA-10.gn.bam >  /Users/baopengfei/Desktop/lulab/tmp/clFTA-10.gn.sort.bam; 
samtools view -bf 0x2 -q 0  /Users/baopengfei/Desktop/lulab/tmp/clFTA-10.gn.sort.bam \
| bedtools bamtobed -i stdin -bedpe -mate1 \
| awk -v s="reverse" 'BEGIN{{FS=OFS="\t"}}{{if($9=="+"){{min=$2;max=$6}} else if($9=="-"){{min=$5;max=$3}} \
            {{if(s=="reverse") print($1,min,max,$7,$8,$10); else if(s=="forward") print($1,min,max,$7,$8,$9)}} }}' \
| pigz -c -p 4 >  /Users/baopengfei/Desktop/lulab/tmp/clFTA-10.gn.bed.gz

#get peak bed
cd /BioII/lulab_b/baopengfei/
  bed="projects/WCHSU-FTC/output/FTC_long/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/domains.bed"
{{ grep -v '^chr' ${bed}  | bash projects/WCHSU-FTC/exSeek-dev/bin/tbed2gbed <(cat projects/WCHSU-FTC/exSeek-dev/genome/hg38/bed/{long_RNA,tRNA}.bed) /dev/stdin /dev/stdout        awk 'BEGIN{{OFS="\t";FS="\t"}}/^chr/{{print $1,$2,$3,$4,$5,$6,0,0,0,1,$3-$2,0}}'  ${bed}  }} | bedtools sort >  projects/WCHSU-FTC/output/FTC_long/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/domains_gn.bed



## get matrix by EnrichedHeatmap in R
#library(EnrichedHeatmap)
# setwd("/Users/baopengfei/Desktop/lulab/tmp/")
setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC")
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggsci)
library(Rmisc)
library(data.table)
conflict_prefer("select", "dplyr")
# signal <- rtracklayer::import.bed(con="elFTA-10_L4.gn.bed.gz") # elFTA-10.gn.bed.gz
# signal2 <- rtracklayer::import.bed(con="clFTA-10.gn.bed.gz") # clFTA-10.gn.bed.gz
# signal3 <- rtracklayer::import.bed(con="projects/WCHSU-FTC/output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/domains_gn.bed") # small domain
# signal4 <- rtracklayer::import.bed(con="projects/WCHSU-FTC/output/FTC_long/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/domains.bed") # small domain
# 
# region <- rtracklayer::import.bed(as.character("shared_reference/RBP/tmp/CDS_8k.bed"))
# region2 <- rtracklayer::import.bed(as.character("shared_reference/RBP/tmp/UTR3_8k.bed"))
# region2 <- rtracklayer::import.bed(as.character("shared_reference/RBP/tmp/start_codon_8k.bed"))


full.len <- 150 # nt
smps <- c("FTA-10") # samples for signals

get.mat <- function(signal,region){
  print(signal)
  print(region)
  signal.gr <- rtracklayer::import.bed(con=as.character(signal)) # elFTA-10.gn.bed.gz
  region.gr <- rtracklayer::import.bed(con=as.character(region))
  #strand info not preserved, not flipped
  mat1 = EnrichedHeatmap::normalizeToMatrix(signal.gr, region.gr, 
                                            #value_column = "cov",  
                                            extend = c(250, 250),
                                            target_ratio = 0.3, # k = 10, # bins of target regions
                                            mean_mode = "coverage",  # c("absolute", "weighted", "w0", "coverage")
                                            w = 10, # bins of extended windows
                                            keep = c(0, 0.99), # Percentiles in the normalized matrix to keep.
                                            background = 0, 
                                            smooth = TRUE # set TRUE may get negative cov 
  )
  mat1 <- as.data.frame(mat1)
  # will normalize whole matrix !
  signal <- basename(as.character(signal))
  signal <- unlist(sapply(strsplit(signal,".",fixed = T),"[",1))
  region <- basename(as.character(region))
  region <- unlist(sapply(strsplit(region,".",fixed = T),"[",1))
  
  res.tmp <- as.data.frame(cbind(bed=region.gr$name,region=rep(region,nrow(mat1)),signal=rep(signal,nrow(mat1)),mat1))
  print(nrow(res.tmp))
  print(ncol(res.tmp))
  return(res.tmp)
}
#"elFTA-10_L4.gn.bed.gz","clFTA-10.gn.bed.gz","projects/WCHSU-FTC/output/FTC_long/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/domains_gn.bed"
cmb <- expand.grid(signal=c("elFTA-10_L4.gn.bed.gz","clFTA-10.gn.bed.gz"),
                   region=c("shared_reference/RBP/tmp/CDS_8k.bed","shared_reference/RBP/tmp/UTR3_8k.bed"))
#one way using sapply
res.list <- lapply(1:nrow(cmb), function(j) {get.mat(signal=as.character(cmb$signal[j]), region=as.character(cmb$region[j]))} )
res.list.df <- do.call(rbind,res.list)
#mat1.df <- as.data.frame(mat1)
#dim(mat1.df)
#mat1.df[1:3,1:4]
#hist(as.matrix(mat1.df))
#table(as.matrix(mat1.df)<0)
#head(res.list[[2]],3)
res.list.df[1:3,1:5]


# plot in R
#(500+500)/20+10=60,再乘以1个样本,再乘以1个区域,为60列
clean_refp <- dplyr::as_tibble(res.list.df)
#clean_refp <- refp #%>% select(c(-1:-3,-5,-6))
colnames(clean_refp)[c(1,3)] <- c('V4','sample')

long_refp <- reshape2::melt(clean_refp,id.vars = c('V4','sample','region'),value.name = 'signal') %>%
  select(-variable)
#add sample name
#long_refp$sample <- rep(c("FTA-10","csFTA-10"),
#                        each = nrow(refp)*22) # (100 + 100 + 20)/10

# add x position
long_refp$pos <- rep(c(1:(ncol(res.list.df)-3)),each = nrow(res.list.df)) # ,times = 2

# calculate means with CI, plot line
filnal_scaler <- long_refp %>%
  # remove na
  tidyr::drop_na() %>%
  dplyr::group_by(pos,sample,region) %>% # sample,
  # mean and 95% interval confidence
  dplyr::summarise(mean_signal = mean(signal), #,trim = 0.05, na.rm=T
                   sd = sd(signal),
                   upper = CI(signal,ci = 0.95)[1],
                   lower = CI(signal,ci = 0.95)[3]
  )


# plot
# p <- ggplot(filnal_scaler,aes(x = pos,y = mean_signal)) +
#   geom_line(size = 1) + #aes(color = sample),
#   theme_classic(base_size = 14) +
#   scale_color_d3(name = '') +
#   # x label
#   scale_x_continuous(breaks = c(0,25,35,ncol(res.list.df)),
#                      labels = c('-500b','Start','End','+500b')) +
#   xlab('') + ylab('Normalized signal') +
#   theme(aspect.ratio = 0.8)
#p
# add CI
ggplot(filnal_scaler,aes(x = pos,y = mean_signal)) +
  # add 0.95 interval
  geom_ribbon(aes(ymin = lower,
                  ymax = upper,
                  fill = sample), # sample
              alpha = 0.5) +
  geom_line(aes(color = sample), size = 1) + # 
  theme_classic(base_size = 16) +
  scale_color_d3(name = 'Peak OR RawReads') +
  scale_fill_d3(name = '') +
  geom_vline(xintercept = c(25,35),color="grey50",linetype="dashed")+
  # x label
  scale_x_continuous(breaks = c(0,25,35,ncol(res.list.df)-1),
                     labels = c('-250nt','Start','End','+250nt')) +
  xlab('') + ylab('Normalized signal') +
  theme(aspect.ratio = 0.4) +  #
  facet_wrap(region~.,scales = 'free_y',nrow = 3) + #
  theme(strip.background = element_rect(color = NA,fill = 'grey'),
        axis.text.x = element_text(angle = 0,hjust = 0.5))



# heatmap (need further optimization)
# plot (slowly)
# long_refp$log10.signal <- log10(long_refp$signal+0.00001)
# long_refp_sort <- long_refp %>% 
#   dplyr::filter(pos %in% 25:35) %>% # only rank/sort by central regions
#   dplyr::group_by(V4) %>%  # sample,
#   #mean of each row/region
#   dplyr::summarise(mean_signal_region = mean(signal,na.rm=T)
#   )
# long_refp_sort <- long_refp_sort[order(long_refp_sort$mean_signal_region,decreasing = F),]
# long_refp2 <- long_refp
# long_refp2$V4 <- factor(long_refp2$V4,levels = unique(long_refp_sort$V4))
# 
# ggplot(long_refp2, # %>% filter(sample == 'H3K27ac')
#        aes(x = pos,y = V4)) + 
#   geom_tile(aes(fill = signal)) +
#   theme_bw() +
#   # coord_cartesian(expand = 0) +
#   scale_x_continuous(breaks = c(0,10,12,22),
#                      labels = c('-100b','Start','End','+100b')) +
#   scale_fill_gradient(low = 'white',high = 'red') +
#   ylab('') + xlab('') +
#   theme(axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         panel.grid = element_blank(),
#         aspect.ratio = 1.5)









# startcodon
get.mat.TSS <- function(signal,region){
  print(signal)
  print(region)
  signal.gr <- rtracklayer::import.bed(con=as.character(signal)) # elFTA-10.gn.bed.gz
  region.gr <- rtracklayer::import.bed(con=as.character(region))
  mat1 = EnrichedHeatmap::normalizeToMatrix(signal.gr, region.gr, 
                                            #value_column = "cov",  
                                            extend = c(250, 250),
                                            target_ratio = 0, # k = 10, # bins of target regions
                                            mean_mode = "coverage",  # c("absolute", "weighted", "w0", "coverage")
                                            w = 10, # bins of extended windows
                                            keep = c(0, 0.99), # Percentiles in the normalized matrix to keep.
                                            background = 0, 
                                            smooth = TRUE # set TRUE may get negative cov 
  )
  mat1 <- as.data.frame(mat1)
  # will normalize whole matrix !
  signal <- basename(as.character(signal))
  signal <- unlist(sapply(strsplit(signal,".",fixed = T),"[",1))
  region <- basename(as.character(region))
  region <- unlist(sapply(strsplit(region,".",fixed = T),"[",1))
  
  res.tmp <- as.data.frame(cbind(bed=region.gr$name,region=rep(region,nrow(mat1)),signal=rep(signal,nrow(mat1)),mat1))
  print(nrow(res.tmp))
  print(ncol(res.tmp))
  return(res.tmp)
} 
cmb <- expand.grid(signal=c("elFTA-10_L4.gn.bed.gz","clFTA-10.gn.bed.gz"), # ,"projects/WCHSU-FTC/output/FTC_long/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/domains_gn.bed"
                   region=c("shared_reference/RBP/tmp/start_codon_8k.bed"))
#one way using sapply
res.list2 <- lapply(1:nrow(cmb), function(j) {get.mat.TSS(signal=as.character(cmb$signal[j]), region=as.character(cmb$region[j]))} )
res.list2.df <- do.call(rbind,res.list2)
#mat1.df <- as.data.frame(mat1)
#dim(mat1.df)
#mat1.df[1:3,1:4]
#hist(as.matrix(mat1.df))
#table(as.matrix(mat1.df)<0)
#head(res.list[[2]],3)
res.list2.df[1:3,1:5]


# plot in R
library(dplyr)
#(500+500)/20+10=60,再乘以1个样本,再乘以1个区域,为60列
clean_refp <- dplyr::as_tibble(res.list2.df)
#clean_refp <- refp #%>% select(c(-1:-3,-5,-6))
colnames(clean_refp)[c(1,3)] <- c('V4','sample')

long_refp <- melt(clean_refp,id.vars = c('V4','sample','region'),value.name = 'signal') %>%
  select(-variable)
#add sample name
#long_refp$sample <- rep(c("FTA-10","csFTA-10"),
#                        each = nrow(refp)*22) # (100 + 100 + 20)/10

# add x position
long_refp$pos <- rep(c(1:(ncol(res.list2.df)-3)),each = nrow(res.list2.df)) # ,times = 2

# calculate means with CI, plot line
filnal_scaler <- long_refp %>%
  # remove na
  tidyr::drop_na() %>%
  dplyr::group_by(pos,sample,region) %>% # sample,
  # mean and 95% interval confidence
  dplyr::summarise(mean_signal = mean(signal), #,trim = 0.05, na.rm=T
                   sd = sd(signal),
                   upper = CI(signal,ci = 0.95)[1],
                   lower = CI(signal,ci = 0.95)[3]
  )

# plot
# p <- ggplot(filnal_scaler,aes(x = pos,y = mean_signal)) +
#   geom_line(size = 1) + #aes(color = sample),
#   theme_classic(base_size = 14) +
#   scale_color_d3(name = '') +
#   # x label
#   scale_x_continuous(breaks = c(0,25,35,ncol(res.list2.df)),
#                      labels = c('-500b','Start','End','+500b')) +
#   xlab('') + ylab('Normalized signal') +
#   theme(aspect.ratio = 0.8)
#p

# add CI
ggplot(filnal_scaler,aes(x = pos,y = mean_signal)) +
  # add 0.95 interval
  geom_ribbon(aes(ymin = lower,
                  ymax = upper,
                  fill = sample), # sample
              alpha = 0.5) +
  geom_line(aes(color = sample), size = 1) + # 
  theme_classic(base_size = 16) +
  scale_color_d3(name = 'Peak OR RawReads') +
  scale_fill_d3(name = '') +
  # x label
  geom_vline(xintercept = c(25),color="grey50",linetype="dashed")+
  scale_x_continuous(breaks = c(0,25,ncol(res.list2.df)-1),
                     labels = c('-250nt','Start_Codon','+250nt')) +
  xlab('') + ylab('Normalized signal') +
  theme(aspect.ratio = 0.4) + 
  facet_wrap(region~.,scales = 'free_y',nrow = 3) + #
  theme(strip.background = element_rect(color = NA,fill = 'grey'),
        axis.text.x = element_text(angle = 0,hjust = .5))



