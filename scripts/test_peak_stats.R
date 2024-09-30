# test localmax

# last 2205 by bpf 
# b.p.f@qq.com
# transmit form local to bioii at 220604

setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/")
#setwd("/Users/baopengfei/Desktop/lulab/tmp/projects/WCHSU-FTC/exSeek-dev")
options(stringsAsFactors = F)
library(ggplot2)
library(conflicted)

# load all func.
source("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/scripts/util.R")
#


# sum alls sample single-base cpm & stat -----------------------------------------------------------------------
library(tidyverse)
library(dplyr)
library(stringr)
library(parallel)
# merge,sort,index bam (each sample)
#samples <- read.table("data/FTC_long/sample_ids_all.txt")$V1
samples <- read.table("data/FTC_small/sample_ids_all.txt")$V1
samples <- samples[sample(1:length(samples),4,replace = F)]
#bedfile <- "../output/FTC_long/call_domain_withRepeats_dedup/domains_localmax/b5_d05_p01.bed"
bedfile <- "../output/FTC_small/call_domain_withRepeats_dedup/domains_localmax/b5_d05_p01.bed"
cores <- 4

# # merge, sort, index bams
# getBam <- function(sample){
#   #sample <- samples[1]
#   print(sample)
#   bamMerge <- paste0("../output/FTC_long/tbam_bowtie2Uniq/",sample,"/bam-deduped/merge16.bam")
#   bamSort <- paste0("../output/FTC_long/tbam_bowtie2Uniq/",sample,"/bam-deduped/merge16_sort") # For asBam asSam, and sortBam this is without the “.bam” file suffix.
#   fs <- paste0("../output/FTC_long/tbam_bowtie2Uniq/",sample,"/bam-deduped/",c(rna,dna),".bam")
#   Rsamtools::mergeBam(files = fs, destination = bamMerge, overwrite=T)
#   Rsamtools::sortBam(file = bamMerge, destination = bamSort, overwrite=T)
#   Rsamtools::indexBam(file = paste0(bamSort,".bam"), overwrite=T)
#   file.remove(list = bamMerge)
# }
# library(parallel)
# mclapply(samples, getBam, mc.cores = cores)

# # bedtools multicov
# bams <- Sys.glob("../output/FTC_long/tbam_bowtie2/*/bam-deduped/merge16_sort.bam")
# #multiCov <- bedtoolsr::bt.multicov(bams = bams, split = T,bed = tmp.tbl, s = T) # too slow in R
# bedtools multicov -s \
#   -bams ../output/FTC_long/tbam_bowtie2/*FTA-17*/bam-deduped/merge16_sort.bam \
#   -bed ../output/FTC_long/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/domains_singleBase.bed \
#   > ../output/FTC_long/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/count_matrix/domains_long_singleBase.txt
# #samtools(bedtools) error when not enough mem: 
# #samtools sort: [E::bgzf_read] Read block operation failed with error 2 after 0 of 4 bytes

# bedtools coverage (seem faster than multicov)
#bedtools coverage -s -split -d \
# -b ../output/FTC_long/tbam_bowtie2/clFTA-17/bam-deduped/merge16_sort.bam \
# -a ../output/FTC_long/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/domains_singleBase.bed | head
getCov <- function(sample){
  #sample <- samples[1]
  print(sample)
  bedtoolsr::bt.coverage(a = bedfile, 
                         # b = paste0("../output/FTC_long/tbam/",sample,"/bam-deduped/merge19_sort.bam"), d = T, s = T, split = T, 
                         b = paste0("../output/FTC_small/tbam/",sample,"/bam-deduped/merge19_sort.bam"), d = T, s = T, split = T,
                         output = paste0("tmp/",sample,"_singleBase.txt"))
}
library(parallel)
mclapply(samples, getCov, mc.cores = cores)

dat <- list()
read.Bed <- function(sample){
  #sample <- samples[1]
  print(sample)
  tmp <- data.table::fread(paste0("tmp/",sample,"_singleBase.txt"),data.table = F,sep = "\t",stringsAsFactors = F,header = F)
  tmp <- dplyr::as_tibble(tmp)
  #rownames(tmp) <- 
  colnames(tmp)[8] <- sample
  # tmp <- tmp %>% 
  #   select() %>% 
  #   rename(sample=s)
  #   mutate(sample_id=sample)
  dat[[sample]] <- tmp
  # return(tmp) 
}
library(parallel)
dat <- mclapply(samples,read.Bed,mc.cores = cores)
#dat <- mclapply(samples,read.Bed)
dat.tbl <- do.call(cbind,dat)
#table(duplicated(colnames(dat.tbl)))
dat.tbl <- dplyr::as_tibble(dat.tbl[,!duplicated(colnames(dat.tbl))])
# dat.tbl[,8:ncol(dat.tbl)] <- dat.tbl[,8:ncol(dat.tbl)]/s
cpm <- edgeR::cpm(y = dat.tbl[,8:ncol(dat.tbl)], prior.count = 0, normalized.lib.sizes = F, log = FALSE)
#cpm[1:30,]

#stats among samples
#cpm like deeptools single base cpm
stats <- dat.tbl[,1:7]
stats$sd.cpm <- apply(cpm,1,sd)
stats$median.cpm <- apply(cpm,1,median)
stats$mean.cpm <- apply(cpm,1,mean)
stats$max.cpm <- apply(cpm,1,max)
stats$min.cpm <- apply(cpm,1,min)
stats$cv.cpm <- stats$sd/stats$mean
stats$freq.cpm <- rowSums(cpm>0)/ncol(cpm)
#hist(stats$freq)
#hist(stats$cv)
stats$V3 <- stats$V2+stats$V7
colnames(stats)[3:4] <- c('end','name')
stats[1:3,]


# convert peak bed to single base bed (only needed for bedtools multicov, not needed for bedtools coverage)
mat <- data.table::fread("./genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt",sep = "\t",stringsAsFactors = F,header = T)
mat[1:3,]
bed <- data.table::fread(bedfile,header = F,sep = "\t",stringsAsFactors = F,data.table = F)
bed <- as_tibble(bed)
colnames(bed) <- c("seqnames","start","end","name","score","strand")
bed$width <- bed$end-bed$start
table( as.character(bed$seqnames) %in% as.character(mat$transcript_id) ) # should be all ture !
bed$tx.type <- mat$transcript_type[match(as.character(bed$seqnames),as.character(mat$transcript_id))]
#bed$start <- bed$start-1
library(parallel)
f1 <- function(i){
  #i <- 1
  # if(i%%500==0){
  #   message(paste0(i,"/",nrow(bed)))
  # }
  numb <- bed$width[i]
  tmp <- tibble(seqnames=rep(bed$seqnames[i],numb),
                start=bed$start[i]:(bed$start[i]+numb-1),
                end=(bed$start[i]+1):(bed$start[i]+numb),
                name=rep(bed$name[i],numb),
                score=rep(bed$score[i],numb),
                strand=rep(bed$strand[i],numb)
  )
  return(tmp)
}
#tmp.list <- lapply(1:nrow(bed),f1)
tmp.list <- mclapply(1:nrow(bed),f1, mc.cores = cores)
tmp.tbl <- as_tibble(do.call("rbind",tmp.list))
sum(bed$width) == nrow(tmp.tbl) # should be TRUE
#data.table::fwrite(tmp.tbl,"../output/FTC_long/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/domains_singleBase.bed",sep = "\t",quote = F,col.names = F,row.names = F)
tmp.tbl$width <- tmp.tbl$end-tmp.tbl$start

all(tmp.tbl$name==stats$name) # all TRUE
all(tmp.tbl$seqnames==stats$V1) # all TRUE

tmp.tbl <- left_join(tmp.tbl,stats)
#tmp.tbl
#   seqnames start   end name   score strand width V1       V2    V5 V6       V7     sd median   mean    cv  freq    max   min
# 1 12985       18    19 peak_1  31.3 +          1 12985    18  31.3 +         1 0.0670 0.0473 0.0473  1.41   0.5 0.0947     0
# 2 12985       19    20 peak_1  31.3 +          1 12985    18  31.3 +         2 0.0670 0.0473 0.0473  1.41   0.5 0.0947     0
# 3 12985       20    21 peak_1  31.3 +          1 12985    18  31.3 +         3 0.0670 0.0473 0.0473  1.41   0.5 0.0947     0

bg.tbl <- tmp.tbl %>% 
  group_by(name) %>% 
  mutate(peak.sum.mean.cpm=sum(width*mean.cpm), peak.min.mean.cpm=min(mean.cpm), peak.max.mean.cpm=max(mean.cpm), peak.mean.mean.cpm=mean(rep(mean.cpm,width)),  peak.median.mean.cov=median(rep(mean.cpm,width)), peak.sd.mean.cpm=sd(rep(mean.cpm,width)), peak.cv.mean.cpm=peak.sd.mean.cpm/peak.mean.mean.cpm ) # mean.cov=total.cov/sum(width), 
bg.tbl$tx.type <- bed$tx.type[match(as.character(bg.tbl$name),as.character(bed$name))]
#table(bg.tbl$tx.type)
#bg.tbl


data.table::fwrite(bg.tbl,paste0(bedfile,".singleBaseCPMstat.txt"),quote = F,sep = '\t',row.names = F,col.names = T)
rownames(cpm) <- paste0(bg.tbl$name,".",bg.tbl$start,".",bg.tbl$end,".",bg.tbl$strand)
data.table::fwrite(cbind(peakID=rownames(cpm),cpm),paste0(bedfile,".singleBaseCPMvalue.txt"),quote = F,sep = '\t',row.names = F,col.names = T)





# plot shape: broad & sharp  ----------------------------------------------
samples <- read.table("data/FTC_long/sample_ids_all.txt")$V1
samples <- samples[sample(1:length(samples),10,replace = F)]
bedfile <- "../output/FTC_long/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/domains.bed"

bg.tbl <- data.table::fread(paste0(bedfile,".singleBaseCPMstat.txt"),sep = '\t',header = T,stringsAsFactors = F,check.names = F)
bg.tbl.peak <- bg.tbl[!duplicated(bg.tbl$name),]
hist(bg.tbl.peak$cv,breaks = 50,col = "grey50")
table(bg.tbl.peak$tx.type)
table(bg.tbl.peak$cv<=mean(bg.tbl.peak$cv))

table(bg.tbl.peak$cv<=0.2)
dat <- as.data.frame(table(bg.tbl.peak$tx.type,bg.tbl.peak$cv>=0.2)) # mean(bg.tbl.peak$cv)
colnames(dat) <- c("RNA","Shape","Number")
dat$Shape <- ifelse(dat$Shape==T,"sharp","broad")
dat.w <- dat %>% 
  tidyr::pivot_wider(names_from = c("Shape"), values_from = "Number") %>% 
  mutate(Ratio=sharp/broad)
dat <- dat.w %>% 
  tidyr::pivot_longer(cols = 2:ncol(dat.w), names_to = "Shape", values_to = "Value") %>% 
  filter(Shape=="Ratio")
dat$RNA <- gsub(".rev","",dat$RNA,fixed = T)
dat$RNA <- gsub(".for","",dat$RNA,fixed = T)
rna <- c("pri_miRNA","lncRNA", "mRNA", "piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA')
#dna <- c("intron.for", "intron.rev", "promoter.for", "promoter.rev", "enhancer.for", "enhancer.rev", "repeats.for", "repeats.rev")
dna <- c("intron","promoter", "enhancer") # repeats
dat$RNA <- factor(dat$RNA,levels = c(rna,dna))
#dat <- as_tibble(dat)
dat <- dat %>% 
  group_by(RNA) %>% 
  summarise(Value=mean(Value))

dat.w$RNA <- gsub(".rev","",dat.w$RNA,fixed = T)
dat.w$RNA <- gsub(".for","",dat.w$RNA,fixed = T)
dat.w <- dat.w %>% 
  select(RNA,broad,sharp) %>% 
  group_by(RNA) %>% 
  summarise(broad=sum(broad),sharp=sum(sharp))

dat <- left_join(dat,dat.w)
dat$lab <- paste0(dat$sharp,"/",dat$broad)
dat$lab.pos <- dat$Value*1.0+0.1
dat$RNA <- factor(dat$RNA,levels = c(rna,dna))
library(ggplot2)
p <- ggplot(dat,aes(x=RNA,y=Value))+ 
  geom_bar(position = "stack",width = 0.8,size=0.3, stat = "sum",color="black")+
  geom_text(aes(x=RNA,y=lab.pos,label=lab),color="black",size=5) +
  xlab("")+
  ylab(paste0("Sharp/broad Ratio"))+ # ," (mean CV:",round(0.2,digits = 2),")"
  #annotate(paste0(" (mean CV: ",round(mean(bg.tbl.peak$cv),digits = 2),")"),x = 0.6,y = 0.8 ) + #, position = c(0.6,0.8)
  geom_vline(xintercept = 0,linetype="dashed",color="grey10")+
  geom_hline(yintercept = 1,linetype="dashed",color="red")+
  #ggsci::scale_fill_d3()+
  scale_fill_manual(name="Species",values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set3"))(length(c(rna,dna))))+ # c(RColorBrewer::brewer.pal(8, "Accent"),RColorBrewer::brewer.pal(9, "Set1"))
  #paletteer::scale_fill_paletteer_d("RColorBrewer::Set1")+
  #scale_y_continuous(trans = "log10")+
  #xlim(c(0,500))+
  #geom_hline(yintercept = c(0))+
  #ggraph::scale_fill_viridis(option = "C") + # name = "value",
  theme_bw() + 
  theme(plot.title = element_text(size = 24,color="black",hjust = 0.5,face="bold"),
        axis.title = element_text(size = 24,color ="black"), 
        axis.text = element_text(size= 20,color = "black"),
        #panel.grid=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey50",linetype = "dashed"), #size= 1,
        #panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5), #  , color = c(rep("#003366",length(rna)),rep("darkred",length(dna)))
        legend.position = "right",#c(.25,.6),
        legend.text = element_text(size= 16),
        legend.title= element_text(size= 16))
p
#ggsave("./output/lulab/FTC_small/localmax_Assign.pdf",width = 7,height = 10)
#ggsave("./output/lulab/FTC_small/localmax_Assign_withPri_miR.pdf",width = 7,height = 10)
ggsave(plot = p, filename = "./localmax_shape.pdf",width = 10,height = 6)




# sum tx reference number/length ----------------------------------------------------
fs <- Sys.glob("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/tbed/*.bed")
# fs <- Sys.glob("/Users/baopengfei/Desktop/lulab/tmp/projects/WCHSU-FTC/exSeek-dev/genome/hg38/tbed/*.bed")
#fs <- fs[-1] # rm all.bed
fs <- fs[!grepl("all|10RNA|11RNA",perl = T,fs)] # 18,19 left
bed.list <- list()
read.txbed <- function(x){
  #x <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/tbed/mRNA.bed"
  print(x)
  type <- gsub(".bed","",basename(x))
  print(type)
  tmp.bed <- data.table::fread(x,sep = "\t",header=F,stringsAsFactors = F,data.table = F)
  #tmp.bed[1:3,]
  tmp.bed$wid <- tmp.bed$V3-tmp.bed$V2
  tmp.bed <- tmp.bed[order(as.numeric(tmp.bed$wid),decreasing = F),]
  #tmp.bed$rank <- 1:nrow(tmp.bed)
  tmp.bed$type <- type
  tmp.bed$total.length <- sum(tmp.bed$wid)
  #tmp.bed$wid.ratio <- tmp.bed$wid/tmp.bed$total.length
  
  tmp.bed$cumsum.wid.ratio <- cumsum(tmp.bed$wid/tmp.bed$total.length)
  tmp.bed$cumsum.rank.ratio <- (1:nrow(tmp.bed))/nrow(tmp.bed)
  
  bed.list[[type]] <- tmp.bed
  #return(bed.list)
}
bed.list <- lapply(fs,read.txbed)
bed.list.df <- do.call("rbind", bed.list)
bed.list.df <- as.data.frame(bed.list.df)
rnas <- c("rRNA","pri_miRNA","lncRNA","mRNA","piRNA","snoRNA","snRNA","srpRNA","tRNA","tucpRNA","Y_RNA","enhancer.for","enhancer.rev","intron.for","intron.rev","promoter.for","promoter.rev","repeats.for","repeats.rev")
# rnas <- c("rRNA","pri_miRNA","lncRNA","mRNA","piRNA","snoRNA","snRNA","srpRNA","tRNA","tucpRNA","Y_RNA","enhancer","intron","promoter","repeats")
#bed.list.df[1:3,]
# bed.list.df$type <- gsub(".rev|.for","",bed.list.df$type,perl = T)
table(bed.list.df$type)

## 1.plot tx num pie
tx.num <- as.data.frame(table(bed.list.df$type))
#tx.num <- tx.num[!grepl("repeats",tx.num$Var1),]
str(tx.num)
tx.num$Freq <- tx.num$Freq/sum(tx.num$Freq)
tx.num$Freq <- round(tx.num$Freq, digits = 4) #format(tx.num$Freq,scientific = 4)

library(ggplot2)
library(dplyr)
#str(df)
df <- tx.num
colnames(df) <- c("group","value")
df$value <- as.numeric(df$value)
df$group <- factor(df$group,levels = rnas)
df <- df[order(df$group,decreasing = F),]
#df <- df[!grepl("repeats",as.character(df$group)),]
# df <- df[!grepl("for|rev",perl = T,as.character(df$group)),]
df2 <- df %>% 
  mutate(csum = rev(cumsum(rev(value))), 
         pos = value/2 + lead(csum, 1),
         pos = if_else(is.na(pos), value/2, pos))
ggplot(df, aes(x = "" , y = value, fill = factor(group))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(name="Group",values = colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(19) ) + # c(RColorBrewer::brewer.pal(8, "Accent"),RColorBrewer::brewer.pal(9, "Set1"))
  ggrepel::geom_label_repel(data = df2,
                            aes(y = pos, label = paste0(100*value, "%")),
                            size = 4.5, nudge_x = 0.75, show.legend = FALSE) +
  guides(fill = guide_legend(title = "Group")) +
  theme_void()
ggsave("./fig/tx_num_pie.pdf")


## 2.plot tx totalLength pie
tx.totalLen <- data.frame(type=bed.list.df$type, total.length=bed.list.df$total.length)
tx.totalLen <- tx.totalLen[!duplicated(paste0(tx.totalLen$type,tx.totalLen$total.length)),]
tx.totalLen
colnames(tx.totalLen) <- c("Var1","Freq")
# table(tx.totalLen$type)

tx.totalLen$Freq <- tx.totalLen$Freq/sum(tx.totalLen$Freq)
tx.totalLen$Freq <- round(tx.totalLen$Freq, digits = 4) #format(tx.num$Freq,scientific = 4)


library(ggplot2)
#str(df)
table(df$group)
df <- tx.totalLen
colnames(df) <- c("group","value")
df$value <- as.numeric(df$value)

df$group <- factor(df$group,levels = rnas)
df <- df[order(df$group,decreasing = F),]
#df <- df[!grepl("repeats",as.character(df$group)),]
#df <- df[!grepl("for|rev",perl = T,as.character(df$group)),]
df2 <- df %>% 
  mutate(csum = rev(cumsum(rev(value))), 
         pos = value/2 + lead(csum, 1),
         pos = if_else(is.na(pos), value/2, pos))
ggplot(df, aes(x = "" , y = value, fill = group)) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(name="Group",values = colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(19)) + # c(RColorBrewer::brewer.pal(8, "Accent"),RColorBrewer::brewer.pal(9, "Set1"))
  ggrepel::geom_label_repel(data = df2,
                            aes(y = pos, label = paste0(100*value, "%")),
                            size = 4.5, nudge_x = 0.75, show.legend = FALSE) +
  guides(fill = guide_legend(title = "Group")) +
  theme_void()
ggsave("./fig/tx_len_pie.pdf")


## 3.plot tx cum length ratio line plot
bed.list.df[1:3,]
bed.list.df$cumsum.wid.ratio1 <- cumsum(bed.list.df$wid/sum(bed.list.df$wid))
bed.list.df$cumsum.rank.ratio1 <- (1:nrow(bed.list.df))/nrow(bed.list.df)


library(ggplot2)
df <- bed.list.df[,c("type","cumsum.wid.ratio","cumsum.rank.ratio")]
df <- df[sample(1:nrow(df),0.05*nrow(df),replace = F),]
df$type <- factor(df$type,levels = rnas)
df <- df[order(df$type,decreasing = F),]
#sample(replace = F)
ggplot(df, aes(x = cumsum.rank.ratio , y = cumsum.wid.ratio, color = factor(type))) +
  geom_line(size=1) +
  scale_color_manual(name="Group",values = colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(19)) + # c(RColorBrewer::brewer.pal(8, "Accent"),RColorBrewer::brewer.pal(9, "Set1"))
  guides(fill = guide_legend(title = "Group")) +
  theme_classic() +
  theme(plot.title = element_text(size=24,color="black",hjust = 0.5),
        axis.title = element_text(size=24,color="black"), 
        axis.text = element_text(size=24,color="black"), # ,face="bold"
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 0, hjust =1 , vjust = 0.5),
        axis.text.y = element_text(angle = 0, hjust = 0.9, vjust = 0.5),
        panel.grid=element_blank(),
        legend.position = "right",
        axis.line = element_line(colour = "black"),
        legend.text = element_text(size= 24),
        legend.title= element_text(size= 24,face="bold"),
        strip.text = element_text(size = 20,face="bold"),
        strip.background = element_blank() )
ggsave("./fig/tx_cumsumlength_line.pdf",width = 10,height = 7)


df <- bed.list.df[,c("type","cumsum.wid.ratio1","cumsum.rank.ratio1")]
df <- df[sample(1:nrow(df),0.05*nrow(df),replace = F),]
df$type <- factor(df$type,levels = rnas)
df <- df[order(df$type,decreasing = F),]
#sample(replace = F)
ggplot(df, aes(x = cumsum.rank.ratio1 , y = cumsum.wid.ratio1, color = factor(type))) +
  geom_line(size=1) +
  scale_color_manual(name="Group",values = colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(19)) + # c(RColorBrewer::brewer.pal(8, "Accent"),RColorBrewer::brewer.pal(9, "Set1"))
  guides(fill = guide_legend(title = "Group")) +
  theme_classic() +
  theme(plot.title = element_text(size=24,color="black",hjust = 0.5),
        axis.title = element_text(size=24,color="black"), 
        axis.text = element_text(size=24,color="black"), # ,face="bold"
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 0, hjust =1 , vjust = 0.5),
        axis.text.y = element_text(angle = 0, hjust = 0.9, vjust = 0.5),
        panel.grid=element_blank(),
        legend.position = "right",
        axis.line = element_line(colour = "black"),
        legend.text = element_text(size= 24),
        legend.title= element_text(size= 24,face="bold"),
        strip.text = element_text(size = 20,face="bold"),
        strip.background = element_blank() )
ggsave("./fig/tx_cumsumlength_line2.pdf",width = 10,height = 7)
# df$cumsum.rank.ratio <- round(df$cumsum.rank.ratio, digits = 2) 
# head(df)
# df <- df[df$cumsum.rank.ratio==0.80,]
# df <- df[order(df$cumsum.wid.ratio),]








# sum tx MAPQ ----------------------------------------------------------------
library(Rsamtools)

#UMI-tools: modeling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy
#they have both identical alignment coordinates and identical UMI sequences (Fig. 1A).

#bamFile是一个bam文件对象，如上图，它有多个属性，如isopen、yieldSize、obeyQname、asMates等等
#opar <- par()
#par(mfrow=c(16,1))
#pdf("test_FTA-17.pdf",width = 5,height = 4)

# dedup reads
mapq.list <- list()
for (i in c("miRNA","lncRNA","mRNA","piRNA","snoRNA","snRNA","srpRNA","tRNA","tucpRNA","Y_RNA","intron.for","intron.rev","promoter.for","promoter.rev","enhancer.for","enhancer.rev") ){
  print(i)
  #i <- "mRNA"
  #x <- paste0("/Users/baopengfei/Desktop/lulab/tmp/projects/WCHSU-FTC/output/FTC_small/tbam/FTA-17_1/bam-deduped/",i,".bam")
  #x <- paste0("/Users/baopengfei/Desktop/lulab/tmp/projects/WCHSU-FTC/output/FTC_small/tbam/FTA-17_1/bam/",i,".bam")
  x <- paste0("/Users/baopengfei/Desktop/lulab/tmp/projects/WCHSU-FTC/test/tbam/clFTA-17/bam-deduped/",i,".bam")
  #x <- paste0("/Users/baopengfei/Desktop/lulab/tmp/projects/WCHSU-FTC/output/FTC_long/tbam/clFTA-17/bam/",i,".bam")
  bamFile <- Rsamtools::BamFile(x)
  #seqinfo(bamFile) # 参考序列信息
  #countBam(bamFile) # bam文件的read数及核酸数等信息
  #idxstatsBam(bamFile) # 序列比对情况
  #scanBamHeader(bamFile) # head信息
  
  #aln是一个列表，其列表元素的名称基本上是bam文件的列名
  aln <- Rsamtools::scanBam(bamFile) # read bam as a whole, not by single read
  #length(aln) # 长度为1，因为sanBam可以读取多个bam文件，并合成一个list文件
  #class(aln) # list
  aln1<-aln[[1]] #取出bam文件的内容，并查看bam文件各个列下的第一个read信息
  #names(aln1) # bam文件的各列名，"qname"  "flag"   "rname"  "strand" "pos" 等共13个
  #lapply(aln1, function(x) x[1])
  #hist(aln1$qwidth, col = "grey50", xlim = c(0,100))
  #summary(aln1$isize)
  mapq.list[[i]] <- data.frame(mapq=aln1$mapq,qwidth=aln1$qwidth,type=i)
  #hist(aln1$mapq,breaks = max(.005*length(aln1$mapq),100), col = "grey50", main = i, xlim = c(0,45))
}
mapq.list.df <- do.call(rbind,mapq.list)
mapq.list.df[1:3,]
#dim(mapq.list.df)
#table(mapq.list.df$type)
mapq.list.df$group <- "dedup" # "all"
#mapq.list.df2 <- mapq.list.df

# all reads
mapq.list2 <- list()
for (i in c("miRNA","lncRNA","mRNA","piRNA","snoRNA","snRNA","srpRNA","tRNA","tucpRNA","Y_RNA","intron.for","intron.rev","promoter.for","promoter.rev","enhancer.for","enhancer.rev") ){
  print(i)
  #x <- paste0("/Users/baopengfei/Desktop/lulab/tmp/projects/WCHSU-FTC/output/FTC_small/tbam/FTA-17_1/bam-deduped/",i,".bam")
  #x <- paste0("/Users/baopengfei/Desktop/lulab/tmp/projects/WCHSU-FTC/output/FTC_small/tbam/FTA-17_1/bam/",i,".bam")
  #x <- paste0("/Users/baopengfei/Desktop/lulab/tmp/projects/WCHSU-FTC/output/FTC_long/tbam/clFTA-17/bam-deduped/",i,".bam")
  x <- paste0("/Users/baopengfei/Desktop/lulab/tmp/projects/WCHSU-FTC/test/tbam/clFTA-17/bam/",i,".bam")
  
  bamFile <- Rsamtools::BamFile(x)
  aln <- Rsamtools::scanBam(bamFile) # read bam as a whole, not by single read
  aln1<-aln[[1]] #取出bam文件的内容，并查看bam文件各个列下的第一个read信息
  mapq.list2[[i]] <- data.frame(mapq=aln1$mapq,qwidth=aln1$qwidth,type=i)
}
mapq.list.df2 <- do.call(rbind,mapq.list2)
mapq.list.df2$group <- "all" # "all"

mapq.list.df3 <- rbind(mapq.list.df,mapq.list.df2)
table(mapq.list.df3$group)
mapq.list.df3$type <- factor(mapq.list.df3$type,levels =  c("miRNA","lncRNA","mRNA","piRNA","snoRNA","snRNA","srpRNA","tRNA","tucpRNA","Y_RNA","intron.for","intron.rev","promoter.for","promoter.rev","enhancer.for","enhancer.rev"))
mapq.list.df3$group <- factor(mapq.list.df3$group,levels = c("all","dedup"))
str(mapq.list.df3)

# plot violin/boxplot
library(ggplot2)
library(gghalves)
#install.packages("gghalves")
#set.seed(1)
#b <- runif(nrow(mapq.list.df3), -0.2, 0.2)
#mapq.list.df4 <- mapq.list.df3[sample(1:nrow(mapq.list.df3),0.1*nrow(mapq.list.df3),replace=F),]
ggplot(mapq.list.df3,aes(x=type,y=qwidth,fill=group))+  # ,label=label
  #gghalves::geom_half_violin() + 
  #geom_split_violin(draw_quantiles = T,scale = T,width=2) +
  # geom_jitter(alpha=0.2,
  #             position=position_jitterdodge(jitter.width = 0.35, 
  #                                           jitter.height = 0, 
  #                                           dodge.width = 0.8))+
  # geom_boxplot(alpha=0.2,width=0.45,
  #              position=position_dodge(width=0.8),
  #              size=0.75,outlier.colour = NA)+
  geom_violin(width=3,
              position=position_dodge(width=0.8),
              size=0.5)+
  #geom_point(aes(x=(as.numeric(group)+b),y=mapq),  
  #           shape=21,size=2,color="black",stroke=0.5,alpha=0.1) +
  ggsci::scale_fill_jco()+
  #scale_x_discrete(label=c("CT","CN","PT","PN","RT","RN"))+
  ylim(c(0,200)) + 
  #labs(title = g)+
  #ylab("ssGSEA score")+
  #ylab("expression level")+
  # ggpubr::stat_compare_means(
  #   label.x.npc = 0.5,hide.ns=F,size =5,step.increase = 0.08,#paired = TRUE,
  #   aes(label = ..p.signif..,group=group),  
  #   symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01,1),symbols = c( "***", "**", "*", "ns")),
  #   method = "wilcox.test") +
  theme_classic(base_size=12) + 
  theme(#axis.ticks.x=element_blank(),  
    #strip.text.y = element_blank(),
    #strip.text.x = element_text(face="bold",family="arial",size=20),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size = 20,hjust = 1,vjust = 0.5,angle = 90),
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size=20),
    legend.position = "right")
#ggsave("./sensitive_small_mapq.pdf",width = 12,height = 7)
#ggsave("./sensitive_long_mapq.pdf",width = 12,height = 7)
#ggsave("./sensitive_small_qwidth.pdf",width = 12,height = 7)
ggsave("./sensitive_long_qwidth.pdf",width = 12,height = 7)
summary(mapq.list.df3$qwidth)


# other func.
#sortBam #排序，具体可以在R中查看帮助（?sortBam）
#indexBam #索引，具体可以在R中查看帮助（?indexBam）
#mergeBam #合并，具体可以在R中查看帮助（?mergeBam）

# summary
#1.
#看起来bowtie2的 --very-sensitive AND --very-fast 对 map.ratio 结果略有影响，sensitive比对率略低一些
#对MAPQ的影响也很有限，mapq的分布差别不大
#以后统一 --very-fast mode





# sum sample-wise tx gini (deprecated, need run single-base cpm & stat: *_singleBase.txt) -----
chr.size <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt",data.table = F,header = T,sep = "\t",stringsAsFactors = F)
chr.size[1:3,]

#bw0 <- data.table::fread("tmp/elFTC-11_L5_singleBase.txt",data.table = F,header = F,sep = "\t",stringsAsFactors = F)
bw0 <- data.table::fread("tmp/csFTC-10_1_singleBase.txt",data.table = F,header = F,sep = "\t",stringsAsFactors = F)
colnames(bw0) <- c("chr","s","e","name","score","strand","pos","count")
bw0[1:30,]

bw0$tx.type <- chr.size$transcript_type[match(as.character(bw0$chr),as.character(chr.size$transcript_id))]
bw0[1:3,]
table(bw0$tx.type)


bw0$cpm <- edgeR::cpm(y = bw0[,8], prior.count = 0, normalized.lib.sizes = F, log = FALSE)
summary(bw0$cpm)
#cor.test(bw0$cpm,bw0$count) # 1

bw0.tbl <- as_tibble(bw0) %>% 
  dplyr::group_by(chr,tx.type) %>% 
  dplyr::summarise(gini.logcount=edgeR::gini(log2(count+1)),gini.count=edgeR::gini(count),gini.cpm=edgeR::gini(cpm)) %>% 
  dplyr::select(chr,gini.logcount,gini.count,gini.cpm,tx.type)
bw0.tbl <- na.omit(bw0.tbl)
#cor.test(bw0.tbl$gini.logcount,bw0.tbl$gini.count) # 0.93
#summary(bw0.tbl$gini.logcount)
bw0.df <- as.data.frame(bw0.tbl)
bw0.df$tx.type <- gsub("\\.for|\\.rev|_for|_rev|\\.ratio","",bw0.df$tx.type,perl=T)
rna <- c("rRNA","pri_miRNA","lncRNA", "mRNA", "piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA')
dna <- c("intron","promoter", "enhancer","repeats") 
bw0.df$tx.type <- factor(bw0.df$tx.type,levels = c(rna,dna))
#hist(bw0.tbl$gini.logcount)
#hist(bw0.tbl$gini.count)
table(bw0.df$tx.type)

library(ggplot2)
my_theme <- theme_minimal() + 
  theme(plot.title = element_text(size = 24,color="black",hjust = 0.5,face="bold"),
        axis.title = element_text(size = 24,color ="black"), 
        axis.text = element_text(size= 20,color = "black"),
        # axis.text.x = element_text(size= 16,color = "black",angle = 90,vjust = 0.5,hjust = 1),
        #panel.grid=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x = element_blank(),
        # panel.grid.major.y = element_line(color = "grey50",linetype = "dashed"), #size= 1,
        #panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        strip.text = element_text(size = 20),
        strip.background = element_rect(fill="white",color="white"),
        #axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, color = c(rep("grey50",length(other)),rep("steelblue",length(rna)),rep("firebrick",length(dna)))), #  
        legend.position = "none",#c(.25,.6),
        legend.text = element_text(size= 16),
        legend.title= element_text(size= 16) 
        )
p1 <- ggplot(bw0.df, aes(x=gini.count))+
  geom_histogram(color="black",bins = 20)+
  # geom_bar(stat="density")+
  xlim(c(0,1))+
  my_theme
p2 <- ggplot(bw0.df, aes(x=gini.count, fill=tx.type))+
  geom_histogram(color="black",bins = 20)+
  # geom_bar(stat="density")+
  xlim(c(0,1))+
  facet_wrap(facets = "tx.type", scales="free", ncol = 3)+
  scale_fill_d3_adaptive()+
  my_theme
ggarrange(p1,p2,nrow = 2,heights = c(1,7),align = "h")
ggsave("gini.density.pdf",width = 12,height = 15)
#not need cpm (only intra-sample calculation, not need cross-sample comparison)
# require(gridExtra)
# pdf("gini.density.pdf",width = 7,height = 14)
# grid.arrange(p1, p2, ncol = 1, heights = c(1, 8))
# dev.off()




# compare potential bias of merge strategy: correlation between sample size AND merged domain length ------------------------------------
s10 <- data.table::fread("/Users/baopengfei/Desktop/lulab/tmp/projects/WCHSU-FTC/output/FTC_long/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax_80samples/domains.bed",data.table = F,sep='\t',stringsAsFactors = F,header = F)
s10[1:3,]
s10$width <- s10$V3-s10$V2
hist(s10$width,breaks = 100,xlim = c(0,2000))
#看起来 >10% 百分比过滤模式本身可以有效减小sample size对merged domain的长度的影响，至少 FTC long RNA-seq 是的，但small RNA-seq 不确定








# sum tx map/dedup/umi ratio (202311 on SLE) --------------------------------------------------------
#setwd("/Users/baopengfei/Desktop/lulab/tmp/projects/WCHSU-FTC/")
#setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/")
setwd("/data2/lulab1/bpf/projects/WCHSU-FTC")
dst <- "SLE"
dedup <- "dedup"

# read.umi <- function(x){
#   #x <- "output/FTC_small/log/dedup/FTA-10_1/Y_RNA.log"
#   #"output/SLE/tbam/20230523-S048-1_BC11/log-dedup/Y_RNA.log"
#   print(x)
#   # sample <- unlist(sapply(strsplit(x,"/",fixed = T),"[",5))
#   rna <- gsub(".log","",basename(x)) 
#   umi <- readLines(x)
#   umi <- umi[grepl("Mean number of unique UMIs per position",umi)]
#   umi <- as.numeric(unlist(sapply(strsplit(umi,": "),"[",2)))
#   umi <- data.frame(sample=sample,type=rna,value=umi)
#   return(umi)
# }
# 
# #umi.df <- lapply(Sys.glob( "output/FTC_small/log/dedup/*/*.log"),read.umi)
# umi.df <- lapply(Sys.glob( "output/SLE/tbam/20230523-S048-1_BC11/log-dedup/*/*.log"),read.umi)
# umi.df <- do.call(rbind,umi.df)
# table(umi.df$type)


# small/SE map ratio (main QC module)
read.SE.map <- function(x){
  #x <- "output/FTC_small/log/tbam/csFTA-10_1/enhancer.for.log"
  #x <- "output/FTC_long/log/clFTA-10/tbam/enhancer.for.log"
  #x <- "output/WSQ_SMARTer_NEB/tbam/NC_ChQ-21_smart_1/log/promoter_for.log"
  # x <- "output/SLE/tbam/20230523-S048-1_BC11/log/Y_RNA.log"
  print(x)
  sample <- unlist(sapply(strsplit(x,"/",fixed = T),"[",4))
  rna <- gsub(".log","",basename(x)) 
  map0 <- readLines(x)
  in.num <- as.numeric(gsub(" reads; of these:","",map0[grepl("reads; of these:",map0)]))
  not.map.num <- as.numeric(unlist(sapply(strsplit(map0[grepl("aligned 0 times",map0)]," (", fixed = T), "[" ,1)))
  uniq.map.num <- as.numeric(unlist(sapply(strsplit(map0[grepl("aligned exactly 1 time",map0)]," (", fixed = T), "[" ,1)))
  multi.map.num <- as.numeric(unlist(sapply(strsplit(map0[grepl("aligned >1 time",map0)]," (", fixed = T), "[" ,1)))
  #map.num <- uniq.map.num + multi.map.num
  #in.num==not.map.num+uniq.map.num+multi.map.num # TURE
  map <- as.numeric(gsub("% overall alignment rate","",map0[6]))
  map <- data.frame(sample=sample,type=rna,overall.align.ratio=map,in.num=in.num, uniq.map.num=uniq.map.num, multi.map.num=multi.map.num,not.map.num=not.map.num)
  return(map)
}

map.df <- lapply(Sys.glob( paste0("output/SLE/tbam/*/log/*.log") ),read.SE.map)
map.df <- do.call(rbind,map.df)
map.df$overall.align.ratio <- map.df$overall.align.ratio/100
#table(map.df$type)
library(tidyr)
#map.df.a <- map.df %>% tidyr::pivot_wider(names_from = type, values_from = value)
#map.df.a2 <- map.df.a
#map.df.a2[,2:ncol(map.df.a2)] <- 1-map.df.a2[,2:ncol(map.df.a2)]
#map.df.a$unmap <- 1-rowSums(map.df.a[,2:ncol(map.df.a)])
#rowSums(map.df.a[,2:ncol(map.df.a)])
#map.df2 <- map.df.a %>% tidyr::pivot_longer(cols = 2:ncol(map.df.a), names_to = "type", values_to = "value")

dat <- as_tibble(map.df)
dat <- dat %>% dplyr::group_by(sample) %>% 
  dplyr::mutate(map.num=uniq.map.num+multi.map.num, 
         all.num=max(in.num),
         unmap.num=min(not.map.num),
         map.ratio=map.num/all.num,
         unmap.ratio=unmap.num/all.num) %>% 
  dplyr::select(sample,type,all.num,map.ratio,unmap.ratio)
dat.wide <- dat %>% tidyr::pivot_wider(names_from = "type", values_from = "map.ratio")
colnames(dat.wide)[1:3] <- c("sample","read_num","unmap")
colnames(dat.wide) <- gsub(".","_",colnames(dat.wide),fixed = T)
colnames(dat.wide) <- gsub("miRNA","pri_miRNA",colnames(dat.wide),fixed = T)
dat.wide <-  dat.wide[,c("sample","read_num","unmap","rRNA","spikein","univec","pri_miRNA","lncRNA", "mRNA", "piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA', "intron_for", "intron_rev", "promoter_for", "promoter_rev", "enhancer_for", "enhancer_rev", "repeats_for", "repeats_rev")]
rowSums(dat.wide[,3:ncol(dat.wide)]) # check: should be all ones !!!
dat <- dat.wide %>% tidyr::pivot_longer(cols = 3:ncol(dat.wide), names_to = "type", values_to = "ratio")

unique(dat$sample)
#table(dat$sample.group)
dat$group <- sample.table0$group[match(dat$sample,sample.table0$sample)]
dat$cell <- sample.table0$cellShort[match(dat$sample,sample.table0$sample)]

tmp <- data.frame(cell=dat$cell,group=dat$group)
tmp <- tmp[!(duplicated(tmp)),]
tmp$cell <- factor(tmp$cell, levels = cellShort)
tmp$group <- factor(tmp$group, levels = groupShort)
tmp <- tmp[order(tmp$group,tmp$cell),]
tmp$lib <- paste(tmp$group,tmp$cell)

dat$lib <- paste(dat$group,dat$cell)
dat$lib <- factor(dat$lib, levels = unique(tmp$lib))
table(dat$lib)
# dat$lib[grepl("NEB",dat$sample)] <- "NEB"
# dat$lib[grepl("NEB_PNK",dat$sample)] <- "NEB_PNK"
# dat$lib[grepl("smart",dat$sample)] <- "SMARTer"
# dat$lib[grepl("smart_PNK",dat$sample)] <- "SMARTer_PNK"

# dat$group <- unlist(sapply(strsplit(dat$sample,"-|_",perl=T),"[",1))
# dat$sample.group <- paste0(dat$group,"|",dat$lib)
# dat$sample.group <- factor(dat$sample.group, levels = c("CRC|NEB","NC|NEB","CRC|NEB_PNK","CRC|SMARTer","NC|SMARTer","CRC|SMARTer_PNK","NC|SMARTer_PNK")) # unique(dat$sample.group)  c("csFTA","csFTC","FTA","FTC")
dat$sample.group <- dat$lib
table(dat$lib)
dat$type <- gsub("miRNA","pri_miRNA",dat$type)
dat$type <- gsub("pri_pri_miRNA","pri_miRNA",dat$type)
dat$RNA <- gsub("\\.for|\\.rev|_for|_rev|\\.ratio","",dat$type,perl=T)
dat$RNA <- gsub("spikein_small|spikein_long","spikein",dat$RNA,perl=T)
other <- c("spikein","univec")
rna <- c("rRNA","pri_miRNA","lncRNA", "mRNA", "piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA')
#dna <- c("intron.for", "intron.rev", "promoter.for", "promoter.rev", "enhancer.for", "enhancer.rev", "repeats.for", "repeats.rev")
dna <- c("intron","promoter", "enhancer","repeats") # repeats
dat$RNA <- factor(dat$RNA,levels = c(other,rna,dna,"unmap"))
dat$source <- ifelse(as.character(dat$RNA) %in% rna,"RNA","DNA")
dat$source[as.character(dat$RNA) %in% other] <- "other"
dat$source[as.character(dat$RNA) %in% "unmap"] <- "unmap"
dat$source <- factor(dat$source,levels = c("other","RNA","DNA","unmap"))
summary(dat$RNA)

### plot stack bar
library(ggplot2)
# order levels by a second numeric variable: unmap ratio, not by sample id/type 
dat.tmp <- dat %>% dplyr::filter(RNA=="unmap") %>% dplyr::select(sample,ratio) %>% dplyr::rename(sample=sample,unmap=ratio)
dat <- left_join(dat,dat.tmp)
dat <- dat[order(dat$unmap,dat$sample,decreasing = F),]
dat$sample <- factor(dat$sample,levels = unique(dat$sample))
# dat %>% 
#   mutate(newsample=fct_reorder(sample, unmap)) %>%
dat.sum <- dat %>% 
  dplyr::group_by(RNA,sample.group) %>% 
  dplyr::summarize(mean.ratio=mean(ratio,trim=0.1,na.rm=T)) # %>% 
  # dplyr::group_by(sample.group) %>% 
  # dplyr::mutate(sum.ratio=sum(mean.ratio),norm.mean.ratio=mean.ratio/sum.ratio) # 
ggplot(dat.sum,aes(x=sample.group,y=mean.ratio,fill=RNA))+ 
  geom_bar(position = "fill",width = 0.8,size=0.3, stat = "identity", color="black") +
  #geom_errorbar()+
  xlab("")+
  ylab(paste0("Reads mapped ratio"))+
  geom_vline(xintercept = 0,linetype="dashed",color="grey10")+
  #ggsci::scale_fill_d3()+
  # scale_fill_nejm_adaptive(name="Species")+
  
  scale_fill_manual(name="Species",values = c(pal_nejm_adaptive()(length(c(other,rna,dna,"unmap"))-1),"grey70") )+ # colorRampPalette(RColorBrewer::brewer.pal(8, "Set3"))(length(c(other,rna,dna,"unmap")))
  #paletteer::scale_fill_paletteer_d("RColorBrewer::Set1")+
  #scale_y_continuous(trans = "log10")+
  #xlim(c(0,500))+
  #geom_hline(yintercept = c(0))+
  #ggraph::scale_fill_viridis(option = "C") + # name = "value",
  theme_bw() + 
  theme(plot.title = element_text(size = 24,color="black",hjust = 0.5,face="bold"),
        axis.title = element_text(size = 24,color ="black"), 
        axis.text = element_text(size= 20,color = "black"),
        axis.text.x = element_text(size= 20,color = "black",angle = 90,vjust = 0.5,hjust = 1),
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
ggsave( paste0("./",dst,"_group_small_map_ratio.pdf"), width = 24,height = 12)




# long/PE map ratio
read.PE.map <- function(x){
  #x <- "output/FTC_small/log/tbam/csFTA-10_1/enhancer.for.log"
  #x <- "output/FTC_long/log/clFTA-10/tbam/enhancer.for.log"
  print(x)
  sample <- unlist(sapply(strsplit(x,"/",fixed = T),"[",4))
  rna <- gsub(".log","",basename(x)) 
  map0 <- readLines(x)
  in.num <- as.numeric(gsub(" reads; of these:","",map0[grepl("reads; of these:",map0)]))
  not.map.num <- as.numeric(unlist(sapply(strsplit(map0[grepl("aligned concordantly 0 times",map0)]," (", fixed = T), "[" ,1)))
  uniq.map.num <- as.numeric(unlist(sapply(strsplit(map0[grepl("aligned concordantly exactly 1 time",map0)]," (", fixed = T), "[" ,1)))
  multi.map.num <- as.numeric(unlist(sapply(strsplit(map0[grepl("aligned concordantly >1 time",map0)]," (", fixed = T), "[" ,1)))
  #map.num <- uniq.map.num + multi.map.num
  map <- as.numeric(gsub("% overall alignment rate","",map0[6]))
  map <- data.frame(sample=sample,type=rna,overall.align.ratio=map,in.num=in.num, uniq.map.num=uniq.map.num, multi.map.num=multi.map.num,not.map.num=not.map.num)
  return(map)
}


map.df <- lapply(Sys.glob( "output/FTC_long/log/*/tbam/*.log"),read.PE.map)
map.df <- do.call(rbind,map.df)
map.df$overall.align.ratio <- map.df$overall.align.ratio/100
#table(map.df$type)
library(tidyr)
#map.df.a <- map.df %>% tidyr::pivot_wider(names_from = type, values_from = value)
#map.df.a2 <- map.df.a
#map.df.a2[,2:ncol(map.df.a2)] <- 1-map.df.a2[,2:ncol(map.df.a2)]
#map.df.a$unmap <- 1-rowSums(map.df.a[,2:ncol(map.df.a)])
#rowSums(map.df.a[,2:ncol(map.df.a)])
#map.df2 <- map.df.a %>% tidyr::pivot_longer(cols = 2:ncol(map.df.a), names_to = "type", values_to = "value")


## plot stacked barplot
dat <- as_tibble(map.df)
dat <- dat %>% group_by(sample) %>% 
  mutate(map.num=uniq.map.num+multi.map.num, 
         all.num=max(in.num),
         unmap.num=min(not.map.num),
         map.ratio=map.num/all.num,
         unmap.ratio=unmap.num/all.num) %>% 
  select(sample,type,map.ratio,unmap.ratio) 
dat <- na.omit(dat) # some NA exist for long RNA-seq
dat.wide <- dat %>% tidyr::pivot_wider(names_from = "type", values_from = "map.ratio")
rowSums(dat.wide[,2:ncol(dat.wide)],na.rm = T) # should be all ones !!!
dat <- dat.wide %>% tidyr::pivot_longer(cols = 2:ncol(dat.wide), names_to = "type", values_to = "ratio")

dat$sample.group <- unlist(sapply(strsplit(dat$sample,"-"),"[",1))
unique(dat$sample.group)
dat$sample.group <- factor(dat$sample.group, levels = c("clFTA","clFTC","elFTA","elFTC"))
#table(dat$sample.group)
dat$type <- gsub("miRNA","pri_miRNA",dat$type)
dat$RNA <- gsub(".for","",dat$type)
dat$RNA <- gsub(".rev","",dat$RNA)
dat$RNA <- gsub(".ratio","",dat$RNA)
other <- c("rRNA","spikein","univec")
rna <- c("pri_miRNA","lncRNA", "mRNA", "piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA')
#dna <- c("intron.for", "intron.rev", "promoter.for", "promoter.rev", "enhancer.for", "enhancer.rev", "repeats.for", "repeats.rev")
dna <- c("intron","promoter", "enhancer","repeats") # repeats
dat$RNA <- factor(dat$RNA,levels = c(other,rna,dna,"unmap"))
dat$source <- ifelse(as.character(dat$RNA) %in% rna,"RNA","DNA")
dat$source[as.character(dat$RNA) %in% other] <- "other"
dat$source[as.character(dat$RNA) %in% "unmap"] <- "unmap"
dat$source <- factor(dat$source,levels = c("other","RNA","DNA","unmap"))
summary(dat$RNA)

head(dat[is.na(dat$ratio),],3)
### plot stack bar
library(ggplot2)
# order levels by a second numeric variable: unmap ratio, not by sample id/type 
dat.tmp <- dat %>% filter(RNA=="unmap") %>% select(sample,ratio) %>% rename(sample=sample,unmap=ratio)
dat <- left_join(dat,dat.tmp)
dat <- dat[order(dat$unmap,dat$sample,decreasing = F),]
dat$sample <- factor(dat$sample,levels = unique(dat$sample))
# dat %>% 
#   mutate(newsample=fct_reorder(sample, unmap)) %>%
ggplot(dat,aes(x=sample,y=ratio,fill=RNA))+ 
  geom_bar(position = "stack",width = 0.8,size=0.3, stat = "identity", color="black") +
  #geom_errorbar()+
  xlab("")+
  ylab(paste0("Reads genome ratio"))+
  geom_vline(xintercept = 0,linetype="dashed",color="grey10")+
  #ggsci::scale_fill_d3()+
  scale_fill_manual(name="Species",values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set3"))(length(c(other,rna,dna,"unmap"))))+ # c(RColorBrewer::brewer.pal(8, "Accent"),RColorBrewer::brewer.pal(9, "Set1"))
  #paletteer::scale_fill_paletteer_d("RColorBrewer::Set1")+
  #scale_y_continuous(trans = "log10")+
  #xlim(c(0,500))+
  #geom_hline(yintercept = c(0))+
  #ggraph::scale_fill_viridis(option = "C") + # name = "value",
  theme_bw() + 
  theme(plot.title = element_text(size = 24,color="black",hjust = 0.5,face="bold"),
        axis.title = element_text(size = 24,color ="black"), 
        axis.text = element_text(size= 20,color = "black"),
        axis.text.x = element_text(size= 16,color = "black",angle = 90,vjust = 0.5,hjust = 1),
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
ggsave("./long_map_ratio.pdf",width = 15,height = 10)






### saturation plot: detected tx,domain numb (UMI data)
cd /BioII/lulab_b/baopengfei/projects/WCHSU-FTC
for j in 10 50 99
do
echo $j
  for i in `cat exSeek-dev/data/FTC_small/sample_ids_test.txt`
  do 
  #echo $i
  cleanfqNum=`cat output/FTC_small/tbam_${j}_1M/${i}/log/rule | grep -v "+" | grep "use fq num" |sed s/"use fq num: "/""/g` #`zcat output/FTC_small/trimmed/${i}.fastq.gz | wc -l`;
  bedNum=`zcat output/FTC_small/call_domain_withRepeats_dedup_${j}_1M/tbed_long_RNA_EM/${i}.bed.gz | wc -l`;
  txNum=`zcat output/FTC_small/call_domain_withRepeats_dedup_${j}_1M/tbed_long_RNA_EM/${i}.bed.gz | cut -f 1 | sort | uniq | wc -l`;
  peakNum=`cat output/FTC_small/call_domain_withRepeats_dedup_${j}_1M/domains_localmax_by_sample_EM/b5_d05_p01/${i}.bed | wc -l`;
  echo -e "$i,$cleanfqNum,$bedNum,$txNum,$peakNum";
  done > small.test.stat.${j}
done
for j in 10 50 # 99
do
echo $j
for i in `cat exSeek-dev/data/FTC_long/sample_ids_test.txt`
do 
#echo $i
cleanfqNum=`cat output/FTC_long/tbam_${j}_5M/${i}/log/rule | grep -v "+" | grep "use fq num" |sed s/"use fq num: "/""/g`; # zcat output/FTC_long/trimGC/${i}_1.fastq.gz | wc -l
bedNum=`zcat output/FTC_long/call_domain_withRepeats_dedup_${j}_5M/tbed_long_RNA_EM/${i}.bed.gz | wc -l`;
txNum=`zcat output/FTC_long/call_domain_withRepeats_dedup_${j}_5M/tbed_long_RNA_EM/${i}.bed.gz | cut -f 1 | sort | uniq | wc -l`;
peakNum=`cat output/FTC_long/call_domain_withRepeats_dedup_${j}_5M/domains_localmax_by_sample_EM/b5_d05_p01/${i}.bed | wc -l`;
echo -e "$i,$cleanfqNum,$bedNum,$txNum,$peakNum";
done > long.test.stat.${j}
done
# for i in `cat exSeek-dev/data/FTC_long/sample_ids_all.txt`;do 
# cleanfqNum=`zcat output/FTC_long/trimGC/${i}_1.fastq.gz | wc -l`;
# bedNum=`zcat output/FTC_long/call_domain_withRepeats_dedup/tbed_long_RNA_EM/${i}.bed.gz | wc -l`;
# txNum=`zcat output/FTC_long/call_domain_withRepeats_dedup/tbed_long_RNA_EM/${i}.bed.gz | cut -f 1 | sort | uniq | wc -l`;
# peakNum=`cat output/FTC_long/call_domain_withRepeats_dedup/domains_localmax_by_sample_EM/b5_d05_p01/${i}.bed | wc -l`;
# echo -e "$i,$cleanfqNum,$bedNum,$txNum,$peakNum";
# done > long.test.stat
# #cat ${echo "sample,bedNum,txNum,peakNum"} test.stat > test.stat2

# fs <- Sys.glob("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/long.test.stat.*")
fs <- Sys.glob("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/small.test.stat.*")
res <- list()
read.stat <- function(x){
  #x <- fs[1]
  #x <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/small.test.stat.10"
  df <- read.table(x,sep = ",",header = F,row.names = 1)
  #df[1:3,]
  colnames(df) <- c("cleanfqNum","bedNum","txNum","peakNum")
  df$lib <- "long_UMI" #"small_UMI"
  df$path <- x
  df$sample <- rownames(df)
  #res[[x]] <- df
  return(df)
}
res <- lapply(fs,read.stat)
res.df <- do.call(rbind,res)

# df2 <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/long.test.stat",sep = ",",header = F,row.names = 1)
# colnames(df2) <- c("cleanfqNum","bedNum","txNum","peakNum")
# df2$lib <- "long_UMI"

# df3 <- rbind(df,df2)
df3 <- res.df
df3$downsample <- as.numeric(paste0("0.",unlist(sapply(strsplit(df3$path,".",fixed=T),"[",4))))
df3$cleanfqNum <- df3$cleanfqNum/1000000 #*df3$downsample
# df3$cleanfqNum <- df3$cleanfqNum/4/1000000*df3$downsample
#df3$txNum <- df3$txNum/1000
df3$downsample <- factor(df3$downsample)

my_theme <- theme(plot.title = element_text(size = 24,color="black",hjust = 0.5,face="bold"),
      axis.title = element_text(size = 24,color ="black"), 
      axis.text = element_text(size= 20,color = "black"),
      strip.text.x = element_text(size = 20, color = "black"), #face="bold",
      strip.text.y = element_text(size = 20, color = "black"),
      #panel.grid=element_blank(),
      panel.grid.major.x=element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey50",linetype = "dashed"), #size= 1,
      #panel.grid.minor.y = element_blank(),
      panel.border = element_blank(),
      axis.text.x = element_text(angle = 0), #  
      legend.position = "none",#c(.25,.6),
      legend.text = element_text(size= 16),
      legend.title= element_text(size= 16))
p1 <- ggplot(df3)+
  # xlim(c(0,25))+
  # ylim(c(0,100000))+
  geom_point(aes(x=cleanfqNum,y=txNum,color=downsample,group=sample))+
  geom_smooth(aes(x=cleanfqNum,y=txNum,group=sample),formula= y~x,color="salmon") +
  theme_bw()+
  geom_vline(xintercept = c(2,5), linetype="dashed")+
  my_theme+
  facet_grid(lib~., scales ="free")
p2 <- ggplot(df3)+
  # xlim(c(0,25))+
  # ylim(c(0,15000))+
  geom_point(aes(x=cleanfqNum,y=peakNum,color=downsample,group=sample))+
  geom_vline(xintercept = c(2,5), linetype="dashed")+
  geom_smooth(aes(x=cleanfqNum,y=peakNum,group=sample), formula= y~x,color="salmon") +
  theme_bw()+
  my_theme+
  facet_grid(lib~., scales ="free")
ggarrange(p1,p2,nrow = 1)



# prepare full-length region -----------------------------------
dst <- "GSE110381_diff"  #"PRJNA540919_diff"
#chr.size <- data.table::fread("/Users/baopengfei/Desktop/lulab/tmp/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/transcriptome_genome_sort_uniq",data.table = F,header = F,sep = "\t",stringsAsFactors = F)
chr.size <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/transcriptome_genome_sort_uniq_newTxID",data.table = F,header = F,sep = "\t",stringsAsFactors = F)

#only select those that can convert to gbed (10 RNAs) ? or might take long time for full length enhancer
domains <- data.table::fread(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_all/cfpeakCNN/b5_d50_p1.bed"),data.table = F,header = F,sep = "\t",stringsAsFactors = F)
#domains <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE110381_diff/call_peak_all/expeak/b5_d50_p1.bed",data.table = F,header = F,sep = "\t",stringsAsFactors = F)

#head(domains,3)
dim(domains)
domains <- domains[!duplicated(domains$V1),]
dim(domains)
domains$V5 <- "."
domains$V4 <- domains$V1
domains$V2 <- 0
#head(chr.size,3)
domains$V3 <- chr.size$V2[match(domains$V1,chr.size$V1)]
#head(domains,3)
#data.table::fwrite(domains,"/Users/baopengfei/Desktop/lulab/tmp/cooperation/yinjianhua/output/TCGA-LIHC_small-call-domain/call_domain/domains_localmax/domains_fullLength.bed",quote = F,sep = "\t",row.names = F,col.names = F)
#data.table::fwrite(domains,"/Users/baopengfei/Desktop/lulab/tmp/projects/WCHSU-FTC/output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/domains_fullLength.bed",quote = F,sep = "\t",row.names = F,col.names = F)
#data.table::fwrite(domains,"/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/domains_fullLength.bed",quote = F,sep = "\t",row.names = F,col.names = F)
data.table::fwrite(domains,paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_all/cfpeakCNN/fullLength.bed"),quote = F,sep = "\t",row.names = F,col.names = F)
#data.table::fwrite(domains,"/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE110381_diff/call_peak_all/expeak/fullLength.bed",quote = F,sep = "\t",row.names = F,col.names = F)


## get shuffle 3 rand regions
cd /BioII/lulab_b/baopengfei/projects/WCHSU-FTC
#cd /BioII/lulab_b/baopengfei/cooperation/yinjianhua # output/TCGA-LIHC_small-call-domain/call_domain/domains_localmax/domains_fullLength.bed
for i in 1 2 3
do
echo $i
bedtools shuffle -chrom -noOverlapping -seed $i \
-i output/TCGA-LIHC_small-call-domain/call_domain/domains_localmax/domains.bed \
-g /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/transcriptome_genome_sort_uniq \
>  output/TCGA-LIHC_small-call-domain/call_domain/domains_localmax/domains_rand${i}.bed
done

## convert tx to gn bed
#cd /Users/baopengfei/Desktop/lulab/tmp/projects/WCHSU-FTC
#i="output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/domains.bed"
i="output/FTC_long/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/domains.bed"
#exSeek-dev/genome/hg38/bed/long_RNA-tRNA-long_DNA_noRepeats.bed
{{
  grep -v '^chr' $i | exSeek-dev/bin/tbed2gbed <(cat exSeek-dev/genome/hg38/bed/{long_RNA,tRNA,pri_miRNA,piRNA}.bed ) /dev/stdin /dev/stdout
  awk 'BEGIN{{OFS="\t";FS="\t"}}/^chr/{{print $1,$2,$3,$4,$5,$6,0,0,0,1,$3-$2,0}}' $i
}} | bedtools sort \
> /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/FTC_long/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/domains_gn.bed

for i in 1 2 3
do
echo $i
in="output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/domains_rand${i}"
#exSeek-dev/genome/hg38/bed/long_RNA-tRNA-long_DNA_noRepeats.bed
{{
  grep -v '^chr' ${in}.bed | exSeek-dev/bin/tbed2gbed <(cat exSeek-dev/genome/hg38/bed/long_RNA.bed exSeek-dev/genome/hg38/bed/tRNA.bed) /dev/stdin /dev/stdout
  awk 'BEGIN{{OFS="\t";FS="\t"}}/^chr/{{print $1,$2,$3,$4,$5,$6,0,0,0,1,$3-$2,0}}' ${in}.bed
}} | bedtools sort \
> ${in}_gn.bed
done


# peak protection ratio ---------------------------------------------------
## POSTAR3 RBP
# cd /BioII/lulab_b/baopengfei/projects/WCHSU-FTC
# bedtools intersect -wao -s \
# -a /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/FTC_long/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/domains_gn.bed \
# -b /BioII/lulab_b/baopengfei/shared_reference/RBP/splitByRBP/all.bed > \
# /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/FTC_long/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/domains_gn_interectRBP.bed
# 
# for i in 1 2 3 
# do 
# echo $i
# in="output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/domains_rand${i}_gn"
# bedtools intersect -wao -s \
# -a ${in}.bed \
# -b /BioII/lulab_b/baopengfei/shared_reference/RBP/splitByRBP/all.bed > \
# ${in}_interectRBP.bed
# done

setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC")
#setwd("/Users/baopengfei/Desktop/lulab/tmp/projects/WCHSU-FTC")
rbp <- data.table::fread("output/FTC_long/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/domains_gn_interectRBP.bed",data.table = F,header = F,sep = "\t",stringsAsFactors = F)
#rbp <- data.table::fread("output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/domains_gn_interectRBP.bed",data.table = F,header = F,sep = "\t",stringsAsFactors = F)
idx <- colnames(rbp)[ncol(rbp)]
rbp$domain.ratio <- as.numeric(rbp[[idx]]/(rbp$V3-rbp$V2))
rbp <- rbp[order(rbp$V4,rbp$domain.ratio,decreasing = T),] # rm dup
rbp <- rbp[!duplicated(rbp$V4),]
#tail(rbp)
hist(rbp[[idx]]/(rbp$V3-rbp$V2))
#rbp$RBP_protect <- rbp$domain.ratio >= 0.1
#table(rbp$domain.ratio >= 0.1)
#FALSE  TRUE 
#5438  2190 
#table(rbp[[idx]] >= 10)
#FALSE  TRUE 
#5634  1994
rbp$RBP_protect <- rbp[[idx]] >= 10
table(rbp$RBP_protect)
#only long_RNA,tRNA gn peak
#FALSE  TRUE 
#5438  2190 
#hist(rbp$V3-rbp$V2,xlim = c(0,100),breaks = 100000)
#table((rbp$V3-rbp$V2)<=100)
rbp <- rbp[,c("V4","RBP_protect")]
rownames(rbp) <- rbp$V4


## G4, i-motif
cd /BioII/lulab_b/baopengfei/projects/WCHSU-FTC
ref="/BioII/lulab_b/baopengfei/shared_reference/structure/G4iMGrinder/bedg4im.bed"
bedtools intersect -wao -s \
-a /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/FTC_long/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/domains_gn.bed \
-b ${ref} > \
/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/FTC_long/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/domains_gn_interectG4im.bed

# for i in 1 2 3 
# do 
# echo $i
# in="output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/domains_rand${i}_gn"
# bedtools intersect -wao -s \
# -a ${in}.bed \
# -b ${ref} > \
# ${in}_interectRBP.bed
# done


setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC")
#setwd("/Users/baopengfei/Desktop/lulab/tmp/projects/WCHSU-FTC")
#G4im <- data.table::fread("output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/domains_gn_interectG4im.bed",data.table = F,header = F,sep = "\t",stringsAsFactors = F)
G4im <- data.table::fread("output/FTC_long/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/domains_gn_interectG4im.bed",data.table = F,header = F,sep = "\t",stringsAsFactors = F)
idx <- colnames(G4im)[ncol(G4im)]
G4im$domain.ratio <- as.numeric(G4im[[idx]]/(G4im$V3-G4im$V2))
G4im <- G4im[order(G4im$V4,G4im$domain.ratio,decreasing = T),] # rm dup
G4im <- G4im[!duplicated(G4im$V4),]
#tail(G4im)
hist(G4im[[idx]]/(G4im$V3-G4im$V2))
#G4im$G4im_protect <- G4im$domain.ratio >= 0.1
G4im$G4im_protect <- G4im[[idx]] >= 10
table(G4im$G4im_protect)
#only long_RNA,tRNA gn peak
#FALSE  TRUE 
#5438  2190 
#hist(G4im$V3-G4im$V2,xlim = c(0,100),breaks = 100000)
#table((G4im$V3-G4im$V2)<=100)
G4im <- G4im[,c("V4","G4im_protect")]
rownames(G4im) <- G4im$V4
G4im[1:3,]


## MFE_v1(RNAfold kcal/mol)
#bedtools slop -g /BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/genome/chrom.size -s -l 20 -r 20 -i /BioII/lulab_b/baopengfei/shared_reference/RBP/tmp/JYF_top500.bed > /BioII/lulab_b/baopengfei/shared_reference/RBP/tmp/JYF_top500_ext.bed
#get genome fasta and MFE (bash)
wkdir="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/FTC_long/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax"
for prefix in domains_gn domains_rand1_gn domains_rand2_gn domains_rand3_gn
do
i="${wkdir}/${prefix}.bed"
echo $i
cat ${i} | awk ' $3-$2 <= 100 {print $0}' > ${i}.tmp  # #filter too large peak (not protected by local MFE), or sometimes too slow
bedtools getfasta -name -s -fi /BioII/lulab_b/shared/genomes/hg38/fasta/genome.fa -bed ${i}.tmp > ${i}.bed.fa
RNAfold  -p -d2 --noLP  ${i}.bed.fa  | grep -E ">|\[" >  ${i}.RNAfold
rm ${i}.tmp
rm ./*.ps
done

conv <- function(x){
  #x <- "JYF_top500"
  #x <- "/Users/baopengfei/Desktop/lulab/tmp/projects/WCHSU-FTC/output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/domains_gn.bed.RNAfold"
  print(x)
  rnafold <- readLines(paste0(x))
  rnafold <- as.data.frame(rnafold)
  #head(rnafold,3)
  colnames(rnafold) <- "V1"
  #rnafold$V3 <- unlist(base::sapply(base::strsplit(rnafold$V1,"[",fixed = T),"[",2))
  V3 <- rnafold$V1[!grepl(">",rnafold$V1)]
  
  rnafold <- rnafold[grepl(">",rnafold$V1),,drop=F]
  mir <- rnafold$V1
  mir <- gsub(">","",mir)
  mir <- unlist(sapply(strsplit(mir,"::",fixed = T),"[",1))
  mir <- unlist(sapply(strsplit(mir,"(",fixed = T),"[",1))
  
  V3 <- unlist(sapply(strsplit(as.character(V3),"[",fixed = T),"[",2))
  rnafold$V3 <- as.numeric(gsub("]","",V3 ))
  
  mfe <- rnafold$V3[rnafold$V3!=""]
  mfe <- as.numeric(mfe[!is.na(mfe)])
  rnafold.out <- data.frame(mir=mir)
  rnafold.out$mfe <- mfe
  #rnafold.out$type <- x
  print(summary(rnafold.out$mfe))
  return(rnafold.out)
}
#outdir <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/"
#outdir <- "/Users/baopengfei/Desktop/lulab/tmp/projects/WCHSU-FTC/output/FTC_long/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/"
outdir <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/FTC_long/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/"
prefix <- c("domains_gn") #,"domains_rand1_gn","domains_rand2_gn","domains_rand3_gn"
options(stringsAsFactors = F)
mfe <- lapply(paste0(outdir,prefix,".bed.RNAfold"),conv)
mfe <- do.call("rbind",mfe)
hist(mfe$mfe)
mfe$Second_structure <- mfe$mfe<=-5
table(mfe$mfe<=-5)
#all long_RNA,tRNA,long_DNA gn peak:
#FALSE  TRUE 
#16941   503
#only long_RNA,tRNA gn peak:
#FALSE  TRUE 
#7114   362 
#mfe <- mfe[mfe$mfe>=-5,] # #only keep MFE >= -5
#table(mfe$group)
head(mfe,4)
#mfe$group <- NULL
#table(duplicated(mfe$mir))
rownames(mfe) <- mfe$mir



## MFE_v2(dinuc shuffle p value)
setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/")
#prefix <- "output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/domains_gn"
prefix <- "output/FTC_long/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/domains_gn"
tmp <- bedtoolsr::bt.slop(i = paste0(prefix,".bed"), b = 20, #l=20, r=20,
                   g = "/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/genome/chrom.size") # output bug: need be "xxxx.bed"
bedtoolsr::bt.getfasta(bed = tmp ,s = T,nameOnly = T,
                       fi = "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/fasta/genome.fa", 
                       fo = paste0(prefix,"_expand.fa") )

fa <- data.table::fread(paste0(prefix,"_expand.fa"),data.table = F,sep = "\t",header = F)
fa <- as.data.frame(fa)
fa$V1 <- gsub("(+)","",fa$V1,fixed = T) # or linux file storage error
fa$V1 <- gsub("(-)","",fa$V1,fixed = T)

dir.create(paste0(dirname(prefix),"/fasta"),showWarnings = F)
for(i in seq(1,(nrow(fa)-1),2)){
  #i <- 1
  if(i %% 500 == 0){
    print(i)
  }
  #fa$V1 <- gsub(">","",fa$V1)
  id <- gsub(">","",fa$V1[i])
  data.table::fwrite(x = fa[i:(i+1),,drop=F], paste0(dirname(prefix),"/fasta/",id,".fa"), sep = "\t",quote = F,col.names = F,row.names = F)
}


PATH=/BioII/lulab_b/jinyunfan/anaconda3/envs/bioinfo_py36/bin:$PATH
#wkdir="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/fasta"
wkdir="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/FTC_long/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/fasta"
cd ${wkdir}
for i in `ls -1 *.fa | tail -n 2007`
do id=`echo $i | sed s/".fa"/""/g`; echo $id ; done \
| parallel -I % -k -j 10 "mkdir -p ${wkdir}/%; cd ${wkdir}/%; mv ${wkdir}/%.fa ${wkdir}/%/%.fa; python3 /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/bin/getStructurePvalueByShuffleDinuc.py \
  %.fa 200 > ${wkdir}/%.txt; cd ${wkdir}; mv ${wkdir}/%/%.fa ${wkdir}/%.fa ; rm -rf ${wkdir}/%" & 
#No more file handles. (too mucl tmp files)
""" # ???


## read in mfe_v2 p
read.MFE.p <- function(x){
  #x <- "output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/fasta/peak_1000.fa.txt"
  tmp <- data.table::fread(x,sep = ":",stringsAsFactors = F,header = F,nThread = 2)
  #colnames(tmp) <- c("id","p")
  tmp <- data.frame(mir=tmp$V2[1],mfe=tmp$V2[2],p=tmp$V2[4])
  return(tmp)
}
fs <- Sys.glob(paste0(dirname(prefix),"/fasta/","*",".txt")) # at most 1547 files ?
library(parallel)
res <- mclapply(fs, read.MFE.p, mc.cores = 4)  # 
res.df <- do.call(rbind,res)
hist(as.numeric(res.df$p))
table(as.numeric(res.df$p)<=0.1)
mfe <- res.df
mfe$mir <- gsub(".fa","",mfe$mir)
rownames(mfe) <- mfe$mir 
colnames(mfe) <- c("mir","originalMFE","p")
mfe$Second_structure <- mfe$p <= 0.1
mfe <- mfe[,c("mir","Second_structure")]


## EV/cf diff 
EV <- read.table("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/lulab/longDomain_diff_EVvsCF_pair.txt",row.names = 1,header = T,sep = "\t",stringsAsFactors = F)
#EV <- read.table("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/lulab/smallDomain_diff_EVvsCF_pair.txt",row.names = 1,header = T,sep = "\t",stringsAsFactors = F)
#EV <- read.table("/Users/baopengfei/Desktop/lulab/tmp/projects/motif-RBP-EDA/output/lulab/longDomain_diff_EVvsCF_pair.txt",row.names = 1,header = T,sep = "\t",stringsAsFactors = F)
EV <- EV[order(EV$pvalue,decreasing = F),]
hist(EV$log2FoldChange)
EV$EV_protect <- EV$padj<=0.1 & EV$log2FoldChange>1
table(EV$EV_protect)
rownames(EV) <- unlist(sapply(strsplit(rownames(EV),"|",fixed = T),"[",4))
EV <- EV[,c("log2FoldChange","EV_protect")]
head(EV)


## combine mat
pl <- Reduce(intersect,list(rownames(rbp),rownames(G4im),rownames(mfe),rownames(EV)))
mat <- as.data.frame(cbind(rbp[pl,"RBP_protect"],G4im[pl,"G4im_protect"],mfe[pl,"Second_structure"],EV[pl,"EV_protect"]))
colnames(mat) <- c("RBP_protect","G4im_protect","Second_structure","EV_protect")
head(mat)

## plot venn
#install.packages("ggvenn") # install via CRAN
library("ggvenn")
#create Venn diagram and display all sets
ggvenn(mat,set_name_size=4,text_size = 5,show_percentage=F, text_color = "white") + ggsci::scale_fill_nejm() + 
  theme(aspect.ratio = 0.8)
#ggsave("tmpVenn.pdf",width = 10,height = 5, limitsize=F)

table(rowSums(mat))
s <- rowSums(mat[,c("RBP_protect","G4im_protect","Second_structure","EV_protect")])

# #mat$RBP_protect_ratio <- as.numeric(mat$RBP_protect)/(s+0.01)
# #mat$RBP_protect_ratio <- mat$RBP_protect_ratio/nrow(mat)
# mat$G4im_protect_ratio <- as.numeric(mat$G4im_protect)/(s+0.01)
# #
# mat$Second_structure_ratio <- as.numeric(mat$Second_structure)/(s+0.01)
# #mat$Second_structure <- mat$Second_structure/nrow(mat)
# mat$EV_protect_ratio <- as.numeric(mat$EV_protect)/(s+0.01)
# mat$Unknown <- s==0
# #mat$EV_protect_ratio <- mat$EV_protect_ratio/nrow(mat)
# #mat <- as.matrix(mat)
# #head(rbp[pl,"RBP_protect"])
# #mat[1:4,]
# mat <- mat[,5:ncol(mat)] # total 6334
# #rowSums(mat)
# mat <- colSums(mat)#/sum(colSums(mat))
mat$Unknown <- s==0
mat[1:3,]


## plot bar
df <- as.data.frame(colSums(mat)/nrow(mat))
colnames(df) <- "value"
df$group <- rownames(df)
ggplot(df, aes(x = group , y = value, fill = group)) +
  geom_bar(width = 1,stat = "summary", color = 1) +
  scale_fill_manual(name="Group",values = colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(18)) + # c(RColorBrewer::brewer.pal(8, "Accent"),RColorBrewer::brewer.pal(9, "Set1"))
  ylim(c(0,1))+
  guides(fill = guide_legend(title = "Group")) +
  theme_bw()+
  theme(plot.title = element_text(size = 24,color="black",hjust = 0.5,face="bold"),
        axis.title = element_text(size = 24,color ="black"), 
        axis.text = element_text(size= 20,color = "black"),
        #panel.grid=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey50",linetype = "dashed"), #size= 1,
        #panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5), #  
        legend.position = "right",#c(.25,.6),
        legend.text = element_text(size= 16),
        legend.title= element_text(size= 16))



## plot pie
library(ggplot2)
df <- as.data.frame(mat)
df$type <- rownames(df)
df <- df[,c(2,1)]
#df <- df[sample(1:nrow(df),0.05*nrow(df),replace = F),]
#str(df)
colnames(df) <- c("group","value")
df$value <- as.integer(df$value)
df$ratio <- df$value/sum(df$value)
#df$group <- factor(df$group)
df$ratio <- as.numeric( round(df$ratio, digits = 4))
#df <- df[order(df$ratio,decreasing = F),]

df2 <- df %>% 
  mutate(csum = rev(cumsum(rev(value))), 
         pos = value/2 + lead(csum, 1),
         pos = if_else(is.na(pos), value/2, pos))

ggplot(df, aes(x = "" , y = value, fill = forcats::fct_inorder(group))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(name="Group",values = colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(18)) + # c(RColorBrewer::brewer.pal(8, "Accent"),RColorBrewer::brewer.pal(9, "Set1"))
  ggrepel::geom_label_repel(data = df2,
                            aes(y = pos, label = paste0(value,": ",100*ratio, "%")),
                            size = 4.5, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "Group")) +
  theme_void()
#ggsave("./domain_pie.pdf",width = 6,height = 5)







# prepare random/full-lenth matrix (deprecated) ----------------------------------------
#setwd("/Users/baopengfei/Desktop/lulab/tmp/projects/motif-RBP-EDA/")
setwd("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/")

#read meta
sample.table <- read.table("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/meta/lulab/FTC/sample_table.txt",sep = "\t",header = T)
sample.table <- reshape2::melt(sample.table,id.var=c("patient_id","group","name"),variable.name="source",value.name = "data_id")
sample.table <- sample.table[sample.table$data_id!="",]
sample.table <- sample.table[sample.table$source %in% c("EV_long"),] # c("EV_small","cf_small")
length(unique(sample.table$patient_id))
positive_samples <- sample.table[sample.table$group=="FTC","data_id"]
negative_samples <- sample.table[sample.table$group=="FTA","data_id"]
negative_samples <- unlist(sapply(strsplit(negative_samples,"_",fixed=T),"[",1))  # trim samples surfix !
positive_samples <- unlist(sapply(strsplit(positive_samples,"_",fixed=T),"[",1))   # trim samples surfix !

#read mat
#count <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/count_matrix/domains_fullLength.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F)
count <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/FTC_long/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/count_matrix/domains_fullLength.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F)
#count <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/count_matrix/domains_long.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F)
#count <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/FTC_long/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/count_matrix/domains_long.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F)
rownames(count) <- count$feature
mat <- count[,2:ncol(count)]
colnames(mat) <- unlist(sapply(strsplit(colnames(mat),"_",fixed=T),"[",1)) 
mat[1:3,1:3]

#diff
samples <- c(positive_samples, negative_samples)
group <- c(rep("positive",length(positive_samples)),rep("negative",length(negative_samples)))
method <- "edger_glmlrt"
norm_method <- "TMM"
res <- diff.v2(mat = mat, samples = samples, group = group,  method = method, norm_method = norm_method, filterType = "small") # patient = patient, small

logcpm <- as.data.frame(res$normMat)
#rownames(logcpm) <- logcpm$gene_id
#colnames(logcpm) <- unlist(sapply(strsplit(colnames(logcpm),"_",fixed=T),"[",1)) 
logcpm[1:3,1:3]
print(table(c(positive_samples,negative_samples) %in% colnames(logcpm))) # should be all TURE !

type.tmp <- "FullLen" #FullLen,localmax
#write.table(cbind(t(logcpm[,positive_samples]),class=positive_samples_class), "/Users/baopengfei/Desktop/lulab/tmp/projects/ML/data_FTC/smallDomain_EV.csv",quote = F,row.names = F,col.names = T,sep = ",")
#write.table(cbind(t(logcpm[,negative_samples]),class=negative_samples_class), "/Users/baopengfei/Desktop/lulab/tmp/projects/ML/data_FTC/smallDomain_cf.csv",quote = F,row.names = F,col.names = T,sep = ",")
#write.table(cbind(gene_id=rownames(logcpm),logcpm[,c(positive_samples,negative_samples)]), paste0("/Share2/home/baopengfei/projects/ML/data_FTC/smallDomain",type.tmp,"_EV.txt"),quote = F,row.names = F,col.names = T,sep = "\t")
write.table(cbind(gene_id=rownames(logcpm),logcpm[,c(positive_samples,negative_samples)]), paste0("/Share2/home/baopengfei/projects/ML/data_FTC/longDomain",type.tmp,"_EV.txt"),quote = F,row.names = F,col.names = T,sep = "\t")

#write.table(positive_samples, paste0("/Share2/home/baopengfei/projects/ML/data_FTC/class_FTC_smallDomain",type.tmp,"_EV.txt"),quote = F,row.names = F,col.names = F,sep = "\t")
#write.table(negative_samples, paste0("/Share2/home/baopengfei/projects/ML/data_FTC/class_FTA_smallDomain",type.tmp,"_EV.txt"),quote = F,row.names = F,col.names = F,sep = "\t")
write.table(positive_samples, paste0("/Share2/home/baopengfei/projects/ML/data_FTC/class_FTC_longDomain",type.tmp,"_EV.txt"),quote = F,row.names = F,col.names = F,sep = "\t")
write.table(negative_samples, paste0("/Share2/home/baopengfei/projects/ML/data_FTC/class_FTA_longDomain",type.tmp,"_EV.txt"),quote = F,row.names = F,col.names = F,sep = "\t")






# evaluate peak/domain pri_miR boundary hist (suppl fig3?) ---------------- ----------------
## get pri-miR & miR hg38 coordinates
##################  miRBase.gff3 和 mirgenedb 注释的pri_RNA区域不同 （mirgenedb更长？），以转录本对齐时出现问题，优先基因组

# #/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008/archive/count_matrix/miRNA_count_matrix.txt
# #filter miR anno. by cf mean abundance 
# tmp <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008/archive/count_matrix/miRNA_count_matrix.txt", header = T, row.names = 1)
# tmp[1:3,1:3]
# 
# miR.anno <- read.table("/BioII/lulab_b/baopengfei/shared_reference/mirbase/hsa.bed")
# miR.anno2 <- miR.anno[miR.anno$V4 %in% rownames(tmp)[rowMedians(as.matrix(tmp))>=5],]
# write.table(miR.anno2,"/BioII/lulab_b/baopengfei/shared_reference/mirbase/hsa_highInPlasma.bed",quote = F,row.names = F,col.names = F,sep = "\t")


#R
setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/")
options(stringsAsFactors = F)
options(bedtools.path = "/BioII/lulab_b/baopengfei/biosoft") # /BioII/lulab_b/baopengfei/anaconda3/envs/exvariance/bin
suppressPackageStartupMessages(library(bedtoolsr))
library(tidyverse)
library(dplyr)

## read overlap bed file
read.mir.intersect <- function(x){
#x <- "output/AGO2_IP/call_domain/domains/5/05_gn_interectMIR.bed"
# x <- "output/GSE133684/call_peak_dedupByPos/domains/b5_p01_gn_interectMIR.bed"
  print(x)
rbp <- data.table::fread(x,data.table = F, header = F,sep = "\t",stringsAsFactors = F)
idx <- colnames(rbp)[ncol(rbp)]

rbp$intersect.ratio <- as.numeric(rbp[[idx]]/(rbp$V3-rbp$V2)) * as.numeric(rbp[[idx]]/(rbp$V15-rbp$V14))

rbp <- rbp[order(rbp$V4,rbp$intersect.ratio,rbp[[idx]],decreasing = T),] # rm dup
rbp <- rbp[!duplicated(rbp$V4),]
rbp <- na.omit(rbp) # na exist in intersect.ratio
# table(rbp[[idx]]==0)
rbp <- rbp[rbp[[idx]]>=10 | rbp$intersect.ratio>=0.01,]

rbp$delta.5 <- rbp$V2 - rbp$V14 # start - start
#hist(rbp$delta.5,breaks = 10,col="grey50",xlim = c(-10,10))
rbp$delta.3 <- rbp$V3 - rbp$V15 # start - start
#hist(rbp$delta.3,breaks = 10,col="grey50",xlim = c(-10,10))
rbp <- rbp[,c(1:6,ncol(rbp)-4,ncol(rbp)-3,ncol(rbp)-2,ncol(rbp)-1,ncol(rbp))]
#plot(x=rbp$delta.5,y=rbp$delta.3)
rbp$path <- x
#rbp$p <- as.numeric(paste0("0.",gsub("_gn_interectMIR.bed","",basename(x),fixed = T)))
return(rbp)
}

#need define best decay and p cutoff statistics/param
#AGO2_IP/call_peak_dedup_bk, GSE71008/call_peak_all
dst <- "GSE71008_NCpool" # GSE71008_NCpool, GSE110381_NCpool, WSQ_SMARTer_NEB_NCpool (NCpool_NEBNext)
smp <- "NCpool"  #"NCpool_NEBNext"
pre <- paste0(dst,"/call_peak_all") # "AGO2_IP/call_peak_dedup_bk", "GSE71008/call_peak_all", "GSE110381/call_peak_all" 
# res1 <- lapply(c(paste0("output/",pre,"/piranha/b5_p01_11RNA_gn_interectMIR.bed"),
#                  paste0("output/",pre,"/clipper/b5_p05_11RNA_gn_interectMIR.bed"),
#                  paste0("output/",pre,"/clam/b5_p005_11RNA_gn_interectMIR.bed"), 
#                  paste0("output/",pre,"/expeak/b5_d50_p1_11RNA_gn_interectMIR.bed") ),
#                read.mir.intersect)
res1 <- lapply(c(paste0("output/",pre,"/piranha_by_sample/b5_p01/intersect/",smp,"_gn_interectMIR.bed"),
                 paste0("output/",pre,"/clipper_by_sample/b5_p05/intersect/",smp,"_gn_interectMIR.bed"),
                 paste0("output/",pre,"/clam_by_sample/b5_p005/intersect/",smp,"_gn_interectMIR.bed"), 
                 paste0("output/",pre,"/expeakCNN_by_sample/b5_d50_p1/intersect/",smp,"_gn_interectMIR.bed") ),
               read.mir.intersect)
#not apply to long seq, no intersect records?

res1.df <- as_tibble(do.call("rbind",res1))
res1.df$method <- "Piranha"
res1.df$method[grepl("clipper",res1.df$path)] <- "CLIPper"
res1.df$method[grepl("clam",res1.df$path)] <- "CLAM"
# res1.df$method[grepl("localmax",res1.df$path)] <- "LocalMax"
res1.df$method[grepl("expeak",res1.df$path)] <- "exPeak"
table(res1.df$method)
# res1.df$method[grepl("EMgini",res1.df$path)] <- "LocalMaxExtGini"
#x <- gsub("_11RNA_gn_interectMIR.bed","",basename(res1.df$path),fixed = T)
#x <- gsub("b5_d50_p|b5_p","",x,perl = T)
x <- "05"
res1.df$p <- as.numeric(paste0("0.",x))
table(res1.df$p )
# 0 0.05  0.1  0.2  0.5  0.9 
# 351  405  405  407  410  414
res1.df2 <- res1.df %>% 
  dplyr::filter(p==0.05) %>% # 0.01
  tidyr::pivot_longer(cols = delta.3:delta.5, names_to = "direction",values_to = "distance")
res1.df2$method <- factor(res1.df2$method,levels = c("Piranha","CLIPper","CLAM","exPeak")) # 

res1.df2$distance[res1.df2$distance< -10] <- -10
res1.df2$distance[res1.df2$distance>10] <- 10
#table(res1.df2$direction=="delta.3")
res1.df2$direction <- ifelse(res1.df2$direction=="delta.3","3'","5'")
res1.df2$direction <- factor(res1.df2$direction,levels = c("5'","3'"))
table(res1.df2$direction)


## plot histplot: boundary postion 
library(plyr); library(dplyr)
mu <- ddply(res1.df2, "method", plyr::summarise, grp.mean=mean(distance))
library(dplyr)
mu <- res1.df2 %>% 
  dplyr::group_by(direction,method) %>% 
  dplyr::summarise(grp.mean=mean(distance), grp.sum=dplyr::n(), grp.zero=sum(distance==0), ratio=grp.zero/grp.sum )
mu$method <- factor(mu$method,levels = c("Piranha","CLIPper","CLAM","exPeak"))
summary(mu$ratio)

#head(res1.df2)
library(ggplot2)
ggplot(res1.df2, aes(x=distance, color=method, fill=method)) +
  geom_histogram( position="identity", alpha=0.9, color="black")+
  # geom_density(alpha=0.1)+
  # geom_vline(data=mu, aes(xintercept=grp.mean, color=method),linetype="dashed")+
  # scale_color_manual(values=c("firebrick", "steelblue4", "grey30"))+
  # scale_fill_manual(values=c("firebrick", "steelblue4", "grey30"))+
  geom_vline(xintercept = 0,color="grey30",linetype="dashed")+
  # geom_label(data=mu, aes(label=ratio, color=method, fill=method )) +
  #annotate(x = 0.5, y=0.5)+
  # geom_text()+
  # ggsci::scale_fill_nejm()+
  # ggsci::scale_color_nejm()+
  scale_fill_d3_adaptive()+
  # scale_color_d3_adaptive()+
  # scale_fill_manual(values = pal_nejm_adaptive()(4))+
  # scale_color_manual(values = pal_nejm_adaptive()(4))+
  labs(title="",x="MIR end relative position (5'/3')", y = "Predicted position frequency")+
  facet_grid(method~direction,scales = "free")+
  # ylim(c(0,1))+xlim(c(-10,10))+
  theme_classic(base_size=12) + 
  theme(#axis.ticks.x=element_blank(),  
    #strip.text.y = element_blank(),
    #strip.text.x = element_text(face="bold",family="arial",size=20),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size = 20,hjust = 0.5,vjust = 0.5), # ,angle = 90
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size=20),
    strip.text = element_text(size = 20),
    legend.position = c(0.5,0.15),# 
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16))
ggsave("mir_boundary.pdf",width = 9,height = 9) # ,height = 14
mu
# Piranha:0.218,0.200
# CLIPper:0.278,0.324
# CLAM:   0.033,0.052
# exPeak: 0.520,0.478



## plot boxplot: domain/mir overlap ratio
library(ggplot2)
library(tidyverse)
library(dplyr)
pre <- paste0(dst,"/call_peak_all") # "AGO2_IP/call_peak_dedup_bk", "GSE71008/call_peak_all", "GSE110381/call_peak_all" 
  
in1 <- Sys.glob( paste0("output/",pre,"/piranha_by_sample/b5_p01/intersect/",smp,"_gn_interectMIR.bed") )
in2 <- Sys.glob( paste0("output/",pre,"/clipper_by_sample/b5_p05/intersect/",smp,"_gn_interectMIR.bed") )
in3 <- Sys.glob( paste0("output/",pre,"/clam_by_sample/b5_p005/intersect/",smp,"_gn_interectMIR.bed") )
# in4 <- Sys.glob(paste0("output/",pre,"/domains_localmax/b5_d05_p*_gn_interectMIR.bed"))
in5 <- Sys.glob( paste0("output/",pre,"/expeak_by_sample/b5_d50_p1/intersect/",smp,"_gn_interectMIR.bed") )
res1 <- lapply(c(in1,in2,in3,in5),
               read.mir.intersect)
#not apply to long seq, need remove file with zero records

res1.df <- as_tibble(do.call("rbind",res1))
res1.df$method <- "Piranha"
res1.df$method[grepl("clipper",res1.df$path)] <- "CLIPper"
res1.df$method[grepl("clam",res1.df$path)] <- "CLAM"
# res1.df$method[grepl("localmax",res1.df$path)] <- "LocalMax"
res1.df$method[grepl("expeak",res1.df$path)] <- "exPeak"
# x <- gsub("_11RNA_gn_interectMIR.bed","",basename(res1.df$path),fixed = T)
# x <- gsub("b5_d50_p||b5_p","",x,perl = T)
x <- "05"
res1.df$p <- as.numeric(paste0("0.",x))
table(res1.df$p )
res1.df <- na.omit(res1.df)

res1.df2 <- res1.df %>% 
  # filter(p==0.01) %>%
  pivot_longer(cols = delta.3:delta.5, names_to = "direction",values_to = "distance")
res1.df2$method <- factor(res1.df2$method,levels =c("Piranha","CLIPper","CLAM","exPeak")) # LocalMax
res1.df2$p <- factor(res1.df2$p)
res1.df2$distance[res1.df2$distance< -10] <- -10
res1.df2$distance[res1.df2$distance>10] <- 10
#table(res1.df2$direction=="delta.3")
res1.df2$direction <- ifelse(res1.df2$direction=="delta.3","3'","5'")
res1.df2$direction <- factor(res1.df2$direction,levels = c("5'","3'"))
table(res1.df2$direction)

res1.df3 <- res1.df2 %>%
  dplyr::filter(direction=="5'") #%>%  #aviod double counting
p1 <- ggplot(res1.df3, aes(x=p, y=intersect.ratio,  fill=method)) +
  geom_boxplot(color="black",outlier.alpha = 0.3)+
  # geom_violin()+
  # geom_histogram(aes(y=..density..), position="identity", alpha=0.5)+
  # geom_density(alpha=0.1)+
  # geom_vline(data=mu, aes(xintercept=grp.mean, color=method),linetype="dashed")+
  # scale_color_manual(values=c("firebrick", "steelblue4", "grey30"))+
  # scale_fill_manual(values=c("firebrick", "steelblue4", "grey30"))+
  # geom_vline(xintercept = 0,color="grey30",linetype="dashed")+
  # geom_label(data=mu, aes(label=ratio, color=method, fill=method )) +
  #annotate(x = 0.5, y=0.5)+
  # geom_text(data=mu, nudge_x = 0.5, nudge_y = 0.5, aes(text=ratio, color=method))+
  scale_fill_d3_adaptive()+
  # scale_color_d3_adaptive()+
  # scale_x_discrete(label=c("Domain Ratio","MIR Ratio"))+
  labs(title="",x="P value", y = "Overlap Ratio")+
  # facet_grid(method~.)+
  #ylim(c(0,1))+
  # xlim(c(0,1))+
  theme_classic(base_size=12) +
  theme(#axis.ticks.x=element_blank(),
    #strip.text.y = element_blank(),
    #strip.text.x = element_text(face="bold",family="arial",size=20),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size = 20,hjust = 0.5,vjust = 0.5), # ,angle = 90
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size=20),
    strip.text = element_text(size = 20),
    legend.position = "right",#
    legend.justification = "top",
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16))


res1.df4 <- res1.df3 %>%
  # pivot_longer(cols = c(intersect.ratio), names_to = "Type",values_to = "Ratio") %>%
  # filter(direction=="5'") %>%  #aviod double counting
  dplyr::group_by(method,p) %>%
  dplyr::mutate(sd=sd(intersect.ratio))
p2 <- ggplot(res1.df4, aes(x=p, y=sd, color=method, fill=method))+
  geom_bar(stat="identity",position = position_dodge(width=0.7),color="black", width=0.5)+
  # geom_violin()+
  # geom_histogram(aes(y=..density..), position="identity", alpha=0.5)+
  # geom_density(alpha=0.1)+
  # geom_vline(data=mu, aes(xintercept=grp.mean, color=method),linetype="dashed")+
  # scale_color_manual(values=c("firebrick", "steelblue4", "grey30"))+
  # scale_fill_manual(values=c("firebrick", "steelblue4", "grey30"))+
  # geom_vline(xintercept = 0,color="grey30",linetype="dashed")+
  # geom_label(data=mu, aes(label=ratio, color=method, fill=method )) +
  #annotate(x = 0.5, y=0.5)+
  # geom_text(data=mu, nudge_x = 0.5, nudge_y = 0.5, aes(text=ratio, color=method))+
  scale_fill_d3_adaptive()+
  # scale_color_d3_adaptive()+
  # scale_x_discrete(label=c("Domain Ratio","MIR Ratio"))+
  labs(title="",x="P value", y = "Overlap Ratio Sd")+ # Overlap Ratio
  # facet_grid(method~.)+
  #ylim(c(0,1))+
  # xlim(c(0,1))+
  theme_classic(base_size=12) +
  theme(#axis.ticks.x=element_blank(),
    #strip.text.y = element_blank(),
    #strip.text.x = element_text(face="bold",family="arial",size=20),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size = 20,hjust = 0.5,vjust = 0.5), # ,angle = 90
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size=20),
    strip.text = element_text(size = 20),
    legend.position = "right",#
    legend.justification = "top",
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16))
cowplot::plot_grid(plotlist = list(p1,p2), nrow = 1, align = "hv", axis = "b")
ggsave(filename = "overlap_miR.pdf", width = 11, height = 8)#p1+p2+patchwork::guide_area+patchwork::plot_layout()



# evaluate peak/domain AGO2-bed miR seed ratio (deprecated) ---------------- ----------------
# preapare AGO2 intersect domain bed
# cd /BioII/lulab_b/baopengfei/projects/WCHSU-FTC
# in="output/AGO2_IP/call_domain/domains/5/50"
# #exSeek-dev/genome/hg38/bed/long_RNA-tRNA-long_DNA_noRepeats.bed
# {{
#   grep -v '^chr' ${in}.bed | exSeek-dev/bin/tbed2gbed <(cat exSeek-dev/genome/hg38/bed/{long_RNA,tRNA,pri_miRNA,piRNA}_newTxID.bed) /dev/stdin /dev/stdout
#   awk 'BEGIN{{OFS="\t";FS="\t"}}/^chr/{{print $1,$2,$3,$4,$5,$6,0,0,0,1,$3-$2,0}}' ${in}.bed
# }} | bedtools sort \
# > ${in}_gn.bed

# bedtools intersect -wao -s \
# -a output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/domains_gn.bed \
# -b /BioII/lulab_b/baopengfei/shared_reference/RBP/splitByRBP/AGO2.bed > \
# output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/domains_gn_interectAGO2.bed
# 
# bedtools intersect -wao -s \
# -a output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains/20/05_gn.bed \
# -b /BioII/lulab_b/baopengfei/shared_reference/RBP/splitByRBP/AGO2.bed > \
# output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains/20/05_gn_interectAGO2.bed

# mkdir -p output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_by_sample/20/05/gn
# mkdir -p output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax_by_sample/gn
# mkdir -p output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_by_sample/20/05/interectAGO2
# mkdir -p output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax_by_sample/interectAGO2
# for i in `cat exSeek-dev/data/FTC_small/sample_ids.txt`
# do
# echo $i
# in="output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_by_sample/20/05/${i}"
# {{
#   grep -v '^chr' ${in}.bed | exSeek-dev/bin/tbed2gbed <(cat exSeek-dev/genome/hg38/bed/{long_RNA,tRNA,pri_miRNA,piRNA}.bed) /dev/stdin /dev/stdout
#   awk 'BEGIN{{OFS="\t";FS="\t"}}/^chr/{{print $1,$2,$3,$4,$5,$6,0,0,0,1,$3-$2,0}}' ${in}.bed
# }} | bedtools sort \
# > ${in}_gn.bed
# bedtools intersect -wao -s \
# -a ${in}_gn.bed \
# -b /BioII/lulab_b/baopengfei/shared_reference/RBP/splitByRBP/AGO2.bed > \
# output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_by_sample/20/05/interectAGO2/${i}.bed
# 
# in="output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax_by_sample/${i}"
# {{
#   grep -v '^chr' ${in}.bed | exSeek-dev/bin/tbed2gbed <(cat exSeek-dev/genome/hg38/bed/{long_RNA,tRNA,pri_miRNA,piRNA}.bed) /dev/stdin /dev/stdout
#   awk 'BEGIN{{OFS="\t";FS="\t"}}/^chr/{{print $1,$2,$3,$4,$5,$6,0,0,0,1,$3-$2,0}}' ${in}.bed
# }} | bedtools sort \
# > ${in}_gn.bed
# bedtools intersect -wao -s \
# -a ${in}_gn.bed \
# -b /BioII/lulab_b/baopengfei/shared_reference/RBP/splitByRBP/AGO2.bed > \
# output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax_by_sample/interectAGO2/${i}.bed
# done

## sub-sample
# dir.create("output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/subset")
# infile="output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/domains_gn.bed"
# tmp <- read.table(infile,header = F)
# i <- seq(0.02,0.98,0.02) # 50 
# n <- unique(as.integer(nrow(tmp)*i))
# for (ni in n){
# print(ni)
# bedtoolsr::bt.sample(output = paste0("output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/subset/domains_gn_",ni,".bed"), i = infile, seed = 1, n = ni)
# bedtoolsr::bt.intersect(a = paste0("output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/subset/domains_gn_",ni,".bed"),
#                         b = "/BioII/lulab_b/baopengfei/shared_reference/RBP/splitByRBP/AGO2.bed",wao = T,s = T,
#                         output = paste0("output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/subset/domains_gn_",ni,"_intersectAGO2.bed") )
# }
# infile <- "output/AGO2_IP/call_domain/domains/5/50_gn"
# bedtoolsr::bt.intersect(a = paste0(infile,".bed"),
#                         b = "/BioII/lulab_b/baopengfei/shared_reference/RBP/splitByRBP/AGO2.bed",wao = T,s = T,
#                         output = paste0(infile,"_intersectAGO2.bed")) 

## get rev complement seq
seq_rev <- function(char) {
  alphabets <- strsplit(char, split = "")[[1]]
  return(rev(alphabets))
}

seq_compl <- function(seq) {
  # Check if there's "T" in the sequence
  RNA <- Reduce(`|`, seq == "U")
  cmplvec <- sapply(seq, function(base) {
    # This makes DNA the default
    # As long as there's no U, the sequence is treated as DNA
    if (RNA) {
      switch(base, "A" = "U", "C" = "G", "G" = "C", "U" = "A")
    } else {
      switch(base, "A" = "T", "C" = "G", "G" = "C", "T" = "A")
    }
  })
  return(paste(cmplvec, collapse = ""))
}

## get seed sequences
top20 <- read.table("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/lulab/miR_diff_FTA_EVvsCF_unpair.txt",header = T,row.names = 1,sep = "\t",stringsAsFactors = F)
top20 <- top20[order(top20$baseMean,decreasing = T),]
top20 <- rownames(top20)[1:nrow(top20)]
seed <- data.table::fread("/BioII/lulab_b/baopengfei/shared_reference/targetScan/miR_Family_Info.txt",header = T,sep = "\t",stringsAsFactors = F)
seed <- as_tibble(seed)
#seed <- seed[grepl("hsa",seed$`Mature sequence`),]
#seed <- seed %>% 
#  filter(grepl("hsa",`MiRBase ID`))
#table(nchar(seed$`Seed+m8`)) # all 7-mer
seed$`Seed+m8` <- gsub("U","T",seed$`Seed+m8`)
#seed2 <- seed[sample(1:nrow(seed),0.1*nrow(seed),replace=F),]
#sequences <- unique(seed$`Seed+m8`)
sequences <- (unique(seed[seed$`MiRBase ID` %in% top20,"Seed+m8"])) #336 only choose top 20 abundant miR seed sequences 
sequences <- sequences$`Seed+m8`
sequences.rev <- sapply(sequences,function(x) seq_compl(seq_rev(x)) )  # add reverse complement
sequences <- unique(c(sequences,sequences.rev)) # sequences.rev
seq.num <- length(sequences)
sequences <- paste((sequences),collapse = "|")
#4^7=16384  
#336/16384=0.02
# (40-7)*366/4^7 # =0.15



## filter domains: AGO2 intersect > 7nt
setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/")
options(stringsAsFactors = F)
options(bedtools.path = "/BioII/lulab_b/baopengfei/biosoft") # /BioII/lulab_b/baopengfei/anaconda3/envs/exvariance/bin
suppressPackageStartupMessages(library(bedtoolsr))
library(tidyverse)
library(dplyr)

sum.match.freq <- function(x){
# x <- "output/AGO2_IP/call_domain/domains/5/50_gn_intersectAGO2.bed"
#x <- "output/AGO2_IP/call_domain/domains/5/50_gn.bed"
  # x <- in1[1]
print(x)
# #x <- "output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax_by_sample/interectAGO2/FTA-12_1.bed"
rbp <- data.table::fread(paste0(x),data.table = F, header = F,sep = "\t",stringsAsFactors = F)

bed <- as_tibble(rbp)
bed <- bed[,1:6]
colnames(bed) <- c("chr","start","end","peak","score","strand")
summary(bed$end-bed$start)
bed <- bed %>% 
  filter(end-start<=200 , end-start>=10)

bedtoolsr::bt.getfasta(fo = paste0("tmp/",basename(x),"tmp.fa"),nameOnly = T, fi = "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/fasta/genome.fa",bed = bed, s = T)
fasta <- rtracklayer::import(paste0("tmp/",basename(x),"tmp.fa"),format = "fasta")
length(fasta)
#fasta2 <- as.data.frame(fasta[sample(1:length(fasta),0.1*length(fasta),replace=F)])
fasta <- as.data.frame(fasta)

tmp <- sum( grepl(sequences,fasta$x,perl=T) )/nrow(fasta)
tmp <- data.frame(ratio=tmp,sample=basename(x),median.len=median(bed$end-bed$start),number=nrow(bed),path=x)
return(tmp)
}

sum.match.freq.intersect <- function(x){
  # x <- "output/GSE133684/call_peak_dedupByPos/domains/b5_p01_gn_interectAGO2.bed"
  print(x)
  # #x <- "output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax_by_sample/interectAGO2/FTA-12_1.bed"
  rbp <- data.table::fread(paste0(x),data.table = F, header = F,sep = "\t",stringsAsFactors = F)
  idx <- colnames(rbp)[ncol(rbp)]
  #rbp$domain.ratio <- as.numeric(rbp[[idx]]/(rbp$V3-rbp$V2))
  rbp$intersect.ratio <- as.numeric(rbp[[idx]]/(rbp$V3-rbp$V2)) * as.numeric(rbp[[idx]]/(rbp$V15-rbp$V14))
  
  rbp <- rbp[order(rbp$V4,rbp$intersect.ratio,rbp[[idx]],decreasing = T),] # rm dup
  rbp <- rbp[!duplicated(rbp$V4),]
  
  bed <- as_tibble(rbp[rbp$intersect.ratio>=0.01 | rbp[[idx]] >=10,])
  bed <- bed[,1:6]
  colnames(bed) <- c("chr","start","end","peak","score","strand")
  #summary(bed$end-bed$start)
  bed <- bed %>% 
    filter(end-start<=500 , end-start>=10)
  
  bedtoolsr::bt.getfasta(fo = paste0("tmp/",basename(x),"tmp.fa"),nameOnly = T, fi = "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/fasta/genome.fa",bed = bed, s = T)
  fasta <- rtracklayer::import(paste0("tmp/",basename(x),"tmp.fa"),format = "fasta")
  length(fasta)
  #fasta2 <- as.data.frame(fasta[sample(1:length(fasta),0.1*length(fasta),replace=F)])
  fasta <- as.data.frame(fasta)
  
  tmp <- sum( grepl(sequences,fasta$x,perl=T) )/nrow(fasta)
  tmp <- data.frame(ratio=tmp,sample=basename(x),median.len=median(bed$end-bed$start),number=nrow(bed),path=x)
  return(tmp)
}

pre <- "GSE71008/call_domain_withRepeats_all" # "AGO2_IP/call_peak_dedup_bk", "GSE71008/call_domain_withRepeats_all", "GSE133684/call_peak_dedupByPos"  
in1 <- Sys.glob(paste0("output/",pre,"/domains/b5_p*_gn_interectAGO2.bed")) # AGO2 only need gn.bed
in2 <- Sys.glob(paste0("output/",pre,"/domains_clipper/b5_p*_gn_interectAGO2.bed"))
in3 <- Sys.glob(paste0("output/",pre,"/domains_localmax/b5_d05_p*_gn_interectAGO2.bed"))
in4 <- Sys.glob(paste0("output/",pre,"/domains_localmax_EM/b5_d05_p*_gn_interectAGO2.bed"))

#seems only small intersect with AGO2 is reasonable, long may produce negative value
#non-CLIP seq need use sum.match.freq.intersect, not sum.match.freq
library(parallel)
res2.df <- do.call(rbind,mclapply(c(in1,in2,in3,in4), sum.match.freq.intersect, mc.cores = 1) )  # sum.match.freq (for AGO2), sum.match.freq.intersect (for small/long)
#seem not applicable to long, get negative ratio

res2.df$method <- "Piranha"
res2.df$method[grepl("localmax",res2.df$path)] <- "LocalMax"
res2.df$method[grepl("EM",res2.df$path)] <- "LocalMaxExt"
res2.df$method[grepl("clipper",res2.df$path)] <- "CLIPper"
res2.df$method <- factor(res2.df$method,levels = c("LocalMaxExt","LocalMax","CLIPper","Piranha"))

#res2.df$p <- unlist(sapply(strsplit(res2.df$sample,"_"),"[",1))
res2.df$p <- gsub("b5_d05_p|b5_p|_gn|_interectAGO2|\\.bed","",res2.df$sample,perl = T)
res2.df$p.value <- as.character(paste0("0.",res2.df$p))
unique(res2.df$p)
table(res2.df$p,res2.df$method)
res2.df$observe.num <- res2.df$number*res2.df$ratio
res2.df$bg.num <- res2.df$observe.num*(res2.df$median.len-7)*seq.num/4^7  # 336, 660 
res2.df$adj.num <- res2.df$observe.num-res2.df$bg.num
res2.df$adj.ratio <- res2.df$adj.num / res2.df$number
#summary(res2.df$true.num )
#res2.df$factor <- (res2.df$median.len-7)*seq.num/4^7 # =0.15
#res2.df$adj.ratio <- res2.df$ratio/res2.df$factor
#table(res2.df$method)
#str(res2.df)
#res2.df <- res2.df[order(res2.df$ratio,de),]
library(ggplot2)
scatter <- ggplot(res2.df, aes(x=number,y=adj.ratio,fill=method,color=method))+
  geom_line(size=1)+
  ggrepel::geom_label_repel(aes(x=number,y=adj.ratio, #fill=method,   
      label=as.numeric(paste0("0.",p))),alpha=1,fill="white", fontface="bold", 
      force = 0.4, force_pull=0.4,
      max.overlaps=10,
      box.padding=unit(0.15, "lines"), point.padding=unit(0.15, "lines"),
      segment.colour = "grey50"
      )+
  # ggrepel::geom_text_repel(aes (label = rownames (mtcars)))+
  scale_fill_d3_adaptive()+
  scale_color_d3_adaptive()+
  labs(title="AGO2 miR seed enrich ratio",x="Significant Domain Number", y = "Adjusted Ratio")+
  ylim(c(0,1))+
  # xlim(c(0,400))+
  theme_classic(base_size=12) + 
  theme(#axis.ticks.x=element_blank(),  
    #strip.text.y = element_blank(),
    #strip.text.x = element_text(face="bold",family="arial",size=20),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size = 20,hjust = 0.5,vjust = 0.5), # ,angle = 90
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size=20),
    legend.position = "right",# c(0.14,0.8),#
    legend.justification = "bottom", #c(0, 0),
    # legend.direction = "horizontal",
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16))
#ggsave(plot = scatter,filename = "seed_ratio_scatter.pdf",width = 9,height = 7)

#bar
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}
bar <- ggplot(res2.df, aes(x=p.value, y=adj.ratio, color=method, fill=method)) +
  geom_bar( stat = "identity", position=position_dodge(width=0.7), width=0.5)+
  scale_fill_d3_adaptive()+
  scale_color_d3_adaptive()+
  scale_x_discrete(label=fancy_scientific)+
  labs(title="",x="", y = "")+
  # facet_grid(method~.)+
  #ylim(c(0,1))+
  # xlim(c(0,1))+
  theme_classic(base_size=12) + 
  theme(#axis.ticks.x=element_blank(),  
    #strip.text.y = element_blank(),
    #strip.text.x = element_text(face="bold",family="arial",size=20),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size = 16,hjust = 1,vjust = 0.5,angle = 90), # 
    axis.text.y = element_text(size = 16),
    plot.title = element_text(size=20),
    strip.text = element_text(size = 20),
    legend.position = "right",# c(0.85,0.7),#
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16))
p <- scatter + patchwork::inset_element(bar,left=0.25,right=1.35,top=1.0,bottom=0.55)
ggsave(plot = p, filename = "seed_ratio.pdf",width = 9,height = 7)









# evaluate peak/domain AGO2-bed miR seed ratio (suppl fig3?) ----------------

## get seed sequences
top20 <- read.table("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/lulab/miR_diff_FTA_EVvsCF_unpair.txt",header = T,row.names = 1,sep = "\t",stringsAsFactors = F)
#top20 <- read.table("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/GSE71008_diff/exPeak_smallDomain_diff_CRCvsNC.diff",header = T,row.names = 1,sep = "\t",stringsAsFactors = F)
top20 <- top20[order(top20$baseMean,decreasing = T),]
top20 <- rownames(top20)[1:100] # nrow(top20) # has great impact in results !!!
seed <- data.table::fread("/BioII/lulab_b/baopengfei/shared_reference/targetScan/miR_Family_Info.txt",header = T,sep = "\t",stringsAsFactors = F)
seed <- as_tibble(seed)
#seed <- seed[grepl("hsa",seed$`Mature sequence`),]
#seed <- seed %>% 
#  filter(grepl("hsa",`MiRBase ID`))
#table(nchar(seed$`Seed+m8`)) # all 7-mer
seed$`Seed+m8` <- gsub("U","T",seed$`Seed+m8`)
#seed2 <- seed[sample(1:nrow(seed),0.1*nrow(seed),replace=F),]
#sequences <- unique(seed$`Seed+m8`)
sequences <- (unique(seed[seed$`MiRBase ID` %in% top20,"Seed+m8"])) #336 only choose top 20 abundant miR seed sequences 
sequences <- sequences$`Seed+m8`
sequences.rev <- sapply(sequences,function(x) seq_compl(seq_rev(x)) )  # add reverse complement, also reverse direction 
sequences <- sequences.rev # unique(c(sequences,sequences.rev)) # sequences.rev
seq.num <- length(sequences)
sequences <- paste((sequences),collapse = "|")
#4^7=16384  
#336/16384=0.02
# (40-7)*366/4^7 # =0.15


## filter domains (AGO2 intersect > 7nt, filter top 1000 ? )
setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/")
options(stringsAsFactors = F)
options(bedtools.path = "/BioII/lulab_b/baopengfei/biosoft") # /BioII/lulab_b/baopengfei/anaconda3/envs/exvariance/bin
suppressPackageStartupMessages(library(bedtoolsr))
library(tidyverse)
library(dplyr)
conflict_prefer("filter", "dplyr")

sum.match.freq.intersect <- function(x){
  # x <- "output/GSE71008_NCpool/call_peak_all/piranha_by_sample/b5_p01/intersect/NCpool_gn_interectAGO2.bed"
  # x <- "output/GSE71008_NCpool/call_peak_all/clam_by_sample/b5_p005/intersect/NCpool_gn_interectAGO2.bed"
  # x <- "output/GSE71008_NCpool/call_peak_all/clipper_by_sample/b5_p05/intersect/NCpool_gn_interectAGO2.bed"
  # x <- "output/GSE71008_NCpool/call_peak_all/expeakCNN_by_sample/b5_d50_p1/intersect/NCpool_gn_interectAGO2.bed"
  x0 <- gsub("intersect/|_gn_interectAGO2","",x)
  print(x)
  # #x <- "output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax_by_sample/interectAGO2/FTA-12_1.bed"
  rbp <- data.table::fread(paste0(x),data.table = F, header = F,sep = "\t",stringsAsFactors = F)
  idx <- colnames(rbp)[ncol(rbp)]
  #rbp$domain.ratio <- as.numeric(rbp[[idx]]/(rbp$V3-rbp$V2))
  rbp$intersect.ratio <- as.numeric(rbp[[idx]]/(rbp$V3-rbp$V2)) * as.numeric(rbp[[idx]]/(rbp$V15-rbp$V14))
  
  rbp <- rbp[order(rbp$V4,rbp$intersect.ratio,rbp[[idx]],decreasing = T),] # rm dup
  rbp <- rbp[!duplicated(rbp$V4),]
  
  bed <- as_tibble(rbp[rbp$intersect.ratio>=0.01 | rbp[[idx]] >= 6,])
  bed <- bed[,1:6]
  colnames(bed) <- c("chr","start","end","peak","score","strand")
  #summary(bed$end-bed$start)
  bed <- bed %>% 
    filter(end-start<=200 , end-start>=10)
  
  bedtoolsr::bt.getfasta(fo = paste0(x,".tmp.fa"),nameOnly = T, fi = "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/fasta/genome.fa",bed = bed, s = T)
  fasta <- rtracklayer::import(paste0(x,".tmp.fa"),format = "fasta")
  #length(fasta)
  #fasta2 <- as.data.frame(fasta[sample(1:length(fasta),0.1*length(fasta),replace=F)])
  fasta <- as.data.frame(fasta)
  fasta$id <- rownames(fasta)
  fasta$id <- gsub("\\(-\\)|\\(\\+\\)","",fasta$id)
  
  if (grepl("piranha",x)){
    method.lab <- "piranha"
    tmp0 <- data.table::fread(paste0(x0),data.table = F, header = F,sep = "\t",stringsAsFactors = F)
    fasta$pval <- tmp0$V7[match(fasta$id,tmp0$V4)]
  } else if (grepl("clipper",x)){
    method.lab <- "clipper"
    tmp0 <- data.table::fread(paste0(x0),data.table = F, header = F,sep = "\t",stringsAsFactors = F)
    fasta$pval <- tmp0$V5[match(fasta$id,tmp0$V4)]
  } else if (grepl("clam",x)){
    method.lab <- "clam"
    tmp0 <- as.data.frame( data.table::fread(paste0(x0),data.table = F, header = F,sep = "\t",stringsAsFactors = F) )
    #12975:19:2.715e-03,12975:22.8040000573:2.673e-03
    tmp0$p1 <- unlist(sapply(strsplit(tmp0$V4,","),"[",1))
    tmp0$p1 <- as.numeric(unlist(sapply(strsplit(tmp0$p1,":"),"[",3)))
    tmp0$p2 <- unlist(sapply(strsplit(tmp0$V4,","),"[",2))
    tmp0$p2 <- as.numeric(unlist(sapply(strsplit(tmp0$p2,":"),"[",3)))
    tmp0$p <- apply(tmp0[,c("p1","p2")],1,min,na.rm=T)
    fasta$pval <- tmp0$p[match(fasta$id,tmp0$V4)]
    #summary(fasta$pval)
  } else if (grepl("expeak",x)){
    method.lab <- "expeak"
    x0 <- gsub("expeakCNN","expeak",x0)
    tmp0 <- data.table::fread(paste0(x0),data.table = F, header = F,sep = "\t",stringsAsFactors = F) # 
    fasta$pval <- tmp0$V10[match(fasta$id,tmp0$V4)]
  }
  
  pval_seq2 <- unique(c(fasta$pval)) #c(-1,pval_seq)
  # pval_seq2 <- c(0.5,0.05,10^seq(0, -20, by = -1)) # NA exist
  if(length(pval_seq2)>=20){
    pval_seq2 <- pval_seq2[sample(1:length(pval_seq2),20,replace = F)]
  }
  pval_seq2 <- unique(c(c(1,0.1,0.05,0.01),pval_seq2)) # max(domain0$pval)+1
  
  tmp.list <- list()
  for (i in 1:length(pval_seq2)){
    #i <- 1
    cutoff <- pval_seq2[i]
    
    fasta.tmp <- fasta[fasta$pval<=cutoff,]
    #bed[1:3,]
    bed.tmp <- bed[bed$peak %in% fasta.tmp$id,] 
    tmp.list[[i]] <- data.frame(ratio=sum( grepl(sequences,fasta.tmp$x,perl=T) )/nrow(fasta.tmp),
                               sample=x,median.len=median(bed.tmp$end-bed.tmp$start),number=nrow(bed.tmp),path=x,pvalue=cutoff)
  }
  tmp.df <- as.data.frame(do.call(rbind,tmp.list))
  file.remove(paste0(x,".tmp.fa"))
  return(tmp.df)
}

pre <- "GSE71008_NCpool/call_peak_all" # "AGO2_IP/call_peak_dedup_bk", "GSE71008/call_domain_withRepeats_all", "GSE133684/call_peak_dedupByPos"  
in1 <- Sys.glob(paste0("output/",pre,"/piranha_by_sample/b5_p01/intersect/*_gn_interectAGO2.bed")) # AGO2 only need gn.bed
in2 <- Sys.glob(paste0("output/",pre,"/clipper_by_sample/b5_p05/intersect/*_gn_interectAGO2.bed"))
in3 <- Sys.glob(paste0("output/",pre,"/clam_by_sample/b5_p005/intersect/*_gn_interectAGO2.bed"))
in4 <- Sys.glob(paste0("output/",pre,"/expeakCNN_by_sample/b5_d50_p1/intersect/*_gn_interectAGO2.bed"))

#seems only small intersect with AGO2 is reasonable, long may produce negative value
#non-CLIP seq need use sum.match.freq.intersect, not sum.match.freq
library(parallel)
res2.df <- do.call(rbind,mclapply(c(in1,in2,in3,in4), sum.match.freq.intersect, mc.cores = 1) )  # sum.match.freq (for AGO2), sum.match.freq.intersect (for small/long)
#seem not applicable to long, get negative ratio
#set cores to 1 !!! or will get error, na exist

res2.df$method <- "Piranha"
res2.df$method[grepl("expeak",res2.df$path)] <- "exPeak"
res2.df$method[grepl("clam",res2.df$path)] <- "CLAM"
res2.df$method[grepl("clipper",res2.df$path)] <- "CLIPper"
res2.df$method <- factor(res2.df$method,levels = c("Piranha","CLIPper","CLAM","exPeak"))
table(res2.df$method)
#res2.df$p <- unlist(sapply(strsplit(res2.df$sample,"_"),"[",1))
#res2.df$p <- gsub("b5_d05_p|b5_p|_gn|_interectAGO2|\\.bed","",res2.df$sample,perl = T)
# table(res2.df$pvalue)
res2.df$p.value <- format(res2.df$pvalue,digits = 3)
res2.df$label <- as.numeric(res2.df$p.value)
res2.df$label[!(res2.df$label %in% c(1,0.1,0.05,0.01))] <- ""
#table(res2.df$label %in% c(1,0.1,0.05,0.01))

#table(res2.df$p,res2.df$method)
res2.df$observe.num <- res2.df$number*res2.df$ratio
res2.df$bg.num <- res2.df$observe.num*(res2.df$median.len-7)*seq.num/4^7  # 336, 660 
res2.df$adj.num <- res2.df$observe.num-res2.df$bg.num
res2.df$adj.ratio <- res2.df$adj.num / res2.df$number
#summary(res2.df$true.num )
#res2.df$factor <- (res2.df$median.len-7)*seq.num/4^7 # =0.15
#res2.df$adj.ratio <- res2.df$ratio/res2.df$factor
#table(res2.df$method)
#str(res2.df)
#res2.df <- res2.df[order(res2.df$ratio,de),]
library(ggplot2)
#res2.df <- na.omit(res2.df)
ggplot(res2.df, aes(x=number,y=adj.ratio,fill=method,color=method))+
  geom_line(size=1)+
  ggrepel::geom_label_repel(aes(x=number,y=adj.ratio, #fill=method,   
                                label=as.character(label)),alpha=1,fill="white", fontface="bold", 
                            force = 0.4, force_pull=0.4, 
                            max.overlaps=10,
                            box.padding=unit(0.15, "lines"), point.padding=unit(0.15, "lines"),
                            segment.colour = "grey50"
  )+
  # ggrepel::geom_text_repel(aes (label = rownames (mtcars)))+
  scale_fill_d3_adaptive()+
  scale_color_d3_adaptive()+
  labs(title="Ratio of peak with seed ",x="Number of filtered peak", y = "Adjusted ratio")+
  ylim(c(0,1))+
  # xlim(c(0,400))+
  theme_classic(base_size=12) + 
  theme(#axis.ticks.x=element_blank(),  
    #strip.text.y = element_blank(),
    #strip.text.x = element_text(face="bold",family="arial",size=20),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size = 20,hjust = 0.5,vjust = 0.5), # ,angle = 90
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size=20),
    legend.position = "right",# c(0.14,0.8),#
    legend.justification = "bottom", #c(0, 0),
    # legend.direction = "horizontal",
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16))
ggsave(filename = "seed_ratio_scatter_rev.pdf",width = 9,height = 5)
# ggsave(filename = "seed_ratio_scatter.pdf",width = 9,height = 5)
# ggplot(res2.df, aes(x=number,y=adj.num,fill=method,color=method))+
#   geom_line(size=1)+
#   ggrepel::geom_label_repel(aes(x=number,y=adj.num, #fill=method,   
#                                 label=as.character(p.value)),alpha=1,fill="white", fontface="bold", 
#                             force = 0.4, force_pull=0.4, 
#                             max.overlaps=10,
#                             box.padding=unit(0.15, "lines"), point.padding=unit(0.15, "lines"),
#                             segment.colour = "white"
#   )+
#   # ggrepel::geom_text_repel(aes (label = rownames (mtcars)))+
#   scale_fill_d3_adaptive()+
#   scale_color_d3_adaptive()+
#   labs(title="Number of peak with seed",x="Number of filtered peak", y = "Adjusted number")+
#   # ylim(c(0,1))+
#   # xlim(c(0,400))+
#   theme_classic(base_size=12) + 
#   theme(#axis.ticks.x=element_blank(),  
#     #strip.text.y = element_blank(),
#     #strip.text.x = element_text(face="bold",family="arial",size=20),
#     axis.title.x = element_text(size=20),
#     axis.title.y = element_text(size=20),
#     axis.text.x = element_text(size = 20,hjust = 0.5,vjust = 0.5), # ,angle = 90
#     axis.text.y = element_text(size = 20),
#     plot.title = element_text(size=20),
#     legend.position = "right",# c(0.14,0.8),#
#     legend.justification = "bottom", #c(0, 0),
#     # legend.direction = "horizontal",
#     legend.text = element_text(size= 16),
#     legend.title= element_text(size= 16))
# ggsave(filename = "seed_num_scatter.pdf",width = 9,height = 7)

# #bar
# fancy_scientific <- function(l) {
#   # turn in to character string in scientific notation
#   l <- format(l, scientific = TRUE)
#   # quote the part before the exponent to keep all the digits
#   l <- gsub("^(.*)e", "'\\1'e", l)
#   # turn the 'e+' into plotmath format
#   l <- gsub("e", "%*%10^", l)
#   # return this as an expression
#   parse(text=l)
# }
# bar <- ggplot(res2.df, aes(x=p.value, y=adj.ratio, color=method, fill=method)) +
#   geom_bar( stat = "identity", position=position_dodge(width=0.7))+
#   scale_fill_d3_adaptive()+
#   scale_color_d3_adaptive()+
#   # scale_x_discrete(label=fancy_scientific)+
#   labs(title="",x="", y = "")+
#   # facet_grid(method~.)+
#   #ylim(c(0,1))+
#   # xlim(c(0,1))+
#   theme_classic(base_size=12) + 
#   theme(#axis.ticks.x=element_blank(),  
#     #strip.text.y = element_blank(),
#     #strip.text.x = element_text(face="bold",family="arial",size=20),
#     axis.title.x = element_text(size=20),
#     axis.title.y = element_text(size=20),
#     axis.text.x = element_text(size = 16,hjust = 1,vjust = 0.5,angle = 90), # 
#     axis.text.y = element_text(size = 16),
#     plot.title = element_text(size=20),
#     strip.text = element_text(size = 20),
#     legend.position = "right",# c(0.85,0.7),#
#     legend.text = element_text(size= 16),
#     legend.title= element_text(size= 16))
# p <- scatter + patchwork::inset_element(bar,left=0.25,right=1.35,top=1.0,bottom=0.55)
# ggsave(plot = p, filename = "seed_ratio.pdf",width = 9,height = 7)








# evaluate peak/domain overlap venn/upset (deprecated, use intervene py pkg now) ---------------------------------
# setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/")
# #setwd("/Users/baopengfei/Desktop/lulab/tmp/projects/WCHSU-FTC/exSeek-dev")
# options(stringsAsFactors = F)
# 
# ## build merged union domain bed 
# cd /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev
# # clipper_pre="../output/AGO2_IP/call_domain_withRepeats_dedup_bk/domains_clipper/b5_p0001"
# # piranha_pre="../output/AGO2_IP/call_domain_withRepeats_dedup_bk/domains/b5_p0001"
# # localmax_pre="../output/AGO2_IP/call_domain_withRepeats_dedup_bk/domains_localmax_primary/b5_d05_p0001"
# # localmaxext_pre="../output/AGO2_IP/call_domain_withRepeats_dedup_bk/domains_localmax_EM/b5_d05_p0001"
# clipper_pre="../output/GSE71008/call_domain_withRepeats_all/domains_clipper/b5_p01"
# piranha_pre="../output/GSE71008/call_domain_withRepeats_all/domains/b5_p01"
# localmax_pre="../output/GSE71008/call_domain_withRepeats_all/domains_localmax/b5_d05_p01"
# localmaxext_pre="../output/GSE71008/call_domain_withRepeats_all/domains_localmax_EM/b5_d05_p01"
# localmaxextgini_pre="../output/GSE71008/call_domain_withRepeats_all/domains_localmax_EMgini/b5_d05_p01"
# # clipper_pre="../output/GSE133684/call_domain_withRepeats_dedupByPos/domains_clipper/b5_p01"
# # piranha_pre="../output/GSE133684/call_domain_withRepeats_dedupByPos/domains/b5_p01"
# # localmax_pre="../output/GSE133684/call_domain_withRepeats_dedupByPos/domains_localmax/b5_d05_p01"
# # localmaxext_pre="../output/GSE133684/call_domain_withRepeats_dedupByPos/domains_localmax_EM/b5_d05_p01"
# mergedBed="./test_merge.bed"
# 
# cat ${clipper_pre}.bed ${piranha_pre}.bed ${localmax_pre}.bed ${localmaxext_pre}.bed | bedtools sort | bedtools merge -i stdin > ${mergedBed}
# bedtools intersect -wao -a ${mergedBed} -b ${clipper_pre}.bed > ${clipper_pre}_intersect3UnionDomain.txt
# 
# 
# #R
# # clipper_pre <- "../output/AGO2_IP/call_domain_withRepeats_dedup_bk/domains_clipper/b5_p0001.bed"
# # piranha_pre <- "../output/AGO2_IP/call_domain_withRepeats_dedup_bk/domains/b5_p0001.bed"
# # localmax_pre <- "../output/AGO2_IP/call_domain_withRepeats_dedup_bk/domains_localmax_primary/b5_d05_p0001.bed"
# # localmaxext_pre <- "../output/AGO2_IP/call_domain_withRepeats_dedup_bk/domains_localmax_EM/b5_d05_p0001.bed"
# clipper_pre <- "../output/GSE71008/call_domain_withRepeats_all/domains_clipper/b5_p01.bed"
# piranha_pre <- "../output/GSE71008/call_domain_withRepeats_all/domains/b5_p01.bed"
# localmax_pre <- "../output/GSE71008/call_domain_withRepeats_all/domains_localmax/b5_d05_p01.bed"
# localmaxext_pre <- "../output/GSE71008/call_domain_withRepeats_all/domains_localmax_EM/b5_d05_p01.bed"
# localmaxextgini_pre <- "../output/GSE71008/call_domain_withRepeats_all/domains_localmax_EMgini/b5_d05_p01.bed"
# # clipper_pre <- "../output/GSE133684/call_domain_withRepeats_dedupByPos/domains_clipper/b5_p01.bed"
# # piranha_pre <- "../output/GSE133684/call_domain_withRepeats_dedupByPos/domains/b5_p01.bed"
# # localmax_pre <- "../output/GSE133684/call_domain_withRepeats_dedupByPos/domains_localmax/b5_d05_p01.bed"
# # localmaxext_pre <- "../output/GSE133684/call_domain_withRepeats_dedupByPos/domains_localmax_EM/b5_d05_p01.bed"
# mergedBed <- "./test_merge.bed"
# 
# CLIPper <- as.data.frame(data.table::fread(clipper_pre,data.table = F,sep = '\t',check.names = F,stringsAsFactors = F))
# CLIPper <- CLIPper[,1:6]
# colnames(CLIPper) <- c("chr","start","end","name","score","strand")
# Piranha <- as.data.frame(data.table::fread(piranha_pre,data.table = F,sep = '\t',check.names = F,stringsAsFactors = F))
# Piranha <- Piranha[,1:6]
# colnames(Piranha) <- c("chr","start","end","name","score","strand")
# LocalMax <- as.data.frame(data.table::fread(localmax_pre,data.table = F,sep = '\t',check.names = F,stringsAsFactors = F))
# LocalMax <- LocalMax[,1:6]
# colnames(LocalMax) <- c("chr","start","end","name","score","strand")
# LocalMaxExt <- as.data.frame(data.table::fread(localmaxext_pre,data.table = F,sep = '\t',check.names = F,stringsAsFactors = F))
# LocalMaxExt <- LocalMaxExt[,1:6]
# colnames(LocalMaxExt) <- c("chr","start","end","name","score","strand")
# LocalMaxExtGini <- as.data.frame(data.table::fread(localmaxextgini_pre,data.table = F,sep = '\t',check.names = F,stringsAsFactors = F))
# LocalMaxExtGini <- LocalMaxExtGini[,1:6]
# colnames(LocalMaxExtGini) <- c("chr","start","end","name","score","strand")
# 
# k <- list(CLIPper,Piranha,LocalMax,LocalMaxExt,LocalMaxExtGini)
# l <- c("CLIPper","Piranha","LocalMax","LocalMaxExt","LocalMaxExtGini")
# 
# domain <- as.data.frame(data.table::fread(mergedBed,data.table = F,sep = '\t',check.names = F,stringsAsFactors = F))
# domain$name <- "X"
# domain$score <- 1
# domain$strand <- "+"
# colnames(domain) <- c("chr","start","end","name","score","strand")
# 
# ref <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)
# 
# res <- list()
# for (i in 1:length(l)){
#   # x <- "CLIPper"
#   x <- l[i]
#   #  i <- 1
#   print(x)
#   
#   domain.ago2 <- bedtoolsr::bt.intersect(a = domain,wao = T,s = T,b = k[i])
#   domain.ago2 <- as_tibble(domain.ago2) %>% 
#     dplyr::mutate(V4=paste0(V1,":",V2,"-",V3),overlap_ratio_product = V13/(V3-V2) * (V13/(V9-V8)), site = x) %>% 
#     dplyr::arrange(desc(overlap_ratio_product)) %>% 
#     dplyr::distinct(V4, .keep_all = TRUE)   # only keep max overlap_ratio_product each peak
#   # dplyr::rename_all(funs(paste0(., "e", "E")))
#   #dplyr::rename()
#   colnames(domain.ago2) <- c( c("chr","start","end","name","score","strand"),paste0("site","_",c("chr","start","end","name","score","strand","overlap_base","overlap_ratio_product") ) ,"site" ) 
#   res[[x]] <- domain.ago2
# }
# res.df <- do.call("rbind",res)
# res.df$RNA <- ref$transcript_type[match(res.df$chr,ref$transcript_id)]
# res.df$RNA <- gsub("_for|_rev|\\.for|\\.rev","",res.df$RNA,perl = T)
# rna <- c("rRNA","pri_miRNA","lncRNA", "mRNA", "piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA')
# dna <- c("intron","promoter", "enhancer","repeats") #
# res.df$RNA <- factor(res.df$RNA,levels = c(rna,dna))
# table(res.df$RNA)
# 
# table(is.na(res.df$site_overlap_ratio_product))
# res.df$site_overlap_ratio_product[is.na(res.df$site_overlap_ratio_product)] <- 0
# 
# num <- nrow(res.df)*0.25
# 
# library(ggplot2)
# res.df$RNA <- gsub("_rev|_for|\\.for|\\.rev","",res.df$RNA)
# rna <- c("rRNA","pri_miRNA","lncRNA", "mRNA", "piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA')
# dna <- c("intron","promoter", "enhancer","repeats") # 
# res.df$RNA <- factor(res.df$RNA,levels = c(rna,dna))
# res.df$site <- factor(res.df$site,levels = c("Piranha","CLIPper","LocalMax","LocalMaxExt","LocalMaxExtGini"))
# res.df$pos <- paste0(res.df$chr,":",res.df$start,"-",res.df$end) 
# table(res.df$site)
# 
# res.tbl <- as_tibble(res.df) %>% 
#   dplyr::mutate(RBP_structure=dplyr::case_when( site=="CLIPper" & (site_overlap_ratio_product>=0.01 | site_overlap_base>=10) ~ T,  
#                                                 site=="Piranha" & (site_overlap_ratio_product>=0.01 | site_overlap_base>=10) ~ T,
#                                                 site=="LocalMax" & (site_overlap_ratio_product>=0.01 | site_overlap_base>=10) ~ T,
#                                                 site=="LocalMaxExt" & (site_overlap_ratio_product>=0.01 | site_overlap_base>=10) ~ T,
#                                                 site=="LocalMaxExtGini" & (site_overlap_ratio_product>=0.01 | site_overlap_base>=10) ~ T,
#                                                 # site=="structure" & (site_overlap_ratio_product<=0.05) ~ T,
#                                                 TRUE ~ F
#   )
#   )
# 
# mat <- res.tbl %>% 
#   dplyr::select(RBP_structure,site,name) %>% 
#   tidyr::pivot_wider(id_cols = name, names_from = site, values_from = RBP_structure)
# 
# mat <- as.data.frame(mat)
# rownames(mat) <- mat$name
# mat$name <- NULL
# mat[1:3,]
# mat <- as.matrix(mat)
# mat[mat==T] <- 1
# mat[mat==F] <- 0
# mat <- as.data.frame(mat)
# 
# setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/")
# pdf("domain_intersect_upset.pdf",width = 4,height = 3)
# UpSetR::upset(data = mat, sets = c("LocalMaxExtGini","LocalMaxExt","LocalMax","CLIPper","Piranha"), nsets = 4, keep.order = T,
#               matrix.color = "salmon",set_size.numbers_size = T,
#               sets.bar.color = pal_d3_adaptive()(5) # colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(5)[1:4]
# )
# dev.off()

if (suppressMessages(!require("UpSetR"))) suppressMessages(install.packages("UpSetR", repos="http://cran.us.r-project.org"))
library("UpSetR")
pdf("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/tmp/GSE71008_NCpool_venn/Intervene_upset.pdf", width=9, height=4, onefile=FALSE, useDingbats=FALSE)
expressionInput <- c('exPeak'=1047,'CLAM'=125,'CLAM&exPeak'=162,'CLIPper'=302,'CLIPper&exPeak'=1217,'CLIPper&CLAM'=4,'CLIPper&CLAM&exPeak'=1616,'Piranha'=0,'Piranha&exPeak'=35,'Piranha&CLAM'=0,'Piranha&CLAM&exPeak'=12,'Piranha&CLIPper'=0,'Piranha&CLIPper&exPeak'=297,'Piranha&CLIPper&CLAM'=0,'Piranha&CLIPper&CLAM&exPeak'=38)
upset(fromExpression(expressionInput), nsets=4, keep.order = F, nintersects=30,text.scale = 1.5, show.numbers="yes", matrix.color = "salmon",set_size.numbers_size = T, main.bar.color="grey30", sets.bar.color = c("#C9121E","#FC6910","#269321","#1B61A5"), empty.intersections=NULL, order.by = "freq", number.angles = 0, mainbar.y.label ="No. of Intersections", sets.x.label ="Set size")
invisible(dev.off())
#colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(4)
#dev.off()



# evaluate peak/domain G4 GG(G) motif freq. (not as expected, deprecated) ----------------------------
# need run in R4 (hub rstudio-server)
# BiocManager::install("monaLisa")
# #meet error: C++14 standard requested but CXX14 is not defined
# #fixed by "https://yuantian1991.github.io/notes/Solve-Bugs:-C++14-standard-requested-but-CXX14-is-not-defined"
library(monaLisa)
dst <- "GSE71008"
pre <- "/BioII/lulab_b/baopengfei/tmp"
fa.path <- Sys.glob(paste0(pre,"/",dst,"*/b5_*/G4.bed.fa"))
seq <- read.table(fa.path[2])
seq <- seq[!grepl(">",seq$V1),]
res <- monaLisa::getKmerFreq(seqs = seq[sample(1:length(seq),2000,replace = F)],  includeRevComp=F, zoops=F, kmerLen = 3) # kmerLen >= 3
names(res)
#>  [1] "freq.obs"    "freq.exp"    "log2enr"     "sqrtDelta"   "z"          
#>  [6] "p"           "padj"        "strata"      "freq.strata" "CpGoe"      
head(res$freq.obs)
#> AAA AAC AAG AAT ACA ACC 
#>   3   0   0   4   0   0 
head(res$freq.exp)
#seem GGG not enriched



# test
# pwms <- TFBSTools::getMatrixSet(JASPAR2020,
#                      opts = list(matrixtype = "PWM",
#                                  tax_group = "vertebrates"))

# db <- file.path(system.file("extdata", package="JASPAR2014"), 
#                 "JASPAR2014.sqlite")
# opts <- list()
# opts[["species"]] <- 9606
# opts[["type"]] <- "SELEX"
# opts[["all_versions"]] <- FALSE
# siteList <- getMatrixSet(db, opts)
db <- file.path("/BioII/lulab_b/baopengfei/shared_reference/meme/motif_databases/RNA/POSTAR3_RBP.meme")
pwms <- TFBSTools::getMatrixSet(db, opts = list(matrixtype = "PWM",tax_group = "vertebrates"))
#failed





# evaluate peak num,length stat (fig3, NCpool box/bar) -------------------------------------------------------------------------
## plot boxplot of Peak num stat
setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/")
options(stringsAsFactors = F)
options(bedtools.path = "/BioII/lulab_b/baopengfei/biosoft") # /BioII/lulab_b/baopengfei/anaconda3/envs/exvariance/bin
suppressPackageStartupMessages(library(bedtoolsr))
library(tidyverse)
library(dplyr)
conflict_prefer("filter", "dplyr")

#read ref
ref <- data.table::fread("../WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt",header = T,sep = "\t", stringsAsFactors = F)
ref <- ref[!is.na(ref$tx.length),]
ref[1:3,]


tmp.list <- list()
sum.num <- function(x){
  #x <- "output/AGO2_IP/call_domain_withRepeats_dedup_bk/domains/b5_p0001_gn.bed"
  print(x)
  
  # #x <- "output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax_by_sample/interectAGO2/FTA-12_1.bed"
  rbp <- data.table::fread(paste0(x),data.table = F, header = F,sep = "\t",stringsAsFactors = F)

  bed <- as_tibble(rbp)
  bed <- bed[,1:6]
  colnames(bed) <- c("chr","start","end","peak","score","strand")
  summary(bed$end-bed$start)
  bed <- bed %>% 
    filter(end-start<=200 , end-start>=10)
  
  # bedtoolsr::bt.getfasta(fo = paste0("tmp/",basename(x),"tmp.fa"),nameOnly = T, fi = "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/fasta/genome.fa",bed = bed, s = T)
  # fasta <- rtracklayer::import(paste0("tmp/",basename(x),"tmp.fa"),format = "fasta")
  # length(fasta)
  # #fasta2 <- as.data.frame(fasta[sample(1:length(fasta),0.1*length(fasta),replace=F)])
  # fasta <- as.data.frame(fasta)
  # ratio <- sum( grepl(sequences,fasta$x,perl=T) )/nrow(fasta)
  tmp <- data.frame(sample=basename(x),median.len=median(bed$end-bed$start),number=nrow(bed),path=x) # ratio=ratio,
  tmp.list[[x]] <- tmp
  # return(tmp)
}
sum.num.single <- function(x){
  #x <- "output/AGO2_IP/call_domain_withRepeats_dedup_bk/domains/b5_p0001_gn.bed"
  print(x)
  
  # #x <- "output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax_by_sample/interectAGO2/FTA-12_1.bed"
  rbp <- data.table::fread(paste0(x),data.table = F, header = F,sep = "\t",stringsAsFactors = F)
  
  bed <- as_tibble(rbp)
  bed <- bed[,1:6]
  colnames(bed) <- c("chr","start","end","peak","score","strand")
  summary(bed$end-bed$start)
  bed <- bed %>% 
    filter(end-start<=200 , end-start>=10)
  
  # bedtoolsr::bt.getfasta(fo = paste0("tmp/",basename(x),"tmp.fa"),nameOnly = T, fi = "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/fasta/genome.fa",bed = bed, s = T)
  # fasta <- rtracklayer::import(paste0("tmp/",basename(x),"tmp.fa"),format = "fasta")
  # length(fasta)
  # #fasta2 <- as.data.frame(fasta[sample(1:length(fasta),0.1*length(fasta),replace=F)])
  # fasta <- as.data.frame(fasta)
  # ratio <- sum( grepl(sequences,fasta$x,perl=T) )/nrow(fasta)
  # tmp <- data.frame(sample=basename(x),median.len=median(bed$end-bed$start),number=nrow(bed),path=x) # ratio=ratio,
  bed$sample=basename(x)
  bed$median.len=median(bed$end-bed$start)
  bed$number=nrow(bed)
  bed$path=x
  tmp.list[[x]] <- bed #tmp
  # return(tmp)
}

pre <- "GSE71008_NCpool/call_peak_all" # "AGO2_IP/call_domain_withRepeats_dedup_bk", "GSE71008/call_domain_withRepeats_all", "GSE133684/call_domain_withRepeats_dedupByPos"  
# in1 <- Sys.glob(paste0("output/",pre,"/piranha_by_sample/b5_p01/intersect/NCpool.bed"))
# in2 <- Sys.glob(paste0("output/",pre,"/clipper_by_sample/b5_p05/intersect/NCpool.bed"))
# in3 <- Sys.glob(paste0("output/",pre,"/clam_by_sample/b5_p005/intersect/NCpool.bed"))
# in5 <- Sys.glob(paste0("output/",pre,"/expeakCNN_by_sample/b5_d50_p1/intersect/NCpool.bed"))
in1 <- Sys.glob(paste0("output/",pre,"/piranha/b5_p01_11RNA.bed"))
in2 <- Sys.glob(paste0("output/",pre,"/clipper/b5_p05_11RNA.bed"))
in3 <- Sys.glob(paste0("output/",pre,"/clam/b5_p005_11RNA.bed"))
in5 <- Sys.glob(paste0("output/",pre,"/expeakCNN/b5_d50_p1_11RNA.bed"))

# in1 <- Sys.glob(paste0("output/",pre,"/piranha/b5_p*_8DNA.bed"))
# in5 <- Sys.glob(paste0("output/",pre,"/expeak/b5_d05_p*_8DNA.bed"))

### num plot
res2.df <- do.call(rbind,mclapply(c(in1,in2,in3,in5), sum.num, mc.cores = 1) )   # 11RNA
# res2.df <- do.call(rbind,mclapply(c(in1,in5), sum.num, mc.cores = 2) )   # 8DNA

res2.df$method <- "Piranha"
res2.df$method[grepl("expeak",res2.df$path)] <- "exPeak"
res2.df$method[grepl("clipper",res2.df$path)] <- "CLIPper"
res2.df$method[grepl("clam",res2.df$path)] <- "CLAM"
res2.df$method <- factor(res2.df$method,levels = c("Piranha","CLIPper","CLAM","exPeak")) # LocalMax
# res2.df$method <- factor(res2.df$method,levels = c("Piranha","exPeak")) # LocalMax

p1 <- ggplot(res2.df, aes(x=method, y=number, color=method, fill=method)) +
  geom_bar( stat = "identity", position=position_dodge(width=0.7),color="black", width=0.5)+
  scale_fill_d3_adaptive()+
  # scale_color_d3()+
  # scale_x_discrete(label=c("Domain Ratio","MIR Ratio"))+
  labs(title="",x="", y = "Peak number")+
  # facet_grid(method~.)+
  #ylim(c(0,1))+
  # xlim(c(0,1))+
  theme_minimal(base_size=12) + 
  theme(#axis.ticks.x=element_blank(),  
    #strip.text.y = element_blank(),
    #strip.text.x = element_text(face="bold",family="arial",size=20),
    aspect.ratio = 1,
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size = 20,hjust = 1,vjust = 0.5,angle = 90), # element_blank(), #
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size=20),
    strip.text = element_text(size = 20),
    legend.position = "right",# c(0.85,0.7),#
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16))
# ggsave("peak_num_stat.pdf")


### num plot (by RNA types)
res2.df <- do.call(rbind,mclapply(c(in1,in2,in3,  in5), sum.num.single, mc.cores = 5) )  # 11RNA
# res2.df <- do.call(rbind,mclapply(c(in1,  in5), sum.num.single, mc.cores = 2) )  # 8DNA

res2.df$method <- "Piranha"
res2.df$method[grepl("expeak",res2.df$path)] <- "exPeak"
res2.df$method[grepl("clipper",res2.df$path)] <- "CLIPper"
res2.df$method[grepl("clam",res2.df$path)] <- "CLAM"
res2.df$method <- factor(res2.df$method,levels = c("Piranha","CLIPper","CLAM","exPeak"))
# res2.df$method <- factor(res2.df$method,levels = c("Piranha","exPeak"))

#add RNA types
res2.df$RNA <- ref$transcript_type[match(res2.df$chr,ref$transcript_id)]
res2.df$RNA <- gsub("_rev|\\.rev|_for|\\.for","",res2.df$RNA,perl = T)
rna <- c("rRNA","pri_miRNA","lncRNA", "mRNA", "piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA')
# dna <- c("intron","promoter", "enhancer","repeats") #
res2.df$RNA <- factor(res2.df$RNA,levels = c(rna))
# res2.df$RNA <- factor(res2.df$RNA,levels = c(dna))
table(res2.df$RNA)
res2.df <- dplyr::as_tibble(res2.df)
res2.df2 <- res2.df %>% 
  dplyr::group_by(sample,method,RNA) %>%  # p.value
  dplyr::summarize(number=dplyr::n_distinct(peak))
res2.df <- as.data.frame(res2.df2)


p11 <- ggplot(res2.df, aes(x=method, y=number, color=RNA, fill=RNA)) +
  geom_bar( stat = "identity", position=position_dodge(width=0.7),color="black", width=0.5)+
  scale_fill_nejm_adaptive()+
  # scale_color_d3()+
  # scale_x_discrete(label=c("Domain Ratio","MIR Ratio"))+
  labs(title="",x="", y = "Peak number")+
  # facet_grid(method~.)+
  #ylim(c(0,1))+
  # xlim(c(0,1))+
  theme_minimal(base_size=12) + 
  theme(#axis.ticks.x=element_blank(),  
    #strip.text.y = element_blank(),
    #strip.text.x = element_text(face="bold",family="arial",size=20),
    aspect.ratio = 0.5,
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size = 20,hjust = 1,vjust = 0.5,angle = 90), # element_blank(), #
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size=20),
    strip.text = element_text(size = 20),
    legend.position = "right",# c(0.85,0.7),#
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16))
# ggsave("peak_num_stat.pdf")


### plot peak length boxplot
bed.list <- list()
sum.len <- function(x){
  # x <- "output/AGO2_IP/call_domain/domains/5/50_gn_intersectAGO2.bed"
  #x <- "output/AGO2_IP/call_domain/domains/5/50_gn.bed"
  print(x)
  
  setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/")
  options(stringsAsFactors = F)
  options(bedtools.path = "/BioII/lulab_b/baopengfei/biosoft") # /BioII/lulab_b/baopengfei/anaconda3/envs/exvariance/bin
  suppressPackageStartupMessages(library(bedtoolsr))
  library(tidyverse)
  library(dplyr)
  
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
res2.df2 <- do.call(rbind,mclapply(c(in1,in2,in3,  in5), sum.len, mc.cores = 5) )  # 11RNA
# res2.df2 <- do.call(rbind,mclapply(c(in1, in5), sum.len, mc.cores = 2) )  # 8DNA

res2.df2$method <- "Piranha"
res2.df2$method[grepl("expeak",res2.df2$path)] <- "exPeak"
res2.df2$method[grepl("clipper",res2.df2$path)] <- "CLIPper"
res2.df2$method[grepl("clam",res2.df2$path)] <- "CLAM"
res2.df2$method <- factor(res2.df2$method,levels = c("Piranha","CLIPper","CLAM","exPeak"))
# res2.df2$method <- factor(res2.df2$method,levels = c("Piranha","exPeak"))

res2.df2$width <- res2.df2$end-res2.df2$start

p2 <- ggplot(res2.df2, aes(x=method, y=width, color=method, fill=method)) +
  # geom_bar(position = "dodge", stat = "identity")+
  geom_boxplot(outlier.size = .05,outlier.colour = "grey", outlier.alpha = 0.9, alpha=0.9, color="black",  position=position_dodge(width=0.7), width=0.5)+
  # geom_jitter()+
  # geom_histogram(aes(y=..density..), position="identity", alpha=0.5)+
  # geom_density(alpha=0.1)+
  # geom_vline(data=mu, aes(xintercept=grp.mean, color=method),linetype="dashed")+
  # scale_color_manual(values=c("firebrick", "steelblue4", "grey30"))+
  # scale_fill_manual(values=c("firebrick", "steelblue4", "grey30"))+
  # geom_vline(xintercept = 0,color="grey30",linetype="dashed")+
  # geom_label(data=mu, aes(label=ratio, color=method, fill=method )) +
  #annotate(x = 0.5, y=0.5)+
  # geom_text(data=mu, nudge_x = 0.5, nudge_y = 0.5, aes(text=ratio, color=method))+
  scale_fill_d3_adaptive()+
  scale_color_d3_adaptive()+
  # scale_x_discrete(label=c("Domain Ratio","MIR Ratio"))+
  labs(title="",x="", y = "Peak length")+
  # facet_grid(method~.)+
  ylim(c(0,50))+
  # xlim(c(0,1))+
  theme_minimal(base_size=12) + 
  theme(#axis.ticks.x=element_blank(),  
    #strip.text.y = element_blank(),
    #strip.text.x = element_text(face="bold",family="arial",size=20),
    aspect.ratio = 1,
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size = 20,hjust = 1,vjust = 0.5,angle = 90), # 
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size=20),
    strip.text = element_text(size = 20),
    legend.position = "right", #c(0.9,0.8),#,# 
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16))




### plot peak length boxplot (by RNA types)
#table(res2.df2$strand)
res2.df2$RNA <- ref$transcript_type[match(res2.df2$chr,ref$transcript_id)]
res2.df2$RNA <- gsub("_rev|\\.rev|_for|\\.for","",res2.df2$RNA,perl = T)
rna <- c("rRNA","pri_miRNA","lncRNA", "mRNA", "piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA')
# dna <- c("intron","promoter", "enhancer","repeats") #
res2.df2$RNA <- factor(res2.df2$RNA,levels = c(rna))
# res2.df2$RNA <- factor(res2.df2$RNA,levels = c(dna))
table(res2.df2$RNA)

p22 <- ggplot(res2.df2, aes(x=method, y=width, color=RNA, fill=RNA)) +
  # geom_bar(position = "dodge", stat = "identity")+
  geom_boxplot(outlier.size = .05,outlier.colour = "grey", outlier.alpha = 0.9, alpha=0.9, color="black",  position=position_dodge(width=0.7), width=0.5)+
  # geom_jitter()+
  # geom_histogram(aes(y=..density..), position="identity", alpha=0.5)+
  # geom_density(alpha=0.1)+
  # geom_vline(data=mu, aes(xintercept=grp.mean, color=method),linetype="dashed")+
  # scale_color_manual(values=c("firebrick", "steelblue4", "grey30"))+
  # scale_fill_manual(values=c("firebrick", "steelblue4", "grey30"))+
  # geom_vline(xintercept = 0,color="grey30",linetype="dashed")+
  # geom_label(data=mu, aes(label=ratio, color=method, fill=method )) +
  #annotate(x = 0.5, y=0.5)+
  # geom_text(data=mu, nudge_x = 0.5, nudge_y = 0.5, aes(text=ratio, color=method))+
  scale_fill_nejm_adaptive()+
  scale_color_nejm_adaptive()+
  # scale_x_discrete(label=c("Domain Ratio","MIR Ratio"))+
  labs(title="",x="", y = "Peak length")+
  # facet_grid(method~.)+
  ylim(c(0,50))+
  # xlim(c(0,1))+
  theme_minimal(base_size=12) + 
  theme(#axis.ticks.x=element_blank(),  
    #strip.text.y = element_blank(),
    aspect.ratio = 0.5,
    #strip.text.x = element_text(face="bold",family="arial",size=20),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size = 20,hjust = 1,vjust = 0.5,angle = 90), # 
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size=20),
    strip.text = element_text(size = 20),
    legend.position = "right", #c(0.9,0.8),#,# 
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16))


cowplot::plot_grid(plotlist = list(p1,p2),rel_heights = c(1,1), nrow = 2, align = "hv", axis = "b")
# ggsave(filename = "domain_num_len.pdf", width = 8, height = 10) # 11RNA
ggsave(filename = "domain_num_len.pdf", width = 7, height = 10) # 8DNA

cowplot::plot_grid(plotlist = list(p11,p22),rel_heights = c(1,1), nrow = 2, align = "hv", axis = "b")
# ggsave(filename = "domain_num_len_pear_peak.pdf", width = 8, height = 10) # 11RNA
ggsave(filename = "domain_num_len_by_RNA.pdf", width = 16, height = 10) # 8DNA






# ## plot EM rescue reads RNA assign 
# #Detected Reads/Tx Num
# setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/")
# options(stringsAsFactors = F)
# options(bedtools.path = "/BioII/lulab_b/baopengfei/biosoft") # /BioII/lulab_b/baopengfei/anaconda3/envs/exvariance/bin
# suppressPackageStartupMessages(library(bedtoolsr))
# suppressPackageStartupMessages(library(Rsamtools))
# library(tidyverse)
# library(dplyr)
# 
# ref <- data.table::fread("exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt",header = T,stringsAsFactors = F)
# ref <- as_tibble(ref) %>% 
#   dplyr::select(transcript_id, transcript_type)
# # ref$transcript_id <- gsub("(-)","_neg",ref$transcript_id,fixed = T)
# # ref$transcript_id <- gsub("(+)","_pos",ref$transcript_id,fixed = T)
# # ref$transcript_id <- gsub("::","__",ref$transcript_id,fixed = T)
# # ref$transcript_id <- gsub(":","___",ref$transcript_id,fixed = T)
# # ref$transcript_id <- gsub("-","____",ref$transcript_id,fixed = T)
# # ref$transcript_id <- gsub(".","_____",ref$transcript_id,fixed = T)
# table(ref$transcript_type)
# #head(ref[ref$transcript_type=="enhancer.for",])
# 
# sum.read.assign <- function(x){
#   #x <- "output/AGO2_IP/tbam/AGO2_IP_SUP_s1/bam-deduped-EM/merge16_sort/realigned.sorted.bam"
#   # x <- in1[1]
#   print(x)
#   
#   bam <- Rsamtools::scanBam(Rsamtools::BamFile(x))
#   bam.df <- as.data.frame(bam)
#   bam.df <- as.data.frame(table(bam.df$rname))
#   bam.df <- bam.df[bam.df$Freq>0,]
#   bam.df$Var1 <- sub("______","__",bam.df$Var1,fixed = T)
#   bam.df[1:3,]
#   ref2 <- ref %>% 
#     filter(ref$transcript_id %in% bam.df$Var1)
#   
#   colnames(bam.df)[1] <- "transcript_id"
#   bam.df <- left_join(bam.df,ref2)
#   # table(bam.df$transcript_type)
#   # bam.df[is.na(bam.df$transcript_type),]
#   # 
#   bam.tbl <- as_tibble(bam.df) %>% 
#     dplyr::group_by(transcript_type) %>% 
#     dplyr::summarise(tx.num=dplyr::n(),read.num=sum(Freq))
#   
#   bam.tbl$path <- x
#   return(bam.tbl)
# }
# 
# in1=Sys.glob("output/AGO2_IP/tbam/*/bam-deduped-EM/merge18_sort/merged.sorted.bam")
# in2=Sys.glob("output/AGO2_IP/tbam/*/bam-deduped/merge10RNA_sort_primary.bam")
# 
# res.df <- do.call(rbind,mclapply(c(in1,in2), sum.read.assign, mc.cores = 1) )  # seem multi-core bug
# res.df.backup <- res.df
# res.df$method <- "Primary"
# res.df$method[grepl("merged",res.df$path)] <- "Merged"
# res.df$method <- factor(res.df$method,levels = c("Primary","Merged"))
# colnames(res.df)[1] <- "RNA"
# res.df[1:3,]
# res.df$RNA <- gsub("_rev|\\.rev","",res.df$RNA,perl = T)
# res.df$RNA <- gsub("_for|\\.for","",res.df$RNA,perl = T)
# rna <- c("pri_miRNA","lncRNA", "mRNA", "piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA')
# #dna <- c("intron.for", "intron.rev", "promoter.for", "promoter.rev", "enhancer.for", "enhancer.rev", "repeats.for", "repeats.rev")
# dna <- c("intron","promoter", "enhancer","repeats") # 
# res.df$RNA <- factor(res.df$RNA,levels = c(rna,dna))
# # res.df$RNA[is.na(res.df$RNA)] <- "repeats"
# # table(is.na(res.df$RNA))
# # res.df$RNA <- as.character(res.df$RNA)
# # str(res.df)
# #table(res2.df2$strand)
# ggplot(res.df, aes(x=RNA, y=log10(tx.num), color=method, fill=method)) +
#   geom_bar(position = "dodge", stat = "identity")+
#   # geom_boxplot()+
#   # geom_histogram(aes(y=..density..), position="identity", alpha=0.5)+
#   # geom_density(alpha=0.1)+
#   # geom_vline(data=mu, aes(xintercept=grp.mean, color=method),linetype="dashed")+
#   # scale_color_manual(values=c("firebrick", "steelblue4", "grey30"))+
#   # scale_fill_manual(values=c("firebrick", "steelblue4", "grey30"))+
#   # geom_vline(xintercept = 0,color="grey30",linetype="dashed")+
#   # geom_label(data=mu, aes(label=ratio, color=method, fill=method )) +
#   #annotate(x = 0.5, y=0.5)+
#   # geom_text(data=mu, nudge_x = 0.5, nudge_y = 0.5, aes(text=ratio, color=method))+
#   ggsci::scale_fill_d3()+
#   ggsci::scale_color_d3()+
#   # scale_x_discrete(label=c("Domain Ratio","MIR Ratio"))+
#   labs(title="Read Num",x="RNA", y = "Detected Tx Num (log10)")+
#   # facet_grid(method~.)+
#   #ylim(c(0,1))+
#   # xlim(c(0,1))+
#   theme_classic(base_size=12) + 
#   theme(#axis.ticks.x=element_blank(),  
#     #strip.text.y = element_blank(),
#     #strip.text.x = element_text(face="bold",family="arial",size=20),
#     axis.title.x = element_text(size=20),
#     axis.title.y = element_text(size=20),
#     axis.text.x = element_text(size = 20,hjust = 1,vjust = 0.5,angle = 90), # 
#     axis.text.y = element_text(size = 20),
#     plot.title = element_text(size=20),
#     strip.text = element_text(size = 20),
#     legend.position = c(0.9,0.8),#"right",# 
#     legend.text = element_text(size= 16),
#     legend.title= element_text(size= 16))
# 
# 
# #Significant Domain Num
# sum.TxPeak.assign <- function(x){
#   #x <- "output/AGO2_IP/call_domain/domains_localmax_EM/001.bed"
#   print(x)
#   
#   bam.df <- data.table::fread(x,stringsAsFactors = F,sep = "\t",header = F)
#   bam.df <- bam.df[,c(1,4)]
#   colnames(bam.df)[1] <- "transcript_id"
#   bam.df$transcript_id <- gsub("______","__",bam.df$transcript_id)
#   ref2 <- ref %>% 
#     filter(ref$transcript_id %in% bam.df$transcript_id)
#   # table(ref$transcript_type)
#   bam.df <- left_join(bam.df,ref2)
#   # res.df[is.na(res.df$transcript_type),]
#   
#   bam.tbl <- as_tibble(bam.df) %>% 
#     dplyr::group_by(transcript_type) %>% 
#     dplyr::summarise(domain.num=dplyr::n())
#   
#   bam.tbl$path <- x
#   return(bam.tbl)
# }
# # in1=Sys.glob("output/AGO2_IP/call_domain_withRepeats_dedup/domains*/5/*.bed")
# # in2=Sys.glob("output/AGO2_IP/call_domain_withRepeats_dedup/domains_localmax_*/*.bed")
# # infile <- c(in2)
# # infile <- infile[!grepl("interect|gn|by_sample|recurrence",infile,perl = T)]
# res <- mclapply(c("output/AGO2_IP/call_domain_withRepeats_dedup/domains/b5_p01.bed",
#                   "output/AGO2_IP/call_domain_withRepeats_dedup/domains_clipper/b5_p01.bed",
#                   "output/AGO2_IP/call_domain_withRepeats_dedup/domains_localmax_primary/b5_d05_p01.bed",
#                   "output/AGO2_IP/call_domain_withRepeats_dedup/domains_localmax_EM/b5_d05_p01.bed"), 
#                 sum.TxPeak.assign, mc.cores = 1) 
# res.df <- do.call(rbind,res)  # seem multi-core bug
# table(res.df$transcript_type)
# res.df$method <- "Primary"
# res.df$method[grepl("EM",res.df$path)] <- "Merged"
# res.df$method <- factor(res.df$method,levels = c("Primary","Merged"))
# colnames(res.df)[1] <- "RNA"
# res.df$RNA <- gsub(".rev","",res.df$RNA,fixed = T)
# res.df$RNA <- gsub(".for","",res.df$RNA,fixed = T)
# res.df$RNA <- gsub("_rev","",res.df$RNA,fixed = T)
# res.df$RNA <- gsub("_for","",res.df$RNA,fixed = T)
# rna <- c("pri_miRNA","lncRNA", "mRNA", "piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA')
# #dna <- c("intron.for", "intron.rev", "promoter.for", "promoter.rev", "enhancer.for", "enhancer.rev", "repeats.for", "repeats.rev")
# dna <- c("intron","promoter", "enhancer", "repeats") # repeats
# res.df$RNA <- factor(res.df$RNA,levels = c(rna,dna))
# table(res.df$RNA)
# unique(res.df$RNA)
# res.df2 <- res.df # res.df[grepl("/0001",res.df$path),]
# 
# ggplot(res.df2, aes(x=RNA, y=log10(domain.num), color=method, fill=method)) +
#   geom_bar(position = "dodge", stat = "identity")+
#   # geom_boxplot()+
#   # geom_histogram(aes(y=..density..), position="identity", alpha=0.5)+
#   # geom_density(alpha=0.1)+
#   # geom_vline(data=mu, aes(xintercept=grp.mean, color=method),linetype="dashed")+
#   # scale_color_manual(values=c("firebrick", "steelblue4", "grey30"))+
#   # scale_fill_manual(values=c("firebrick", "steelblue4", "grey30"))+
#   # geom_vline(xintercept = 0,color="grey30",linetype="dashed")+
#   # geom_label(data=mu, aes(label=ratio, color=method, fill=method )) +
#   #annotate(x = 0.5, y=0.5)+
#   # geom_text(data=mu, nudge_x = 0.5, nudge_y = 0.5, aes(text=ratio, color=method))+
#   ggsci::scale_fill_d3()+
#   ggsci::scale_color_d3()+
#   # scale_x_discrete(label=c("Domain Ratio","MIR Ratio"))+
#   labs(title="Domain Num",x="RNA", y = "Domain Num (log10)")+
#   # facet_grid(method~.)+
#   #ylim(c(0,1))+
#   # xlim(c(0,1))+
#   theme_classic(base_size=12) + 
#   theme(#axis.ticks.x=element_blank(),  
#     #strip.text.y = element_blank(),
#     #strip.text.x = element_text(face="bold",family="arial",size=20),
#     axis.title.x = element_text(size=20),
#     axis.title.y = element_text(size=20),
#     axis.text.x = element_text(size = 20,hjust = 1,vjust = 0.5,angle = 90), # 
#     axis.text.y = element_text(size = 20),
#     plot.title = element_text(size=20),
#     strip.text = element_text(size = 20),
#     legend.position = c(0.9,0.8),#"right",# 
#     legend.text = element_text(size= 16),
#     legend.title= element_text(size= 16))










# evaluate peak CPM/RPPM/depth stat (fig3/suppl, NCpool sample-wise/consensus half violin plot) ----------------
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"

## read ref
ref <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)
table(ref$transcript_type)

dst <- "GSE71008_NCpool"
dedup <- "all"
gold <- paste0(pre,"/exSeek-dev/tmp/4dst_filter2_11RNA.bed")
smp <- "NCpool"

#need pre-run norma matrix in bash !!!!
l <- list()
## consensus peak
for (method in c("piranha_b5_p01","clipper_b5_p05","clam_b5_p005","expeakCNN_b5_d50_p1") ){
  #method <- "clipper_b5_p05"
  print(method)
  y.cpm <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_",dedup,"/count_matrix/",method,"_CPM.txt"),header = T)
  
  y.cpm$method <- method
  y.cpm$chr <- unlist(sapply(strsplit(y.cpm$gene_id, "|", fixed = T), "[", 3))
  y.cpm$start <- as.numeric(unlist(sapply(strsplit(y.cpm$gene_id, "|", fixed = T), "[", 7)))
  y.cpm$end <- as.numeric(unlist(sapply(strsplit(y.cpm$gene_id, "|", fixed = T), "[", 6)))
  y.cpm$id <- (unlist(sapply(strsplit(y.cpm$gene_id, "|", fixed = T), "[", 4)))
  y.cpm$width <- y.cpm$start - y.cpm$end 
  # hist(y.cpm$width,breaks = 100,xlim = c(0,60))
  
  l[[method]] <- y.cpm
}
## sample-wise peak
# for (method in c("piranha_by_sample/b5_p01","clipper_by_sample/b5_p05","clam_by_sample/b5_p005","expeakCNN_by_sample/b5_d50_p1") ){
#   #method <- "clipper_by_sample/b5_p05"
#   print(method)
#   count <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_",dedup,"/",method,"/intersect/",smp,".bed6.count"),header = F)
#   if (grepl("expeak|clam",method,perl = T)){
#     lib <- "EM"
#   }else if(grepl("piranha|clipper",method,perl = T)){
#     lib <- "primary"
#   }
#   count$width <- count$V3-count$V2
#   #summary(count$width)
#   #hist(count$width,breaks = 100)
#   count <- count[count$width>=10 & count$V5>=3,]
#   
#   print(lib)
#   libSize <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_",dedup,"/count_matrix/",lib,".libSize"),header = F)$V1
#   print(sum(count$V5)) # should be TRUE !!!!
#   print(libSize)
#   
#   mat <- data.frame(row.names = paste0(count$V4,"|",count$V1,"|",count$V2), smp=count$V5)
#   colnames(mat) <- smp
#   y <- edgeR::DGEList(counts=mat,lib.size=libSize)
#   # y <- edgeR::calcNormFactors(y, method="none") # not needed for 1 smp
#   y.cpm <- as.data.frame(edgeR::cpm(y,normalized.lib.sizes = F,
#                                     log = F, 
#                                     prior.count = 1) 
#                          # group=y$sample$group
#   )
#   y.cpm$method <- method
#   y.cpm$chr <- count$V1
#   y.cpm$start <- count$V2
#   y.cpm$end <- count$V3
#   y.cpm$id <- count$V4
#   
#   l[[method]] <- y.cpm
# }

df <- do.call(rbind,l)
df$path <- df$method
df$method <- ""
df$method[grepl("piranha",df$path)] <- "Piranha"
df$method[grepl("clipper",df$path)] <- "CLIPper"
df$method[grepl("expeak",df$path)] <- "exPeak"
df$method[grepl("clam",df$path)] <- "CLAM"
df$method <- factor(df$method, levels = c("Piranha", "CLIPper","CLAM","exPeak")) #  c("exPeak", "CLAM","CLIPper","Piranha")
df[1:3,1:3]
table(df$method)

df$RNA <- ref$transcript_type[match(df$chr,ref$transcript_id)]
rna <- c("rRNA","pri_miRNA","lncRNA", "mRNA", "piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA')
# dna <- c("intron","promoter", "enhancer","repeats") #
df <- df[df$RNA %in% rna,]
df$RNA <- factor(df$RNA,levels = c(rna))
table(df$RNA)

library(ggridges)
#table(is.na(df2$Sample1))
summary(df$NCpool)
df$NCpool.log <- log10(df$NCpool+1)
ggplot(df, aes(x=NCpool.log, y=method, fill=method))+ #
  geom_density_ridges_gradient(scale = 0.6,na.rm = T, rel_min_height = 0.01, gradient_lwd = 3) + 
  # theme_ridges() + 
  # theme(legend.position = "none")
  coord_flip() +
  scale_fill_manual(values = pal_nejm_adaptive()(4)[c(2,3,4,1)]) +
  # scale_fill_nejm_adaptive() + # manual(name = '', values = alpha(colour = c("#2C68A9","#843692","#FD6905"),alpha = 0.8)) +
  xlab('') + ylab('') + # Peak count per million
  scale_x_continuous(breaks = c(0,1,2,3,4),limits = c(0,max(df$NCpool.log)), labels = c("0","10",expression(paste("10"^{2})),expression(paste("10"^{3})),expression(paste("10"^{4})))) +
  # theme_bw() + 
  # xlim(c(0,2))+
  # ylim(c(1,4))+
  # facet_wrap(facets = "RNA") +
  theme_minimal() +  # base_size=12
  theme(#axis.ticks.x=element_blank(),
    #strip.text.y = element_blank(),
    aspect.ratio = 1,
    strip.text = element_text(size=20),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size = 20, angle = 45,vjust = 0.5,hjust = 1),
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size=20),
    # strip.text = element_blank(),
    legend.position = "none", #c(0.9,0.8),#,#
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16))
ggsave("cpm_half_violin.pdf",width = 4,height = 5) # Fig3

ggplot(df, aes(x=method, y=NCpool.log, color=RNA, fill=RNA))+ #
  # geom_density_ridges_gradient(scale = 0.6,na.rm = T, rel_min_height = 0.01, gradient_lwd = 3) + 
  geom_boxplot(outlier.size = .05,outlier.colour = "grey", outlier.alpha = 0.9, alpha=0.9, color="black",  position=position_dodge(width=0.7), width=0.5)+
  # coord_flip() +
  scale_fill_nejm_adaptive()+
  scale_color_nejm_adaptive()+
  # scale_fill_nejm_adaptive() + # manual(name = '', values = alpha(colour = c("#2C68A9","#843692","#FD6905"),alpha = 0.8)) +
  xlab('') + ylab('') + # Peak count per million
  scale_y_continuous(breaks = c(0,1,2,3,4),limits = c(0,max(df$NCpool.log)), labels = c("0","10",expression(paste("10"^{2})),expression(paste("10"^{3})),expression(paste("10"^{4})))) +
  theme_minimal(base_size=12) + 
  theme(#axis.ticks.x=element_blank(),  
    #strip.text.y = element_blank(),
    #strip.text.x = element_text(face="bold",family="arial",size=20),
    aspect.ratio = 0.5,
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size = 20,hjust = 1,vjust = 0.5,angle = 90), # 
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size=20),
    strip.text = element_text(size = 20),
    legend.position = "none", #c(0.9,0.8),#,# 
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16))
ggsave("cpm_half_violin_by_RNA.pdf",width = 8, height = 5) # sup fig3?



# cmp peak depth (from bigwig)
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
dst <- "GSE71008_NCpool" # GSE71008_NCpool, "GSE71008"
smp <- "NCpool" # "NCpool", "SAMN03863396"

library("GenomicRanges")
getMaxDep <- function(bw,chr,start,width){
# bw <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008_NCpool/call_peak_all/tbigwig_RNA_EM/NCpool.transcriptome.bigWig"
# chr <- "NR_023363_____1"
# start <- 83
# width <- 37
  
  #bg <- rtracklayer::import.bw(bw, which=GRanges(c(chr), IRanges(start = start+1, width = width)))
bg <- subsetByOverlaps(bw, GRanges(seqnames=chr,ranges=IRanges(start=start, end=start+width)) )
return(max(bg$score))
}

res.list <- list()
## sample-wise peak
# readSmpPeak <- function(method){
#   # method <- "expeak_by_sample/b5_d50_p1"
#   print(method)
#   peakIn <- paste0(pre,"/output/",dst,"/call_peak_all/",method,"/intersect/",smp,".bed")
#   bed <- data.table::fread(peakIn,header = F,sep = "\t",check.names = F)
#   bed <- as.data.frame(bed[,1:6])
#   colnames(bed) <- c("seqnames","start","end","width","score","strand")
#   bed$width <- bed$end-bed$start
#   bed$method <- method
#   bed[1:3,]
#   # bed <- bed[1:10,]
#   
#   if(grepl("piranha|clipper",method,perl = T)){
#     bwType <- "tbigwig_11RNA_primary"
#   }else if(grepl("clam|expeak",method,perl = T)){
#     bwType <- "tbigwig_RNA_EM"
#   }
#   
#   # library(parallel)
#   # bed$maxDep <- do.call(c,mclapply(1:nrow(bed), 
#   #                        function(j) {getMaxDep(bw=paste0(pre,"/output/",dst,"/call_peak_all/",bwType,"/",smp,".transcriptome.bigWig"), 
#   #                                                          chr=bed$seqnames[j],start=bed$start[j],width=bed$width[j])},
#   #                        mc.cores = 10)
#   # )
#   bed$maxDep <- do.call(c,lapply(1:nrow(bed), function(j) {getMaxDep(bw=paste0(pre,"/output/",dst,"/call_peak_all/",bwType,"/",smp,".transcriptome.bigWig"),
#                                                            chr=bed$seqnames[j],start=bed$start[j],width=bed$width[j])} )
#                         )
#   return(bed)
# }
# res.list <- lapply(c("piranha_by_sample/b5_p01","clipper_by_sample/b5_p05","clam_by_sample/b5_p005","expeak_by_sample/b5_d50_p1"), readSmpPeak)
## consensus peak
readSmpPeak <- function(method){
  # method <- "expeakCNN/b5_d50_p1"
  print(method)
  peakIn <- paste0(pre,"/output/",dst,"/call_peak_all/",method,"_11RNA.bed")
  bed <- data.table::fread(peakIn,header = F,sep = "\t",check.names = F)
  bed <- as.data.frame(bed[,1:6])
  colnames(bed) <- c("seqnames","start","end","width","score","strand")
  bed$width <- bed$end-bed$start
  bed$method <- method
  # bed <- bed[1:10,]
  
  if(grepl("piranha|clipper",method,perl = T)){
    bwType <- "tbigwig_11RNA_primary"
  }else if(grepl("clam|expeak",method,perl = T)){
    bwType <- "tbigwig_RNA_EM"
  }
  
  # library(parallel)
  # bed$maxDep <- do.call(c,mclapply(1:nrow(bed), 
  #                        function(j) {getMaxDep(bw=paste0(pre,"/output/",dst,"/call_peak_all/",bwType,"/",smp,".transcriptome.bigWig"), 
  #                                                          chr=bed$seqnames[j],start=bed$start[j],width=bed$width[j])},
  #                        mc.cores = 10)
  # )
  bg <- rtracklayer::import.bw(paste0(pre,"/output/",dst,"/call_peak_all/",bwType,"/",smp,".transcriptome.bigWig") )
  bed$maxDep <- do.call(c,lapply(1:nrow(bed), function(j) {getMaxDep(bw=bg, chr=bed$seqnames[j],start=bed$start[j],width=bed$width[j])} )
  )
  return(bed)
}
res.list <- lapply(c("piranha/b5_p01","clipper/b5_p05","clam/b5_p005","expeak/b5_d50_p1"), readSmpPeak)
df <- as.data.frame(do.call(rbind,res.list))
# tmp <- (res.list[[2]])
# summary(tmp$maxDep)
# table(tmp$maxDep==4)



#plot
df$path <- df$method
df$method[grepl("piranha",df$path)] <- "Piranha"
df$method[grepl("clipper",df$path)] <- "CLIPper"
df$method[grepl("expeak",df$path)] <- "exPeak"
df$method[grepl("clam",df$path)] <- "CLAM"
df$method <- factor(df$method, levels = c("Piranha", "CLIPper","CLAM","exPeak"))
table(df$method)
library(ggridges)
#table(is.na(df2$Sample1))
df$maxDep2 <- as.numeric(unlist(df$maxDep))

df$NCpool.log <- log10(as.numeric(df$maxDep2)+1)
df$RNA <- ref$transcript_type[match(df$seqnames,ref$transcript_id)]
df$RNA <- gsub("_rev|_for","",df$RNA)
table(df$RNA,df$method)
# ggplot(df, aes(x=NCpool.log, y=method, fill=method))+ #
#   geom_density_ridges_gradient(scale = 0.6,na.rm = T, rel_min_height = 0.01, gradient_lwd = 3) + 
#   # theme_ridges() + 
#   # theme(legend.position = "none")
#   coord_flip() +
#   scale_fill_manual(values = pal_nejm_adaptive()(4)[c(2,3,4,1)]) +
#   # scale_fill_nejm_adaptive() + # manual(name = '', values = alpha(colour = c("#2C68A9","#843692","#FD6905"),alpha = 0.8)) +
#   xlab('') + ylab('') + # Peak count per million
#   # scale_x_continuous(breaks = c(0,1,2,3,4),limits = c(0,max(df$NCpool.log)), labels = c("0","10",expression(paste("10"^{2})),expression(paste("10"^{3})),expression(paste("10"^{4})))) +
#   # theme_bw() + 
#   # xlim(c(0,2))+
#   # ylim(c(1,4))+
#   # facet_wrap(facets = "RNA") +
#   theme_minimal() +  # base_size=12
#   theme(#axis.ticks.x=element_blank(),
#     #strip.text.y = element_blank(),
#     aspect.ratio = 1,
#     strip.text = element_text(size=20),
#     axis.title.x = element_text(size=20),
#     axis.title.y = element_text(size=20),
#     axis.text.x = element_text(size = 20, angle = 45,vjust = 0.5,hjust = 1),
#     axis.text.y = element_text(size = 20),
#     plot.title = element_text(size=20),
#     # strip.text = element_blank(),
#     legend.position = "none", #c(0.9,0.8),#,#
#     legend.text = element_text(size= 16),
#     legend.title= element_text(size= 16)) +
#   facet_grid(RNA~.)
# #ggsave("depth_half_violin.pdf",width = 4,height = 5)
# ggsave("depth_half_violin_by_RNA.pdf",width = 4,height = 25)

summary(unlist(df$maxDep2))
df2.list <- list()
for(method in c("Piranha", "CLIPper","CLAM","exPeak")){
  #method <- "exPeak"
  df.tmp <- df[as.character(df$method)==method,]
  #table(duplicated(df.tmp$seqnames))
  df2.tmp <- df.tmp[df.tmp$seqnames %in% df.tmp$seqnames[duplicated(df.tmp$seqnames)],] # only keep tx with >= 2 peaks
  df2.tmp <- df2.tmp[order(df2.tmp$method,df2.tmp$seqnames,df2.tmp$start,df2.tmp$end),]
  
  df2.tmp$min.dist <- -1 # add distance to nearest peak
  df2.tmp$foldchange <- -1 # add fold change of depth
  chr.idx <- base::which(!duplicated(df2.tmp$seqnames))
  chr.idx <- c(chr.idx,nrow(df2.tmp)+1)
  
  df2.tmp.list <- list()
  for(i in 1:(length(chr.idx)-1)){
    j <- chr.idx[i]
    k <- chr.idx[i+1]
    df2.tmp2 <- df2.tmp[j:(k-1),]
    # df2.tmp2.list <- list()
    for (m in 1:(nrow(df2.tmp2))){
      df2.tmp2$min.dist[m] <- min(c(df2.tmp2$start[m]-df2.tmp2$end[m-1],df2.tmp2$start[m+1]-df2.tmp2$end[m]),na.rm = T)
      df2.tmp2$foldchange[m] <- max(c(df2.tmp2$maxDep2[m]/df2.tmp2$maxDep2[m-1],df2.tmp2$maxDep2[m-1]/df2.tmp2$maxDep2[m],df2.tmp2$maxDep2[m+1]/df2.tmp2$maxDep2[m],df2.tmp2$maxDep2[m]/df2.tmp2$maxDep2[m+1]),na.rm = T)
    }
    # m <- 1
    # df2.tmp2$min.dist[m] <- df2.tmp2$start[m+1]-df2.tmp2$end[m])
    # m <- nrow(df2.tmp2)
    # df2.tmp2$min.dist[m] <- df2.tmp2$start[m]-df2.tmp2$end[m-1])
    df2.tmp.list[[i]] <- as.data.frame(df2.tmp2)
  }
  df2.tmp <- as.data.frame(do.call(rbind,df2.tmp.list))
  table(df2.tmp$min.dist==-1)

  df2.list[[method]] <- df2.tmp
}
df2 <- do.call(rbind,df2.list)
df2 <- df2[order(df2$method,df2$seqnames,df2$start,df2$end),]
table(is.na(df2$min.dist))
summary(df2$min.dist)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1       8      22     769     142   41956
summary(df2$min.dist[df2$method=="exPeak"])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1       2      15     446      34   41954
summary(df2$min.dist[df2$method=="CLIPper"])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2      20      37     827     196   41953
hist(df2$min.dist[df2$method=="exPeak"],xlim=c(0,50), breaks = 50000)
table(df2$min.dist>=5 & df2$min.dist<=200) # filter min distance
#df2 <- df2[df2$min.dist>=5 & df2$min.dist<=200,]


summary(df2$foldchange)
df2$log10.fc <- log10(as.numeric(df2$foldchange))
ggplot(df2, aes(x=log10.fc, y=method, fill=method))+ #
  geom_density_ridges_gradient(scale = 0.6,na.rm = T, rel_min_height = 0.01, gradient_lwd = 3) + 
  # theme_ridges() + 
  # theme(legend.position = "none")
  coord_flip() +
  scale_fill_manual(values = pal_nejm_adaptive()(4)[c(2,3,4,1)]) +
  # scale_fill_nejm_adaptive() + # manual(name = '', values = alpha(colour = c("#2C68A9","#843692","#FD6905"),alpha = 0.8)) +
  xlab('') + ylab('') + # Peak count per million
  scale_x_continuous(breaks = c(0,1,2,3,4),limits = c(0,max(df2$log10.fc)), labels = c("0","10",expression(paste("10"^{2})),expression(paste("10"^{3})),expression(paste("10"^{4})))) +
  # theme_bw() + 
  # xlim(c(0,2))+
  # ylim(c(1,4))+
  # facet_wrap(facets = "RNA") +
  theme_minimal() +  # base_size=12
  theme(#axis.ticks.x=element_blank(),
    #strip.text.y = element_blank(),
    aspect.ratio = 1,
    strip.text = element_text(size=20),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size = 20, angle = 45,vjust = 0.5,hjust = 1),
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size=20),
    # strip.text = element_blank(),
    legend.position = "none", #c(0.9,0.8),#,#
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16))
ggsave("filtered_depthFoldChange_half_violin.pdf",width = 4,height = 5)
dev.off()
#




# 
# # CPM rank bar plot
# #cpm <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_",dedup,"/count_matrix/expeak_b5_d50_p1_CPM.txt"),header = T,row.names = 1)
# l2 <- list()
# for(method in as.character(unique(df$method))){
#   #method <- "exPeak"
#   print(method)
#   cpm <- df[df$method==method,]
#   cpm <- cpm[order(cpm[[smp]],decreasing = T),,drop=F]
#   #cpm[1:3,]
#   cpm$rank <- 1:nrow(cpm)
#   # cpm$rank <- factor(cpm$rank,levels = unique(cpm$id))
#   cpm$cumsumCPM <- cumsum(cpm[[smp]])
#   l2[[method]] <- cpm
# }
# df2 <- do.call(rbind,l2)
# #str(df2)
# ggplot(df2, aes(x=rank, y=cumsumCPM))+ #
#   geom_bar(aes(fill=method,color=method),stat="identity",) + # ,color="grey50", fill="grey50", linewidth=0.001
#   geom_hline(yintercept = 1000000,color="salmon",linetype="dashed") +
#   # geom_vline(xintercept = (1:max(df2$rank))+0.5,color="salmon",linetype="dashed",size=0.0001) +
#   scale_fill_manual(values = pal_nejm_adaptive()(4)[c(2,3,4,1)]) +
#   scale_color_manual(values = pal_nejm_adaptive()(4)[c(2,3,4,1)]) +
#   # scale_fill_nejm_adaptive() + # manual(name = '', values = alpha(colour = c("#2C68A9","#843692","#FD6905"),alpha = 0.8)) +
#   # scale_color_nejm_adaptive() +
#   xlab('peak rank by CPM') + ylab('cumsum CPM') + # Peak count per million
#   scale_y_continuous(breaks = c(0,1000000), limits = c(0,1000000), labels = c("0M","1M")) +
#   facet_grid(method~.) +
#   theme_minimal() +  # base_size=12
#   # facet_grid(method~.,rows = 4) +
#   theme(#axis.ticks.x=element_blank(),
#     #strip.text.y = element_blank(),
#     # aspect.ratio = 0.1,
#     # strip.background = element_rect(color="grey",fill=alpha(colour = pal_nejm_adaptive()(4),alpha = 0.6)),
#     strip.text = element_text(size=20),
#     axis.title.x = element_text(size=20),
#     axis.title.y = element_text(size=20),
#     axis.text.x = element_blank(), # element_text(size = 20, angle = 45,vjust = 0.5,hjust = 1),
#     axis.text.y = element_text(size = 20),
#     plot.title = element_text(size=20),
#     # strip.text = element_blank(),
#     legend.position = "none", #c(0.9,0.8),#,#
#     legend.text = element_text(size= 16),
#     legend.title= element_text(size= 16))
# ggsave("CPM_rank_bar.pdf")
# 
# 
# ggplot(df2, aes(x=RNA, y=rank, fill=RNA))+ #
#   geom_boxplot() +
#   # geom_hline(yintercept = 1000000,color="salmon",linetype="dashed") +
#   # geom_vline(xintercept = (1:max(df2$rank))+0.5,color="salmon",linetype="dashed",size=0.0001) +
#   # scale_fill_manual(values = pal_nejm_adaptive()(4)[c(2,3,4,1)]) +
#   # scale_fill_nejm_adaptive() + # manual(name = '', values = alpha(colour = c("#2C68A9","#843692","#FD6905"),alpha = 0.8)) +
#   # xlab('peak rank by CPM') + ylab('cumsum CPM') + # Peak count per million
#   # scale_y_continuous(breaks = c(0,1000000), limits = c(0,1000000), labels = c("0M","1M")) +
#   facet_grid(method~.) +
#   theme_minimal() +  # base_size=12
#   # facet_grid(method~.,rows = 4) +
#   theme(#axis.ticks.x=element_blank(),
#     #strip.text.y = element_blank(),
#     # aspect.ratio = 0.1,
#     strip.text = element_text(size=20),
#     axis.title.x = element_text(size=20),
#     axis.title.y = element_text(size=20),
#     axis.text.x = element_blank(), # element_text(size = 20, angle = 45,vjust = 0.5,hjust = 1),
#     axis.text.y = element_text(size = 20),
#     plot.title = element_text(size=20),
#     # strip.text = element_blank(),
#     legend.position = "right", #c(0.9,0.8),#,#
#     legend.text = element_text(size= 16),
#     legend.title= element_text(size= 16))
# ggsave("CPM_rank_box_byRNA.pdf")




# evaluate peak num,length stat (fig4 deprecated?, box/bar, 8DNA) -------------------------------------------------------------------------


#Fig3: standard-overlapped peak number see plot_AUC_AUPR.R !!!!


## plot boxplot of Peak num stat
setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/")
options(stringsAsFactors = F)
options(bedtools.path = "/BioII/lulab_b/baopengfei/biosoft") # /BioII/lulab_b/baopengfei/anaconda3/envs/exvariance/bin
suppressPackageStartupMessages(library(bedtoolsr))
library(tidyverse)
library(dplyr)
conflict_prefer("filter", "dplyr")

#read ref
ref <- data.table::fread("../WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt",header = T,sep = "\t", stringsAsFactors = F)
ref <- ref[!is.na(ref$tx.length),]
ref[1:3,]


tmp.list <- list()
sum.num <- function(x){
  #x <- "output/AGO2_IP/call_domain_withRepeats_dedup_bk/domains/b5_p0001_gn.bed"
  print(x)
  
  # #x <- "output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax_by_sample/interectAGO2/FTA-12_1.bed"
  rbp <- data.table::fread(paste0(x),data.table = F, header = F,sep = "\t",stringsAsFactors = F)
  
  bed <- as_tibble(rbp)
  bed <- bed[,1:6]
  colnames(bed) <- c("chr","start","end","peak","score","strand")
  summary(bed$end-bed$start)
  bed <- bed %>% 
    filter(end-start<=200 , end-start>=10)
  
  # bedtoolsr::bt.getfasta(fo = paste0("tmp/",basename(x),"tmp.fa"),nameOnly = T, fi = "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/fasta/genome.fa",bed = bed, s = T)
  # fasta <- rtracklayer::import(paste0("tmp/",basename(x),"tmp.fa"),format = "fasta")
  # length(fasta)
  # #fasta2 <- as.data.frame(fasta[sample(1:length(fasta),0.1*length(fasta),replace=F)])
  # fasta <- as.data.frame(fasta)
  # ratio <- sum( grepl(sequences,fasta$x,perl=T) )/nrow(fasta)
  tmp <- data.frame(sample=basename(x),median.len=median(bed$end-bed$start),number=nrow(bed),path=x) # ratio=ratio,
  tmp.list[[x]] <- tmp
  # return(tmp)
}
sum.num.single <- function(x){
  #x <- "output/AGO2_IP/call_domain_withRepeats_dedup_bk/domains/b5_p0001_gn.bed"
  print(x)
  
  # #x <- "output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax_by_sample/interectAGO2/FTA-12_1.bed"
  rbp <- data.table::fread(paste0(x),data.table = F, header = F,sep = "\t",stringsAsFactors = F)
  
  bed <- as_tibble(rbp)
  bed <- bed[,1:6]
  colnames(bed) <- c("chr","start","end","peak","score","strand")
  summary(bed$end-bed$start)
  bed <- bed %>% 
    filter(end-start<=200 , end-start>=10)
  
  # bedtoolsr::bt.getfasta(fo = paste0("tmp/",basename(x),"tmp.fa"),nameOnly = T, fi = "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/fasta/genome.fa",bed = bed, s = T)
  # fasta <- rtracklayer::import(paste0("tmp/",basename(x),"tmp.fa"),format = "fasta")
  # length(fasta)
  # #fasta2 <- as.data.frame(fasta[sample(1:length(fasta),0.1*length(fasta),replace=F)])
  # fasta <- as.data.frame(fasta)
  # ratio <- sum( grepl(sequences,fasta$x,perl=T) )/nrow(fasta)
  # tmp <- data.frame(sample=basename(x),median.len=median(bed$end-bed$start),number=nrow(bed),path=x) # ratio=ratio,
  bed$sample=basename(x)
  bed$median.len=median(bed$end-bed$start)
  bed$number=nrow(bed)
  bed$path=x
  tmp.list[[x]] <- bed #tmp
  # return(tmp)
}

pre <- "GSE71008_NCpool/call_domain_withRepeats_all" # "AGO2_IP/call_domain_withRepeats_dedup_bk", "GSE71008/call_domain_withRepeats_all", "GSE133684/call_domain_withRepeats_dedupByPos"  
# in1 <- Sys.glob(paste0("output/",pre,"/domains/b5_p*_11RNA.bed"))
# in2 <- Sys.glob(paste0("output/",pre,"/domains_clipper/b5_p*_11RNA.bed"))
# in3 <- Sys.glob(paste0("output/",pre,"/domains_clam/b5_p*_11RNA.bed"))
# # in4 <- Sys.glob(paste0("output/",pre,"/domains_localmax/b5_d05_p*_11RNA.bed")) # 
# in5 <- Sys.glob(paste0("output/",pre,"/domains_localmax_EM2/b5_d05_p*_11RNA.bed"))

in1 <- Sys.glob(paste0("output/",pre,"/domains/b5_p*_8DNA.bed"))
in5 <- Sys.glob(paste0("output/",pre,"/domains_localmax_EM2/b5_d05_p*_8DNA.bed"))

### num plot

### num plot (by RNA types)
# res2.df <- do.call(rbind,mclapply(c(in1,in2,in3,  in5), sum.num.single, mc.cores = 5) )  # 11RNA
res2.df <- do.call(rbind,mclapply(c(in1,  in5), sum.num.single, mc.cores = 2) )  # 8DNA

res2.df$method <- "Piranha"
# res2.df$method[grepl("localmax",res2.df$path)] <- "LocalMax"
res2.df$method[grepl("EM2",res2.df$path)] <- "exPeak"
# res2.df$method[grepl("clipper",res2.df$path)] <- "CLIPper"
# res2.df$method[grepl("clam",res2.df$path)] <- "CLAM"
# res2.df$method <- factor(res2.df$method,levels = c("Piranha","CLIPper","CLAM","exPeak"))
res2.df$method <- factor(res2.df$method,levels = c("Piranha","exPeak"))
res2.df <- res2.df[res2.df$method!="Piranha",]

#add RNA types
res2.df$RNA <- ref$transcript_type[match(res2.df$chr,ref$transcript_id)]
res2.df$RNA <- gsub("_rev|\\.rev|_for|\\.for","",res2.df$RNA,perl = T)
# rna <- c("rRNA","pri_miRNA","lncRNA", "mRNA", "piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA')
dna <- c("intron","promoter", "enhancer","repeats") #
# res2.df$RNA <- factor(res2.df$RNA,levels = c(rna))
res2.df$RNA <- factor(res2.df$RNA,levels = c(dna))
table(res2.df$RNA)
res2.df <- dplyr::as_tibble(res2.df)
res2.df2 <- res2.df %>% 
  dplyr::group_by(sample,method,RNA) %>%  # p.value
  dplyr::summarize(number=dplyr::n_distinct(peak))
res2.df <- as.data.frame(res2.df2)


p11 <- ggplot(res2.df, aes(x=method, y=number, fill=RNA)) +
  geom_bar( stat = "identity", position=position_dodge(width=0.7),color="black", width=0.5)+
  # scale_fill_nejm_adaptive()+
  scale_fill_manual(values = RColorBrewer::brewer.pal(8,"Greys")[2:5]) +
  # scale_color_d3()+
  # scale_x_discrete(label=c("Domain Ratio","MIR Ratio"))+
  labs(title="",x="", y = "Peak number")+
  # facet_grid(method~.)+
  #ylim(c(0,1))+
  # xlim(c(0,1))+
  theme_classic(base_size=12) + 
  theme(#aspect.ratio = 1, #axis.ticks.x=element_blank(),  
    #strip.text.y = element_blank(),
    #strip.text.x = element_text(face="bold",family="arial",size=20),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size = 20,hjust = 1,vjust = 0.5,angle = 90), # element_blank(), #
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size=20),
    strip.text = element_text(size = 20),
    legend.position = "right",# c(0.85,0.7),#
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16))
# ggsave("peak_num_stat.pdf")




### plot peak length boxplot
bed.list <- list()
sum.len <- function(x){
  # x <- "output/AGO2_IP/call_domain/domains/5/50_gn_intersectAGO2.bed"
  #x <- "output/AGO2_IP/call_domain/domains/5/50_gn.bed"
  print(x)
  
  setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/")
  options(stringsAsFactors = F)
  options(bedtools.path = "/BioII/lulab_b/baopengfei/biosoft") # /BioII/lulab_b/baopengfei/anaconda3/envs/exvariance/bin
  suppressPackageStartupMessages(library(bedtoolsr))
  library(tidyverse)
  library(dplyr)
  
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
res2.df2 <- do.call(rbind,mclapply(c(in1, in5), sum.len, mc.cores = 2) )  # 8DNA

res2.df2$method <- "Piranha"
# res2.df2$method[grepl("localmax",res2.df2$path)] <- "LocalMax"
res2.df2$method[grepl("EM2",res2.df2$path)] <- "exPeak"
# res2.df2$method[grepl("clipper",res2.df2$path)] <- "CLIPper"
# res2.df2$method[grepl("clam",res2.df2$path)] <- "CLAM"
# res2.df2$method <- factor(res2.df2$method,levels = c("Piranha","CLIPper","CLAM","exPeak"))
res2.df2$method <- factor(res2.df2$method,levels = c("Piranha","exPeak"))
res2.df2$width <- res2.df2$end-res2.df2$start
res2.df2 <- res2.df2[res2.df2$method!="Piranha",]

### plot peak length boxplot (by RNA types)
#table(res2.df2$strand)
res2.df2$RNA <- ref$transcript_type[match(res2.df2$chr,ref$transcript_id)]
res2.df2$RNA <- gsub("_rev|\\.rev|_for|\\.for","",res2.df2$RNA,perl = T)
# rna <- c("rRNA","pri_miRNA","lncRNA", "mRNA", "piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA')
dna <- c("intron","promoter", "enhancer","repeats") #
# res2.df2$RNA <- factor(res2.df2$RNA,levels = c(rna))
res2.df2$RNA <- factor(res2.df2$RNA,levels = c(dna))
table(res2.df2$RNA)

p22 <- ggplot(res2.df2, aes(x=method, y=width, fill=RNA)) +
  # geom_bar(position = "dodge", stat = "identity")+
  geom_boxplot(outlier.size = .05,outlier.colour = "grey", outlier.alpha = 0.9, alpha=0.9, color="black",  position=position_dodge(width=0.7), width=0.5)+
  # geom_jitter()+
  # geom_histogram(aes(y=..density..), position="identity", alpha=0.5)+
  # geom_density(alpha=0.1)+
  # geom_vline(data=mu, aes(xintercept=grp.mean, color=method),linetype="dashed")+
  # scale_color_manual(values=c("firebrick", "steelblue4", "grey30"))+
  # scale_fill_manual(values=c("firebrick", "steelblue4", "grey30"))+
  # geom_vline(xintercept = 0,color="grey30",linetype="dashed")+
  # geom_label(data=mu, aes(label=ratio, color=method, fill=method )) +
  #annotate(x = 0.5, y=0.5)+
  # geom_text(data=mu, nudge_x = 0.5, nudge_y = 0.5, aes(text=ratio, color=method))+
  # scale_fill_nejm_adaptive()+
  scale_fill_manual(values = RColorBrewer::brewer.pal(8,"Greys")[2:5]) +
  
  # scale_color_nejm_adaptive()+
  # scale_x_discrete(label=c("Domain Ratio","MIR Ratio"))+
  labs(title="",x="", y = "Peak length")+
  # facet_grid(method~.)+
  ylim(c(0,50))+
  # xlim(c(0,1))+
  theme_classic(base_size=12) + 
  theme(#aspect.ratio = 1, #axis.ticks.x=element_blank(),  
    #strip.text.y = element_blank(),
    #strip.text.x = element_text(face="bold",family="arial",size=20),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size = 20,hjust = 1,vjust = 0.5,angle = 90), # 
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size=20),
    strip.text = element_text(size = 20),
    legend.position = "right", #c(0.9,0.8),#,# 
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16))


cowplot::plot_grid(plotlist = list(p11,p22),rel_widths = c(1,1), nrow = 1, align = "hv", axis = "b")
# ggsave(filename = "domain_num_len_pear_peak.pdf", width = 8, height = 10) # 11RNA
ggsave(filename = "8DNA_domain_num_len_pear_peak.pdf", width = 10, height = 6) # 8DNA


# evaluate peak num bar/pie (Fig4 RNA+DNA) -------------------------------------------------------------------------
## plot boxplot of Peak num stat
setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/")
options(stringsAsFactors = F)
options(bedtools.path = "/BioII/lulab_b/baopengfei/biosoft") # /BioII/lulab_b/baopengfei/anaconda3/envs/exvariance/bin
suppressPackageStartupMessages(library(bedtoolsr))
library(tidyverse)
library(dplyr)
library(conflicted)
conflict_prefer("filter", "dplyr")

#read ref
ref <- data.table::fread("../WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt",header = T,sep = "\t", stringsAsFactors = F)
ref <- ref[!is.na(ref$tx.length),]
ref[1:3,]
table(ref$transcript_type)


tmp.list <- list()
sum.num.single <- function(x){
  #x <- "output/AGO2_IP/call_domain_withRepeats_dedup_bk/domains/b5_p0001_gn.bed"
  print(x)
  
  # #x <- "output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax_by_sample/interectAGO2/FTA-12_1.bed"
  rbp <- data.table::fread(paste0(x),data.table = F, header = F,sep = "\t",stringsAsFactors = F)
  
  bed <- as_tibble(rbp)
  bed <- bed[,1:6]
  colnames(bed) <- c("chr","start","end","peak","score","strand")
  summary(bed$end-bed$start)
  bed <- bed %>% 
    filter(end-start<=200 , end-start>=10)
  
  # bedtoolsr::bt.getfasta(fo = paste0("tmp/",basename(x),"tmp.fa"),nameOnly = T, fi = "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/fasta/genome.fa",bed = bed, s = T)
  # fasta <- rtracklayer::import(paste0("tmp/",basename(x),"tmp.fa"),format = "fasta")
  # length(fasta)
  # #fasta2 <- as.data.frame(fasta[sample(1:length(fasta),0.1*length(fasta),replace=F)])
  # fasta <- as.data.frame(fasta)
  # ratio <- sum( grepl(sequences,fasta$x,perl=T) )/nrow(fasta)
  # tmp <- data.frame(sample=basename(x),median.len=median(bed$end-bed$start),number=nrow(bed),path=x) # ratio=ratio,
  bed$sample=basename(x)
  bed$median.len=median(bed$end-bed$start)
  bed$number=nrow(bed)
  bed$path=x
  tmp.list[[x]] <- bed #tmp
  # return(tmp)
}

### num plot (by RNA types)
pre <- "GSE71008_NCpool/call_peak_all"
in1 <- Sys.glob(paste0("output/",pre,"/expeakCNN_by_sample/b5_d50_p1/NCpool.bed"))
res2.df <- do.call(rbind,mclapply(c(in1), sum.num.single, mc.cores = 5) )
res2.df$method <- ""
res2.df$method[grepl("expeak",res2.df$path)] <- "exPeak"

#add RNA types
res2.df$RNA <- ref$transcript_type[match(res2.df$chr,ref$transcript_id)]
res2.df$RNA <- gsub("_rev|\\.rev|_for|\\.for","",res2.df$RNA,perl = T)
#mir <- "pri_miRNA"
rna <- c("rRNA","pri_miRNA","lncRNA", "mRNA", "piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA')
dna <- c("intron","promoter", "enhancer","repeats") #
res2.df$RNA <- factor(res2.df$RNA,levels = c(rna,dna)) # mir
# res2.df$RNA <- factor(res2.df$RNA,levels = c(dna))
table(res2.df$RNA)
res2.df$RNA.group <- "" # primary miRNA
res2.df$RNA.group[res2.df$RNA %in% rna] <- "annotated RNA species"
res2.df$RNA.group[res2.df$RNA %in% dna] <- "unannotated regions"
res2.df$RNA.group <- factor(res2.df$RNA.group, levels = c("annotated RNA species","unannotated regions"))
table(res2.df$RNA.group)

res2.df <- dplyr::as_tibble(res2.df)
res2.df <- res2.df %>% 
  dplyr::group_by(sample,method,RNA,RNA.group) %>%  # p.value
  dplyr::summarize(number=dplyr::n_distinct(peak))
res2.df <- as.data.frame(res2.df)



# p11 <- ggplot(res2.df, aes(x=method, y=number, color=RNA, fill=RNA)) +
#   geom_bar( stat = "identity", position = 'fill', color="black", width=0.5) + # position=position_dodge(width=0.7)
#   # scale_fill_nejm_adaptive() +
#   scale_fill_manual(values = c(RColorBrewer::brewer.pal(5,"Oranges")[4], 
#                     colorRampPalette(RColorBrewer::brewer.pal(8, "Blues"))(11)[1:10],
#                     RColorBrewer::brewer.pal(8,"Greys")[2:5])
#                     )+ 
#   # scale_color_d3()+
#   # scale_x_discrete(label=c("Domain Ratio","MIR Ratio"))+
#   labs(title="",x="", y = "Peak number ")+
#   # facet_grid(method~.)+
#   #ylim(c(0,1))+
#   # xlim(c(0,1))+
#   theme_void(base_size=12) + 
#   theme(
#     # aspect.ratio = 1,#axis.ticks.x=element_blank(),  
#     #strip.text.y = element_blank(),
#     #strip.text.x = element_text(face="bold",family="arial",size=20),
#     axis.title.x = element_text(size=20),
#     axis.title.y = element_text(size=20),
#     axis.text.x = element_text(size = 20,hjust = 1,vjust = 0.5,angle = 90), # element_blank(), #
#     axis.text.y = element_text(size = 20),
#     plot.title = element_text(size=20),
#     strip.text = element_text(size = 20),
#     legend.position = "right",# c(0.85,0.7),#
#     legend.text = element_text(size= 16),
#     legend.title= element_text(size= 16))
# 
# res2.df2 <- res2.df %>% 
#   dplyr::group_by(sample,method,RNA.group) %>%  # p.value
#   dplyr::summarize(number=sum(number))
# 
# p12 <- ggplot(res2.df2, aes(x=method, y=number, color=RNA.group, fill=RNA.group)) +
#   geom_bar( stat = "identity", position = 'fill', color="black", width=0.5) + # position=position_dodge(width=0.7)
#   # scale_fill_npg_adaptive() +
#   scale_fill_manual(values = c(RColorBrewer::brewer.pal(5,"Oranges")[4], 
#                     colorRampPalette(RColorBrewer::brewer.pal(8, "Blues"))(11)[11],
#                     RColorBrewer::brewer.pal(8,"Greys")[6])
#                     )+ 
#   # scale_color_d3()+
#   # scale_x_discrete(label=c("Domain Ratio","MIR Ratio"))+
#   labs(title="",x="", y = "Peak number ")+
#   # facet_grid(method~.)+
#   #ylim(c(0,1))+
#   # xlim(c(0,1))+
#   theme_void(base_size=12) + 
#   theme(
#     # aspect.ratio = 1,#axis.ticks.x=element_blank(),  
#         #strip.text.y = element_blank(),
#         #strip.text.x = element_text(face="bold",family="arial",size=20),
#         axis.title.x = element_text(size=20),
#         axis.title.y = element_text(size=20),
#         axis.text.x = element_text(size = 20,hjust = 1,vjust = 0.5,angle = 90), # element_blank(), #
#         axis.text.y = element_text(size = 20),
#         plot.title = element_text(size=20),
#         strip.text = element_text(size = 20),
#         legend.position = "right",# c(0.85,0.7),#
#         legend.text = element_text(size= 16),
#         legend.title= element_text(size= 16))
# # ggsave("peak_num_stat.pdf")
# # gridExtra::grid.arrange(p11,p12,ncol=2)
# cowplot::plot_grid(plotlist = list(p11,p12),rel_heights = c(1,1), nrow = 1, align = "hv", axis = "b")
# ggsave(filename = "domain_num_len_fig5.pdf", width = 16, height = 10) # 8DNA


#pie
df2 <- as.data.frame(res2.df[,c("RNA","number")])   # as.data.frame(table(peak$type))
colnames(df2) <- c("Var1","Freq")
sum(df2$Freq[1:11]) # 4171
sum(df2$Freq[12:15]) # 4377
df2$Var1 <- factor(df2$Var1,levels = c("lncRNA", "mRNA", "piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA', 'intron','promoter', 'enhancer','repeats', "rRNA","pri_miRNA"))
df2 <- df2[order(df2$Var1),]
df2$lab <- round(df2$Freq/sum(df2$Freq),digits = 3)
#peak
df2 <- df2 %>% 
  dplyr::mutate(csum = rev(cumsum(rev(Freq))), 
         pos = Freq/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Freq/2, pos))

ggplot(df2, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width=1, color="black") +
  coord_polar("y", start=0) +
  # geom_text(aes(label = lab), position = position_stack(vjust=0.5), size = 5, colour="white") +
  ggrepel::geom_label_repel(data = df2, col="white",
                            aes(y = pos, label = paste0(100*lab, "%")),
                            size = 8, nudge_x = 0.75, show.legend = FALSE) +
  labs(x = NULL, y = NULL, title = "") + # paste0("Total repeats peak: ",nrow(peak))
  # scale_color_manual("black") +
  scale_fill_manual(values = c(pal_nejm_adaptive()(15)[3:14],"#11838D",pal_nejm_adaptive()(15)[1:2])) +
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
    legend.position = "none",#c(.25,.6),
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16))
ggsave("all_peak_by_RNA_pie.pdf" ) # ,width = 9,height = 9



#



# phastCons100way conservation analysis (8DNA,seem not needed) ---------------------------
#BiocManager::install("phastCons100way.UCSC.hg38")
library(GenomicRanges)
library(phastCons100way.UCSC.hg38)
phast <- phastCons100way.UCSC.hg38
# ##计算指定区域的平均保守性得分
# gscores(phast, GRanges(seqnames="chr7", IRanges(start=117232380, width=5)))
# ##计算指定区域的保守性得分
# gscores(phast, GRanges(seqnames="chr7", IRanges(start=117232380:117232384, width=1)))
# ##计算多个区域的平均保守性得分
# gscores(phast, GRanges(c("chr7:117232380-117232384","chr2:115262390-115262395","chr3:19597000-19597005")))

peak.gn <- rtracklayer::import.bed("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008_NCpool/call_peak_all/expeak/b5_d50_p1_8DNA_gn.bed")
# #shuf 8DNA and convert to gn
# pre="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
# /usr/bin/Rscript $pre/exSeek-dev/bin/shufPeakBed.R --coord gn --RNAtype RNA \
#  -i /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008_NCpool/call_peak_all/expeak/b5_d50_p1_8DNA.bed
# in=/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008_NCpool/call_peak_all/expeak/b5_d50_p1_8DNA.bed.shuffle
# {{
#   grep -v '^chr' ${in} | $pre/exSeek-dev/bin/tbed2gbed <(cat $pre/exSeek-dev/genome/hg38/bed/{long_RNA,tRNA,pri_miRNA,piRNA,rRNA}_newTxID.bed) /dev/stdin /dev/stdout
#   awk 'BEGIN{{OFS="\t";FS="\t"}}/^chr/{{print $1,$2,$3,$4,$5,$6,0,0,0,1,$3-$2,0}}' ${in}
# }} > ${in}_gn.bed

# shuf.gn <- rtracklayer::import.bed("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008_NCpool/call_peak_all/expeak/b5_d50_p1_8DNA.bed.shuffle_gn.bed")
# flank.gn <- rtracklayer::import.bed("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008_NCpool/call_peak_all/expeak/b5_d50_p1_8DNA_gn.bed.flank")
#peak.gn
# summary(flank.gn@ranges@width)

# # filter length
# keep <- which(shuf.gn@ranges@width>=10 & shuf.gn@ranges@width<=30 & flank.gn@ranges@width>=10 & flank.gn@ranges@width<=30)
# peak.gn <- peak.gn[keep]
# shuf.gn <- shuf.gn[keep]
# flank.gn <- flank.gn[keep]


## line plot
#fixed length from peak center
exp <- 30 # nt
dwsmp <- 200
get.exp.GR <- function(GR){
  exp.gn <- GRanges("seqnames"=GR@seqnames,
                    ranges = IRanges(start = GR@ranges@start + as.integer(GR@ranges@width * 0.5)-exp,width = exp*2, names =GR$name ), 
                    strand = GR@strand)
  return(exp.gn)
}
exp.gn <- get.exp.GR(peak.gn) # shuf.gn, peak.gn
exp.gn <- exp.gn[sample(1:length(exp.gn),dwsmp,replace = F)]
exp.gn <- exp.gn[exp.gn@ranges@width==exp*2]

get.singleBase.phast <- function(GR) {
  #single GR obj as input
  chr <- GR@seqnames
  start <- GR@ranges@start
  end <- GR@ranges@start+GR@ranges@width
  s <- gscores(phast, GRanges(seqnames=chr, IRanges(start=start:(end-1), width=1)))
  return(s$default)
}
exp.gn.l <- list()
for (i in 1:length(exp.gn)){
  exp.gn.l[[i]] <- exp.gn[i]
}
#length(exp.gn.l)
l.df <- lapply(exp.gn.l, FUN = get.singleBase.phast) 
#summary(shuf.gn@ranges@width)
res <- do.call("cbind",l.df)
res <- na.omit(t(res))
# res <- t(maxmin.normalize((res))) # 0<phcons<1 no need to norm
long_refp <- reshape2::melt(res,value.name = 'signal') 
# add x position
colnames(long_refp)[1:2] <- c("region","pos")

# calculate means with CI, plot line
filnal_scaler <- long_refp %>%
  # remove na
  tidyr::drop_na() %>%
  dplyr::group_by(pos) %>% # sample,
  # mean and 95% interval confidence
  dplyr::summarise(mean_signal = mean(signal), #,trim = 0.05, na.rm=T
                   sd = sd(signal),
                   upper = Rmisc::CI(signal,ci = 0.95)[1],
                   lower = Rmisc::CI(signal,ci = 0.95)[3]
  )

# add CI
ggplot(filnal_scaler,aes(x = pos,y = mean_signal)) +
  # add 0.95 interval
  geom_ribbon(aes(ymin = lower,
                  ymax = upper), # sample
              alpha = 0.5) +
  geom_line(size = 1) + # 
  theme_classic(base_size = 16) +
  # ylim(c(0,1.2)) +
  # scale_color_nejm(name = 'Data type') +
  scale_fill_manual(name = '', values = c("#FD6905","#843692","#2C68A9")) +
  scale_color_manual(name = 'Data type', values = c("#FD6905","#843692","#2C68A9")) +
  # geom_vline(xintercept = c(up_bins,up_bins+target_bins),color="grey50",linetype="dashed")+
  # geom_vline(xintercept = c((1+up_bins+target_bins+down_bins)*0.5),color="grey50",linetype="dashed")+ # need +1 for 1-based coord
  # # x label
  # scale_x_continuous(breaks = c(1,up_bins+1,up_bins+target_bins,up_bins+target_bins+down_bins), # need +1 for 1-based coord
  #                    labels = c(paste0("-",up,"bp"),'Start','End',paste0("+",down,"bp")) ) +
  # scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1.00), limits = c(0,1),labels = c("0","0.25","0.50","0.75","1.00")) +
  xlab('') + ylab('Min-max scaled phastCons100way') +
  ylim(c(0,1))+
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
ggsave("cov_conv.pdf",width = 9,height = 6)






## barplot
dwsmp <- 500
keep <- sample(1:length(peak.gn),dwsmp,replace = F)
peak.gn.phast <- gscores(phast, peak.gn[keep])
shuf.gn.phast <- gscores(phast, shuf.gn[keep])
flank.gn.phast <- gscores(phast, flank.gn[keep])

###get RNA type
ref <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)
peak.tx <- rtracklayer::import.bed("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008_NCpool/call_peak_all/expeak/b5_d50_p1_8DNA.bed")
peak.tx$RNA <- ref$transcript_type[match(as.character(peak.tx@seqnames),ref$transcript_id)]
peak.tx$RNA <- gsub("_for|_rev","",peak.tx$RNA,perl=T)
table(peak.tx$RNA)


df <- data.frame(phast=c(peak.gn.phast$default,shuf.gn.phast$default),
                 RNA=c(peak.tx$RNA[match(peak.gn.phast$name,peak.tx$name)], peak.tx$RNA[match(unlist(sapply(strsplit(shuf.gn.phast$name,"--"),"[",1)),peak.tx$name)]),
                 type=c(rep("Peak",dwsmp),rep("Background",dwsmp)))
df$RNA <- factor(df$RNA, levels = c("intron","promoter","enhancer","repeats"))
df$type <- factor(df$type, levels = c("Peak","Background"))
table(df$RNA,df$type)
ggbarplot(data = df, x = "type", y = "phast", fill = "type", add="mean_se",position = position_dodge(),  short.panel.labs = T) + #  facet.by=c("RNA"),
  stat_compare_means(
    aes(x=type, y=phast, group=type), # ref.group="Flank",  #
    # comparisons = list(c("Peak","Flank")),
    label="p.format",  # p.signif / p.format
    # symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01,1),symbols = c( "***", "**", "*", "ns")),
    method = "wilcox.test", # method.args = list(alternative = "greater"),  # greater means ref.group less, 
    # label.x.npc = 0.8, label.y.npc = 0.8, # 11RNA 
    # step.increase = 0.08, vjust = -0.6, # 11RNA 
    label.x.npc = 0.8, label.y.npc = 0.8, # 8DNA
    bracket.size = 0.5, step.increase = 0.08, vjust = 6.6, hjust = 1, # 8DNA
    hide.ns=T,size =6, paired = F # F?
  ) + 
  # geom_bracket(
  #   xmin = "0.5", xmax = "1", y.position = 30,
  #   label = "t-test, p < 0.05"
  # ) + 
  # facet_grid(site~.,scales = "free")+ # facet_grid seem not work well with ggbarplot
  scale_fill_manual(values = c('firebrick','salmon'))+ #c("red","grey90")
  # scale_x_discrete(label=c("Peak Ratio","MIR Ratio"))+
  labs(title="",x="", y = "")+
  # xlim(c(0,1))+   
  # ylim(c(-20,0))+
  theme_minimal() +  # base_size=12
  theme(aspect.ratio = 2,
        #axis.ticks.x=element_blank(),
        #strip.text.y = element_blank(),
        strip.text = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size = 20,hjust = 1,vjust = 0.5,angle = 90), #
        axis.text.y = element_text(size = 20),
        plot.title = element_text(size=20),
        # strip.text = element_blank(),
        legend.position = "none", #c(0.9,0.8),#,#
        legend.text = element_text(size= 16),
        legend.title= element_text(size= 16))
ggsave(filename = "domain_phacons_sum.pdf", width = 3, height = 6) # 8DNA







# repeats/rmsk peak num/type (Fig4) -----------------------------------
#rmsk <- read.table("/lulabdata/baopengfei/shared_reference/hg38/backup/rmsk.txt",header = T,sep = "\t")
rmsk <- data.table::fread("/lulabdata/baopengfei/shared_reference/hg38/backup/rmsk.txt",header = F,sep = "\t")
rmsk <- rmsk[,c("V11","V12","V13")]
rmsk$type <- rmsk$V12
rmsk$type[!(rmsk$type %in% c("LINE","SINE","LTR","Retroposon","Satellite","Simple_repeat","Low_complexity"))] <- "Other"
# table(rmsk$V11) # 14572
# table(rmsk$V12) # 30
# table(rmsk$V13) # 200
table(rmsk$type) # 8

dst <- "GSE71008"
l <- list()
smps <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/",dst,"/sample_ids_NCpool.txt"))$V1
for(smp in smps){
  print(smp)
  peak <- data.table::fread(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_all/expeak_by_sample/b5_d50_p1/",smp,".bed"),header = F,sep = "\t")
  peak$id <- unlist(sapply(strsplit(peak$V1,"__chr"),"[",1))
  peak <- peak[((peak$id) %in% unique(rmsk$V11)),]
  #table(unique(peak$id) %in% unique(rmsk$V13))
  peak$type <- rmsk$type[match(peak$id,rmsk$V11)]
  peak$type <- factor(peak$type,levels = c("LINE","SINE","LTR","Satellite","Simple_repeat","Other"))
  
  ## plot pie
  df2 <- as.data.frame(table(peak$type))
  df2$smp <- smp
  l[[smp]] <- df2
}
df2 <-dplyr:: as_tibble(do.call(rbind,l)) %>% 
  dplyr::group_by(Var1) %>% 
  dplyr::summarize(Freq=mean(Freq))
df2$lab <- round(df2$Freq/sum(df2$Freq),digits = 2)
df2 <- df2 %>%
  dplyr::mutate(csum = rev(cumsum(rev(Freq))),
         pos = Freq/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Freq/2, pos))
ggplot(df2, aes(x="", y=Freq, fill=Var1)) + 
  geom_bar(stat="identity", width=1, color="black") + 
  coord_polar("y", start=0) + 
  # geom_text(aes(label = lab), position = position_stack(vjust=0.5), size = 5, colour="white") + 
  ggrepel::geom_label_repel(data = df2, col="white", 
                            aes(y = pos, label = paste0(100*lab, "%")),
                            size = 8, nudge_x = 0.75, show.legend = FALSE) +
  labs(x = NULL, y = NULL, title =paste0("Total repeats peak: ",nrow(peak))) +
  scale_color_manual("black") +
  scale_fill_manual(name="repeats",values = c("#1AC0AD",pal_d3_adaptive()(25)[21:25]) ) + 
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
ggsave("repeats_pie.pdf")







# predict snoRNA/tRNA/k-turn (Fig4)  ------------------------------------------------------------
# #exp100.bed
# #G004418__chr1___145982791____146044553_pos      52795   53032   peak_4250       39      +
# #G004418__chr1___145982791____146044553_pos      52838   53074   peak_4251       37      +
# G004418__chr1___145982791____146044553_pos  52795 53074
# 
# #pred.fa
# >peak_4250.trna1 peak_4250:106-177 (+) Glu (CTC) 72 bp Sc: 73.2
# TCCCTGGTGGTCTAGTGGTTAGGATTCGGCGCTCTCACCGCCGCGGCCCGGGTTCGATTC
# CCGGTCAGGGAA
# >peak_4251.trna1 peak_4251:63-134 (+) Glu (CTC) 72 bp Sc: 73.2
# TCCCTGGTGGTCTAGTGGTTAGGATTCGGCGCTCTCACCGCCGCGGCCCGGGTTCGATTC
# CCGGTCAGGGAA

## input peak bed (merged)
pre="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
dst <- "GSE71008_NCpool"
peak0 <- read.table(paste0(pre,"/output/",dst,"/call_peak_all/expeak_by_sample/b5_d50_p1/intersect_8DNA/NCpool.bed6.exp100.merge"),header = F)
peak <- peak0
peak$V1 <- paste0(peak$V4,"::",peak$V1,":",peak$V2,"-",peak$V3) # #peak_7460::enhancer__chr11___235400____236600_neg:25-242
peak$e <- peak$V3-peak$V2
peak$s <- 0
#rm adjacent multi record
table(grepl(",",peak$V4))
peak <- peak[!grepl(",",peak$V4),]
rownames(peak) <- peak$V4

rawpeak <- read.table(paste0(pre,"/output/",dst,"/call_peak_all/expeak_by_sample/b5_d50_p1/intersect_8DNA/NCpool.bed6"),header = F)
rawpeak <- rawpeak[rawpeak$V4 %in% peak$V4,]
rownames(rawpeak) <- rawpeak$V4
rawpeak <- rawpeak[peak$V4,]
rawpeak$V1 <- peak$V1
rawpeak$V2 <- rawpeak$V2-peak$V2
rawpeak$V3 <- rawpeak$V3-peak$V2
rawpeak[1:2,]
#rawpeak <- na.omit(rawpeak)
write.table(rawpeak,paste0(pre,"/output/",dst,"/call_peak_all/expeak_by_sample/b5_d50_p1/intersect_8DNA/NCpool.bed6.exp100.merge.rawpeak"),quote = F,sep = "\t",row.names = F,col.names = F)

peak <- peak[,c(1,8,7,4:6)]
peak[1:2,]
write.table(peak,paste0(pre,"/output/",dst,"/call_peak_all/expeak_by_sample/b5_d50_p1/intersect_8DNA/NCpool.bed6.exp100.merge.newpeak"),quote = F,sep = "\t",row.names = F,col.names = F)


##define tx2gn function
mergeTxCoor2GnCoor <- function(chr,start,end){
  # chr <- "peak_5249::G078490__chr7___52897807____52917741_neg:2670:2885"
  # chr <- "peak_8189::ERVL____E____int__chr4___93815152____93815645_pos:152:369"
  # start <- 88
  # end <- 126
  raw.chr <- unlist(sapply(strsplit(as.character(chr),":",fixed = T),"[",3))
  raw.start <- as.numeric(unlist(sapply(strsplit(chr,":",fixed = T),"[",4)))
  raw.end <- as.numeric(unlist(sapply(strsplit(chr,":",fixed = T),"[",5)))
  gn.chr <- unlist(sapply(strsplit(raw.chr,"__chr",fixed = T),"[",2))
  gn.chr <- paste0("chr",unlist(sapply(strsplit(gn.chr,"___",fixed = T),"[",1)))
  gn.start <- unlist(sapply(strsplit(raw.chr,"__chr",fixed = T),"[",2))
  gn.start <- unlist(sapply(strsplit(gn.start,"___",fixed = T),"[",2))
  gn.start <- unlist(sapply(strsplit(gn.start,"____",fixed = T),"[",1))
  gn.end <- unlist(sapply(strsplit(raw.chr,"__chr",fixed = T),"[",2))
  gn.end <- unlist(sapply(strsplit(gn.end,"____",fixed = T),"[",2))
  gn.end <- unlist(sapply(strsplit(gn.end,"_",fixed = T),"[",1))
  gn.start <- as.numeric(gn.start)
  gn.end <- as.numeric(gn.end)
  gn.strand <- ifelse(grepl("_pos",chr),"+","-")
  start <- as.numeric(start)
  end <- as.numeric(end)
  if(gn.strand=="+"){
    s <- gn.start+raw.start+start
    e <- gn.start+raw.start+end
  }else if(gn.strand=="-"){
    s <- gn.end-raw.start-end
    e <- gn.end-raw.start-start
  }
  return(paste0(gn.chr,":",s,"-",e))
}



## snoscan (lowelab.ucsc.edu/snoscan/snoscanReadme.html)
snoscan <- readLines("/BioII/lulab_b/baopengfei/gitsoft/snoscan-1.0/output_peak_merge")
head(snoscan)

snoscan <- snoscan[grepl(">> peak_",snoscan)]
snoscan.list <- strsplit(snoscan,"\t")
snoscan.df <- as.data.frame(do.call(rbind,snoscan.list))

peak.id <- unlist(sapply(strsplit(snoscan," "),"[",2))
peak.id <- peak.id[!duplicated(peak.id)] # 142 unique expanded peak
table(grepl("Gs-DpBox",snoscan.df$V1))
#2512 28333
#Gs-D box,  Gs-DpBox
snoscan.df <- snoscan.df$V1
#head(strsplit(snoscan.df," "))
snoscan.bed <- data.frame(chr=unlist(sapply(strsplit(snoscan.df," "),"[",2)),
                          s=unlist(sapply(strsplit(snoscan.df," "),"[",6)),
                          # e=unlist(sapply(strsplit(snoscan.df," "),"[",4)),
                          name="X",
                          score=as.numeric(unlist(sapply(strsplit(snoscan.df," "),"[",4))),
                          strand="+",
                          meta=snoscan.df
                          )
snoscan.bed$s <- gsub("\\(|\\)","",snoscan.bed$s,perl = T)
snoscan.bed$e <- as.numeric(unlist(sapply(strsplit(snoscan.bed$s,"-"),"[",2)))
snoscan.bed$s <- as.numeric(unlist(sapply(strsplit(snoscan.bed$s,"-"),"[",1)))
snoscan.bed <- snoscan.bed[,c("chr","s","e","name","score","strand","meta")]
#snoscan.bed6[1:3,]
snoscan.bed6 <- snoscan.bed[,1:6]
snoscan.bed6$name <- paste0(snoscan.bed6$chr,"_",snoscan.bed6$s,"_",snoscan.bed6$e)
#rm adjacent multi record
table(grepl(",",snoscan.bed6$name))
snoscan.bed6 <- snoscan.bed6[!grepl(",",snoscan.bed6$name),]
write.table(snoscan.bed6,paste0(pre,"/output/",dst,"/call_peak_all/expeak_by_sample/b5_d50_p1/intersect_8DNA/NCpool.bed6.exp100.merge.snoscan"),quote = F,sep = "\t",row.names = F,col.names = F)
#str(snoscan.bed)
table(as.numeric(snoscan.bed$s)<as.numeric(snoscan.bed$e)) # filter reverse strand ?
#why not ballence ???????????????
#FALSE  TRUE 
#27152  2623 
snoscan.bed <- snoscan.bed[snoscan.bed$s<snoscan.bed$e,]
summary(as.numeric(snoscan.bed$e)-as.numeric(snoscan.bed$s))
snoscan.bed <- snoscan.bed[order(snoscan.bed$chr,as.numeric(snoscan.bed$score),as.numeric(snoscan.bed$e)-as.numeric(snoscan.bed$s),decreasing = T),] # only keep high score & widest region for each exp peak region
snoscan.bed <- snoscan.bed[!duplicated(paste(snoscan.bed$chr)),] # 59 

#head(strsplit(snoscan.bed$meta," "))
snoscan.bed$type <- unlist(sapply(strsplit(snoscan.bed$meta," "),"[",15))
table(snoscan.bed$type) # Whether the guide region is adjacent to the D' box ("Gs-DpBox") or D box ("Gs-D box")
#peak_7473::enhancer__chr11___62853176____62855776_neg:2319-2547

#highest socre: 
#>> peak_7460::enhancer__chr11___235400____236600_neg:25-242  23.49  (9-189)  Cmpl: NR_146144_____1-Gm8698 (-)  25/3 bp  Gs-DpBox: 88 (80)  Len: 181  TS
#
#peak_7460::enhancer__chr11___235400____236600_neg:25-242: (raw peak: 125-142)
#CCTCCGGGGCGAGAGATGAGACACCAGACTAGGGAACTTCCTCTACTCGCCTCTACGTCTACCTGGGCGCGCCGGAGGCGTCGGCAAGGAggcgggggcggggcgcggaggcgggaccggggcgccgggggcgggggcggtaggcggaggcggggctggggcgccgggggcgggggtgggaggcggaggcggggccggggcgccgggggcggggCGA
#

# >> peak_7460::enhancer__chr11___235400____236600_neg:25-242  23.49  (9-189)  Cmpl: NR_146144_____1-Gm8698 (-)  25/3 bp  Gs-DpBox: 88 (80)  Len: 181 TS
# 
# C Box:  GAGAUGA   Sc: 5.51    (13-19)             C-D box dist: 164 bp
# D Box:  CGGA      Sc: 3.77
# D'Box:  CCGG      Sc: 1.30              D'box Guide Transit Sc: -0.67
# 
# No known meth site found        Guide Seq Sc: 22.52  (47.25 -5.64 -18.10 -1.00)
# 
#                        *
#   Db seq:  5'-      CCCGUCUCCGCCCCCCGGCCCCGCG -3'  NR_146144_____1,45S pre-ribosomal N2        (8695-8718)
#                  ||||:|||||| ||||| |||||| |||
#   Qry seq: 3'- GGCCAGGGCGGAGGCGCGGGGCGGGGGCGG -5'  peak_7460::enhancer__chr11___235400____236600_neg:25-242       (115-88):27
# 
# 
# C Box-> Guide Seq Gap Sc: -8.08 (68 bp)         Guide Seq-> D Box Gap Sc: -2.16 (63 bp)
# 
# Terminal stem:            +-[C Box] -N-GCGGGGCC - 5'            Stem Sc: 2.42 (4 bp)
#                           |             | |  ||
#                           +---[D Box] - GGCGGGGC - 3'           Stem Transit Sc: -1.11
# 
# 
# >Summary      [ C Box ] --         -- [ Cmpl/ Mism ]  X [D'Bx] --       -- [D Bx]  Length
# >Meth Gm8698  [GAGAUGA] --  68 bp  -- [  25 / 3    ]  1 [CCGG] -- 63 bp -- [CGGA]  181 bp
# >Sc    23.49  [  5.51 ] --  -8.08  -- [ 47.25 bits ]    [1.30]    -2.16    [3.77]
# 
# Seq: >peak_7460::enhancer__chr11___235400____236600_neg:25-242  23.49  (9-189)  Cmpl: NR_146144_____1-Gm8698  Len: 181
# Seq: GCGAGAGAUGAGACACCAGACUAGGGAACUUCCUCUACUCGCCUCUACGUCUACCUGGGCGCGCCGGAGGCGUCGGCAAGGAGGCGGGGGCGGGGCGCGGAGGCGGGACCGGGGCGCCGGGGGCGGGGGCGGUAGGCGGAGGCGGGGCUGGGGCGCCGGGGGCGGGGGUGGGAGGCGGAGG
# 125, 142
# 125-9=116
# 142-125=17
#1-116:grey 117-134:salmon 135-181:grey
#guide region: 81-108


#final best eg（by search UGAUGA box and select pos strand candidate that spans exp region center ）
# >> peak_7473::enhancer__chr11___62853176____62855776_neg:2319-2547  14.33  (66-129)  Cmpl: NR_146151_____1-Cm2494 (-)  10/1 bp  Gs-DpBox: 85 (20)  L
# 
# C Box:  GUGAUGA   Sc: 10.76    (70-76)            C-D box dist: 47 bp
# D Box:  CUGA      Sc: 8.05
# D'Box:  CUGA      Sc: 7.34              D'box Guide Transit Sc: -0.67
# 
# No known meth site found        Guide Seq Sc: -1.90  (19.67 -1.12 -18.10 -2.35)
# 
#                                        *
#   Db seq:  5'-                    GACUGCGGCGG -3'  NR_146151_____1,45S pre-ribosomal N3        (2489-2501)
#                                   ||||| |||||
#   Qry seq: 3'-                 AGUCUGACCCCGCC -5'  peak_7473::enhancer__chr11___62853176____62855776_neg:2319-2547        (95-85)
# 
# 
# C Box-> Guide Seq Gap Sc: -1.59 (8 bp)          Guide Seq-> D Box Gap Sc: -8.08 (25 bp)
# 
# Terminal stem:            +-[C Box] -N- CCACUCGU - 5'           Stem Sc: 1.53 (4 bp)
#                           |             | || | :
#                           +---[D Box] - GCUGUGAG - 3'           Stem Transit Sc: -1.11
# 
# 
# Seq: >peak_7473::enhancer__chr11___62853176____62855776_neg:2319-2547  14.33  (66-129)  Cmpl: NR_146151_____1-Cm2494  Len: 64
# Seq: ACCAGUGAUGAGUUGAAUACCGCCCCAGUCUGAUCAAUGUGUGACUGAAAGGUAUUUUCUGAGC
# 
# summary(snoscan.bed$e-snoscan.bed$s)
# #Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# #44.0    81.5    94.0    99.3   117.0   180.0
# 2419    2447
# 100-66
# peak:34-62
# chr11:62853176-62855776
# 2319+66=2385
# 2319+129=2448
# 62855776-2448
# snoRNA gn: chr11:62853328-62853391
# snoRNA gn: chr11:62853263-62853455 (1/3 in center is snoRNA)
# snoRNA tx: enhancer__chr11___62853176____62855776_neg:2321-2512


### get random region snoscan
get.snoscan <- function(x){
  # x <- "/BioII/lulab_b/baopengfei/gitsoft/snoscan-1.0/output_peak_merge"
  print(x)
  snoscan <- readLines(x)
  # head(snoscan)
  
  snoscan <- snoscan[grepl(">> peak_",snoscan)]
  snoscan.list <- strsplit(snoscan,"\t")
  snoscan.df <- as.data.frame(do.call(rbind,snoscan.list))
  
  peak.id <- unlist(sapply(strsplit(snoscan," "),"[",2))
  peak.id <- peak.id[!duplicated(peak.id)] # 142 unique expanded peak

  snoscan.df <- snoscan.df$V1
  #head(strsplit(snoscan.df," "))
  snoscan.bed <- data.frame(chr=unlist(sapply(strsplit(snoscan.df," "),"[",2)),
                            s=unlist(sapply(strsplit(snoscan.df," "),"[",6)),
                            # e=unlist(sapply(strsplit(snoscan.df," "),"[",4)),
                            name="X",
                            score=as.numeric(unlist(sapply(strsplit(snoscan.df," "),"[",4))),
                            strand="+",
                            meta=snoscan.df
  )
  snoscan.bed$s <- gsub("\\(|\\)","",snoscan.bed$s,perl = T)
  snoscan.bed$e <- as.numeric(unlist(sapply(strsplit(snoscan.bed$s,"-"),"[",2)))
  snoscan.bed$s <- as.numeric(unlist(sapply(strsplit(snoscan.bed$s,"-"),"[",1)))
  snoscan.bed <- snoscan.bed[,c("chr","s","e","name","score","strand","meta")]
  # table(as.numeric(snoscan.bed$s)<as.numeric(snoscan.bed$e)) # filter reverse strand ?
  #why not ballence ???????????????
  #FALSE  TRUE 
  #27152  2623 
  snoscan.bed <- snoscan.bed[snoscan.bed$s<snoscan.bed$e,]
  # summary(as.numeric(snoscan.bed$e)-as.numeric(snoscan.bed$s))
  snoscan.bed <- snoscan.bed[order(snoscan.bed$chr,as.numeric(snoscan.bed$score),as.numeric(snoscan.bed$e)-as.numeric(snoscan.bed$s),decreasing = T),] # only keep high score & widest region for each exp peak region
  snoscan.bed <- snoscan.bed[!duplicated(paste(snoscan.bed$chr)),] # 59 
  return(nrow(snoscan.bed))
}

fs <- Sys.glob("/BioII/lulab_b/baopengfei/gitsoft/snoscan-1.0/shuf/output_peak_merge_*")
obser <- get.snoscan("/BioII/lulab_b/baopengfei/gitsoft/snoscan-1.0/output_peak_merge")
snoscan.l <- lapply(fs,get.snoscan)
snoscan.v <- as.data.frame(do.call(rbind,snoscan.l))
ggplot(snoscan.v,aes(x=V1))+
  # geom_histogram(binwidth = 0.1) +
  geom_density( color="grey30", fill="grey50") + # stat = "count",  # xlim(c(0,obser*1.1)) + 
  labs(x="number",title = "" ) +
  # xlim(c(0,obser*1.1)) + 
  geom_vline(xintercept = 59, linetype="dashed", color="salmon") +
  labs(x="number") +
  scale_x_continuous(breaks = c(0,59),labels=c(0,59), limits = c(0,59)) +
  theme_void() +  # base_size=12
  theme(#axis.ticks.x=element_blank(),
    #strip.text.y = element_blank(),
    aspect.ratio = 0.3,
    strip.text = element_text(size=20),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size = 20, color = c("black","salmon")),
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size=20),
    # strip.text = element_blank(),
    legend.position = "none", #c(0.9,0.8),#,#
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16))
ggsave("snoscan_bg.pdf",width = 6,height = 4)






## snoGPS (lowelab.ucsc.edu/snoGPS/snogpsReadme.html)
snoGPS <- readLines("/BioII/lulab_b/baopengfei/gitsoft/snoGPS-0.2/snoGPS_hits_merge.fa")
snoGPS <- snoGPS[grepl(">peak_",snoGPS)]
head(snoGPS)

snoGPS.list <- strsplit(snoGPS,"\t")
snoGPS.df <- as.data.frame(do.call(rbind,snoGPS.list))
str(snoGPS.df)
snoGPS.df$V2 <- as.numeric(snoGPS.df$V2)
snoGPS.df[1:3,]
table(duplicated(snoGPS.df$V1))

snoGPS.df$V3 <- gsub("\\(|\\)","",snoGPS.df$V3,perl = T)
snoGPS.df$e <- as.numeric(unlist(sapply(strsplit(snoGPS.df$V3,"-"),"[",2)))
snoGPS.df$s <- as.numeric(unlist(sapply(strsplit(snoGPS.df$V3,"-"),"[",1)))

snoGPS.df$chr <- unlist(sapply(strsplit(snoGPS.df$V1,".",fixed=T),"[",1))
table(duplicated(snoGPS.df$chr))
#FALSE  TRUE 
#164    63
summary(snoGPS.df$e- snoGPS.df$s) #

snoGPS.df <- snoGPS.df[order(snoGPS.df$chr,as.numeric(snoGPS.df$V2),as.numeric(snoGPS.df$e)-as.numeric(snoGPS.df$s),decreasing = T),] # only keep high score & widest region for each exp peak region
snoGPS.df <- snoGPS.df[!duplicated(paste(snoGPS.df$chr)),] # 164 
snoGPS.df$chr <- gsub(">","",snoGPS.df$chr)
snoGPS.bed6 <- snoGPS.df[,c("chr","s","e","chr","V2","V4")]

snoGPS.bed6$V4 <- "+"
snoGPS.bed6$chr.1 <- paste0(snoGPS.bed6$chr.1,"_",snoGPS.bed6$s,"_",snoGPS.bed6$e)
#rm adjacent multi record
table(grepl(",",snoGPS.bed6$chr.1))
snoGPS.bed6 <- snoGPS.bed6[!grepl(",",snoGPS.bed6$chr.1),]
snoGPS.bed6[1:3,]
write.table(snoGPS.bed6,paste0(pre,"/output/",dst,"/call_peak_all/expeak_by_sample/b5_d50_p1/intersect_8DNA/NCpool.bed6.exp100.merge.snoGPS"),quote = F,sep = "\t",row.names = F,col.names = F)
summary(snoGPS.bed6$e-snoGPS.bed6$s)
table(duplicated(snoGPS.bed6$chr.1))
#highest socre: 
# >peak_6747::enhancer__chr17___78129929____78136129_pos:2570-2786.1       37.86   (7-177)         Cmpl: NR_023363_____1.U12       Nd     Pairs: 10/2/
#   7/3/H        171 bp (W)
#AGUGUGA
# GAGTGAGTGTGAGGCCCACCAGGGGAAGTGCTCCGGTGCCCACCTGCGCCATGGGGGGGCTGGCCCAGGGCCCAGAGCGTGGCGACAACGCTGGGCGTGGTCCTGCCGTGCAGGCCCCGGGGCTCTCTCTCCCTGACCTCGCCTTGTGTGGGGCACGCCTTTGGCACATCC
# #  Gap Len: 48   Sc: -2.56       Compl Len Sc: -2.74     IntStem Len Sc: -2.21   IntStem St Sc: -0.47
# #  IntStem Dither: 0     Sc: 0.00 ExtStem Dither: 1      Sc: 0.00        HACA-RCstart Int: 13    Sc: -3.12
# # HACA: AGAGCG   (74) Sc: 1.78
# #    H: AGAGCG   (74) Sc: 1.78   XSTEM2-H Int: 16        Sc: -1.27
# #   HACA-ISTEMRend Int: 13       POST_ACA: TCAGCCCCCTGC  Sc: -2.69
# #
# #  LC Sc: 4.46 (6)       RC Sc: 11.06 (61)       IntStem Sc: 9.81 (13,52)        ExtStem Sc: 3.01 (1,70)
# #  LCompl: AGUGUGA-3'    RCompl: UGGCC-3'        IntStem: 5'GGCCCACCA--+\        ExtStem: 5'GAGU
# #    rRNA: UCCCACC-5'          : ACCGG-5'                   CGGGGGGGU--+/        ExtStem: 3'CCCG
# #stem1 seq
# GAGUGAGUGUGAGGCCCACCAGGGGAAGUGCUCCGGUGCCCACCUGCGCCAUGGGGGGGCUGGCCCAGGGCCCAGAGCG
# XXXX LLLLLLLIIIIIIIII                              IIIIIIIIIRRRRR    XXXXHCHCHC
# 
# #    H: AGAGCG   (74) Sc: 1.78   Stem1-Score: 19.01      Stem2 Gap Len: 46       Sc: -2.56
# #  ACA: ACAUCC   (166) Sc: 4.64  HACA2-ISTEM2Rend Int: 14        Sc: -1.12
# #
# # IntStem2 Sc: 11.82 (106,143)   ExtStem2 Sc: 10.02 (90,156)
# # IntStem2: 5'CCGUGCAGG--+\      ExtStem2: 5'GCUGGGCGUG
# #             GGUGUGUUC--+/      ExtStem2: 3'CGGUUUCCGC
# #stem2 seq
# GCUGGGCGUGGUCCUGCCGUGCAGGCCCCGGGGCUCUCUCUCCCUGACCUCGCCUUGUGUGGGGCACGCCUUUGGCACAUCC
# XXXXXXXXXX      IIIIIIIII                            IIIIIIIII    XXXXXXXXXXCCCCCC
# 
#
#


#final best eg
# >peak_4537::G044704__chr2___88897784____88968045_pos:40801-41017.1       19.25   (32-174)        Cmpl: NR_023363_____1.U12       Nd     Pairs: 6/0/4/3/H         143 bp (W)
# TCCTGGAGCCTCAGTGGGGTGGCCAACAGTGAGCCCCAAAGGCCTGGTTGGGGCCACAGAACAATGTGACTCTGGCTTCTTAGAGTCAACAGAACAGGAGAATTCTTACATCCAAGTGTCCTGCCTTAAACAATTGCACAAAC
# #  Gap Len: 38   Sc: -3.22       Compl Len Sc: -4.97     IntStem Len Sc: -3.22   IntStem St Sc: -0.47
# #  IntStem Dither: 0     Sc: 0.00 ExtStem Dither: 7      Sc: 0.00        HACA-RCstart Int: 13    Sc: -3.12
# # HACA: AGAACA   (58) Sc: 2.96
# #    H: AGAACA   (58) Sc: 2.96   XSTEM2-H Int: 17        Sc: -1.27
# #   HACA-ISTEMRend Int: 14       POST_ACA: ATAATTAGTAGG  Sc: -1.30
# #
# #  LC Sc: 6.38 (4)       RC Sc: 6.38 (45)        IntStem Sc: 7.40 (8,40)         ExtStem Sc: 4.62 (1,51)
# #  LCompl: UGG-3'        RCompl: UGG-3'  IntStem: 5'GCCU--+\     ExtStem: 5'UCC
# #    rRNA: ACC-5'              : ACC-5'             CGGA--+/     ExtStem: 3'GGG
# #stem1 seq
# UCCUGGAGCCUCAGUGGGGUGGCCAACAGUGAGCCCCAAAGGCCUGGUUGGGGCCACAGAACA
# XXXLLL IIII                            IIII RRR   XXX    HCHCHC
# 
# #    H: AGAACA   (58) Sc: 2.96   Stem1-Score: 12.74      Stem2 Gap Len: 34       Sc: -3.22
# #  ACA: ACAAAC   (138) Sc: 0.86  HACA2-ISTEM2Rend Int: 14        Sc: -1.12
# #
# # IntStem2 Sc: 8.54 (90,117)     ExtStem2 Sc: 4.02 (75,136)
# # IntStem2: 5'CAGAACA--+\        ExtStem2: 5'GC
# #             GUCCUGU--+/        ExtStem2: 3'CG
# #stem2 seq
# GCUUCUUAGAGUCAACAGAACAGGAGAAUUCUUACAUCCAAGUGUCCUGCCUUAAACAAUUGCACAAAC
# XX             IIIIIII                    IIIIIII            XXCCCCCC
# 

### get random region snoGPS 
get.snoGPS <- function(x){
  # x <- "/BioII/lulab_b/baopengfei/gitsoft/snoGPS-0.2/snoGPS_hits_merge.fa"
  print(x)
  snoGPS <- readLines(x)
  snoGPS <- snoGPS[grepl(">peak_",snoGPS)]
  # head(snoGPS)
  
  snoGPS.list <- strsplit(snoGPS,"\t")
  snoGPS.df <- as.data.frame(do.call(rbind,snoGPS.list))
  # str(snoGPS.df)
  snoGPS.df$V2 <- as.numeric(snoGPS.df$V2)
  # snoGPS.df[1:3,]
  # table(duplicated(snoGPS.df$V1))
  
  snoGPS.df$V3 <- gsub("\\(|\\)","",snoGPS.df$V3,perl = T)
  snoGPS.df$e <- as.numeric(unlist(sapply(strsplit(snoGPS.df$V3,"-"),"[",2)))
  snoGPS.df$s <- as.numeric(unlist(sapply(strsplit(snoGPS.df$V3,"-"),"[",1)))
  
  snoGPS.df$chr <- unlist(sapply(strsplit(snoGPS.df$V1,".",fixed=T),"[",1))
  # table(duplicated(snoGPS.df$chr))
  #FALSE  TRUE 
  #164    63
  # summary(snoGPS.df$e- snoGPS.df$s) #
  
  snoGPS.df <- snoGPS.df[order(snoGPS.df$chr,as.numeric(snoGPS.df$V2),as.numeric(snoGPS.df$e)-as.numeric(snoGPS.df$s),decreasing = T),] # only keep high score & widest region for each exp peak region
  snoGPS.df <- snoGPS.df[!duplicated(paste(snoGPS.df$chr)),] # 164 
  return(nrow(snoGPS.df))
}

obser <- get.snoGPS("/BioII/lulab_b/baopengfei/gitsoft/snoGPS-0.2/snoGPS_hits_merge.fa")
fs <- Sys.glob("/BioII/lulab_b/baopengfei/gitsoft/snoGPS-0.2/shuf/snoGPS_hits_merge_*.fa")
#get.snoscan("/BioII/lulab_b/baopengfei/gitsoft/snoscan-1.0/output_peak_merge")
snoGPS.l <- lapply(fs,get.snoGPS)
snoGPS.v <- as.data.frame(do.call(rbind,snoGPS.l))

# ggplot(snoGPS.v,aes(x=V1))+
#   # geom_histogram(binwidth = 0.1) +
#   geom_density( color="grey30", fill="grey50") + # stat = "count",  # xlim(c(0,obser*1.1)) + 
#   labs(x="number",title = obser ) +
#   scale_x_continuous(breaks = c(100,max(snoGPS.v$V1))) + 
#   theme_void() +  # base_size=12
#   theme(#axis.ticks.x=element_blank(),
#     #strip.text.y = element_blank(),
#     aspect.ratio = 0.3,
#     strip.text = element_text(size=20),
#     axis.title.x = element_text(size=20),
#     axis.title.y = element_text(size=20),
#     axis.text.x = element_text(size = 20),
#     axis.text.y = element_text(size = 20),
#     plot.title = element_text(size=20),
#     # strip.text = element_blank(),
#     legend.position = "none", #c(0.9,0.8),#,#
#     legend.text = element_text(size= 16),
#     legend.title= element_text(size= 16))
ggplot(snoGPS.v,aes(x=V1)) +
  # geom_histogram(binwidth = 0.1) +
  geom_density(color="grey30", fill="grey50") + # stat = "count",
  labs(x="number",y="Background freq.") +
  # geom_label(label="test",position = c(92,40)) +
  scale_x_continuous(breaks = c(0,164),labels = c(0,164), limits=c(0,164)) +  #
  geom_vline(xintercept = 164, linetype="dashed", color="salmon") +
  theme_void() +  # base_size=12
  theme(#axis.ticks.x=element_blank(),
    #strip.text.y = element_blank(),
    aspect.ratio = 0.3,
    strip.text = element_text(size=20),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20, angle=90),
    axis.text.x = element_text(size = 20, color=c("black","salmon")),
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size=20),
    # strip.text = element_blank(),
    legend.position = "none", #c(0.9,0.8),#,#
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16))
ggsave("snoGPS_bg.pdf",width = 6,height = 4)




## tRNA-SE
tRNA <- read.table("/BioII/lulab_b/baopengfei/gitsoft/tRNAscan-SE-2.0.12/peak_merge_bed",header = F)
tRNA$peak <- unlist(sapply(strsplit(tRNA$V4,".",fixed = T),"[",1))
tRNA$codon <- unlist(sapply(strsplit(tRNA$V4,"-"),"[",2))
table(duplicated(tRNA$peak)) # all False, no need for rmdup, 93 left
table(tRNA$V6) # 1 minus， pseudogene: peak_7691::enhancer__chr2___202619187____202620387_neg:321-544
tRNA <- tRNA[tRNA$V6=="+",]
tRNA.bed6 <- tRNA[,1:6]
tRNA.bed6 <- tRNA.bed6[order(tRNA.bed6$V5),]
#rm adjacent multi record
table(grepl(",",tRNA.bed6$V4))
tRNA.bed6 <- tRNA.bed6[!grepl(",",tRNA.bed6$V4),]
write.table(tRNA.bed6,paste0(pre,"/output/",dst,"/call_peak_all/expeak_by_sample/b5_d50_p1/intersect_8DNA/NCpool.bed6.exp100.merge.tRNAscan"),quote = F,sep = "\t",row.names = F,col.names = F)
#stats:
# tRNAs decoding Standard 20 AA:              88
# Selenocysteine tRNAs (TCA):                 0
# Possible suppressor tRNAs (CTA,TTA,TCA):    0
# tRNAs with undetermined/unknown isotypes:   0
# Predicted pseudogenes:                      5
#   Total tRNAs:                                93

#best standard eg:
#>peak_7856.trna1 peak_7856:105-177 (+) Val (AAC) 73 bp Sc: 79.4 ， more like tRNA 
#peak_7856::enhancer__chr6___27752463____27755463_neg:1887-2108,104,177
#rawpeak: peak_7856::enhancer__chr6___27752463____27755463_neg:1887-2108,100,121,peak_7856
#exppeak: peak_7856::enhancer__chr6___27752463____27755463_neg:1887-2108,0,221,peak_7856
#1-17:salmon 18-73:grey
# chr6:27752463-27755463
# tRNA gn: chr6:27753399-27753472
# tRNA exp gn: chr6:27753326-27753545 (1/3 center is tRNA)
# tRNA exp tx: enhancer__chr6___27752463____27755463_neg:1918-2137
# 1887+104=1991
# 1887+177=2064
# 27755463-2064
# 27755463-1991
# 1991-73=1918
# 2064+73=2137
# 27753472+73

#best pseudo eg:
# >peak_7691::enhancer__chr2___202619187____202620387_neg:321-544.trna1
# peak_7691::enhancer__chr2___202619187____202620387_neg:321-544:91-156 (-) Tyr (GTA) 66 bp Sc: 24.1 Possible pseudogene
# GGTAAAATGGCTGAGCAAGCATTAGACTGTAAATCTAAAGACAGAGGTTAAGGCCTCTTTTTACCA
# chr2:202619187-202620387
#gn: chr2:202619910-202619976
#exp gn: chr2:202619844-202620042
#exp tx: enhancer__chr2___202619187____202620387_neg:346-544
# 321+91=412
# 321+156=477
# 202620387-477=202619910
# 202620387-412=202619975
# 202619910-66=202619844
# 202619976+66=202620042
# 412-66=346
# 477+66=543


### get random tRNA-SE
get.tRNA <- function(x){
  # x <- "/BioII/lulab_b/baopengfei/gitsoft/tRNAscan-SE-2.0.12/peak_merge_bed"
  #x <- "/BioII/lulab_b/baopengfei/gitsoft/tRNAscan-SE-2.0.12/shuf/peak_merge_bed_19"
  # tryCatch(
  #   expr = {
  #     tRNA <- read.table(x,header = F)
  #     tRNA$peak <- unlist(sapply(strsplit(tRNA$V4,".",fixed = T),"[",1))
  #     tRNA$codon <- unlist(sapply(strsplit(tRNA$V4,"-"),"[",2))
  #     # table(duplicated(tRNA$peak)) # all False, no need for rmdup, 93 left
  #     # table(tRNA$V6) # 1 minus， pseudogene: peak_7691::enhancer__chr2___202619187____202620387_neg:321-544
  #     tRNA <- tRNA[tRNA$V6=="+",]
  #     # return(nrow(tRNA))
  #   },
  #   error = {
  #     message("null records")
  #     tRNA <- data.frame()
  #     # return(0)
  #   }
  # )
  if (file.info(x)$size==0){
    return(0)
  }else{
    tRNA <- read.table(x,header = F)
    tRNA$peak <- unlist(sapply(strsplit(tRNA$V4,".",fixed = T),"[",1))
    tRNA$codon <- unlist(sapply(strsplit(tRNA$V4,"-"),"[",2))
    # table(duplicated(tRNA$peak)) # all False, no need for rmdup, 93 left
    # table(tRNA$V6) # 1 minus， pseudogene: peak_7691::enhancer__chr2___202619187____202620387_neg:321-544
    tRNA <- tRNA[tRNA$V6=="+",]
    return(nrow(tRNA))
  }
}

obser <- get.tRNA("/BioII/lulab_b/baopengfei/gitsoft/tRNAscan-SE-2.0.12/peak_merge_bed")
fs <- Sys.glob("/BioII/lulab_b/baopengfei/gitsoft/tRNAscan-SE-2.0.12/shuf/peak_merge_bed_*")
tRNA.l <- lapply(fs,get.tRNA)
tRNA.v <- as.data.frame(do.call(rbind,tRNA.l))
#hist(tRNA.v$V1,breaks = 100)

# ggplot(tRNA.v,aes(x=V1))+
#   geom_den(stat = "identity")+
#   xlim(c(0,obser*1.1))
# hist(tRNA.v$V1,breaks = 3, 
ggplot(tRNA.v,aes(x=V1)) +
  # geom_histogram(binwidth = 0.1) +
  geom_density( color="grey30", fill="grey50") + # stat = "count",
  labs(x="number",y="Background freq.") +
  # geom_label(label="test",position = c(92,40)) +
  scale_x_continuous(breaks = c(0,92), labels = c(0,92), limits=c(0,92)) +  # 
  geom_vline(xintercept = 92, linetype="dashed", color="salmon") +
  theme_void() +  # base_size=12
  theme(#axis.ticks.x=element_blank(),
    #strip.text.y = element_blank(),
    aspect.ratio = 0.3,
    strip.text = element_text(size=20),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20, angle=90),
    axis.text.x = element_text(size = 20, color=c("black","salmon")),
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size=20),
    # strip.text = element_blank(),
    legend.position = "none", #c(0.9,0.8),#,#
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16))
ggsave("tRNA_bg.pdf",width = 6,height = 4)




## pie plot
df2 <- data.frame(Var1=c("H/ACA box","C/D box"),Freq=c(164,59))
df2$Var1 <- factor(df2$Var1, levels = c("H/ACA box","C/D box")) # need factor, or label wrong position
df2$lab <- round(df2$Freq/sum(df2$Freq),digits = 2)
df2 <- df2 %>%
  dplyr::mutate(csum = rev(cumsum(rev(Freq))),
                pos = Freq/2 + lead(csum, 1),
                pos = if_else(is.na(pos), Freq/2, pos))
ggplot(df2, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  # geom_text(aes(label = lab), position = position_stack(vjust=0.5), size = 5, colour="white") +
  ggrepel::geom_label_repel(data = df2, col="white",
                            aes(y = pos, label = paste0(100*lab, "%")),
                            size = 8, nudge_x = 0.75, show.legend = FALSE) +
  labs(x = NULL, y = NULL, title =paste0("Total:",sum(df2$Freq))) +
  scale_color_manual("black") +
  scale_fill_manual(name="snoRNA",values = c(pal_nejm_adaptive()(10)[c(2,6)]) ) +  # "#11838D",
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
ggsave("sno_pie.pdf",width = 6,height = 4)



df2 <- data.frame(Var1=c("standard","pseudo"),Freq=c(88,4)) 
df2$Var1 <- factor(df2$Var1, levels = c("standard","pseudo")) # need factor, or label wrong position
df2$lab <- round(df2$Freq/sum(df2$Freq),digits = 2)
df2 <- df2 %>%
  dplyr::mutate(csum = rev(cumsum(rev(Freq))),
                pos = Freq/2 + lead(csum, 1),
                pos = if_else(is.na(pos), Freq/2, pos))
ggplot(df2, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  # geom_text(aes(label = lab), position = position_stack(vjust=0.5), size = 5, colour="white") +
  ggrepel::geom_label_repel(data = df2, col="white",
                            aes(y = pos, label = paste0(100*lab, "%")),
                            size = 8, nudge_x = 0.75, show.legend = FALSE) +
  labs(x = NULL, y = NULL, title =paste0("Total:",sum(df2$Freq))) +
  scale_color_manual("black") +
  scale_fill_manual(name="tRNA",values = c(pal_d3_adaptive()(25)[c(14,16)]) ) + 
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
ggsave("tRNA_pie.pdf",width = 6,height = 4)



## kink-turn
kt <- read.table("/BioII/lulab_b/baopengfei/gitsoft/kturnSeeker-1.1/expeak.txt",sep = "\t",header = F)
colnames(kt) <- c("seqnames","start","end","name","score","strand","type","stemLen","seqLen","5pDist","3pDist","kturnSeq","orgSeq","struc")
#https://github.com/sysu-software/kturnSeeker
kt$s1 <- toupper(unlist(sapply(strsplit(kt$struc,"stem1|stem2|pairs",perl = T),"[",2)))
kt$s1 <- gsub("5-|\\:","",kt$s1,perl = T)
kt$pair <- toupper(unlist(sapply(strsplit(kt$struc,"stem1|stem2|pairs",perl = T),"[",3)))
kt$pair <- gsub("\\||\\:|X|\\.","",kt$pair,perl = T)
kt$s2 <- toupper(unlist(sapply(strsplit(kt$struc,"stem1|stem2|pairs",perl = T),"[",4)))
kt$s2 <- gsub("3-|\\:","",kt$s2,perl = T)
#minus coord exist !?
# table(kt$start>=0)
# table(kt$end>=0)
# table(kt$end-kt$start>=0)
# kt <- kt[kt$end-kt$start>=0,]
#filter too many
kt <- kt[grepl("::G",kt$seqnames,fixed = T) & !grepl(",",kt$seqnames) & kt$start<120 & kt$end>100 ,]  # 221
hist(kt$end-kt$start)
table(kt$type)
# backward  forward 
# 113      108 

### k-turn pie
df2 <- data.frame(Var1=c("backward","forward"),Freq=c(113,108)) 
df2$Var1 <- factor(df2$Var1, levels = c("backward","forward")) # need factor, or label wrong position
df2$lab <- round(df2$Freq/sum(df2$Freq),digits = 2)
df2 <- df2 %>%
  dplyr::mutate(csum = rev(cumsum(rev(Freq))),
                pos = Freq/2 + lead(csum, 1),
                pos = if_else(is.na(pos), Freq/2, pos))
ggplot(df2, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  # geom_text(aes(label = lab), position = position_stack(vjust=0.5), size = 5, colour="white") +
  ggrepel::geom_label_repel(data = df2, col="white",
                            aes(y = pos, label = paste0(100*lab, "%")),
                            size = 8, nudge_x = 0.75, show.legend = FALSE) +
  labs(x = NULL, y = NULL, title =paste0("Total:",sum(df2$Freq))) +
  scale_color_manual("black") +
  scale_fill_manual(name="Kink-turn",values = c(pal_d3_adaptive()(25)[c(20,21)]) ) + 
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
ggsave("ktRNA_pie.pdf",width = 6,height = 4)


### get random kink-turn
get.kturn <- function(x){
  # x <- "oII/lulab_b/baopengfei/gitsoft/kturnSeeker-1.1/shuf/expeak_10.txt"
  if (file.info(x)$size==0){
    return(0)
  }else{
    kt <- read.table(x,sep = "\t",header = F)
    colnames(kt) <- c("seqnames","start","end","name","score","strand","type","stemLen","seqLen","5pDist","3pDist","kturnSeq","orgSeq","struc")
    kt <- kt[grepl("::G",kt$seqnames,fixed = T) & !grepl(",",kt$seqnames) & kt$start<120 & kt$end>100 ,] # filter
    return(nrow(kt))
  }
}

obser <- get.kturn("/BioII/lulab_b/baopengfei/gitsoft/kturnSeeker-1.1/expeak.txt")
fs <- Sys.glob("/BioII/lulab_b/baopengfei/gitsoft/kturnSeeker-1.1/shuf/expeak_*.txt")
kturn.l <- lapply(fs,get.kturn)
kturn.v <- as.data.frame(do.call(rbind,kturn.l))
#hist(kturn.v$V1,breaks = 100)

ggplot(kturn.v,aes(x=V1)) +
  # geom_histogram(binwidth = 0.1) +
  geom_density( color="grey30", fill="grey50") + # stat = "count",
  labs(x="number",y="") +
  # xlim(c(0,obser*1.1)) + 
  # geom_label(label="test",position = c(92,40)) +
  scale_x_continuous(breaks = c(0,obser), labels = c(0,obser), limits=c(0,obser*1.1)) +  # 
  geom_vline(xintercept = obser, linetype="dashed", color="salmon") +
  theme_void() +  # base_size=12
  theme(#axis.ticks.x=element_blank(),
    #strip.text.y = element_blank(),
    aspect.ratio = 0.3,
    strip.text = element_text(size=20),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20, angle=90),
    axis.text.x = element_text(size = 20, color=c("black","salmon")),
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size=20),
    # strip.text = element_blank(),
    legend.position = "none", #c(0.9,0.8),#,#
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16))
ggsave("kturn_bg.pdf",width = 6,height = 4)




#select eg:
#intron
#one peak per region
#median len: 40 nt
#span within peak
#motif:CTGA, GTAG(GATG)

kt2 <- kt[grepl("::G",kt$seqnames,fixed = T) & !grepl(",",kt$seqnames) & kt$start>80 & kt$end<140 & (kt$end-kt$start)<60 & grepl("CTGA",kt$s1) ,] #  & grepl("GTAG",kt$s2)

# #peak_4462::G033992__chr17___22520709____22523110_pos:2068:2300 (88,105) 17 nt
# #tx coord: G033992__chr17___22520709____22523110_pos:2156-2173 
# #tx expand coord: G033992__chr17___22520709____22523110_pos:2139-2190
# #gn expand coord: chr17:22522848-22522899
# 
# backward:
# stem1: 5-CTGA---GG\
# pairs:   |:xx   || GT
# stem2: 3-GGAGGAACC/
# #not correct exRNA signal!
  
# #peak_5249::G078490__chr7___52897807____52917741_neg:2670:2885 (88,126) 38 nt
# #tx coord: G078490__chr7___52897807____52917741_neg:2758-2796
# #tx expand coord: G078490__chr7___52897807____52917741_neg:2720-2834 (3fold)
# #gn expand coord: chr7:52914907-52915021 (3fold)
# forward:
# stem1: 5-tctgctgagc\ 
# pairs:   |||   xxx| aaaagacaatgtcaggtgaat 
# stem2: 3-aga---agag/
# #has ex signal !
#   



ktRNA.ref <- data.table::fread("/lulabdata/baopengfei/shared_reference/2023_natbiotech_ktRNA/ktRNA.txt",header = T,sep = "\t",check.names = F,stringsAsFactors = F)
table(ktRNA.ref$geneType,ktRNA.ref$RNAType)
#seem some derived from snoRNA (novel sncRNA ?!)

## change to gn coord
kt <- read.table("/BioII/lulab_b/baopengfei/gitsoft/kturnSeeker-1.1/expeak.txt",sep = "\t",header = F)
colnames(kt) <- c("seqnames","start","end","name","score","strand","type","stemLen","seqLen","5pDist","3pDist","kturnSeq","orgSeq","struc")
kt <- kt[,1:6]
kt <- kt[kt$end-kt$start>0,] # rm illegal records
# kt <- kt[grepl("::G",kt$seqnames,fixed = T) & !grepl(",",kt$seqnames) & kt$start<120 & kt$end>100 ,]
kt <- kt[!grepl(",",kt$seqnames),] # rm merged peaks
kt$gn <- ""
for (i in 1:nrow(kt)){
  kt$gn[i] <- mergeTxCoor2GnCoor( chr = kt$seqnames[i], start = kt$start[i], end = kt$end[i])
}

kt$gn.chr <- unlist(sapply(strsplit(kt$gn,":",fixed = T),"[",1))
kt$gn.start <- unlist(sapply(strsplit(kt$gn,"\\-|:",perl = T),"[",2))
kt$gn.end <- unlist(sapply(strsplit(kt$gn,"\\-|:",perl = T),"[",3))
kt.gn <- kt[,c("gn.chr","gn.start","gn.end","name","score")]
kt.gn$strand <- ifelse(grepl("_pos",kt$seqnames),"+","-")
tmp <- bedtoolsr::bt.intersect(a=kt.gn,b=ktRNA.ref[,1:6],wao = T) 
# only 1 record: chr11 62853345 62853367 forward_kturn-880  8  -  (SNORD31, gencode.v27 only has SNORD31B in chr13)

#best eg in fig4:
#peak_7473::enhancer__chr11___62853176____62855776_neg:2319:2547 (90-112) 22nt
# tx coord: enhancer__chr11___62853176____62855776_neg:2409-2431
# tx exp coord (5fold): enhancer__chr11___62853176____62855776_neg:2365-2475
# gn exp coord (5fold): chr11:62853301-62853411
# stem1: 5-CAGTCTGATC\ 
# pairs:   |||   xx.| AATGT 
# stem2: 3-GTC---AGTG/
#mergeTxCoor2GnCoor("peak_7473::enhancer__chr11___62853176____62855776_neg:2319:2547",90,112)
#chr11:62853345-62853367
#seem repeated with detection by snoscan (SNORD321)






# ## snoRNA/tRNA & peak cov meta across exp peak region
# library(ggplot2)
# library(tidyverse)
# library(reshape2)
# library(ggsci)
# library(Rmisc)
# #library(data.table)
# conflict_prefer("select", "dplyr")
# 
# #tx
# up=0
# down=0
# bin=1 # 10 for test
# ratio=1
# dst <- "GSE71008_NCpool" # GSE71008,GSE110381,GSE123972
# cmb <- expand.grid(signal=c(paste0(pre,"/output/",dst,"/call_peak_all/expeak_by_sample/b5_d50_p1/intersect_8DNA/NCpool.bed6.exp100.merge.tRNAscan"),
#                             paste0(pre,"/output/",dst,"/call_peak_all/expeak_by_sample/b5_d50_p1/intersect_8DNA/NCpool.bed6.exp100.merge.snoscan"),
#                             paste0(pre,"/output/",dst,"/call_peak_all/expeak_by_sample/b5_d50_p1/intersect_8DNA/NCpool.bed6.exp100.merge.snoGPS"),
#                             paste0(pre,"/output/",dst,"/call_peak_all/expeak_by_sample/b5_d50_p1/intersect_8DNA/NCpool.bed6.exp100.merge.rawpeak")
#                             ),
#                    region=c(paste0(pre,"/output/",dst,"/call_peak_all/expeak_by_sample/b5_d50_p1/intersect_8DNA/NCpool.bed6.exp100.merge.newpeak")
#                    )
# )
# cmb$signal.label <- ""
# cmb$signal.label[grepl("snoscan",cmb$signal)] <- c("snoscan") 
# cmb$signal.label[grepl("snoGPS",cmb$signal)] <- c("snoGPS") 
# cmb$signal.label[grepl("tRNAscan",cmb$signal)] <- c("tRNAscan") 
# cmb$signal.label[grepl("rawpeak",cmb$signal)] <- c("rawpeak") 
# # table(cmb$signal.label)
# cmb$region.label <- "exp"
# 
# #one way using sapply
# res.list <- lapply(1:nrow(cmb), function(j) {get.mat(signal=as.character(cmb$signal[j]), region=as.character(cmb$region[j]),
#                                                      signal.label=as.character(cmb$signal.label[j]),region.label=as.character(cmb$region.label[j]),
#                                                      up = up, down = down, bin = bin, ratio = ratio, )} )
# res.list.df <- do.call(rbind,res.list)
# 
# # plot in R
# table(res.list.df$signal)
# 
# clean_refp <- dplyr::as_tibble(res.list.df)
# #clean_refp <- refp #%>% select(c(-1:-3,-5,-6))
# colnames(clean_refp)[c(1,3)] <- c('V4','sample')
# #plot(x=1:(ncol(clean_refp)-3), y=as.numeric(clean_refp[1,4:ncol(clean_refp)]))
# table(clean_refp$region,clean_refp$sample)
# 
# 
# clean_refp <- clean_refp[base::which( apply(clean_refp[,4:ncol(clean_refp)] ,1, sd)>=0),] # filter low sd
# #~ half num of peaks have no ref overlap: sd==0
# 
# #max-min scale
# clean_refp[,4:ncol(clean_refp)] <- t(maxmin.normalize(t(clean_refp[,4:ncol(clean_refp)])))
# long_refp <- reshape2::melt(clean_refp,id.vars = c('V4','sample','region'),value.name = 'signal') %>%
#   select(-variable)
# 
# # add x position
# long_refp$pos <- rep(c(1:(ncol(clean_refp)-3)),each = nrow(clean_refp)) # ,times = nrow(cmb)
# # dim(long_refp)
# 
# #lowess model
# correct.loess <- function(long_refp){
#   lo = stats::loess(signal ~ pos, data = long_refp)
#   long_refp$signal = mean(long_refp$signal, na.rm = TRUE) * long_refp$bc/stats::predict(lo, newdata = long_refp[,c("signal","pos")])
#   long_refp$signal = round(long_refp$signal, digits = 2)
#   if (any(long_refp$signal < 0, na.rm = TRUE))
#     long_refp$signal[long_refp$signal < 0] = 0
#   return(long_refp)
# }
# get.loess <- function(long_refp){
#   lo = stats::loess(signal ~ pos, data = long_refp, span = 0.5) # default: span=0.75, higher is smoother
#   long_refp$signal = stats::predict(lo, newdata = long_refp_tmp$pos)
# 
#   long_refp$signal = round(long_refp$signal, digits = 2)
#   if (any(long_refp$signal < 0, na.rm = TRUE))
#     long_refp$signal[long_refp$signal < 0] = 0
#   return(long_refp)
# }
# 
# l <- list()
# for(smp in as.character(unique(long_refp$sample))){
#   #smp <- "RBPs"
#   print(smp)
#   for (region in as.character(unique(long_refp$region))){
#     #region <- "exPeak"
#     print(region)
#     long_refp_tmp <- long_refp[long_refp$sample==smp & long_refp$region==region,]
#     #plot(long_refp_tmp$pos, long_refp_tmp$signal)
#     long_refp_tmp2 <- get.loess(long_refp_tmp[,c("pos","signal")])
#     long_refp_tmp$signal <- long_refp_tmp2$signal
#     colnames(long_refp_tmp)[colnames(long_refp_tmp)=="signal"] <- "mean_signal"
#     l[[paste0(smp,region)]] <- long_refp_tmp
#     #table(long_refp_tmp$pos==long_refp_tmp2$pos)
#   }
# }
# filnal_scaler <- do.call(rbind,l)
# summary(filnal_scaler$mean_signal)
# 
# up_bins <- as.integer(up/bin) # ceiling
# down_bins <- as.integer(down/bin) # ceiling
# target_bins <- 20 #round((up_bins + down_bins)*(ratio/(1-ratio)),digits = 0) #as.integer((up_bins + down_bins)*(ratio/(1-ratio)))
# table(filnal_scaler$sample)
# table(filnal_scaler$region)
# #filnal_scaler$region <- factor(filnal_scaler$region,levels = c("snoscan","snoGPS","tRNAscan","rawpeak"))
# 
# #2nd max-min scale (only display mean line, no ribbon)
# filnal_scaler <- filnal_scaler %>%
#   dplyr::group_by(sample,region) %>%
#   dplyr::mutate(mean_signal_norm= (mean_signal - min(mean_signal))/(max(mean_signal)-min(mean_signal))^as.logical(sd(mean_signal)) )
# 
# ggplot(filnal_scaler,aes(x = pos,y = mean_signal_norm)) + # mean_signal,mean_signal_norm
#   geom_line(aes(color = region), size = 2) + # 
#   theme_classic(base_size = 16) +
#   scale_color_d3(name = 'Data type') +
#   scale_fill_d3(name = '') +
#   geom_vline(xintercept = c((target_bins+1)*0.5),color="grey50",linetype="dashed")+
#   # geom_vline(xintercept = c((1+up_bins+target_bins+down_bins)*0.5),color="grey50",linetype="dashed")+ # need +1 for 1-based coord
#   # x label
#   scale_x_continuous(breaks = c(1,up_bins+1,up_bins+target_bins,up_bins+target_bins+down_bins), # need +1 for 1-based coord
#                      labels = c(paste0("-",up,"nt"),'Start','End',paste0("+",down,"nt")) ) +
#   scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1.0), limits = c(0,1), labels = c("0","0.25","0.50","0.75","1.00")) +
#   xlab('') + ylab('Scaled depth') +
#   theme(aspect.ratio = 0.6,
#         strip.background = element_rect(color = NA,fill = 'grey'),
#         strip.text = element_text(size=30),
#         axis.title.x = element_text(size=20),
#         axis.title.y = element_text(size=30),
#         axis.text.x = element_blank(), #  element_text(size=20,angle = 0,hjust = 0.5)
#         axis.text.y = element_text(size=30),
#         legend.position =  "right",#c(0.9,0.8),
#         legend.text = element_text(size= 24),
#         legend.title= element_text(size= 30)) + 
#   facet_grid(#.~sample, 
#     cols = vars(sample), 
#     scales = "free_y")
# ggsave("predRNA_cov_meta.pdf",width = 20,height = 8)
# 
# filnal_scaler2 <- filnal_scaler[filnal_scaler$sample %in% c("rawpeak","snoscan","snoGPS"),]
# ggplot(filnal_scaler2,aes(x = pos,y = mean_signal_norm)) + # mean_signal,mean_signal_norm
#   geom_line(aes(color = sample), size = 2) + # 
#   theme_classic(base_size = 16) +
#   scale_color_manual(name="snoRNA",values = c("grey50","#11838D",pal_d3_adaptive()(25)[21:25]) ) + 
#   # scale_color_d3(name = 'Data type') +
#   # scale_fill_d3(name = '') +
#   geom_vline(xintercept = c((target_bins+1)*0.5),color="grey50",linetype="dashed")+
#   # geom_vline(xintercept = c((1+up_bins+target_bins+down_bins)*0.5),color="grey50",linetype="dashed")+ # need +1 for 1-based coord
#   # x label
#   scale_x_continuous(breaks = c(1,up_bins+1,up_bins+target_bins,up_bins+target_bins+down_bins), # need +1 for 1-based coord
#                      labels = c(paste0("-",up,"nt"),'Start','End',paste0("+",down,"nt")) ) +
#   scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1.0), limits = c(0,1), labels = c("0","0.25","0.50","0.75","1.00")) +
#   xlab('') + ylab('') +
#   theme(aspect.ratio = 0.2,
#         strip.background = element_rect(color = NA,fill = 'grey'),
#         strip.text = element_text(size=30),
#         axis.title.x = element_text(size=20),
#         axis.title.y = element_text(size=30),
#         axis.text.x = element_blank(), #  element_text(size=20,angle = 0,hjust = 0.5)
#         axis.text.y = element_text(size=30),
#         legend.position =  "right",#c(0.9,0.8),
#         legend.text = element_text(size= 24),
#         legend.title= element_text(size= 30)) 
# ggsave("sno_cov_meta.pdf",width = 20,height = 8)
# #
# 
# table(filnal_scaler$sample)
# filnal_scaler2 <- filnal_scaler[filnal_scaler$sample %in% c("rawpeak","tRNAscan"),]
# ggplot(filnal_scaler2,aes(x = pos,y = mean_signal_norm)) + # mean_signal,mean_signal_norm
#   geom_line(aes(color = sample), size = 2) + # 
#   theme_classic(base_size = 16) +
#   scale_color_manual(name="snoRNA",values = c("grey50","#11838D",pal_d3_adaptive()(25)[21:25]) ) + 
#   # scale_color_d3(name = 'Data type') +
#   # scale_fill_d3(name = '') +
#   geom_vline(xintercept = c((target_bins+1)*0.5),color="grey50",linetype="dashed")+
#   # geom_vline(xintercept = c((1+up_bins+target_bins+down_bins)*0.5),color="grey50",linetype="dashed")+ # need +1 for 1-based coord
#   # x label
#   scale_x_continuous(breaks = c(1,up_bins+1,up_bins+target_bins,up_bins+target_bins+down_bins), # need +1 for 1-based coord
#                      labels = c(paste0("-",up,"nt"),'Start','End',paste0("+",down,"nt")) ) +
#   scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1.0), limits = c(0,1), labels = c("0","0.25","0.50","0.75","1.00")) +
#   xlab('') + ylab('') +
#   theme(aspect.ratio = 0.2,
#         strip.background = element_rect(color = NA,fill = 'grey'),
#         strip.text = element_text(size=30),
#         axis.title.x = element_text(size=20),
#         axis.title.y = element_text(size=30),
#         axis.text.x = element_blank(), #  element_text(size=20,angle = 0,hjust = 0.5)
#         axis.text.y = element_text(size=30),
#         legend.position =  "right",#c(0.9,0.8),
#         legend.text = element_text(size= 24),
#         legend.title= element_text(size= 30)) 
# ggsave("sno_cov_meta.pdf",width = 20,height = 8)

#need further prove enrichment relative to shuffle 







# GSE110381   CPM vs. LASSO coef scatter + DCB filtering/overlap (Fig6, suppl Fig6) ----------------------------
## get tx ref
ref <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt") #,header=T
ref[1:3,1:3]
#ref2 <- NULL
#ref2 <- ref[ref$transcript_type=="tRNA",]
#summary(ref2$tx.length)


## get GSE110381_diff ML LASSO coef
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
dst <- "GSE110381_diff"
read.coef <- function(dst, type, k){ # seed, 
  # dst <- "GSE110381_diff"
  # type <- "mRNA"
  # bootstrap <- 1
  # print(paste0(type,"_",bootstrap))
  #_auroc.txt
  #pre <- "/BioII/lulab_b/wangtaiwei/peak_calling/output/GSE110381_diff/CV-20fold"
  l <- list()
  for(seed in c(10086:10090)){ # 000,111,222,333,444,555,666,777,888,999;   100,101,102,103,104,105,106,107,108,109
      # seed <- ""
      # print(seed)
      for(k in 0:9){
        #"/BioII/lulab_b/wangtaiwei/peak_calling/output/GSE110381_diff/stCV10_LR2_C0.01/enhancer_5_coef.txt"
        #k <- 0
        # print(k)
        tmp <- read.table( paste0("/BioII/lulab_b/wangtaiwei/peak_calling/output/",dst,"/stCV10_LR2_C0.01/",type,"_",k,"_seed",seed,"_coef.txt"), header = T)
        # tmp <- read.table(paste0("/BioII/lulab_b/wangtaiwei/peak_calling/output/",dst,"/CV-20fold/Round_",seed,"_",type,"_fold",k,"_coef.txt"),header = T)
        # tmp$dst <- dst
        tmp$type <- type
        tmp$k <- k
        tmp$seed <- seed
        return(tmp)   
      }
    }
}
rna <- c("rRNA","pri_miRNA","lncRNA", "mRNA", "piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA')
dna <- c("intron","promoter", "enhancer", "repeats") 
cmb <- expand.grid(type=c(rna,dna), seed=c("","_seed667","_seed668","_seed669","_seed700"),k=0:9)
res.list <- lapply(1:nrow(cmb), function(j) {read.coef(dst=dst, type=as.character(cmb$type[j]), k=as.character(cmb$k[j]) )} ) # seed=as.character(cmb$seed[j]), 
df <- do.call(rbind,res.list)
colnames(df)[1] <- "feature"
#df$bootstrap
df[1:3,]
#hist(df$coef,breaks = 10000,xlim = c(-1,1))
df.mean <- dplyr::as_tibble(df) %>% 
  dplyr::group_by(type,feature) %>% 
  dplyr::summarise(mean.coef=mean(coef,trim=0.05))
df.mean <- as.data.frame(df.mean)
rownames(df.mean) <- df.mean$feature
hist(df$coef) # 
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.071107 -0.001795  0.000053 -0.000076  0.001827  0.074676
table(R.utils::isZero(df$coef)) # near 
#seem all peak kept and with coef not equal to 0 !!!



cpm.list <- list()
diff.list <- list()



## get GSE110381 diff/CPM/RPKM stats
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
dst <- "GSE110381_diff"

### read smp table
sample.table <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/GSE110381_diff/sample_table.txt",check.names = F,header = T)
#QC
passQC <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/GSE110381_diff/sample_ids_diff_passQC2_rmBatch.txt",check.names = F,sep = "\t",header = F)$V1
sample.table <- sample.table[sample.table$sample %in% passQC,]
rownames(sample.table) <- sample.table$sample

### read count
count <- read.table(paste0(pre,"/output/",dst,"/call_peak_all/count_matrix/expeakCNN_b5_d50_p1.txt"),check.names = F,header = T)
rownames(count) <- count$feature
count <- count[,2:ncol(count)]
sample.table <- sample.table[sample.table$sample %in% colnames(count),]
sample.table$group <- gsub("adenoma","CRC",sample.table$group)
#sample.table <- sample.table[order(sample.table$group),]
positive_samples <- sample.table[sample.table$group=="CRC","sample"]
negative_samples <- sample.table[sample.table$group=="NC","sample"]
samples <- c(positive_samples, negative_samples)
sample.table <- sample.table[match(samples,sample.table$sample),]
group <- c(rep("positive",length(positive_samples)),rep("negative",length(negative_samples)))
method <- "edger_glmlrt"
norm_method <- "TMM"

# diff.res <- diff.v2.dcb(mat = count, samples = samples, group = group, method = method, norm_method = norm_method, filterType = "NULL", featureType = "domain")
# data.table::fwrite(as.data.frame(cbind("feature"=rownames(diff.res$diffTable),diff.res[["diffTable"]])),paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_all/count_matrix/expeak_b5_d50_p1.txt.diff"),quote = F,sep = "\t",row.names = F,col.names = T)
# data.table::fwrite(as.data.frame(cbind("feature"=rownames(diff.res$diffTable),diff.res[["normMat"]])),paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_all/count_matrix/expeak_b5_d50_p1.txt.rowsumLogCPM"),quote = F,sep = "\t",row.names = F,col.names = T)
# diff.res <- NULL
cpm <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",dst,"/filerSmp_exPeakCNN_smallDomain_diff_CRCvsNC.cpm"),check.names = F,header = T)
#cpm <- read.table(paste0(pre,"/output/",dst,"/call_peak_all/count_matrix/expeak_b5_d50_p1.txt.rowsumLogCPM"),check.names = F,header = T)
#cpm <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE110381_diff/call_peak_all/count_matrix/expeak_b5_d50_p1_CPM.txt",check.names = F,header = T)
cpm$feature <- rownames(cpm)
cpm <- cpm[,c(ncol(cpm),1:(ncol(cpm)-1))]
#colnames(cpm)[1] <- "feature"
#cpm[,2:ncol(cpm)] <- log2(cpm[,2:ncol(cpm)] + 1) # log2
#log2(5+1)==2.585
cpm.list[[dst]] <- cpm
#diff.list[[dst]] <- read.table(paste0(pre,"/output/",dst,"/call_peak_all/count_matrix/expeak_b5_d50_p1.txt.diff"),check.names = F,header = T)
diff.list[[dst]] <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",dst,"/filerSmp_exPeakCNN_smallDomain_diff_CRCvsNC.diff"),check.names = F,header = T)

diff.list[[dst]]$feature <- rownames(diff.list[[dst]])
diff.list[[dst]] <- merge(diff.list[[dst]],df.mean,by="feature")
# tmp <- cpm.list[[dst]]
# tmp2 <- diff.list[[dst]]
# length()

diff.list[[dst]]$cpm.median <- apply(cpm.list[[dst]][,2:ncol(cpm.list[[dst]])], 1, median ) # re-cal cpm in all group
diff.list[[dst]]$cpm.mean <- apply(cpm.list[[dst]][,2:ncol(cpm.list[[dst]])], 1, function(x) mean(x, trim=0.05) ) # re-cal cpm in all group
diff.list[[dst]]$log.abs.mean.coef <- log2(abs(diff.list[[dst]]$mean.coef))
#diff.list[[dst]]$RNA <- unlist(sapply(strsplit(diff.list[[dst]]$feature,"|",fixed=T),"[",2))
diff.list[[dst]]$peak <- unlist(sapply(strsplit(diff.list[[dst]]$feature,"|",fixed=T),"[",3))
diff.list[[dst]]$peak.id <- unlist(sapply(strsplit(diff.list[[dst]]$feature,"|",fixed=T),"[",4))
diff.list[[dst]]$tx.length <- ref$tx.length[match(diff.list[[dst]]$peak,ref$transcript_id)]
diff.list[[dst]]$peak.width <- as.numeric(unlist(sapply(strsplit(diff.list[[dst]]$feature,"|",fixed=T),"[",7)))-as.numeric(unlist(sapply(strsplit(diff.list[[dst]]$feature,"|",fixed=T),"[",6)))
hist(diff.list[[dst]]$peak.width,xlim = c(0,200),breaks = 1000)
table(diff.list[[dst]]$peak.width >=10 & diff.list[[dst]]$peak.width<=200)
diff.list[[dst]] <- diff.list[[dst]][diff.list[[dst]]$peak.width >=10 & diff.list[[dst]]$peak.width<=200,] # filter peak length
diff.list[[dst]]$type2 <- ""
diff.list[[dst]]$type2[diff.list[[dst]]$type %in% rna] <- "peak in annotated tx"
diff.list[[dst]]$type2[diff.list[[dst]]$type == "pri_miRNA"] <- "peak in primary miRNA"
diff.list[[dst]]$type2[diff.list[[dst]]$type %in% dna] <- "peak in unannotated region"
diff.list[[dst]]$type2 <- factor(diff.list[[dst]]$type2,levels = c("peak in primary miRNA","peak in annotated tx","peak in unannotated region"))
table(diff.list[[dst]]$type2)




### plot scatter 
library(ggplot2)
p0 <- ggplot(diff.list[[dst]], aes(x=log.abs.mean.coef, y=cpm.mean) ) +
  geom_bin2d(bins = 200) +
  # geom_point(aes(fill=type2),shape=21,color="white",alpha=0.1) +
  scale_fill_continuous(type = "viridis") +
  # scale_fill_manual(values = c("#ED782F","#B52B1A","#3E8C9D")) +
  # geom_smooth(method='lm', formula= y~x,color="black") +
  # ggpubr::stat_cor(label.x.npc = 0.1,size=7,label.y.npc = 0.72)+ #this means at 35th unit in the y axis, the r squared and p value will be shown
  # ggpubr::stat_regline_equation(label.x.npc = 0.1,size=7,label.y.npc = 0.79)+ #this means at 30th unit regresion line equation will be shown
  # # geom_abline(slope = 1,intercept = 0,linetype="dashed")+
  
  # geom_vline(xintercept = quantile(res$CPM.NC.EV, probs = 0.9),linetype="dashed")+
  # geom_hline(yintercept = quantile(res$CPM.NC.cf, probs = 0.9),linetype="dashed")+
  # theme_classic() +
  # theme(axis.text = element_text(size = 18,colour = "black"),
  #       axis.title = element_text(size = 24,colour = "black"),
  #       axis.line = element_blank(),
#       legend.position = "right",  legend.box = "vertical",legend.box.just = "left",legend.direction = "vertical",
#       legend.text = element_text(size = 18,colour = "black"),
#       legend.title = element_text(size = 16,face = "bold",colour = "black"))
theme_minimal() +  # base_size=12
  theme(#axis.ticks.x=element_blank(),
    #strip.text.y = element_blank(),
    aspect.ratio = 1,
    strip.text = element_text(size=20),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size = 20), #,hjust = 0,vjust = 1
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size=20),
    # strip.text = element_blank(),
    legend.position = c(.2,.6), #c(0.9,0.8),#,#
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16))
p0
# smoothScatter(nbin = 200, x=diff.list[[dst]]$log.abs.mean.coef, y=diff.list[[dst]]$cpm.mean, nrpoints = 0)

# p0 <- ggplot(diff.list[[dst]], aes(x=log.abs.mean.coef, y=cpm.mean) ) +
#   geom_bin2d(bins = 200,) +
#   scale_fill_continuous(type = "viridis") +
#   # geom_smooth(method='lm', formula= y~x,color="black") +
#   # ggpubr::stat_cor(label.x.npc = 0.1,size=7,label.y.npc = 0.72)+ #this means at 35th unit in the y axis, the r squared and p value will be shown
#   # ggpubr::stat_regline_equation(label.x.npc = 0.1,size=7,label.y.npc = 0.79)+ #this means at 30th unit regresion line equation will be shown
#   # # geom_abline(slope = 1,intercept = 0,linetype="dashed")+
#   
#   # geom_vline(xintercept = quantile(res$CPM.NC.EV, probs = 0.9),linetype="dashed")+
#   # geom_hline(yintercept = quantile(res$CPM.NC.cf, probs = 0.9),linetype="dashed")+
#   # theme_classic() +
#   # theme(axis.text = element_text(size = 18,colour = "black"),
#   #       axis.title = element_text(size = 24,colour = "black"),
#   #       axis.line = element_blank(),
#   #       legend.position = "right",  legend.box = "vertical",legend.box.just = "left",legend.direction = "vertical",
#   #       legend.text = element_text(size = 18,colour = "black"),
#   #       legend.title = element_text(size = 16,face = "bold",colour = "black"))
#   theme_minimal() +  # base_size=12
#   theme(#axis.ticks.x=element_blank(),
#     #strip.text.y = element_blank(),
#     aspect.ratio = 1,
#     strip.text = element_text(size=20),
#     axis.title.x = element_text(size=20),
#     axis.title.y = element_text(size=20),
#     axis.text.x = element_text(size = 20), #,hjust = 0,vjust = 1
#     axis.text.y = element_text(size = 20),
#     plot.title = element_text(size=20),
#     # strip.text = element_blank(),
#     legend.position = "right", #c(0.9,0.8),#,#
#     legend.text = element_text(size= 16),
#     legend.title= element_text(size= 16))
# p0
ggsave('coef_com.pdf',width = 5,height = 4)


dst <- "GSE110381_diff"
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/"
# gold <- read.table(paste0(pre,"/output/",dst,"/call_peak_all/expeak/b5_d50_p1_11RNA_expeakGold.bed"),check.names = F,header = F)
# goldPeak <- diff.list[[dst]]
# #goldPeak <- goldPeak[goldPeak$peak.id %in% gold$V4,]
# goldPeak$type3 <- "all"
# goldPeak$type3[goldPeak$peak.id %in% gold$V4] <- "gold"

#peak.id <- read.table(paste0(pre,"/output/TCGA_small_diff/call_peak_all/count_matrix/expeakCNN_b5_d50_p1_GSE110381.txt.diff"),check.names = F,header = T)
diff.tcga <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/","TCGA_small_diff","/exPeakCNN_smallDomain_diff_COADvsLAML.diff"), header = T, sep="\t",check.names = F)
diff.gse110381 <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/","GSE110381_diff","/filerSmp_exPeakCNN_smallDomain_diff_CRCvsNC.diff"), header = T, sep="\t",check.names = F)
#diff.gse110381 <- read.table(paste0(pre,"/WCHSU-FTC/output/WSQ_SMARTer_NEB_diff/call_peak_all/count_matrix/expeak_b5_d50_p1_GSE110381.txt.diff.NEB"), header = T, sep="\t",check.names = F)
diff.tcga$feature <- rownames(diff.tcga)
diff.gse110381$feature <- rownames(diff.gse110381)
table(diff.tcga$feature==diff.gse110381$feature) # all TURE !!!


# op1: filter CRC diff high expressed peaks
highTissue <- diff.tcga$feature[diff.tcga$pvalue<=0.05 & diff.tcga$log2FoldChange>0]
length(highTissue) # p=0.05:1542

# # op2: filter Tissue/blood specific peak in TCGA
# cutoff.crc <- 0.1
# cutoff.blood <- 0.1
# table(diff.tcga$posMeanCPM>=cutoff.crc)
# table(diff.tcga$negMeanCPM<=cutoff.blood)
# table(diff.tcga$posMeanCPM>=cutoff.crc & diff.tcga$negMeanCPM<=cutoff.blood)
# highTissue <- diff.tcga$feature[diff.tcga$negMeanCPM<=cutoff.blood & diff.tcga$posMeanCPM>=cutoff.crc]
# length(highTissue) # 0.01: 1591; 0.1:3123

#peak.id <- read.table(paste0(pre,"/output/TCGA_small_diff/call_peak_all/count_matrix/expeak_b5_d50_p1_GSE110381.txt.diffDcb2Union"),check.names = F,header = F)$V1
CRCindexPeak <- diff.list[[dst]]
CRCindexPeak$type3 <- "Other peak"
CRCindexPeak$type3[CRCindexPeak$feature %in% highTissue] <- "TissueHigh peak"
table(CRCindexPeak$type3)
CRCindexPeak$type3 <- factor(CRCindexPeak$type3,levels = c("Other peak","TissueHigh peak"))
wilcox.test(CRCindexPeak$log.abs.mean.coef[CRCindexPeak$type3=="Other peak"],CRCindexPeak$log.abs.mean.coef[CRCindexPeak$type3=="TissueHigh peak"],alternative = "less")
p1 <- ggplot(CRCindexPeak, aes(x=log.abs.mean.coef, y=cpm.mean, group="type3") ) + # type2, type3
  geom_point(aes(fill=type3,color=type3,alpha=cpm.mean*cpm.mean/100), shape=21) + # fill="grey",color="grey",
  scale_fill_manual(values = c("grey50","#ED782F")) + #"#ED782F","#B52B1A","#3E8C9D"
  scale_color_manual(values = c("grey50","#ED782F")) +
  theme_void() +  # base_size=12
  theme(#axis.ticks.x=element_blank(),
    #strip.text.y = element_blank(),
    aspect.ratio = 0.7,
    strip.text = element_text(size=20),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size = 20), #,hjust = 0,vjust = 1
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size=20),
    # strip.text = element_blank(),
    legend.position = "none", #c(0.9,0.8),#,#
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16))
p <- ggExtra::ggMarginal(p1,groupFill = TRUE, groupColour = TRUE, type = "density",size=4,alpha=0.6)
# p
#ggsave(plot = p,filename = 'coef_marg_gold.pdf',width = 9,height = 4)
#ggsave(plot = p,filename = 'coef_marg_index.pdf',width = 9,height = 4)
# ggsave(plot = p,filename = 'coef_marg_CRChighBloodlow.pdf',width = 9,height = 4)
ggsave(plot = p,filename = 'coef_marg_CRChighdiff.pdf',width = 9,height = 4)




### GSE110381 heatmap 
### filter selected features
#coef.cutoff <- -5  # -2
#cpm.cutoff <- .05 #.01
cpm.cutoff <- 50 #no filter
#table(df.mean$log.abs.mean.coef>=coef.cutoff)
table(diff.list[[dst]]$cpm.mean<=cpm.cutoff )
#-3
#FALSE  TRUE 
#9067   833 
#-2
#FALSE  TRUE 
#9570   330
#-1
#FALSE  TRUE 
#9756   144
#df.mean.filter <- as.data.frame(df.mean[df.mean$log.abs.mean.coef>=coef.cutoff,])
df.mean.filter <- as.data.frame(diff.list[[dst]][diff.list[[dst]]$cpm.mean<=cpm.cutoff,]) #  & diff.list[[dst]]$log.abs.mean.coef>=coef.cutoff
rownames(df.mean.filter) <- df.mean.filter$feature


### filter logcpm mat
cpm.filter <- cpm.list[[dst]]
rownames(cpm.filter) <- cpm.filter$feature
cpm.filter$feature <- NULL
cpm.filter <- cpm.filter[df.mean.filter$feature, sample.table$sample] #cpm.filter[rownames(cpm.filter) %in% df.mean.filter$feature, colnames(cpm.filter) %in% sample.table$sample]


## order anno tbls
annotation_col_gse <- sample.table[colnames(cpm.filter),]
annotation_col_gse$sample <- NULL
annotation_row_gse <- data.frame( RNA = df.mean.filter[rownames(cpm.filter),"type"],
                              row.names = rownames(cpm.filter)) # , diff=genes$diff
# ann_colors = list(
#   group =  c(pal_nejm_adaptive()(2)), # c("white", "firebrick"),
#   RNA = c(pal_npg_adaptive()(15))
)
# pheatmap::pheatmap(mat = cpm.filter,
#                    annotation_col = annotation_col_gse, #data frame that specifies the annotations shown on left side of the heatmap.
#                    annotation_row = annotation_row_gse,
#                    # annotation_colors = ann_colors, # RColorBrewer::brewer.pal(6, "Set3")
#                    border_color="white",
#                    scale = "row",
#                    #labels_col = 3, labels_row = 6,
#                    cluster_cols = F,cluster_rows = T,
#                    #cutree_cols = 2,cutree_rows = 3,
#                    show_colnames=F, show_rownames=F,
#                    #fontsize = 8,
#                    height = 5,width = 8,
#                    colorRampPalette(c("blue","white","red"))(100),  # steelblue/dodgerblue4-firebrick2
#                    #fontsize_row = 5,
#                    filename ="./GSE_heatmap.pdf")



### get DCB candidates peaks
dst <- "GSE110381_diff"
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/"
# trimMeanCPM.NC <- apply(cpm.filter[,rownames(annotation_col_gse)[annotation_col_gse$group=="NC"] ], 1, function(x) mean(x,trim=0.1))
# trimMeanCPM.CRC <- apply(cpm.filter[,rownames(annotation_col_gse)[annotation_col_gse$group=="CRC"] ], 1, function(x) mean(x,trim=0.1))
# p90.NC <- apply(cpm.filter[,rownames(annotation_col_gse)[annotation_col_gse$group=="NC"] ], 1, function(x) quantile(x, probs = 0.9))
# p90.CRC <- apply(cpm.filter[,rownames(annotation_col_gse)[annotation_col_gse$group=="CRC"] ], 1, function(x) quantile(x, probs = 0.9))
# p10.NC <- apply(cpm.filter[,rownames(annotation_col_gse)[annotation_col_gse$group=="NC"] ], 1, function(x) quantile(x, probs = 0.1))
# p10.CRC <- apply(cpm.filter[,rownames(annotation_col_gse)[annotation_col_gse$group=="CRC"] ], 1, function(x) quantile(x, probs = 0.1))
# sdCPM.NC <- apply(cpm.filter[,rownames(annotation_col_gse)[annotation_col_gse$group=="NC"] ], 1, function(x) sd(x))
# sdCPM.CRC <- apply(cpm.filter[,rownames(annotation_col_gse)[annotation_col_gse$group=="CRC"] ], 1, function(x) sd(x))
# dcb.gse <- rownames(cpm.filter)[(trimMeanCPM.CRC >= trimMeanCPM.NC) & (p90.CRC > p90.NC) & (p10.CRC >= p10.NC) & (sdCPM.CRC > sdCPM.NC)]  #
tmp.diff <- diff.list[[dst]][match(rownames(cpm.filter),diff.list[[dst]]$feature),]

dcb.gse <- rownames(cpm.filter)[(tmp.diff$pvalue <= 0.1) & (tmp.diff$log2FoldChange > 0) & (tmp.diff$posGT1RatioCount >= max(0.05,2/50)) & (tmp.diff$negGT1RatioCount <= 0.1) & (tmp.diff$negCent90CPM <= 0.1)]  # (tmp.diff$pvalue <= 0.5) &  & (tmp.diff$negCent90CPM <= 10)  & (tmp.diff$posGT1RatioCount >= 0.2) & (tmp.diff$negGT1RatioCount <= 0.2)
#50 CRC vs. 34 NC, require minimum 2 samples in postive samples: 2/50
#table(tmp.diff$negGT1RatioCount <= 0.2)  # 44844
#table(tmp.diff$negMeanCPM <= 1)
length(dcb.gse) # 5030



# ## get inhouse diff/CPM/RPKM stats
# dst <- "WSQ_SMARTer_NEB_diff"
# lib <- "NEB"
# pnk <- "PNK-"
# ### read smp table
# sample.table.inhouse <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/WSQ_SMARTer_NEB_diff/sample_table.txt",check.names = F, header = T)
# sample.table.inhouse$sample <- paste0(sample.table.inhouse$sample, "_1")
# rownames(sample.table.inhouse) <- sample.table.inhouse$sample
# #sample.table.inhouse <- sample.table.inhouse[sample.table.inhouse$sample %in% colnames(cpm.inhouse),]
# sample.table.inhouse$lib <- "NEB"
# sample.table.inhouse$lib[grepl("smart",sample.table.inhouse$sample)] <- "SMARTer"
# sample.table.inhouse$PNK <- "PNK-"
# sample.table.inhouse$PNK[grepl("PNK",sample.table.inhouse$sample)] <- "PNK+"
# sample.table.inhouse$group <- "NC"
# sample.table.inhouse$group[grepl("CRC",sample.table.inhouse$sample)] <- "CRC"
# sample.table.inhouse$group <- factor(sample.table.inhouse$group,levels = c("CRC","NC"))
# sample.table.inhouse <- sample.table.inhouse[order(sample.table.inhouse$lib,sample.table.inhouse$group,sample.table.inhouse$PNK),]
# sample.table.inhouse <- sample.table.inhouse[sample.table.inhouse$lib==lib & sample.table.inhouse$PNK==pnk,] # filter lib:SMARTer,NEB
# table(sample.table.inhouse$lib,sample.table.inhouse$PNK,sample.table.inhouse$group)
# 
# # ### read count mat
# count <- read.table(paste0(pre,"/output/",dst,"/call_peak_all/count_matrix/expeak_b5_d50_p1_GSE110381.txt"),check.names = F,header = T)
# rownames(count) <- count$feature
# count <- count[,2:ncol(count)]
# sample.table.inhouse <- sample.table.inhouse[sample.table.inhouse$sample %in% colnames(count),]
# # sample.table.inhouse$group <- gsub("adenoma","CRC",sample.table.inhouse$group)
# # #sample.table.inhouse <- sample.table.inhouse[order(sample.table.inhouse$group),]
# # positive_samples <- sample.table.inhouse[sample.table.inhouse$group=="CRC","sample"]
# # negative_samples <- sample.table.inhouse[sample.table.inhouse$group=="NC","sample"]
# # samples <- c(positive_samples, negative_samples)
# # sample.table.inhouse <- sample.table.inhouse[match(samples,sample.table.inhouse$sample),]
# # group <- c(rep("positive",length(positive_samples)),rep("negative",length(negative_samples)))
# # method <- "edger_glmlrt"
# # norm_method <- "TMM"
# 
# # ????? ????? ????? ????? ????? ?????
# # ???? not match mat samples ?????
# # ????? ????? ????? ????? ?????
# #   
# # diff.res <- diff.v2.dcb(mat = count, samples = samples, group = group, method = method, norm_method = norm_method, filterType = "NULL", featureType = "domain")
# # data.table::fwrite(as.data.frame(cbind("feature"=rownames(diff.res$diffTable),diff.res[["diffTable"]])),paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_all/count_matrix/expeak_b5_d50_p1_GSE110381.txt.diff.",lib),quote = F,sep = "\t",row.names = F,col.names = T) # .SMARTer,NEB
# # data.table::fwrite(as.data.frame(cbind("feature"=rownames(diff.res$diffTable),diff.res[["normMat"]])),paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_all/count_matrix/expeak_b5_d50_p1_GSE110381.txt.rowsumLogCPM.",lib),quote = F,sep = "\t",row.names = F,col.names = T) # .SMARTer,NEB
# # diff.res <- NULL
# 
# cpm.inhouse <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",dst,"/exPeak_smallDomain_diff_CRCvsNC_",lib,"_",pnk,".cpm"),row.names = 1,check.names = F,header = T)
# #cpm.inhouse <- read.table(paste0(pre,"/output/",dst,"/call_peak_all/count_matrix/expeak_b5_d50_p1_GSE110381.txt.rowsumLogCPM"),row.names = 1,check.names = F,header = T)
# #cpm.inhouse <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/WSQ_SMARTer_NEB_diff/call_peak_all/count_matrix/expeak_b5_d50_p1_GSE110381_CPM.txt",check.names = F, header = T)
# #cpm.inhouse <- log2(cpm.inhouse + 1) # log2colnames(cpm)[1] <- "feature"
# cpm.list[[dst]] <- cpm.inhouse
# diff.list[[dst]] <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",dst,"/exPeak_smallDomain_diff_CRCvsNC_",lib,"_",pnk,".diff"),check.names = F,header = T)
# #diff.list[[dst]] <- read.table(paste0(pre,"/output/",dst,"/call_peak_all/count_matrix/expeak_b5_d50_p1_GSE110381.txt.diff"),check.names = F,header = T)
# diff.list[[dst]]$feature <- rownames(diff.list[[dst]])
# 
# ### filter logcpm mat
# cpm.inhouse.filter <- cpm.inhouse[df.mean.filter$feature, sample.table.inhouse$sample] 
# cpm.inhouse.filter <- cpm.inhouse.filter[apply(cpm.inhouse.filter,1,sum)>0,] # rm 0 sd
# #table(apply(cpm.inhouse.filter,1,sum)==0)
# 
# ### order anno tbls
# annotation_col_inhouse <- sample.table.inhouse[colnames(cpm.inhouse.filter),]
# annotation_col_inhouse$sample <- NULL
# annotation_col_inhouse <- annotation_col_inhouse[,c("lib","group","PNK")]
# annotation_row_inhouse <- data.frame( RNA = df.mean.filter[rownames(cpm.inhouse.filter),"type"],
#                               row.names = rownames(cpm.inhouse.filter)) # , diff=genes$diff
# # 
# # pheatmap::pheatmap(mat = cpm.inhouse.filter,
# #                    annotation_col = annotation_col_inhouse, #data frame that specifies the annotations shown on left side of the heatmap.
# #                    annotation_row = annotation_row_inhouse,
# #                    # annotation_colors = ann_colors, # RColorBrewer::brewer.pal(6, "Set3")
# #                    border_color="white",
# #                    scale = "row",
# #                    #labels_col = 3, labels_row = 6,
# #                    cluster_cols = F,cluster_rows = T,
# #                    #cutree_cols = 2,cutree_rows = 3,
# #                    show_colnames=F, show_rownames=F,
# #                    #fontsize = 8,
# #                    height = 5,width = 8,
# #                    color = colorRampPalette(c("blue","white","red"))(100),  # steelblue/dodgerblue4-firebrick2
# #                    
# #                    #fontsize_row = 5,
# #                    filename ="./inhouse_heatmap.pdf")
# 
# ### get DCB candidates peaks
# # trimMeanCPM.inhouse.NC <- apply(cpm.inhouse.filter[,rownames(annotation_col_inhouse)[annotation_col_inhouse$group=="NC"] ], 1, function(x) mean(x,trim=0.1))
# # trimMeanCPM.inhouse.CRC <- apply(cpm.inhouse.filter[,rownames(annotation_col_inhouse)[annotation_col_inhouse$group=="CRC"] ], 1, function(x) mean(x,trim=0.1))
# # p90.inhouse.NC <- apply(cpm.inhouse.filter[,rownames(annotation_col_inhouse)[annotation_col_inhouse$group=="NC"] ], 1, function(x) quantile(x, probs = 0.9))
# # p90.inhouse.CRC <- apply(cpm.inhouse.filter[,rownames(annotation_col_inhouse)[annotation_col_inhouse$group=="CRC"] ], 1, function(x) quantile(x, probs = 0.9))
# # p10.inhouse.NC <- apply(cpm.inhouse.filter[,rownames(annotation_col_inhouse)[annotation_col_inhouse$group=="NC"] ], 1, function(x) quantile(x, probs = 0.1))
# # p10.inhouse.CRC <- apply(cpm.inhouse.filter[,rownames(annotation_col_inhouse)[annotation_col_inhouse$group=="CRC"] ], 1, function(x) quantile(x, probs = 0.1))
# # sdCPM.inhouse.NC <- apply(cpm.inhouse.filter[,rownames(annotation_col_inhouse)[annotation_col_inhouse$group=="NC"] ], 1, function(x) sd(x))
# # sdCPM.inhouse.CRC <- apply(cpm.inhouse.filter[,rownames(annotation_col_inhouse)[annotation_col_inhouse$group=="CRC"] ], 1, function(x) sd(x))
# # dcb.inhouse <- rownames(cpm.inhouse.filter)[(trimMeanCPM.inhouse.CRC >= trimMeanCPM.inhouse.NC) & (p90.inhouse.CRC > p90.inhouse.NC) & (p10.inhouse.CRC >= p10.inhouse.NC) & (sdCPM.inhouse.CRC > sdCPM.inhouse.NC)]  #  
# # dcb.inhouse <- rownames(cpm.inhouse.filter)[(trimMeanCPM.inhouse.CRC >= trimMeanCPM.inhouse.NC) & (p90.inhouse.CRC > p90.inhouse.NC) & (p10.inhouse.CRC >= p10.inhouse.NC) & (sdCPM.inhouse.CRC > sdCPM.inhouse.NC)]  #
# tmp.diff <- diff.list[[dst]][match(rownames(cpm.inhouse.filter),diff.list[[dst]]$feature),]
# dcb.inhouse <- rownames(cpm.inhouse.filter)[(tmp.diff$pvalue <= 0.1) & (tmp.diff$log2FoldChange > 0) & (tmp.diff$posGT1RatioCount >= 0.05) & (tmp.diff$negGT1RatioCount <= 0.1) & (tmp.diff$negCent90CPM <= 0.1)]  # this version not use ex libsize
# #table(tmp.diff$negGT1RatioCount <= 0.2) 
# length(dcb.inhouse) # 1946,699,2073
# nrow(cpm.inhouse.filter)



## tcga heatmap 
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
dst <- "TCGA_small_diff"
### read smp table
sample.table.tcga <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/TCGA_small_diff/sample_table.txt",check.names = F, sep = "\t", header = T)
sample.table.tcga <- sample.table.tcga[,1:2]
colnames(sample.table.tcga) <- c("sample","group")
sample.table.tcga$sample <- gsub(".bam","",sample.table.tcga$sample)
rownames(sample.table.tcga) <- sample.table.tcga$sample
#table(sample.table.tcga$sample %in% colnames(cpm.tcga))
#sample.table.tcga <- sample.table.tcga[sample.table.tcga$sample %in% colnames(cpm.tcga),]
sample.table.tcga$group <- gsub("TCGA-","",sample.table.tcga$group)
table(sample.table.tcga$group)
sample.table.tcga[1:3,]
smpid <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/TCGA_small_diff/sample_ids.txt")$V1
#write.table(sample.table.tcga[sample.table.tcga$sample %in% smpid,],"/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/TCGA_small_diff/call_peak_all/too/smps.csv",sep = ",",quote = F,row.names = F,col.names = T)
sample.table.tcga$group <- gsub("COAD","CRC",sample.table.tcga$group)
sample.table.tcga$group <- gsub("LAML","Blood",sample.table.tcga$group)
sample.table.tcga <- sample.table.tcga[sample.table.tcga$group %in% c("CRC", "Blood"),]
sample.table.tcga$group <- factor(sample.table.tcga$group,levels = c("CRC", "Blood")) # , "PAAD", "PRAD"
sample.table.tcga <- sample.table.tcga[order(sample.table.tcga$group),]


### read count mat
count <- read.table(paste0(pre,"/output/",dst,"/call_peak_all/count_matrix/expeakCNN_b5_d50_p1_GSE110381.txt"),check.names = F,header = T)
rownames(count) <- count$feature
count <- count[,2:ncol(count)]
sample.table.tcga <- sample.table.tcga[sample.table.tcga$sample %in% colnames(count),]
sample.table.tcga$group <- gsub("adenoma","CRC",sample.table.tcga$group)
# #sample.table.tcga <- sample.table.tcga[order(sample.table.tcga$group),]
# positive_samples <- sample.table.tcga[sample.table.tcga$group=="CRC","sample"]
# negative_samples <- sample.table.tcga[sample.table.tcga$group=="Blood","sample"] # Blood !
# samples <- c(positive_samples, negative_samples)
# sample.table.tcga <- sample.table.tcga[match(samples,sample.table.tcga$sample),]
# group <- c(rep("positive",length(positive_samples)),rep("negative",length(negative_samples)))
# method <- "edger_glmlrt"
# norm_method <- "TMM"

# diff.res <- diff.v2.dcb(mat = count, samples = samples, group = group, method = method, norm_method = norm_method, filterType = "NULL", featureType = "domain")
# data.table::fwrite(as.data.frame(cbind("feature"=rownames(diff.res$diffTable),diff.res[["diffTable"]])),paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_all/count_matrix/expeak_b5_d50_p1_GSE110381.txt.diff"),quote = F,sep = "\t",row.names = F,col.names = T)
# data.table::fwrite(as.data.frame(cbind("feature"=rownames(diff.res$diffTable),diff.res[["normMat"]])),paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_all/count_matrix/expeak_b5_d50_p1_GSE110381.txt.rowsumLogCPM"),quote = F,sep = "\t",row.names = F,col.names = T)
# diff.res <- NULL

#diff.list <- list()
cpm.tcga <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",dst,"/exPeakCNN_smallDomain_diff_COADvsLAML.cpm"),row.names = 1,check.names = F,header = T)
#cpm.tcga <- read.table(paste0(pre,"/output/",dst,"/call_peak_all/count_matrix/expeak_b5_d50_p1_GSE110381.txt.rowsumLogCPM"),row.names = 1,check.names = F,header = T)
#cpm.tcga <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/TCGA_small_diff/call_peak_all/count_matrix/expeak_b5_d50_p1_GSE110381_CPM.txt",row.names = 1,check.names = F, header = T)
cpm.list[[dst]] <- cpm.tcga
diff.list[[dst]] <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",dst,"/exPeakCNN_smallDomain_diff_COADvsLAML.diff"),check.names = F,header = T)
#diff.list[[dst]] <- read.table(paste0(pre,"/output/",dst,"/call_peak_all/count_matrix/expeak_b5_d50_p1_GSE110381.txt.diff"),check.names = F,header = T)
diff.list[[dst]]$feature <- rownames(diff.list[[dst]])
dim(cpm.tcga)

### filter logcpm mat
cpm.tcga.filter <- cpm.tcga[df.mean.filter$feature, sample.table.tcga$sample] 
#cpm.tcga.filter <- cpm.tcga
cpm.tcga.filter <- cpm.tcga.filter[apply(cpm.tcga.filter,1,sum)>0,] # rm 0 sum
cpm.tcga.filter <- cpm.tcga.filter[apply(cpm.tcga.filter,1,sd)>0,] # rm 0 std
#table(apply(cpm.tcga.filter,1,sum)==0)

### order anno tbls
annotation_col_tcga <- sample.table.tcga[colnames(cpm.tcga.filter),]
annotation_col_tcga$sample <- NULL
annotation_row_tcga <- data.frame( RNA = df.mean.filter[rownames(cpm.tcga.filter),"type"],
                              row.names = rownames(cpm.tcga.filter)) # , diff=genes$diff

# pheatmap::pheatmap(mat = cpm.tcga.filter,
#                    annotation_col = annotation_col_tcga, #data frame that specifies the annotations shown on left side of the heatmap.
#                    annotation_row = annotation_row_tcga,
#                    # annotation_colors = ann_colors, # RColorBrewer::brewer.pal(6, "Set3")
#                    border_color="white",
#                    scale = "row",
#                    #labels_col = 3, labels_row = 6,
#                    cluster_cols = F,cluster_rows = T,
#                    #cutree_cols = 2,cutree_rows = 3,
#                    show_colnames=F, show_rownames=F,
#                    #fontsize = 8,
#                    height = 5,width = 8,
#                    color = colorRampPalette(c("blue","white","red"))(100),  # steelblue/dodgerblue4-firebrick2
                   # 
                   # #fontsize_row = 5,
                   # filename ="./tcga_heatmap.pdf")

### get DCB candidates peaks
# trimMeanCPM.tcga.NC <- apply(cpm.tcga.filter[,rownames(annotation_col_tcga)[annotation_col_tcga$group=="Blood"] ], 1, function(x) mean(x,trim=0.1))
# trimMeanCPM.tcga.CRC <- apply(cpm.tcga.filter[,rownames(annotation_col_tcga)[annotation_col_tcga$group=="CRC"] ], 1, function(x) mean(x,trim=0.1))
# p90.tcga.NC <- apply(cpm.tcga.filter[,rownames(annotation_col_tcga)[annotation_col_tcga$group=="Blood"] ], 1, function(x) quantile(x, probs = 0.9))
# p90.tcga.CRC <- apply(cpm.tcga.filter[,rownames(annotation_col_tcga)[annotation_col_tcga$group=="CRC"] ], 1, function(x) quantile(x, probs = 0.9))
# p10.tcga.NC <- apply(cpm.tcga.filter[,rownames(annotation_col_tcga)[annotation_col_tcga$group=="Blood"] ], 1, function(x) quantile(x, probs = 0.1))
# p10.tcga.CRC <- apply(cpm.tcga.filter[,rownames(annotation_col_tcga)[annotation_col_tcga$group=="CRC"] ], 1, function(x) quantile(x, probs = 0.1))
#sdCPM.tcga.NC <- apply(cpm.tcga.filter[,rownames(annotation_col_tcga)[annotation_col_tcga$group=="Blood"] ], 1, function(x) sd(x))
#sdCPM.tcga.CRC <- apply(cpm.tcga.filter[,rownames(annotation_col_tcga)[annotation_col_tcga$group=="CRC"] ], 1, function(x) sd(x))
#dcb.tcga <- rownames(cpm.tcga.filter)[(trimMeanCPM.tcga.CRC >= trimMeanCPM.tcga.NC) & (p90.tcga.CRC > p90.tcga.NC) & (p10.tcga.CRC >= p10.tcga.NC) ]  #  & (sdCPM.tcga.CRC > sdCPM.tcga.NC) (LAML  v.s. COAD)
tmp.diff <- diff.list[[dst]][match(rownames(cpm.tcga.filter),diff.list[[dst]]$feature),]
#table(tmp.diff$negCent90CPM<=0.5)
dcb.tcga <- rownames(cpm.tcga.filter)[(tmp.diff$pvalue <= 0.2) & (tmp.diff$log2FoldChange > 0) & (tmp.diff$posGT1RatioCount >= max(0.1,1/8))  & (tmp.diff$negCent90CPM <= 0.1)]  # this version not use ex libsize
#8COAD vs. 30Blood, required min 1 sampels in positive sample: tmp.diff$posGT1RatioCount >= max(0.1,1/8)
#table(tmp.diff$negGT1RatioCount <= 0.2) 
# data.table::fwrite(x = as.data.frame(dcb.tcga),file = paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_all/count_matrix/expeakCNN_b5_d50_p1_GSE110381.txt.diff.dcb"),quote = F,sep = "\t",row.names = F,col.names = F)
length(dcb.tcga) # pval0.2:2012; pval0.1:1851; pval0.05:762; pval0.01:328; 
#nrow(cpm.tcga.filter)



## merge  DCB candidates peaks, plot venn (2 dst)
length(dcb.tcga)
df <- data.frame(peak=unique(c(dcb.gse,dcb.tcga)),gse=0,tcga=0)
rownames(df) <- df$peak
df[dcb.gse,"gse"] <- 1
df[dcb.tcga,"tcga"] <- 1
df$peak <- NULL
dcb <- rownames(df)[rowSums(df)==2]
#dcb <- dcb[1:20]
dcb.peak <- unlist(sapply(strsplit(dcb,"|",fixed=T),"[",3))
dcb.RNA <- unlist(sapply(strsplit(dcb,"|",fixed=T),"[",2))
dcb.df <- data.frame(dcb=dcb,dcb.RNA=dcb.RNA,dcb.peak=dcb.peak)
dcb.df$dcb.RNA <- gsub("_rev|_for","",dcb.df$dcb.RNA,perl = T)
dcb.df$dcb.RNA <- factor(dcb.df$dcb.RNA, levels = c(rna,dna))
dcb.df <- dcb.df[order(dcb.df$dcb.RNA,decreasing = F),]
dcb <- dcb.df$dcb
dcb.peak <- dcb.df$dcb.peak
#dcb.peak[4:5]
length(dcb.peak)
dcb[duplicated(dcb.peak)]

library(ggvenn)
x <- list(dcb.gse,dcb.tcga)
names(x) <- c("GSE110381 plasma", "TCGA tissue")
ggvenn(data = x, stroke_color = "white",text_color = "black", text_size = 10,
       fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF"), # , "#868686FF"
       stroke_size = 0.5, set_name_size = 8
)
ggsave("dcb_venn.pdf",width = 9,height = 9)

# ## merge  DCB candidates peaks, plot venn (3 dst)
# df <- data.frame(peak=unique(c(dcb.gse,dcb.inhouse,dcb.tcga)),gse=0,inhouse=0,tcga=0)
# rownames(df) <- df$peak
# df[dcb.gse,"gse"] <- 1
# df[dcb.inhouse,"inhouse"] <- 1
# df[dcb.tcga,"tcga"] <- 1
# df$peak <- NULL
# dcb <- rownames(df)[rowSums(df)==3]
# #dcb <- dcb[1:20]
# dcb.peak <- unlist(sapply(strsplit(dcb,"|",fixed=T),"[",3))
# dcb.RNA <- unlist(sapply(strsplit(dcb,"|",fixed=T),"[",2))
# dcb.df <- data.frame(dcb=dcb,dcb.RNA=dcb.RNA,dcb.peak=dcb.peak)
# dcb.df$dcb.RNA <- gsub("_rev|_for","",dcb.df$dcb.RNA,perl = T)
# dcb.df$dcb.RNA <- factor(dcb.df$dcb.RNA, levels = c(rna,dna))
# dcb.df <- dcb.df[order(dcb.df$dcb.RNA,decreasing = F),]
# dcb <- dcb.df$dcb
# dcb.peak <- dcb.df$dcb.peak
# #dcb.peak[4:5]
# length(dcb.peak)
# dcb[duplicated(dcb.peak)]
# 
# library(ggvenn)
# # ggvenn(df,set_name_size=4,text_size = 5,show_percentage=F, text_color = "white") + ggsci::scale_fill_nejm() +
# #   theme(aspect.ratio = 0.8)
# x <- list(dcb.gse,dcb.inhouse,dcb.tcga)
# names(x) <- c("GSE110381 plasma","in-house plasma", "TCGA tissue")
# ggvenn(data = x, stroke_color = "white",text_color = "black", text_size = 10,
#   fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF"), # , "#868686FF"
#   stroke_size = 0.5, set_name_size = 8
# )
# ggsave("dcb_venn3.pdf",width = 9,height = 9)
# #dev.off()



## plot merge heatmap for overlapped feature
### pre-scale each matrix to prevent 3-dst-scaling together)
cpm.filter.norm <- t(scale(t(cpm.filter))) # default: zscore each column
# cpm.inhouse.filter.norm <- t(scale(t(cpm.inhouse.filter)))
cpm.tcga.filter.norm <- t(scale(t(cpm.tcga.filter)))
#colMeans(cpm.filter.norm)
#colSds(cpm.filter.norm)
#rowMeans(cpm.filter.norm)
#rowSds(cpm.filter.norm)

RNA_colors <- list()
for(i in 1:length(c(rna,dna))){
  j <- c(rna,dna)[i]
  RNA_colors[[j]] <- c(pal_nejm_adaptive()(15)[1:14],"#11838D")[i] #pal_d3_adaptive()(15)[i]
}
RNA_colors <- do.call("c",RNA_colors)

# # #1dst: GSE
# cpm.merge <- cbind(cpm.filter.norm[dcb,])
# annotation_col_merge <- rbind(annotation_col_gse[,"group",drop=F])
# annotation_col_merge$dataset <- c(rep("GSE110381",nrow(annotation_col_gse)))
# annotation_col_merge$group <- gsub("Blood","NC",annotation_col_merge$group)
# annotation_row_merge <- annotation_row_gse[dcb,,drop=F]
# 
# ann_colors <- list(group=c(CRC="firebrick",NC="grey50"),RNA=RNA_colors) # ,"salmon",seagreen,orange,chocolate
# table(annotation_col_merge$dataset )
# pheatmap::pheatmap(mat = cpm.merge,
#                    annotation_col = annotation_col_merge, #data frame that specifies the annotations shown on left side of the heatmap.
#                    annotation_row = annotation_row_merge,
#                    annotation_colors = ann_colors, # RColorBrewer::brewer.pal(6, "Set3")
#                    border_color="grey30",
#                    scale = "none", # row
#                    #labels_col = 3, labels_row = 6,
#                    cluster_cols = F,cluster_rows = F,
#                    # gaps_col = 24,
#                    #cutree_cols = 2,cutree_rows = 3,
#                    show_colnames=F, show_rownames=F,
#                    fontsize = 12,
#                    height = 7,width =10,
#                    color = viridis::viridis_pal()(100), #colorRampPalette(c("blue","white","red"))(100),  # steelblue/dodgerblue4-firebrick2
#                    #fontsize_row = 5,
#                    filename ="./merge1GSE_heatmap.pdf")

#2dst (GSE110381+TCGA)
cpm.merge <- cbind(cpm.filter.norm[dcb,],cpm.tcga.filter.norm[dcb,]) # cpm.inhouse.filter.norm, cpm.filter.norm
annotation_col_merge <- rbind(annotation_col_gse[,"group",drop=F],annotation_col_tcga[,"group",drop=F])
annotation_col_merge$dataset <- c(rep("GSE110381 plasma",nrow(annotation_col_gse)),rep("TCGA tissue",nrow(annotation_col_tcga)))
annotation_col_merge$group[annotation_col_merge$dataset=="GSE110381 plasma"] <- gsub("CRC","CRC plasma",annotation_col_merge$group[annotation_col_merge$dataset=="GSE110381 plasma"] )
annotation_col_merge$group[annotation_col_merge$dataset=="GSE110381 plasma"] <- gsub("NC","NC plasma",annotation_col_merge$group[annotation_col_merge$dataset=="GSE110381 plasma"] )
annotation_col_merge$group[annotation_col_merge$dataset=="TCGA tissue"] <- gsub("CRC","CRC tissue",annotation_col_merge$group[annotation_col_merge$dataset=="TCGA tissue"] )
annotation_row_merge <- annotation_row_gse[dcb,,drop=F]
ann_colors <- list(group=c("CRC plasma"="firebrick","NC plasma"="grey70","CRC tissue"="firebrick4","Blood"="grey50"),dataset=c(`GSE110381 plasma`="chocolate",`TCGA tissue`="steelblue"),RNA=RNA_colors) # ,"salmon",seagreen,orange,chocolate
table(annotation_col_merge$dataset )
# annotation_col_merge$group <- factor(annotation_col_merge$group,levels = c("CRC plasma","NC plasma","CRC tissue","Blood"))
pheatmap::pheatmap(mat = cpm.merge,
                   annotation_col = annotation_col_merge, #data frame that specifies the annotations shown on left side of the heatmap.
                   annotation_row = annotation_row_merge,
                   annotation_colors = ann_colors, # RColorBrewer::brewer.pal(6, "Set3")
                   border_color="grey20",
                   scale = "none", # row
                   #labels_col = 3, labels_row = 6,
                   cluster_cols = F,cluster_rows = F,
                   gaps_col = 84,
                   #cutree_cols = 2,cutree_rows = 3,
                   show_colnames=F, show_rownames=F,
                   fontsize = 12,
                   # height = 7,width =10,
                   height = 7,width =15,
                   color = viridis::viridis_pal()(100), #colorRampPalette(c("blue","white","red"))(100),  # steelblue/dodgerblue4-firebrick2
                   #fontsize_row = 5,
                   filename ="./merge_gse-tcga_heatmap.pdf") # merge_gse-tcga_heatmap_2pos.pdf
pheatmap::pheatmap(mat = cpm.merge,
                   annotation_col = annotation_col_merge, #data frame that specifies the annotations shown on left side of the heatmap.
                   annotation_row = annotation_row_merge,
                   annotation_colors = ann_colors, # RColorBrewer::brewer.pal(6, "Set3")
                   border_color="grey20",
                   scale = "none", # row
                   #labels_col = 3, labels_row = 6,
                   cluster_cols = F,cluster_rows = F,
                   gaps_col = 84,
                   #cutree_cols = 2,cutree_rows = 3,
                   show_colnames=T, show_rownames=T,
                   fontsize = 12,
                   # height = 7,width =10,
                   height = 20,width =45,
                   color = viridis::viridis_pal()(100), #colorRampPalette(c("blue","white","red"))(100),  # steelblue/dodgerblue4-firebrick2
                   #fontsize_row = 5,
                   filename ="./merge_gse-tcga_heatmap2.pdf")

# #2dst (inhouse+TCGA)
# cpm.merge <- cbind(cpm.inhouse.filter.norm[dcb,],cpm.tcga.filter.norm[dcb,]) # cpm.inhouse.filter.norm, cpm.filter.norm
# annotation_col_merge <- rbind(annotation_col_inhouse[,"group",drop=F],annotation_col_tcga[,"group",drop=F])
# annotation_col_merge$dataset <- c(rep("in-house plasma",nrow(annotation_col_inhouse)),rep("TCGA tissue",nrow(annotation_col_tcga)))
# table(annotation_col_merge$dataset=="in-house plasma")
# annotation_col_merge$group <- as.character(annotation_col_merge$group)
# annotation_col_merge$group[annotation_col_merge$dataset=="in-house plasma"] <- gsub("CRC","CRC plasma",annotation_col_merge$group[annotation_col_merge$dataset=="in-house plasma"] )
# annotation_col_merge$group[annotation_col_merge$dataset=="in-house plasma"] <- gsub("NC","NC plasma",annotation_col_merge$group[annotation_col_merge$dataset=="in-house plasma"] )
# annotation_col_merge$group[annotation_col_merge$dataset=="TCGA tissue"] <- gsub("CRC","CRC tissue",annotation_col_merge$group[annotation_col_merge$dataset=="TCGA tissue"] )
# # annotation_col_merge$group <- gsub("Blood","NC",annotation_col_merge$group)
# annotation_row_merge <- annotation_row_gse[dcb,,drop=F]
# ann_colors <- list(group=c("CRC plasma"="firebrick","NC plasma"="grey70","CRC tissue"="firebrick4","Blood"="grey50"),dataset=c(`in-house plasma`="chocolate",`TCGA tissue`="steelblue"),RNA=RNA_colors) # ,"salmon",seagreen,orange,chocolate
# table(annotation_col_merge$dataset )
# pheatmap::pheatmap(mat = cpm.merge,
#                    annotation_col = annotation_col_merge, #data frame that specifies the annotations shown on left side of the heatmap.
#                    annotation_row = annotation_row_merge,
#                    annotation_colors = ann_colors, # RColorBrewer::brewer.pal(6, "Set3")
#                    border_color="grey30",
#                    scale = "none", # row
#                    #labels_col = 3, labels_row = 6,
#                    cluster_cols = F,cluster_rows = F,
#                    # gaps_col = 24,
#                    #cutree_cols = 2,cutree_rows = 3,
#                    show_colnames=F, show_rownames=F,
#                    fontsize = 12,
#                    # height = 7,width =10,
#                    height =8,width =16,
# 
#                    color = viridis::viridis_pal()(100), #colorRampPalette(c("blue","white","red"))(100),  # steelblue/dodgerblue4-firebrick2
# 
#                    #fontsize_row = 5,
#                    filename ="./merge_inhouse-tcga_heatmap.pdf")

# #3dst
# cpm.merge <- cbind(cpm.filter.norm[dcb,],cpm.inhouse.filter.norm[dcb,],cpm.tcga.filter.norm[dcb,])
# annotation_col_merge <- rbind(annotation_col_gse[,"group",drop=F],annotation_col_inhouse[,"group",drop=F],annotation_col_tcga[,"group",drop=F])
# annotation_col_merge$dataset <- c(rep("GSE110381 plasma",nrow(annotation_col_gse)),rep("in-house plasma",nrow(annotation_col_inhouse)),rep("TCGA tissue",nrow(annotation_col_tcga)))
# annotation_col_merge$group[annotation_col_merge$dataset!="TCGA tissue"] <- gsub("CRC","CRC plasma",annotation_col_merge$group[annotation_col_merge$dataset!="TCGA tissue"] )
# annotation_col_merge$group[annotation_col_merge$dataset!="TCGA tissue"] <- gsub("NC","NC plasma",annotation_col_merge$group[annotation_col_merge$dataset!="TCGA tissue"] )
# annotation_col_merge$group[annotation_col_merge$dataset=="TCGA tissue"] <- gsub("CRC","CRC tissue",annotation_col_merge$group[annotation_col_merge$dataset=="TCGA tissue"] )
# annotation_row_merge <- annotation_row_gse[dcb,,drop=F]
# ann_colors <- list(group=c("CRC plasma"="firebrick","NC plasma"="grey70","CRC tissue"="firebrick4","Blood"="grey50"),dataset=c(`TCGA tissue`="steelblue", `GSE110381 plasma`="chocolate",`in-house plasma`="orange"),RNA=RNA_colors) # ,"salmon",seagreen,orange,chocolate
# table(annotation_col_merge$dataset )
# # annotation_col_merge$group <- factor(annotation_col_merge$group,levels = c("CRC plasma","NC plasma","CRC tissue","Blood"))
# pheatmap::pheatmap(mat = cpm.merge,
#                    annotation_col = annotation_col_merge, #data frame that specifies the annotations shown on left side of the heatmap.
#                    annotation_row = annotation_row_merge,
#                    annotation_colors = ann_colors, # RColorBrewer::brewer.pal(6, "Set3")
#                    border_color="grey30",
#                    scale = "none", # row
#                    #labels_col = 3, labels_row = 6,
#                    cluster_cols = F,cluster_rows = F,
#                    # gaps_col = 24,
#                    #cutree_cols = 2,cutree_rows = 3,
#                    show_colnames=F, show_rownames=F,
#                    fontsize = 12,
#                    # height = 7,width =10,
#                    height = 4,width =20,
#                    
#                    color = viridis::viridis_pal()(100), #colorRampPalette(c("blue","white","red"))(100),  # steelblue/dodgerblue4-firebrick2
#                    
#                    #fontsize_row = 5,
#                    filename ="./merge_gse-tcga-inhouse_heatmap.pdf")



# #peak <- "T295851_13316_13332_+|tucpRNA|T295851|peak_57431|T295851|13316|13332"
# for (dataset in unique(annotation_col_merge$dataset)){
#   #dataset <- "in-house plasma"
#   for (group in unique(tmp1.annotation_col_merge$group)){
#     #group <- "CRC"
#     for(i in 1:3){
#       annotation_row_merge[peak,paste0(dataset,"_",group,"_",i)] <- ""
#     }
#   }
# }

for (peak in rownames(annotation_row_merge)){
  #peak <- "ENST00000385214_____1_30_49_+|pri_miRNA|ENST00000385214_____1|peak_6181|ENST00000385214_____1|30|49"  # "T295851_13316_13332_+|tucpRNA|T295851|peak_57431|T295851|13316|13332"
  # peak <- "L2a__chr16___3509264____3509549_neg_84_99_+|repeats_rev|L2a__chr16___3509264____3509549_neg|peak_48841|L2a__chr16___3509264____3509549_neg|84|99"
  for (dataset in unique(annotation_col_merge$dataset)){
    #dataset <- "GSE110381 plasma"
    tmp1.annotation_col_merge <- annotation_col_merge[annotation_col_merge$dataset==dataset,]
    tmp1 <- as.data.frame(cpm.merge[peak,annotation_col_merge$dataset==dataset,drop=F])
    for (group in unique(tmp1.annotation_col_merge$group)){
      #group <- "CRC plasma"
      # annotation_row_merge[peak,paste0(dataset,"_",group)] <- ""
      
      tmp2.annotation_col_merge <- tmp1.annotation_col_merge[tmp1.annotation_col_merge$group==group,]
      tmp2 <- as.data.frame(t(tmp1[,tmp1.annotation_col_merge$group==group,drop=F]))
      if(grepl("CRC",group)){
        tmp2 <- tmp2[order(tmp2[,1],decreasing = T),,drop=F]
      } else {
        tmp2 <- tmp2[order(tmp2[,1],decreasing = F),,drop=F]
      }
      top3 <- rownames(tmp2)[1:3]
      annotation_row_merge[peak,paste0(dataset,"_",group,"_",1:3)] <- top3
    }
  }
}
annotation_row_merge$peak <- unlist(sapply(strsplit(rownames(annotation_row_merge),"|",fixed=T),"[",4))
annotation_row_merge$chr <- unlist(sapply(strsplit(rownames(annotation_row_merge),"|",fixed=T),"[",5))
annotation_row_merge$start <- unlist(sapply(strsplit(rownames(annotation_row_merge),"|",fixed=T),"[",6))
annotation_row_merge$end <- unlist(sapply(strsplit(rownames(annotation_row_merge),"|",fixed=T),"[",7))
annotation_row_merge$width <- as.numeric(annotation_row_merge$end) - as.numeric(annotation_row_merge$start)
annotation_row_merge$idx <- 1:nrow(annotation_row_merge)
annotation_row_merge <- annotation_row_merge[,c((ncol(annotation_row_merge)-5):ncol(annotation_row_merge),1:(ncol(annotation_row_merge)-6))]
#use annotation_row_merge for IGV eg.


write.table(x = annotation_row_merge, file = "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/tmp/annotation_row_merge.txt",quote = F,sep = "\t",row.names = F,col.names = T)




# cpm.merge.long[1:3,]
# cpm.merge.long <- dplyr::as_tibble(cpm.merge) %>% 
#   tidyr::pivot_longer(cols = 1:ncol(cpm.merge), names_to = "sample")
# cpm.merge.long$dataset <- annotation_col_merge$dataset[match(cpm.merge.long$samples,rownames(annotation_col_merge))]
# cpm.merge.long$group <- annotation_col_merge$group[match(cpm.merge.long$samples,rownames(annotation_col_merge))]
# table(cpm.merge.long$group )
# cpm.merge.long.top3 <- cpm.merge.long %>% 
#   # dplyr::arrange(desc(site),desc(site_overlap_ratio_product),desc(site_overlap_base)) %>% 
#   # dplyr::top_n(name, site, .keep_all = TRUE) 




# ## complex heatmap
# ## try complexheatmap
# library(ComplexHeatmap)
# library(circlize)
# library(grid)
# RColorBrewer::brewer.pal(11,"Spectral")[4]
# color1 <- circlize::colorRamp2(c(-3,0,3), c("blue","white","red"))
# color2 <- rev(RColorBrewer::brewer.pal(11, "RdBu")) # Spectral
# #color3 <- circlize::colorRamp2(seq(min(mat), max(mat), length=3), c("blue","#EEEEEE", "red"), space = "RGB")
# #color4 <- colorRampPalette(RColorBrewer::brewer.pal(4, "Paired"))(4)
# color4 <- RColorBrewer::brewer.pal(6, "Paired")
# color4 <- color4[c(1,5,2,6)]
# #RColorBrewer::display.brewer.all()
# 
# cormat[1:3,1:3]
# all(colnames(cormat)==rownames(cormat))
# #cormat
# 
# ge <- rownames(cormat)
# ge[grepl("^FTC",ge)] <- "FTC-EV"
# ge[grepl("^FTA",ge)] <- "FTA-EV"
# ge[grepl("csFTC",ge)] <- "FTC-cf"
# ge[grepl("csFTA",ge)] <- "FTA-cf"
# table(ge)
# 
# ht_row = ComplexHeatmap::rowAnnotation(
#   name = "Type",
#   "Tissue Type" = ge, 
#   col = list("Tissue Type"=c("FTA-cf"=color4[1],
#                              "FTA-EV"=color4[2],
#                              "FTC-cf"=color4[3],
#                              "FTC-EV"=color4[4])
#   ),
#   annotation_legend_param = list("Tissue Type"=list(title="Tissue Types",
#                                                     color_bar = "discrete",
#                                                     #direction = "horizontal",
#                                                     #nrow = 2,
#                                                     at=c("FTA-cf",
#                                                          "FTA-EV",
#                                                          "FTC-cf",
#                                                          "FTC-EV"),
#                                                     labels=c("FTA-cf",
#                                                              "FTA-EV",
#                                                              "FTC-cf",
#                                                              "FTC-EV"),
#                                                     legend_gp = gpar(fontsize = 50)
#   )
#   )
# ) #, col = col_letters
# 
# ht_col = ComplexHeatmap::HeatmapAnnotation(
#   name = "Type",
#   "Tissue Type" = ge, 
#   col = list("Tissue Type"=c("FTA-cf"=color4[1],
#                              "FTA-EV"=color4[2],
#                              "FTC-cf"=color4[3],
#                              "FTC-EV"=color4[4])
#   ),
#   annotation_legend_param = list("Tissue Type"=list(title="Tissue Types",
#                                                     color_bar = "discrete",
#                                                     #direction = "horizontal",
#                                                     #nrow = 2,
#                                                     at=c("FTA-cf",
#                                                          "FTA-EV",
#                                                          "FTC-cf",
#                                                          "FTC-EV"),
#                                                     labels=c("FTA-cf",
#                                                              "FTA-EV",
#                                                              "FTC-cf",
#                                                              "FTC-EV"),
#                                                     legend_gp = gpar(fontsize = 50)
#   )
#   )
# ) #, col = col_letters
# 
# #summary(cormat)
# ht_heat <- ComplexHeatmap::Heatmap(cormat, # not zscore default  
#                                    col = color2, na_col = "black", #ComplexHeatmap允许数据中含有NA,只需要通过参数na_col来控制NA的颜色
#                                    # name = "legend", column_title = "Column", column_title_side = "bottom", # change title
#                                    #row_title_gp = grid::gpar(fontsize=3, fontface="bold"), 
#                                    #column_title_gp = grid::gpar(fontsize=3, fontface="bold"),  # change title font
#                                    #row_title_rot = 0, # change row title rotation
#                                    show_row_names = FALSE, show_column_names = FALSE,
#                                    
#                                    
#                                    ## annotation
#                                    show_heatmap_legend = TRUE,
#                                    heatmap_legend_param = list(title="Pearson Cor",
#                                                                legend_gp = gpar(fontsize = 48)
#                                                                #,at = c(-1, 0, 1) 
#                                                                #,labels = c(-1, 0, 1)
#                                    ),
#                                    top_annotation = ht_col, 
#                                    left_annotation = ht_row,
#                                    
#                                    ## cluster
#                                    cluster_rows = FALSE, cluster_columns = FALSE,  
#                                    # cluster_rows/columns ：是否进行聚类
#                                    # show_column/row_dend ：是否显示聚类树
#                                    # column/row_dend_side ：聚类图绘制的位置
#                                    # column_dend_height/row_dend_widht ：聚类树的高度 和 宽度
#                                    
#                                    ## heat box
#                                    rect_gp = grid::gpar(col = "white", lty = 1, lwd = 0.1)
# )
# ht_heat
# 
# 
# pdf("./output/lulab/FTC_small/miR_cpm_cor_complexheatmap.pdf",height = 7,width = 8)
# ComplexHeatmap::draw(
#   ht_heat, 
#   #heatmap_legend_side = "left",
#   row_title = "", 
#   row_title_gp = gpar(col = "black",fontsize = 16),
#   column_title = "Sample Pearson Cor",
#   column_title_gp = gpar(col = "black",fontsize = 8)
# )
# dev.off()










# PRJNA540919 CPM vs. LASSO coef scatter + DCB filtering/overlap (Fig6, suppl Fig6) ----------------------------
## get tx ref
ref <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt") #,header=T
ref[1:3,1:3]


## get GSE110381_diff ML LASSO coef
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
dst <- "PRJNA540919_diff"
read.coef <- function(dst, type){ # seed, 
  # dst <- "GSE110381_diff"
  # type <- "mRNA"
  # bootstrap <- 1
  # print(paste0(type,"_",bootstrap))
  #_auroc.txt
  #pre <- "/BioII/lulab_b/wangtaiwei/peak_calling/output/GSE110381_diff/CV-20fold"
  l <- list()
  for(seed in c(1:99)){ # 000,111,222,333,444,555,666,777,888,999;   100,101,102,103,104,105,106,107,108,109
    # seed <- ""
    # print(seed)
    # for(k in 0:9){
      #"/BioII/lulab_b/wangtaiwei/peak_calling/output/GSE110381_diff/stCV10_LR2_C0.01/enhancer_5_coef.txt"
      #k <- 0
      # print(k)
      tmp <- read.table( paste0("/BioII/lulab_b/wangtaiwei/peak_calling/output/",dst,"/call_peak_all/BOO_l2_0.0001/",type,"_",seed,"_coef.txt"), header = T)
      tmp$type <- type
      # tmp$k <- k
      tmp$seed <- seed
      return(tmp)   
    # }
  }
}
rna <- c("rRNA","pri_miRNA","lncRNA", "mRNA", "piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA')
dna <- c("intron","promoter", "enhancer", "repeats") 
cmb <- expand.grid(type=c(rna,dna), seed=1:99)
res.list <- lapply(1:nrow(cmb), function(j) {read.coef(dst=dst, type=as.character(cmb$type[j]) )} ) # seed=as.character(cmb$seed[j]), 
df <- do.call(rbind,res.list)
colnames(df)[1] <- "feature"
#df$bootstrap
df[1:3,]
#hist(df$coef,breaks = 10000,xlim = c(-1,1))
df.mean <- dplyr::as_tibble(df) %>% 
  dplyr::group_by(type,feature) %>% 
  dplyr::summarise(mean.coef=mean(coef,trim=0.05))
df.mean <- as.data.frame(df.mean)
rownames(df.mean) <- df.mean$feature
hist(df$coef) # 
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.071107 -0.001795  0.000053 -0.000076  0.001827  0.074676
table(R.utils::isZero(df$coef)) # near 
table(R.utils::isZero(df.mean$mean.coef)) # near 
dim(df.mean)
# df.mean <- df.mean[!R.utils::isZero(df.mean$mean.coef),]
# dim(df.mean)
#ENST00000362125_____3: hsa-miR-340-5p	logfc<0

cpm.list <- list()
diff.list <- list()



## get PRJNA540919 diff/CPM/RPKM stats
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
dst <- "PRJNA540919_diff"

### read smp table
sample.table <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/",dst,"/sample_table.txt"),check.names = F,header = T)
rownames(sample.table) <- sample.table$sample
table(sample.table$group)

### read count
count <- read.table(paste0(pre,"/output/",dst,"/call_peak_all/count_matrix/cfpeakCNN_b5_d50_p1_GSE110381.txt"),check.names = F,header = T)
rownames(count) <- count$feature
count <- count[,2:ncol(count)]
sample.table <- sample.table[sample.table$sample %in% colnames(count),]
# positive_samples <- sample.table[sample.table$group=="CRC","sample"]
# negative_samples <- sample.table[sample.table$group=="NC","sample"]
# samples <- c(positive_samples, negative_samples)
# sample.table <- sample.table[match(samples,sample.table$sample),]
# group <- c(rep("positive",length(positive_samples)),rep("negative",length(negative_samples)))
# method <- "edger_glmlrt"
# norm_method <- "TMM"

cpm <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",dst,"/cfPeakCNN_smallDomain_diff_CRCvsNC.cpm"),check.names = F,header = T)
cpm$feature <- rownames(cpm)
cpm <- cpm[,c(ncol(cpm),1:(ncol(cpm)-1))]
cpm.list[[dst]] <- cpm
diff.list[[dst]] <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",dst,"/cfPeakCNN_smallDomain_diff_CRCvsNC.diff"),check.names = F,header = T)

diff.list[[dst]]$feature <- rownames(diff.list[[dst]])
diff.list[[dst]] <- merge(diff.list[[dst]],df.mean,by="feature")

diff.list[[dst]]$cpm.median <- apply(cpm.list[[dst]][,2:ncol(cpm.list[[dst]])], 1, median ) # re-cal cpm in all group
diff.list[[dst]]$cpm.mean <- apply(cpm.list[[dst]][,2:ncol(cpm.list[[dst]])], 1, function(x) mean(x, trim=0.05) ) # re-cal cpm in all group
diff.list[[dst]]$log.abs.mean.coef <- log2(abs(diff.list[[dst]]$mean.coef))
#diff.list[[dst]]$RNA <- unlist(sapply(strsplit(diff.list[[dst]]$feature,"|",fixed=T),"[",2))
diff.list[[dst]]$peak <- unlist(sapply(strsplit(diff.list[[dst]]$feature,"|",fixed=T),"[",3))
diff.list[[dst]]$peak.id <- unlist(sapply(strsplit(diff.list[[dst]]$feature,"|",fixed=T),"[",4))
diff.list[[dst]]$tx.length <- ref$tx.length[match(diff.list[[dst]]$peak,ref$transcript_id)]
diff.list[[dst]]$peak.width <- as.numeric(unlist(sapply(strsplit(diff.list[[dst]]$feature,"|",fixed=T),"[",7)))-as.numeric(unlist(sapply(strsplit(diff.list[[dst]]$feature,"|",fixed=T),"[",6)))
hist(diff.list[[dst]]$peak.width,xlim = c(0,200),breaks = 1000)
table(diff.list[[dst]]$peak.width >=10 & diff.list[[dst]]$peak.width<=200)
diff.list[[dst]] <- diff.list[[dst]][diff.list[[dst]]$peak.width >=10 & diff.list[[dst]]$peak.width<=200,] # filter peak length
diff.list[[dst]]$type2 <- ""
diff.list[[dst]]$type2[diff.list[[dst]]$type %in% rna] <- "peak in annotated tx"
diff.list[[dst]]$type2[diff.list[[dst]]$type == "pri_miRNA"] <- "peak in primary miRNA"
diff.list[[dst]]$type2[diff.list[[dst]]$type %in% dna] <- "peak in unannotated region"
diff.list[[dst]]$type2 <- factor(diff.list[[dst]]$type2,levels = c("peak in primary miRNA","peak in annotated tx","peak in unannotated region"))
table(diff.list[[dst]]$type2)
#diff.list[[dst]][grepl("ENST00000362125_____3",diff.list[[dst]]$feature),] # seem top importance


### plot scatter 
library(ggplot2)
dst <- "PRJNA540919_diff"
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/"

diff.tcga <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/","TCGA_small_diff3","/cfPeakCNN_GSE110381_smallDomain_diff_COADvsLAML.diff"), header = T, sep="\t",check.names = F)
diff.gse110381 <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",dst,"/filterSmp_exPeakCNN_smallDomain_diff_CRCvsNC.diff"), header = T, sep="\t",check.names = F)
diff.tcga$feature <- rownames(diff.tcga)
diff.gse110381$feature <- rownames(diff.gse110381)
table(diff.tcga$feature==diff.gse110381$feature) # all TURE !!!


# op1: filter CRC diff high expressed peaks
highTissue <- diff.tcga$feature[diff.tcga$pvalue<=0.05 & diff.tcga$log2FoldChange>0]
length(highTissue) # p=0.05:4858

# # op2: filter Tissue/blood specific peak in TCGA
# cutoff.crc <- 0.1
# cutoff.blood <- 0.1
# table(diff.tcga$posMeanCPM>=cutoff.crc)
# table(diff.tcga$negMeanCPM<=cutoff.blood)
# table(diff.tcga$posMeanCPM>=cutoff.crc & diff.tcga$negMeanCPM<=cutoff.blood)
# highTissue <- diff.tcga$feature[diff.tcga$negMeanCPM<=cutoff.blood & diff.tcga$posMeanCPM>=cutoff.crc]
# length(highTissue) # 0.01: 1591; 0.1:3123

#peak.id <- read.table(paste0(pre,"/output/TCGA_small_diff/call_peak_all/count_matrix/expeak_b5_d50_p1_GSE110381.txt.diffDcb2Union"),check.names = F,header = F)$V1
CRCindexPeak <- diff.list[[dst]]
CRCindexPeak <- CRCindexPeak[CRCindexPeak$mean.coef!=0 ,]
#table(R.utils::isZero(CRCindexPeak$mean.coef))
#table(CRCindexPeak$mean.coef==0)

CRCindexPeak$type3 <- "Other peak"
CRCindexPeak$type3[CRCindexPeak$feature %in% highTissue] <- "TissueHigh peak"
table(CRCindexPeak$type3)
CRCindexPeak$type3 <- factor(CRCindexPeak$type3,levels = c("Other peak","TissueHigh peak"))
CRCindexPeak <- CRCindexPeak[order(CRCindexPeak$type3),]
wilcox.test(CRCindexPeak$log.abs.mean.coef[CRCindexPeak$type3=="Other peak"],CRCindexPeak$log.abs.mean.coef[CRCindexPeak$type3=="TissueHigh peak"],alternative = "less")
p1 <- ggplot(CRCindexPeak, aes(x=log.abs.mean.coef, y=cpm.mean, group="type3") ) + # type2, type3
  geom_point(aes(fill=type3,color=type3,alpha=cpm.mean*cpm.mean/100), shape=21) + # fill="grey",color="grey",
  scale_fill_manual(values = c("grey50","#ED782F")) + #"#ED782F","#B52B1A","#3E8C9D"
  scale_color_manual(values = c("grey50","#ED782F")) +
  theme_void() +  # base_size=12
  theme(#axis.ticks.x=element_blank(),
    #strip.text.y = element_blank(),
    aspect.ratio = 0.7,
    strip.text = element_text(size=20),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size = 20), #,hjust = 0,vjust = 1
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size=20),
    # strip.text = element_blank(),
    legend.position = "none", #c(0.9,0.8),#,#
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16))
p <- ggExtra::ggMarginal(p1,groupFill = TRUE, groupColour = TRUE, type = "density",size=4,alpha=0.6)
p
#ggsave(plot = p,filename = 'coef_marg_gold.pdf',width = 9,height = 4)
#ggsave(plot = p,filename = 'coef_marg_index.pdf',width = 9,height = 4)
# ggsave(plot = p,filename = 'coef_marg_CRChighBloodlow.pdf',width = 9,height = 4)
ggsave(plot = p,filename = 'coef_marg_CRChighdiff.pdf',width = 9,height = 4)




### PRJNA540919 heatmap 
### filter selected features
#coef.cutoff <- -5  # -2
#cpm.cutoff <- .05 #.01
cpm.cutoff <- 50 #no filter
#table(df.mean$log.abs.mean.coef>=coef.cutoff)
table(diff.list[[dst]]$cpm.mean<=cpm.cutoff )
df.mean.filter <- as.data.frame(diff.list[[dst]][diff.list[[dst]]$cpm.mean<=cpm.cutoff,]) #  & diff.list[[dst]]$log.abs.mean.coef>=coef.cutoff
rownames(df.mean.filter) <- df.mean.filter$feature


### filter logcpm mat
cpm.filter <- cpm.list[[dst]]
rownames(cpm.filter) <- cpm.filter$feature
cpm.filter$feature <- NULL
cpm.filter <- cpm.filter[df.mean.filter$feature, sample.table$sample] #cpm.filter[rownames(cpm.filter) %in% df.mean.filter$feature, colnames(cpm.filter) %in% sample.table$sample]


## order anno tbls
annotation_col_gse <- sample.table[colnames(cpm.filter),]
annotation_col_gse$sample <- NULL
annotation_row_gse <- data.frame( RNA = df.mean.filter[rownames(cpm.filter),"type"],
                                  row.names = rownames(cpm.filter)) # , diff=genes$diff



### get DCB candidates peaks
dst <- "PRJNA540919_diff"
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/"
tmp.diff <- diff.list[[dst]][match(rownames(cpm.filter),diff.list[[dst]]$feature),]
dcb.gse <- rownames(cpm.filter)[(tmp.diff$pvalue <= 0.1) & (tmp.diff$log2FoldChange > 0) & (tmp.diff$posGT1RatioCount >= max(0.1,2/15)) & (tmp.diff$negGT1RatioCount <= 0.1) & (tmp.diff$negCent90CPM <= 0.1)]  # (tmp.diff$pvalue <= 0.5) &  & (tmp.diff$negCent90CPM <= 10)  & (tmp.diff$posGT1RatioCount >= 0.2) & (tmp.diff$negGT1RatioCount <= 0.2)
#15 CRC vs. 10 NC, require minimum 2 samples in postive samples: 2/15
length(dcb.gse) # 3894




## tcga heatmap 
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
dst <- "TCGA_small_diff3"
### read smp table
sample.table.tcga <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/archive/",dst,"/sample_table.txt"), check.names = F, sep = "\t", header = T)
sample.table.tcga <- sample.table.tcga[,c("File Name","Project ID")]
colnames(sample.table.tcga) <- c("sample","group")
sample.table.tcga$sample <- gsub(".bam","",sample.table.tcga$sample)
rownames(sample.table.tcga) <- sample.table.tcga$sample
sample.table.tcga$group <- gsub("TCGA-","",sample.table.tcga$group)
table(sample.table.tcga$group)
sample.table.tcga[1:3,]
smpid <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/archive/",dst,"/sample_ids.txt"))$V1
sample.table.tcga$group <- gsub("COAD","CRC",sample.table.tcga$group)
sample.table.tcga$group <- gsub("LAML","Blood",sample.table.tcga$group)
sample.table.tcga <- sample.table.tcga[sample.table.tcga$group %in% c("CRC", "Blood"),]
sample.table.tcga$group <- factor(sample.table.tcga$group,levels = c("CRC", "Blood")) # , "PAAD", "PRAD"
sample.table.tcga <- sample.table.tcga[order(sample.table.tcga$group),]


### read count mat
count <- read.table(paste0(pre,"/output/",dst,"/call_peak_all/count_matrix/cfpeakCNN_b5_d50_p1_GSE110381.txt"),check.names = F,header = T)
rownames(count) <- count$feature
count <- count[,2:ncol(count)]
sample.table.tcga <- sample.table.tcga[sample.table.tcga$sample %in% colnames(count),]
sample.table.tcga$group <- gsub("adenoma","CRC",sample.table.tcga$group)

cpm.tcga <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",dst,"/cfPeakCNN_GSE110381_smallDomain_diff_COADvsLAML.cpm"),row.names = 1,check.names = F,header = T)
cpm.list[[dst]] <- cpm.tcga
diff.list[[dst]] <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",dst,"/cfPeakCNN_GSE110381_smallDomain_diff_COADvsLAML.diff"),check.names = F,header = T)
diff.list[[dst]]$feature <- rownames(diff.list[[dst]])
dim(cpm.tcga)

### filter logcpm mat
cpm.tcga.filter <- cpm.tcga[df.mean.filter$feature, sample.table.tcga$sample] 
#cpm.tcga.filter <- cpm.tcga
cpm.tcga.filter <- cpm.tcga.filter[apply(cpm.tcga.filter,1,sum)>0,] # rm 0 sum
cpm.tcga.filter <- cpm.tcga.filter[apply(cpm.tcga.filter,1,sd)>0,] # rm 0 std
#table(apply(cpm.tcga.filter,1,sum)==0)

### order anno tbls
annotation_col_tcga <- sample.table.tcga[colnames(cpm.tcga.filter),]
annotation_col_tcga$sample <- NULL
annotation_row_tcga <- data.frame( RNA = df.mean.filter[rownames(cpm.tcga.filter),"type"],
                                   row.names = rownames(cpm.tcga.filter)) # , diff=genes$diff
table(annotation_col_tcga$group)

### get DCB candidates peaks
tmp.diff <- diff.list[[dst]][match(rownames(cpm.tcga.filter),diff.list[[dst]]$feature),]
dcb.tcga <- rownames(cpm.tcga.filter)[(tmp.diff$pvalue <= 0.1) & (tmp.diff$log2FoldChange > 0) & (tmp.diff$posGT1RatioCount >= max(0.1,2/60))  & (tmp.diff$negGT1RatioCount <= 0.1) & (tmp.diff$negCent90CPM <= 0.1)]  # this version not use ex libsize
#60COAD vs. 60Blood, required min 2 sampels in positive sample: tmp.diff$posGT1RatioCount >= max(0.1,2/60)
length(dcb.tcga) # pval0.2:1963; pval0.1:1745; pval0.05:; pval0.01:; 
#nrow(cpm.tcga.filter)



## merge  DCB candidates peaks, plot venn (2 dst)
length(dcb.tcga)
df <- data.frame(peak=unique(c(dcb.gse,dcb.tcga)),gse=0,tcga=0)
rownames(df) <- df$peak
df[dcb.gse,"gse"] <- 1
df[dcb.tcga,"tcga"] <- 1
df$peak <- NULL
dcb <- rownames(df)[rowSums(df)==2]
#dcb <- dcb[1:20]
dcb.peak <- unlist(sapply(strsplit(dcb,"|",fixed=T),"[",3))
dcb.RNA <- unlist(sapply(strsplit(dcb,"|",fixed=T),"[",2))
dcb.df <- data.frame(dcb=dcb,dcb.RNA=dcb.RNA,dcb.peak=dcb.peak)
dcb.df$dcb.RNA <- gsub("_rev|_for","",dcb.df$dcb.RNA,perl = T)
dcb.df$dcb.RNA <- factor(dcb.df$dcb.RNA, levels = c(rna,dna))
dcb.df <- dcb.df[order(dcb.df$dcb.RNA,decreasing = F),]
dcb <- dcb.df$dcb
dcb.peak <- dcb.df$dcb.peak
#dcb.peak[4:5]
length(dcb.peak)
dcb[duplicated(dcb.peak)]

library(ggvenn)
x <- list(dcb.gse,dcb.tcga)
names(x) <- c("PRJNA540919 plasma", "TCGA tissue")
ggvenn(data = x, stroke_color = "white",text_color = "black", text_size = 10,
       fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF"), # , "#868686FF"
       stroke_size = 0.5, set_name_size = 8
)
ggsave("dcb_venn.pdf",width = 9,height = 9)



## plot merge heatmap for overlapped feature
### pre-scale each matrix to prevent 3-dst-scaling together)
cpm.filter.norm <- t(scale(t(cpm.filter))) # default: zscore each column
# cpm.inhouse.filter.norm <- t(scale(t(cpm.inhouse.filter)))
cpm.tcga.filter.norm <- t(scale(t(cpm.tcga.filter)))

RNA_colors <- list()
for(i in 1:length(c(rna,dna))){
  j <- c(rna,dna)[i]
  RNA_colors[[j]] <- c(pal_nejm_adaptive()(15)[1:14],"#11838D")[i] #pal_d3_adaptive()(15)[i]
}
RNA_colors <- do.call("c",RNA_colors)

# # #1dst: GSE
# cpm.merge <- cbind(cpm.filter.norm[dcb,])
# annotation_col_merge <- rbind(annotation_col_gse[,"group",drop=F])
# annotation_col_merge$dataset <- c(rep("GSE110381",nrow(annotation_col_gse)))
# annotation_col_merge$group <- gsub("Blood","NC",annotation_col_merge$group)
# annotation_row_merge <- annotation_row_gse[dcb,,drop=F]
# 
# ann_colors <- list(group=c(CRC="firebrick",NC="grey50"),RNA=RNA_colors) # ,"salmon",seagreen,orange,chocolate
# table(annotation_col_merge$dataset )
# pheatmap::pheatmap(mat = cpm.merge,
#                    annotation_col = annotation_col_merge, #data frame that specifies the annotations shown on left side of the heatmap.
#                    annotation_row = annotation_row_merge,
#                    annotation_colors = ann_colors, # RColorBrewer::brewer.pal(6, "Set3")
#                    border_color="grey30",
#                    scale = "none", # row
#                    #labels_col = 3, labels_row = 6,
#                    cluster_cols = F,cluster_rows = F,
#                    # gaps_col = 24,
#                    #cutree_cols = 2,cutree_rows = 3,
#                    show_colnames=F, show_rownames=F,
#                    fontsize = 12,
#                    height = 7,width =10,
#                    color = viridis::viridis_pal()(100), #colorRampPalette(c("blue","white","red"))(100),  # steelblue/dodgerblue4-firebrick2
#                    #fontsize_row = 5,
#                    filename ="./merge1GSE_heatmap.pdf")

#2dst (PRJNA540919+TCGA)
cpm.merge <- cbind(cpm.filter.norm[dcb,],cpm.tcga.filter.norm[dcb,]) # cpm.inhouse.filter.norm, cpm.filter.norm
annotation_col_merge <- rbind(annotation_col_gse[,"group",drop=F],annotation_col_tcga[,"group",drop=F])
annotation_col_merge$dataset <- c(rep("PRJNA540919 plasma",nrow(annotation_col_gse)),rep("TCGA tissue",nrow(annotation_col_tcga)))
annotation_col_merge$group[annotation_col_merge$dataset=="PRJNA540919 plasma"] <- gsub("CRC","CRC plasma",annotation_col_merge$group[annotation_col_merge$dataset=="PRJNA540919 plasma"] )
annotation_col_merge$group[annotation_col_merge$dataset=="PRJNA540919 plasma"] <- gsub("NC","NC plasma",annotation_col_merge$group[annotation_col_merge$dataset=="PRJNA540919 plasma"] )
annotation_col_merge$group[annotation_col_merge$dataset=="TCGA tissue"] <- gsub("CRC","CRC tissue",annotation_col_merge$group[annotation_col_merge$dataset=="TCGA tissue"] )
annotation_row_merge <- annotation_row_gse[dcb,,drop=F]
ann_colors <- list(group=c("CRC plasma"="firebrick1","NC plasma"="grey80","CRC tissue"="firebrick4","Blood"="grey30"),dataset=c(`PRJNA540919 plasma`="chocolate",`TCGA tissue`="steelblue"),RNA=RNA_colors) # ,"salmon",seagreen,orange,chocolate
table(annotation_col_merge$dataset )
# annotation_col_merge$group <- factor(annotation_col_merge$group,levels = c("CRC plasma","NC plasma","CRC tissue","Blood"))
pheatmap::pheatmap(mat = cpm.merge,
                   annotation_col = annotation_col_merge, #data frame that specifies the annotations shown on left side of the heatmap.
                   annotation_row = annotation_row_merge,
                   annotation_colors = ann_colors, # RColorBrewer::brewer.pal(6, "Set3")
                   border_color="grey20",
                   scale = "none", # row
                   #labels_col = 3, labels_row = 6,
                   cluster_cols = F,cluster_rows = F,
                   gaps_col = c(15,25,85),
                   #cutree_cols = 2,cutree_rows = 3,
                   show_colnames=F, show_rownames=F,
                   fontsize = 12,
                   # height = 7,width =10,
                   height = 7,width =15,
                   color = viridis::viridis_pal()(100), # colorRampPalette(RColorBrewer::brewer.pal(n = 10, name = "RdBu"))(100),  
                   #fontsize_row = 5,
                   filename ="./merge_gse-tcga_heatmap.pdf") # merge_gse-tcga_heatmap_2pos.pdf
pheatmap::pheatmap(mat = cpm.merge,
                   annotation_col = annotation_col_merge, #data frame that specifies the annotations shown on left side of the heatmap.
                   annotation_row = annotation_row_merge,
                   annotation_colors = ann_colors, # RColorBrewer::brewer.pal(6, "Set3")
                   border_color="grey20",
                   scale = "none", # row
                   #labels_col = 3, labels_row = 6,
                   cluster_cols = F,cluster_rows = F,
                   gaps_col = 84,
                   #cutree_cols = 2,cutree_rows = 3,
                   show_colnames=T, show_rownames=T,
                   fontsize = 12,
                   # height = 7,width =10,
                   height = 20,width =45,
                   color = viridis::viridis_pal()(100), #colorRampPalette(c("blue","white","red"))(100),  # steelblue/dodgerblue4-firebrick2
                   #fontsize_row = 5,
                   filename ="./merge_gse-tcga_heatmap2.pdf")

# #peak <- "T295851_13316_13332_+|tucpRNA|T295851|peak_57431|T295851|13316|13332"
# for (dataset in unique(annotation_col_merge$dataset)){
#   #dataset <- "in-house plasma"
#   for (group in unique(tmp1.annotation_col_merge$group)){
#     #group <- "CRC"
#     for(i in 1:3){
#       annotation_row_merge[peak,paste0(dataset,"_",group,"_",i)] <- ""
#     }
#   }
# }

for (peak in rownames(annotation_row_merge)){
  #peak <- "ENST00000385214_____1_30_49_+|pri_miRNA|ENST00000385214_____1|peak_6181|ENST00000385214_____1|30|49"  # "T295851_13316_13332_+|tucpRNA|T295851|peak_57431|T295851|13316|13332"
  # peak <- "L2a__chr16___3509264____3509549_neg_84_99_+|repeats_rev|L2a__chr16___3509264____3509549_neg|peak_48841|L2a__chr16___3509264____3509549_neg|84|99"
  for (dataset in unique(annotation_col_merge$dataset)){
    #dataset <- "GSE110381 plasma"
    tmp1.annotation_col_merge <- annotation_col_merge[annotation_col_merge$dataset==dataset,]
    tmp1 <- as.data.frame(cpm.merge[peak,annotation_col_merge$dataset==dataset,drop=F])
    for (group in unique(tmp1.annotation_col_merge$group)){
      #group <- "CRC plasma"
      # annotation_row_merge[peak,paste0(dataset,"_",group)] <- ""
      
      tmp2.annotation_col_merge <- tmp1.annotation_col_merge[tmp1.annotation_col_merge$group==group,]
      tmp2 <- as.data.frame(t(tmp1[,tmp1.annotation_col_merge$group==group,drop=F]))
      if(grepl("CRC",group)){
        tmp2 <- tmp2[order(tmp2[,1],decreasing = T),,drop=F]
      } else {
        tmp2 <- tmp2[order(tmp2[,1],decreasing = F),,drop=F]
      }
      top3 <- rownames(tmp2)[1:3]
      annotation_row_merge[peak,paste0(dataset,"_",group,"_",1:3)] <- top3
    }
  }
}
annotation_row_merge$peak <- unlist(sapply(strsplit(rownames(annotation_row_merge),"|",fixed=T),"[",4))
annotation_row_merge$chr <- unlist(sapply(strsplit(rownames(annotation_row_merge),"|",fixed=T),"[",5))
annotation_row_merge$start <- unlist(sapply(strsplit(rownames(annotation_row_merge),"|",fixed=T),"[",6))
annotation_row_merge$end <- unlist(sapply(strsplit(rownames(annotation_row_merge),"|",fixed=T),"[",7))
annotation_row_merge$width <- as.numeric(annotation_row_merge$end) - as.numeric(annotation_row_merge$start)
annotation_row_merge$idx <- 1:nrow(annotation_row_merge)
annotation_row_merge <- annotation_row_merge[,c((ncol(annotation_row_merge)-5):ncol(annotation_row_merge),1:(ncol(annotation_row_merge)-6))]
#use annotation_row_merge for IGV eg.
annotation_row_merge$enst <- ""
enst.idx <- which(grepl("ENST",annotation_row_merge$chr))
annotation_row_merge$enst[enst.idx] <- unlist(sapply(strsplit(annotation_row_merge$chr[enst.idx],"_"),"[",1))
tmp <- as.data.frame(mygene::queryMany(annotation_row_merge$enst[enst.idx],scopes="ensembl.transcript",fields=c("ensembl.gene","symbol"),species="human"))
dim(tmp)
table(duplicated(tmp$ensembl.gene))
table(duplicated(tmp$query))
library(dplyr)
library(biomaRt)
mart <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")
ann <- biomaRt::getBM(useCache = FALSE,c("hgnc_symbol","description","chromosome_name","band","strand","start_position","end_position","ensembl_gene_id"),"ensembl_gene_id", tmp$ensembl.gene, mart)
dim(ann)
table(duplicated(ann$hgnc_symbol))
table(duplicated(ann$ensembl_gene_id))
tmp$func <- ann$description[match(tmp$ensembl.gene,ann$ensembl_gene_id)]
length(enst.idx)==nrow(tmp) # T !!!
annotation_row_merge$func <- ""
annotation_row_merge$func[enst.idx] <- tmp$func
annotation_row_merge$symbol <- ""
annotation_row_merge$symbol[enst.idx] <- tmp$symbol
write.table(x = annotation_row_merge, file = "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/tmp/PRNNA-TCGA3_annotation_row_merge.txt",quote = T,sep = "\t",row.names = F,col.names = T)


#



# CRC-peak-index TCGA_small+GSE110381+PRJNA+inhouse (GSE110381 peak) (cfpeak suppl Fig10)  ------------------------------------

## TCGA 
dst <- "TCGA_small_diff3"
#p <- 0.01
#feature <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_all/count_matrix/expeakCNN_b5_d50_p1_GSE110381.txt.diff.tissueHighP",p))$V1 # use CRC diff-high seem poor
feature <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_all/count_matrix/cfpeakCNN_b5_d50_p1_GSE110381.txt.diff.dcb"))$V1
#head(feature)

df <- data.frame(peak.id = unlist(sapply(strsplit(feature,"|",fixed=T),"[",4)),
                 RNA = unlist(sapply(strsplit(feature,"|",fixed=T),"[",2)))
df$RNA <- gsub("_rev|_for","",df$RNA,perl=T)
df$RNA <- factor(df$RNA,levels = c(rna,dna))
table(df$RNA)



## plot diff peak RNA type (pie plot)
#df <- df[df$peak.id %in% peak.ids[['top1000']],] # 50
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
  my_theme_pie
ggsave("TDCPs_pie.pdf",width = 7,height = 5)





## diff volcano plot of TDCP
dst <- "TCGA_small_diff3"
disease.label <- "COAD"
normal.label <- "LAML"

diff.plot <- list()
# tmp <- read.table(paste0("./output/",dst,"/exPeak_smallDomain_diff_",disease.label,"vs",normal.label,".diff"),sep = "\t",header = T,check.names = F)
tmp <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",dst,"/cfPeakCNN_GSE110381_smallDomain_diff_",disease.label,"vs",normal.label,".diff"),sep = "\t",header = T,check.names = F,stringsAsFactors = F)

tmp$minus.log10.pvalue <- -log10(tmp$pvalue)
tmp$log10.baseMean <- log10(tmp$baseMean)
table(feature %in% rownames(tmp)) # all True
tmp$TDCP <- "No"
tmp$TDCP[rownames(tmp) %in% feature] <- "Yes"
table(tmp$TDCP)
tmp <- tmp[order(tmp$TDCP),]
tmp$alpha <- maxmin.normalize.vec(tmp$minus.log10.pvalue)
diff.plot[[disease.label]] <- ggplot(tmp, aes(x=log10.baseMean,y=log2FoldChange)) + 
  geom_point(aes(fill=TDCP),shape=21,color="white",alpha=0.9) + 
  labs(title = paste0(disease.label," vs. ",normal.label)) +
  ylab("log2.FoldChange") +
  xlab("log10.Mean.CPM") +
  # ggraph::scale_fill_viridis(direction = -1) +
  scale_fill_manual(values = c("grey80","firebrick")) +
  geom_hline(yintercept = 0,linetype="dashed",color="salmon") + 
  theme_bw()+
  my_theme_volcano

# p <- ggpubr::ggarrange(plotlist = diff.plot, rel_widths = c(1,1,1), nrow = 1,  common.legend = T,align = "hv", axis = "b", legend = "right") #  labels = c("RBPs", "", "EV", "", "G4"),
ggsave(plot = diff.plot[[disease.label]], paste0("./TDCPs_log2FC_basemean.pdf"), width=6, height=6) # "_",sample






# TCGA_small_diff CRCvsNC AUROC
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
dst <- "TCGA_small_diff3"
NC.label <- "Blood"
CRC.label <- "COAD"

### read smp table
sample.table.tcga <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/archive/",dst,"/sample_table.txt"),check.names = F, sep = "\t", header = T)
#sample.table.tcga <- sample.table.tcga[,1:2]
sample.table.tcga <- sample.table.tcga[,c(2,5)]
colnames(sample.table.tcga) <- c("sample","group")
sample.table.tcga$sample <- gsub(".bam","",sample.table.tcga$sample)
rownames(sample.table.tcga) <- sample.table.tcga$sample
#table(sample.table.tcga$sample %in% colnames(cpm.tcga))
#sample.table.tcga <- sample.table.tcga[sample.table.tcga$sample %in% colnames(cpm.tcga),]
sample.table.tcga$group <- gsub("TCGA-","",sample.table.tcga$group)
table(sample.table.tcga$group)
sample.table.tcga[1:3,]
# sample.table.tcga$group <- gsub("COAD","CRC",sample.table.tcga$group)
sample.table.tcga$group <- gsub("LAML","Blood",sample.table.tcga$group)
sample.table.tcga <- sample.table.tcga[sample.table.tcga$group %in% c("COAD", "Blood"),] # this peak index can not classify other cancer types !
sample.table.tcga$group <- factor(sample.table.tcga$group,levels = c("Blood","COAD")) # , "PAAD", "PRAD"
sample.table.tcga <- sample.table.tcga[order(sample.table.tcga$group),]
sample.table.tcga <- sample.table.tcga[order(sample.table.tcga$group),]


#logcpm <- read.table("count_matrix",check.names = F,header = T) # cpm
#logcpm <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",dst,"/cfPeakCNN_GSE110381_smallDomain_diff_COADvsLAML.cpm"),check.names = F,header = T)  # log2cpm+1
logcpm <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_all/count_matrix/cfpeakCNN_b5_d50_p1_GSE110381_CPMrowsum.txt"),check.names = F,header = T)  # log2cpm+1
rownames(logcpm) <- logcpm$gene_id # feature
logcpm <- logcpm[,2:ncol(logcpm)] ## cpm, not log2cpm !!!
max(logcpm[,1])
cpm <- logcpm

# logcpm <- logcpm[feature,]
sample.table.tcga <- sample.table.tcga[sample.table.tcga$sample %in% colnames(logcpm),]
table(sample.table.tcga$group)
# Blood  COAD 
# 60    60 
#sample.table.tcga$group <- gsub("adenoma","CRC",sample.table.tcga$group)
#sample.table.tcga <- sample.table.tcga[order(sample.table.tcga$group),]
positive_samples <- sample.table.tcga[sample.table.tcga$group==CRC.label,"sample"]
negative_samples <- sample.table.tcga[sample.table.tcga$group==NC.label,"sample"]
samples <- c(positive_samples, negative_samples)
# sample.table.tcga <- sample.table.tcga[match(samples,sample.table.tcga$sample),]
# group <- c(rep("positive",length(positive_samples)),rep("negative",length(negative_samples)))
# logcpm <- logcpm[,sample.table.tcga$sample]



## build ref db
nrow(df)==length(feature) # TRUE !!!
rownames(df) <- feature
feature.lab.sort <- c(rna,dna)
feature.lab.sort <- feature.lab.sort[feature.lab.sort %in% as.character(unique(df$RNA))]
feature.lab.sort <- c(feature.lab.sort, "RNA", "DNA","all","highROC")

highROC <- c("pri_miRNA","lncRNA","mRNA","tucpRNA","intron","promoter","enhancer","repeats") # #only keep RNA species with roc >0.9 in TCGA

feature.list <- list()
mean.list <- list()
sd.list <- list()
for(i in feature.lab.sort){
  print(i)
  if(i=="all"){
    feature.list[[i]] <- rownames(df)
  }else if(i=="RNA"){
    feature.list[[i]] <- rownames(df)[df$RNA %in% rna]
  }else if(i=="DNA"){
    feature.list[[i]] <- rownames(df)[df$RNA %in% dna]
  }else if(i=="highROC"){
    feature.list[[i]] <- rownames(df)[df$RNA %in% highROC]
  }else{
    feature.list[[i]] <- rownames(df)[df$RNA==i]
  }
  
  #extract fixed coad&blood group
  blood.mat <- cpm[feature.list[[i]], negative_samples] # cpm, not log2cpm !!!
  coad.mat <- cpm[feature.list[[i]], positive_samples]
  mean.list[[i]] <- apply(blood.mat,1,mean) # length(mean.vec)
  sd.list[[i]] <- apply(cbind(coad.mat,blood.mat),1,sd) # blood mat has many 0 : Zhu - 2021 - NC - Tissue-specific cell-free DNA degradation...
  # sd.list <- (rep(sd(as.matrix(blood.mat)),length(mean.vec)))
  #print(table(sd.list==0))
  print( table(sd.list[[i]]==0) ) # need all FALSE !!! or else enlarge cohort when calculating sd (eg: add other groups)
}




### read count mat
#x <- paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",dst,"/cfPeakCNN_GSE110381_smallDomain_diff_COADvsLAML.cpm")
#logcpm <- read.table(x,check.names = F,header = T)  # log2cpm+1
#dim(logcpm)
logcpm <- log2(logcpm+1)
#max(logcpm[,1])


logcpm.sum.list <- list()
for(i in feature.lab.sort ){
  print(i)
  if(i=="all"){
    feature.list[[i]] <- rownames(df)
  }else if(i=="RNA"){
    feature.list[[i]] <- rownames(df)[df$RNA %in% rna]
  }else if(i=="DNA"){
    feature.list[[i]] <- rownames(df)[df$RNA %in% dna]
  }else if(i=="highROC"){
    feature.list[[i]] <- rownames(df)[df$RNA %in% highROC]
  }else{
    feature.list[[i]] <- rownames(df)[df$RNA==i]
  }
  logcpm.sum.list[[i]] <- getPeakIndex(feature.lab = i,logcpm = logcpm, sample.table = sample.table.tcga, prescale = T, mean.list = mean.list, sd.list = sd.list, zscore.cutoff = 100, onlyCountZscoreOutlier = F)
}


### plot boxplot
for(y in names(logcpm.sum.list)) {
  logcpm.sum.list[[y]]$feature.lab <- y
  logcpm.sum.list[[y]]$value <- logcpm.sum.list[[y]]$mean.log2cpm.scale # sum.log2cpm.scale, mean.log2cpm
  logcpm.sum.list[[y]]$group <- factor(logcpm.sum.list[[y]]$group,levels = c(NC.label,CRC.label))
} 
logcpm.sum.df <- as.data.frame(do.call(rbind,logcpm.sum.list))
logcpm.sum.df$feature.lab <- factor(logcpm.sum.df$feature.lab,levels = feature.lab.sort)

plotViolin(logcpm.sum.df,sig.size = 4) + # 2e-16
  facet_grid(~feature.lab)
ggsave(filename = paste0("./",dst,"_CRC-idx_box_all19.pdf"),width = 48,height = 8)
plotViolin(logcpm.sum.df[logcpm.sum.df$feature.lab=="all",],sig.size = 4) # 2e-16
ggsave(filename = paste0("./",dst,"_CRC-idx_box_all.pdf"),width = 8,height = 8)

plotMultiROC(logcpm.sum.list, direction = "<", RNA_colors = RNA_colors)
ggsave(filename = paste0("./",dst,"_CRC-idx_multiroc_all19.pdf"),width = 9,height = 9)

# #only keep RNA species with roc >0.9 in TCGA
# highROC <- c("pri_miRNA","lncRNA","mRNA","tucpRNA","intron","promoter","enhancer","repeats")






## GSE110381 CRCvsNC 
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
dst <- "GSE110381_diff"
NC.label <- "NC"
CRC.label <- "CRC"
### read smp table
sample.table <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/GSE110381_diff/sample_table.txt",check.names = F,header = T) 
# passQC <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/GSE110381_diff/sample_ids_diff_passQC2_rmBatch.txt",check.names = F,header = F)$V1
# sample.table <- sample.table[sample.table$sample %in% passQC,]
rownames(sample.table) <- sample.table$sample
sample.table$group <- gsub("adenoma","CRC",sample.table$group)

x <- paste0(pre,"/output/",dst,"/call_peak_all/count_matrix/expeakCNN_b5_d50_p1_CPMrowsum.txt") # use all samples: better auroc for this dst
#x <- paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",dst,"/filerSmp_exPeakCNN_smallDomain_diff_CRCvsNC.cpm")
logcpm <- read.table(x,check.names = F,header = T)  # not log2cpm+1
max(logcpm[,2])
rownames(logcpm) <- logcpm$gene_id
logcpm$gene_id <- NULL
max(as.matrix(logcpm))
logcpm <- log2(logcpm+1)


logcpm.sum.list <- list()
for(i in feature.lab.sort ){
  print(i)
  if(i=="all"){
    feature.list[[i]] <- rownames(df)
  }else if(i=="RNA"){
    feature.list[[i]] <- rownames(df)[df$RNA %in% rna]
  }else if(i=="DNA"){
    feature.list[[i]] <- rownames(df)[df$RNA %in% dna]
  }else if(i=="highROC"){
    feature.list[[i]] <- rownames(df)[df$RNA %in% highROC]
  }else{
    feature.list[[i]] <- rownames(df)[df$RNA==i]
  }
  logcpm.sum.list[[i]] <- getPeakIndex(feature.lab = i,logcpm = logcpm, sample.table = sample.table, prescale = T, mean.list = mean.list, sd.list = sd.list, zscore.cutoff = 100, onlyCountZscoreOutlier = F)
}


### plot boxplot
for(y in names(logcpm.sum.list)) {
  logcpm.sum.list[[y]]$feature.lab <- y
  logcpm.sum.list[[y]]$value <- logcpm.sum.list[[y]]$mean.log2cpm.scale # sum.log2cpm.scale, mean.log2cpm
  logcpm.sum.list[[y]]$group <- factor(logcpm.sum.list[[y]]$group,levels = c(NC.label,CRC.label))
} 
logcpm.sum.df <- as.data.frame(do.call(rbind,logcpm.sum.list))
logcpm.sum.df$feature.lab <- factor(logcpm.sum.df$feature.lab,levels = feature.lab.sort)

plotViolin(logcpm.sum.df,sig.size = 4) + # 2e-16
  facet_grid(~feature.lab)
ggsave(filename = paste0("./",dst,"_CRC-idx_box_all19.pdf"),width = 48,height = 8)
plotViolin(logcpm.sum.df[logcpm.sum.df$feature.lab=="all",],sig.size = 4) # 2e-16
ggsave(filename = paste0("./",dst,"_CRC-idx_box_all.pdf"),width = 8,height = 8)

plotMultiROC(logcpm.sum.list, direction = "<", RNA_colors = RNA_colors)
ggsave(filename = paste0("./",dst,"_CRC-idx_multiroc_all19.pdf"),width = 9,height = 9)





## PRJNA540919_diff CRCvsNC
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
dst <- "PRJNA540919_diff"
NC.label <- "NC"
CRC.label <- "CRC"
### read smp table
sample.table <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/PRJNA540919_diff/sample_table.txt",check.names = F,header = T)
colnames(sample.table) <- c("sample","group")
table(sample.table$group)

### read cpm
logcpm <- read.table(paste0(pre,"/output/",dst,"/call_peak_all/count_matrix/cfpeakCNN_b5_d50_p1_GSE110381_CPMrowsum.txt"),check.names = F,header = T)
# x <- paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",dst,"/filterSmp_exPeakCNN_smallDomain_diff_CRCvsNC.cpm")
# logcpm <- read.table(x,check.names = F,header = T)  # log2cpm+1
rownames(logcpm) <- logcpm$gene_id
logcpm$gene_id <- NULL
max(as.matrix(logcpm))
logcpm <- log2(logcpm+1)

logcpm.sum.list <- list()
for(i in feature.lab.sort ){
  print(i)
  if(i=="all"){
    feature.list[[i]] <- rownames(df)
  }else if(i=="RNA"){
    feature.list[[i]] <- rownames(df)[df$RNA %in% rna]
  }else if(i=="DNA"){
    feature.list[[i]] <- rownames(df)[df$RNA %in% dna]
  }else if(i=="highROC"){
    feature.list[[i]] <- rownames(df)[df$RNA %in% highROC]
  }else{
    feature.list[[i]] <- rownames(df)[df$RNA==i]
  }
  logcpm.sum.list[[i]] <- getPeakIndex(feature.lab = i,logcpm = logcpm, sample.table = sample.table, prescale = T, mean.list = mean.list, sd.list = sd.list, zscore.cutoff = 100, onlyCountZscoreOutlier = F)
}


### plot boxplot
for(y in names(logcpm.sum.list)) {
  logcpm.sum.list[[y]]$feature.lab <- y
  logcpm.sum.list[[y]]$value <- logcpm.sum.list[[y]]$mean.log2cpm.scale # sum.log2cpm.scale, mean.log2cpm
  logcpm.sum.list[[y]]$group <- factor(logcpm.sum.list[[y]]$group,levels = c(NC.label,CRC.label))
} 
logcpm.sum.df <- as.data.frame(do.call(rbind,logcpm.sum.list))
logcpm.sum.df$feature.lab <- factor(logcpm.sum.df$feature.lab,levels = feature.lab.sort)

plotViolin(logcpm.sum.df,sig.size = 4) + # 2e-16
  facet_grid(~feature.lab)
ggsave(filename = paste0("./",dst,"_CRC-idx_box_all19.pdf"),width = 48,height = 8)
plotViolin(logcpm.sum.df[logcpm.sum.df$feature.lab=="all",],sig.size = 4) # 2e-16
ggsave(filename = paste0("./",dst,"_CRC-idx_box_all.pdf"),width = 8,height = 8)
plotViolin(logcpm.sum.df[logcpm.sum.df$feature.lab %in% c("all","RNA","DNA"),],sig.size = 0) + # 2e-16
  facet_grid(~feature.lab)
ggsave(filename = paste0("./",dst,"_CRC-idx_box_all3.pdf"),width = 8,height = 8)

plotMultiROC(logcpm.sum.list, direction = "<", RNA_colors = RNA_colors) # poor if use onlyCountZscoreOutlier=3 !!!!
ggsave(filename = paste0("./",dst,"_CRC-idx_multiroc_all19.pdf"),width = 9,height = 9)
plotMultiROC(list("all"=logcpm.sum.list[["all"]],"RNA"=logcpm.sum.list[["RNA"]],"DNA"=logcpm.sum.list[["DNA"]]), direction = "<", RNA_colors = RNA_colors) # poor if use onlyCountZscoreOutlier=3 !!!!
ggsave(filename = paste0("./",dst,"_CRC-idx_multiroc_all3.pdf"),width = 9,height = 9)
# plotMultiROC(list("all"=logcpm.sum.list[["all"]],"pri_miRNA"=logcpm.sum.list[["pri_miRNA"]])) 
# ggsave(filename = paste0("./",dst,"_CRC-idx_multi2roc_2.pdf"),width = 9,height = 9)

#







# ## GSE71008 CRCvsNC (poor auc~0.5, even rm batch, only use this dst for multi-classifier tSNE )
# pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
# dst <- "GSE71008_diff"
# NC.label <- "NC"
# CRC.label <- "CRC"
# ### read smp table
# sample.table <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/GSE71008/sample_table.txt",check.names = F,header = T)
# table(sample.table$disease_status)
# sample.table <- sample.table[sample.table$disease_status %in% c("Colorectal","Healthy"),]
# sample.table$sample <- rownames(sample.table)
# #sample.table <- sample.table[,c("source","sample")]
# sample.table <- sample.table[,c("disease_type","sample")]
# colnames(sample.table) <- c("group","sample")
# sample.table$group <- gsub("Cancer","CRC",sample.table$group)
# sample.table$group <- gsub("Control","NC",sample.table$group)
# table(sample.table$group)
# 
# ### read cpm (all fig6 use rowsum cpm)
# # logcpm <- read.table(paste0(pre,"/output/",dst,"/call_peak_all/count_matrix/expeakCNN_b5_d50_p1_GSE110381_CPMrowsum.txt"),check.names = F,header = T)
# # rownames(logcpm) <- logcpm$gene_id
# # logcpm$gene_id <- NULL
# # max(as.matrix(logcpm))
# # logcpm <- log2(logcpm+1)
# logcpm <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/GSE71008_diff/allSmp_rmBatch_exPeakCNN_GSE110381_smallDomain_diff_CRCvsNC.cpm"),check.names = F,header = T)
# max(as.matrix(logcpm))
# 
# logcpm.sum.list <- list()
# for(i in feature.lab.sort ){
#   print(i)
#   if(i=="all"){
#     feature.list[[i]] <- rownames(df)
#   }else if(i=="RNA"){
#     feature.list[[i]] <- rownames(df)[df$RNA %in% rna]
#   }else if(i=="DNA"){
#     feature.list[[i]] <- rownames(df)[df$RNA %in% dna]
#   }else if(i=="highROC"){
#     feature.list[[i]] <- rownames(df)[df$RNA %in% highROC]
#   }else{
#     feature.list[[i]] <- rownames(df)[df$RNA==i]
#   }
#   logcpm.sum.list[[i]] <- getPeakIndex(feature.lab = i,logcpm = logcpm, sample.table = sample.table, prescale = T, mean.list = mean.list, sd.list = sd.list, zscore.cutoff = 100, onlyCountZscoreOutlier = F)
# }
# 
# 
# ### plot boxplot
# for(y in names(logcpm.sum.list)) {
#   logcpm.sum.list[[y]]$feature.lab <- y
#   logcpm.sum.list[[y]]$value <- logcpm.sum.list[[y]]$mean.log2cpm.scale # sum.log2cpm.scale, mean.log2cpm
#   logcpm.sum.list[[y]]$group <- factor(logcpm.sum.list[[y]]$group,levels = c(NC.label,CRC.label))
# } 
# logcpm.sum.df <- as.data.frame(do.call(rbind,logcpm.sum.list))
# logcpm.sum.df$feature.lab <- factor(logcpm.sum.df$feature.lab,levels = feature.lab.sort)
# 
# plotViolin(logcpm.sum.df,sig.size = 4) + # 2e-16
#   facet_grid(~feature.lab)
# ggsave(filename = paste0("./",dst,"_CRC-idx_box_all19.pdf"),width = 48,height = 8)
# plotViolin(logcpm.sum.df[logcpm.sum.df$feature.lab=="all",],sig.size = 4) # 2e-16
# ggsave(filename = paste0("./",dst,"_CRC-idx_box_all.pdf"),width = 8,height = 8)
# 
# plotMultiROC(logcpm.sum.list, direction = "<", RNA_colors = RNA_colors)
# ggsave(filename = paste0("./",dst,"_CRC-idx_multiroc_all19.pdf"),width = 9,height = 9)







## inhouse CRCvsNC
dst <- "WSQ_SMARTer_NEB_diff"
lib <- "SMARTer"  #"NEB", SMARTer
pnk <- "PNK+"  #"PNK-", "PNK+"
NC.label <- "NC"
CRC.label <- "CRC"
### read smp table
sample.table.inhouse <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/WSQ_SMARTer_NEB_diff/sample_table.txt",check.names = F, header = T)
sample.table.inhouse$sample <- paste0(sample.table.inhouse$sample, "_1")
rownames(sample.table.inhouse) <- sample.table.inhouse$sample
#sample.table.inhouse <- sample.table.inhouse[sample.table.inhouse$sample %in% colnames(cpm.inhouse),]
sample.table.inhouse$lib <- "NEB"
sample.table.inhouse$lib[grepl("smart",sample.table.inhouse$sample)] <- "SMARTer"
sample.table.inhouse$PNK <- "PNK-"
sample.table.inhouse$PNK[grepl("PNK",sample.table.inhouse$sample)] <- "PNK+"
sample.table.inhouse$group <- "NC"
sample.table.inhouse$group[grepl("CRC",sample.table.inhouse$sample)] <- "CRC"
sample.table.inhouse$group <- factor(sample.table.inhouse$group,levels = c("CRC","NC"))
sample.table.inhouse <- sample.table.inhouse[order(sample.table.inhouse$lib,sample.table.inhouse$group,sample.table.inhouse$PNK),]
table(sample.table.inhouse$lib,sample.table.inhouse$group,sample.table.inhouse$PNK)
table(sample.table.inhouse$group)
# sample.table.inhouse <- sample.table.inhouse[sample.table.inhouse$lib==lib,] # filter lib: SMARTer, NEB
# sample.table.inhouse <- sample.table.inhouse[sample.table.inhouse$PNK==pnk,] # filter PNK: PNK-,PNK+
# table(sample.table.inhouse$group)
# #table(sample.table.inhouse$sample %in% colnames(logcpm))

### read cpm (all fig6 use rowsum cpm)
#logcpm <- read.table(paste0(pre,"/output/",dst,"/call_peak_all/count_matrix/expeak_b5_d50_p1_GSE110381.txt.rowsumLogCPM.",lib),check.names = F,header = T) #log
logcpm <- read.table(paste0(pre,"/output/",dst,"/call_peak_all/count_matrix/cfpeakCNN_b5_d50_p1_GSE110381_CPMrowsum.txt"),check.names = F,header = T) #log
dim(logcpm)
max(logcpm[,2])
rownames(logcpm) <- logcpm$gene_id # gene_id,feature
logcpm$gene_id <- NULL # gene_id,feature
max(as.matrix(logcpm))
logcpm <- log2(logcpm+1)
sample.table <- sample.table.inhouse[sample.table.inhouse$sample %in% colnames(logcpm),]
logcpm <- logcpm[,sample.table$sample]


logcpm.sum.list <- list()
for(i in feature.lab.sort ){
  print(i)
  if(i=="all"){
    feature.list[[i]] <- rownames(df)
  }else if(i=="RNA"){
    feature.list[[i]] <- rownames(df)[df$RNA %in% rna]
  }else if(i=="DNA"){
    feature.list[[i]] <- rownames(df)[df$RNA %in% dna]
  }else if(i=="highROC"){
    feature.list[[i]] <- rownames(df)[df$RNA %in% highROC]
  }else{
    feature.list[[i]] <- rownames(df)[df$RNA==i]
  }
  logcpm.sum.list[[i]] <- getPeakIndex(feature.lab = i,logcpm = logcpm, sample.table = sample.table, prescale = T, mean.list = mean.list, sd.list = sd.list, zscore.cutoff = 100, onlyCountZscoreOutlier = T)
}


### plot boxplot
for(y in names(logcpm.sum.list)) {
  logcpm.sum.list[[y]]$feature.lab <- y
  logcpm.sum.list[[y]]$value <- logcpm.sum.list[[y]]$mean.log2cpm.scale # sum.log2cpm.scale, mean.log2cpm
  logcpm.sum.list[[y]]$group <- factor(logcpm.sum.list[[y]]$group,levels = c(NC.label,CRC.label))
} 
logcpm.sum.df <- as.data.frame(do.call(rbind,logcpm.sum.list))
logcpm.sum.df$feature.lab <- factor(logcpm.sum.df$feature.lab,levels = feature.lab.sort)

# plotViolin(logcpm.sum.df,sig.size = 4) + # 2e-16
#   facet_grid(~feature.lab)
# ggsave(filename = paste0("./",dst,"_CRC-idx_box_all19.pdf"),width = 48,height = 8)
# plotViolin(logcpm.sum.df[logcpm.sum.df$feature.lab=="all",],sig.size = 4) # 2e-16
# ggsave(filename = paste0("./",dst,"_CRC-idx_box_all.pdf"),width = 8,height = 8)
# 
# plotMultiROC(logcpm.sum.list, direction = "<",RNA_colors = RNA_colors)
# ggsave(filename = paste0("./",dst,"_CRC-idx_multiroc_all19.pdf"),width = 9,height = 9)
# 
# # plotMultiROC(list("all"=logcpm.sum.list[["all"]],"pri_miRNA"=logcpm.sum.list[["pri_miRNA"]])) 
# # ggsave(filename = paste0("./",dst,"_CRC-idx_multi2roc_2.pdf"),width = 9,height = 9)
# 
# #






## finds relation with clinical phenotype (tumor stage et al.)
dst <- "WSQ_SMARTer_NEB_diff"
NC.label <- "NC"
CRC.label <- "CRC"
### read smp table
sample.table.inhouse <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/WSQ_SMARTer_NEB_diff/sample_table.txt",check.names = F, header = T)
sample.table.inhouse$sample <- paste0(sample.table.inhouse$sample, "_1")
rownames(sample.table.inhouse) <- sample.table.inhouse$sample
#sample.table.inhouse <- sample.table.inhouse[sample.table.inhouse$sample %in% colnames(cpm.inhouse),]
sample.table.inhouse$lib <- "NEB"
sample.table.inhouse$lib[grepl("smart",sample.table.inhouse$sample)] <- "SMARTer"
sample.table.inhouse$PNK <- "PNK-"
sample.table.inhouse$PNK[grepl("PNK",sample.table.inhouse$sample)] <- "PNK+"
sample.table.inhouse$group <- "NC"
sample.table.inhouse$group[grepl("CRC",sample.table.inhouse$sample)] <- "CRC"
sample.table.inhouse$group <- factor(sample.table.inhouse$group,levels = c("CRC","NC"))
sample.table.inhouse <- sample.table.inhouse[order(sample.table.inhouse$lib,sample.table.inhouse$group,sample.table.inhouse$PNK),]
table(sample.table.inhouse$lib,sample.table.inhouse$group,sample.table.inhouse$PNK)

sample.table.inhouse$pid <- unlist(sapply(strsplit(sample.table.inhouse$sample,"_"),"[",1))
pid <- gsub("CRC-","",sample.table.inhouse$pid)
pid <- pid[pid != "NC"]
#clinical <- read.table("/BioII/lulab_b/baopengfei/projects/multi-omics-explore/meta/clinical_table.txt",header = T,sep = "\t")
#clinical$pid <- gsub("CRC-","",clinical$Lulab.ID)
#clinical <- clinical[clinical$pid %in% pid,]
clinical <- read.table("/lulabdata/baopengfei/shared_data/lulab_plasma/meta/all_CRC_patient_clinical.txt",header = T,sep = "\t")
clinical <- clinical[(clinical$PatientID %in% pid) & clinical$Cancer_type=="CRC",]
clinical$sample <- paste0(clinical$Cancer_type,"-",clinical$PatientID)

sample.table.inhouse$group2 <- clinical$Stage[match(sample.table.inhouse$pid,clinical$sample)]
sample.table.inhouse$group2[grepl(1,sample.table.inhouse$group2)] <- "I"
sample.table.inhouse$group2[grepl(2,sample.table.inhouse$group2)] <- "II"
sample.table.inhouse$group2[grepl(3,sample.table.inhouse$group2)] <- "III"
#sample.table.inhouse$group2[grepl(4,sample.table.inhouse$group2)] <- 4
sample.table.inhouse$group2[is.na(sample.table.inhouse$group2)] <- "NC"
sample.table.inhouse <- sample.table.inhouse[sample.table.inhouse$sample %in% colnames(logcpm),]
table(sample.table.inhouse$group2)
# I  II III  NC 
# 7  20  15  25

#cfpeakCNN_b5_d50_p1_GSE110381_CPMrowsum.txt
#archive/cfpeak_b5_d50_p1_GSE110381.txt.rowsumCPM
# logcpm <- read.table(paste0(pre,"/output/",dst,"/call_peak_all/count_matrix/cfpeakCNN_b5_d50_p1_GSE110381_CPMrowsum.txt"),check.names = F,header = T) #log
# dim(logcpm)
# rownames(logcpm) <- logcpm$gene_id # gene_id,feature
# logcpm$gene_id <- NULL # gene_id,feature
# max(as.matrix(logcpm))
# table(unlist(sapply(strsplit(rownames(logcpm),"|", fixed = T),"[",2)))
# # enhancer_for enhancer_rev   intron_for   intron_rev       lncRNA         mRNA        piRNA    pri_miRNA promoter_for promoter_rev 
# # 10359         7001         7272         8216        17530         9779          172         1767         4316         3272 
# # repeats_for  repeats_rev         rRNA       snoRNA        snRNA       srpRNA         tRNA      tucpRNA        Y_RNA 
# # 3461         2474          427          198          234           94          672         6820          193 
# 
# # logcpm2 <- read.table(paste0(pre,"/output/",dst,"/call_peak_all/count_matrix/cfpeakCNN_b5_d50_p1_GSE110381_CPMrowsum.txt"),check.names = F,header = T) #log
# # dim(logcpm2)
# # rownames(logcpm2) <- logcpm2$gene_id # gene_id,feature
# # logcpm2$gene_id <- NULL # gene_id,feature
# # max(as.matrix(logcpm2))
# # f2 <- unlist(sapply(strsplit(rownames(logcpm2),"|", fixed = T),"[",1))
# 
# # fea <- intersect(rownames(logcpm),rownames(logcpm2))
# # # mRNA repeats_for repeats_rev        tRNA 
# # # 81         592         248         672 
# # table(unlist(sapply(strsplit(fea,"|", fixed = T),"[",2)))
# # all(colnames(logcpm2)==colnames(logcpm))
# # cor.test(x = logcpm2[f2 %in% fea,10], y=logcpm[f1 %in% fea,10])
# # table(f2 %in% fea)

# logcpm <- log2(logcpm+1)
# table(sample.table.inhouse$sample %in% colnames(logcpm))
# sample.table.inhouse <- sample.table.inhouse[sample.table.inhouse$sample %in% colnames(logcpm),]
# logcpm <- logcpm[,sample.table.inhouse$sample]

# logcpm <- getPeakIndex(feature.lab = "all", logcpm = logcpm, sample.table = sample.table.inhouse, prescale = T, mean.list = mean.list, sd.list = sd.list, zscore.cutoff = 100, onlyCountZscoreOutlier = F)
logcpm.sum.df[1:3,]
logcpm.sum.df$group2 <- sample.table.inhouse$group2[match(logcpm.sum.df$sample,sample.table.inhouse$sample)]
logcpm.sum.df <- logcpm.sum.df[logcpm.sum.df$group2!="NC",]
logcpm.sum.df$group2 <- factor(logcpm.sum.df$group2,levels = c(c('I','II','III'))) # c('1','2','3') # NC.label
logcpm.sum.df[1:3,]
# logcpm.sum <- logcpm %>%
#   group_by(sample,group2) %>%
#   summarise(sum.log2cpm=sum(sum.log2cpm.scale), mean.log2cpm=sum(mean.log2cpm.scale))
# # summarise(sum.log2cpm=mean(log2cpm,trim=0.05))
#NC.label <- "NC"
str(logcpm.sum.df)
tmp <- logcpm.sum.df[logcpm.sum.df$feature.lab %in% c(rna,dna),]
ggplot(tmp)+
  geom_boxplot(aes(x=group2,y=mean.log2cpm.scale,fill=group2))+
  scale_fill_manual(name="Stage",values = c("firebrick2","salmon","firebrick4"))+#"grey50",
  facet_grid(.~feature.lab,scales = "free")+
  ggpubr::stat_compare_means(
    # label.x.npc = 0.5,size =12,step.increase = 0.08,#paired = TRUE,
    aes(x=group2,y=mean.log2cpm.scale,label = ..p.signif..), # p.signif,..p.format.. group = Group,
    symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05,1),
                     symbols = c( "***", "**", "*", "ns")),
    ref.group = "I", # NC.label
    hide.ns=T, size=0,
    # method.args = list(alternative="greater"),
    method = "wilcox.test"
  ) +
  theme_bw() +
  ylab("Scaled peak-index") +
  theme(
    aspect.ratio = 2,
        plot.title = element_text(size = 24,color="black",hjust = 0.5),
        axis.title = element_text(size = 24,color ="black"),
        axis.text = element_text(size= 24,color = "black"),
        axis.text.x = element_text(size= 24,color = "black"),
        #panel.grid=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(), #element_line(color = "grey50",linetype = "dashed"), #size= 1,
        panel.grid.minor.y = element_blank(),
        #panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        legend.position = "right",#c(.25,.6),
        legend.text = element_text(size= 24),
        legend.title= element_text(size= 24),
        strip.text.y = element_blank(),
        strip.text.x = element_text(size=24)
  )
#ggsave(filename = paste0(dst,"_bystage_CRC-idx_box_2.pdf"),width = 8,height = 8)
ggsave(filename = paste0(dst,"_bystage_CRC-idx_box_2.pdf"),width = 28,height = 5)
# still opposite using all samples
# #


# #finds relation with clinical phenotype (tumor stage et al.)
# #to finish !!!!!
# sample.table.tcga$pid <- substr(sample.table.tcga$sample,1,15)
# table(sample.table.tcga$pid %in% clinical$data_id)
# pid <- gsub("CRC-","",sample.table.tcga$pid)
# pid <- pid[pid != "NC"]
# 
# clinical <- rio::import("/BioII/lulab_b/baopengfei/projects/xena/TCGA-LIHC/TCGA-LIHC_all.xlsx")
# clinical$submitter_id <- gsub(".bam","",clinical$submitter_id)
# clinical <- clinical[,c("patient_id","submitter_id","DFI","DFI.time","gender","race","tumor_stage","type")]
# colnames(clinical) <- c("patient_id","data_id","DFI","DFI_time","gender","race","tumor_stage","group")
# clinical$data_id <- substr(clinical$data_id,1,15)
# 
# clinical <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/TCGA_small/COAD_meta.txt",header = T,sep = "\t")
# clinical <- clinical[(clinical$PatientID %in% pid) & clinical$Cancer_type=="CRC",]
# clinical$sample <- paste0(clinical$Cancer_type,"-",clinical$PatientID)
# 
# 
# sample.table.tcga$group2 <- clinical$Stage[match(sample.table.tcga$pid,clinical$sample)]
# sample.table.tcga$group2[grepl(2,sample.table.tcga$group2)] <- 2
# sample.table.tcga$group2[grepl(3,sample.table.tcga$group2)] <- 3
# #sample.table.tcga$group2[grepl(4,sample.table.tcga$group2)] <- 4
# sample.table.tcga$group2[is.na(sample.table.tcga$group2)] <- "NC"
# 
# logcpm$group2 <- sample.table.tcga$group2[match(logcpm$sample,sample.table.tcga$sample)]
# logcpm.sum <- logcpm %>% 
#   group_by(sample,group2) %>% 
#   summarise(sum.log2cpm=sum(log2cpm))
# # summarise(sum.log2cpm=mean(log2cpm,trim=0.05))
# NC.label <- "NC"
# table(logcpm.sum$group2)
# logcpm.sum$group2 <- factor(logcpm.sum$group2,levels = c(NC.label,c('1','2','3')))
# ggplot(logcpm.sum)+
#   geom_boxplot(aes(x=group2,y=sum.log2cpm,fill=group2))+
#   scale_fill_manual(name="Group",values = c("grey50","firebrick","firebrick","firebrick","firebrick"))+#colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(unique(logcpm.sum$group))))+)
#   # facet_grid(.~name,scales = "free")+
#   ggpubr::stat_compare_means(
#     # label.x.npc = 0.5,size =12,step.increase = 0.08,#paired = TRUE,
#     aes(x=group2,y=sum.log2cpm,label = ..p.signif..), # p.signif, group = Group,
#     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01,1),
#                      symbols = c( "***", "**", "*", "ns")),
#     ref.group = NC.label,hide.ns=T, size=12,
#     # method.args = list(alternative="greater"),
#     method = "wilcox.test"
#   ) +
#   theme_bw() + 
#   theme(aspect.ratio = 2,
#         plot.title = element_text(size = 24,color="black",hjust = 0.5),
#         axis.title = element_text(size = 24,color ="black"), 
#         axis.text = element_text(size= 24,color = "black"),
#         axis.text.x = element_text(size= 24,color = "black",angle = 90,vjust = 0.5,hjust = 1),
#         #panel.grid=element_blank(),
#         panel.grid.major.x=element_blank(),
#         panel.grid.minor.x = element_blank(),
#         panel.grid.major.y = element_blank(), #element_line(color = "grey50",linetype = "dashed"), #size= 1,
#         panel.grid.minor.y = element_blank(),
#         #panel.grid.minor.y = element_blank(),
#         panel.border = element_blank(),
#         legend.position = "right",#c(.25,.6),
#         legend.text = element_text(size= 24),
#         legend.title= element_text(size= 24),
#         strip.text.y = element_blank(),
#         strip.text.x = element_text(size=24)
#   )
# ggsave(filename = paste0("./TCGA_small_diff_",lib,"_",pnk,"_CRC-idx_box_bystage.pdf"),width = 8,height = 12)


#








# LIHC-peak-index TCGA-LIHC_small+GSE123972 (GSE110381 peak) (cfpeak suppl Fig? fail)  ------------------------------------

## TCGA 
dst <- "TCGA-LIHC_small_diff"
#p <- 0.01
#feature <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_all/count_matrix/expeakCNN_b5_d50_p1_GSE110381.txt.diff.tissueHighP",p))$V1 # use CRC diff-high seem poor
feature <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_all/count_matrix/cfpeakCNN_b5_d50_p1_GSE110381.txt.diff.dcb"))$V1
#head(feature)

df <- data.frame(peak.id = unlist(sapply(strsplit(feature,"|",fixed=T),"[",4)),
                 RNA = unlist(sapply(strsplit(feature,"|",fixed=T),"[",2)))
df$RNA <- gsub("_rev|_for","",df$RNA,perl=T)
df$RNA <- factor(df$RNA,levels = c(rna,dna))
table(df$RNA)



## plot diff peak RNA type (pie plot)
#df <- df[df$peak.id %in% peak.ids[['top1000']],] # 50
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
  my_theme_pie
ggsave("TDCPs_pie.pdf",width = 7,height = 5)





## diff volcano plot of TDCP
dst <- "TCGA-LIHC_small_diff"
disease.label <- "LIHC"
normal.label <- "LAML"

diff.plot <- list()
# tmp <- read.table(paste0("./output/",dst,"/exPeak_smallDomain_diff_",disease.label,"vs",normal.label,".diff"),sep = "\t",header = T,check.names = F)
tmp <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",dst,"/cfPeakCNN_GSE110381_smallDomain_diff_",disease.label,"vs",normal.label,".diff"),sep = "\t",header = T,check.names = F,stringsAsFactors = F)

tmp$minus.log10.pvalue <- -log10(tmp$pvalue)
tmp$log10.baseMean <- log10(tmp$baseMean)
table(feature %in% rownames(tmp)) # all True
tmp$TDCP <- "No"
tmp$TDCP[rownames(tmp) %in% feature] <- "Yes"
table(tmp$TDCP)
tmp <- tmp[order(tmp$TDCP),]
tmp$alpha <- maxmin.normalize.vec(tmp$minus.log10.pvalue)
diff.plot[[disease.label]] <- ggplot(tmp, aes(x=log10.baseMean,y=log2FoldChange)) + 
  geom_point(aes(fill=TDCP),shape=21,color="white",alpha=0.9) + 
  labs(title = paste0(disease.label," vs. ",normal.label)) +
  ylab("log2.FoldChange") +
  xlab("log10.Mean.CPM") +
  # ggraph::scale_fill_viridis(direction = -1) +
  scale_fill_manual(values = c("grey80","firebrick")) +
  geom_hline(yintercept = 0,linetype="dashed",color="salmon") + 
  theme_bw()+
  my_theme_volcano

# p <- ggpubr::ggarrange(plotlist = diff.plot, rel_widths = c(1,1,1), nrow = 1,  common.legend = T,align = "hv", axis = "b", legend = "right") #  labels = c("RBPs", "", "EV", "", "G4"),
ggsave(plot = diff.plot[[disease.label]], paste0("./TDCPs_log2FC_basemean.pdf"), width=6, height=6) # "_",sample






# TCGA-LIHC_small_diff LIHCvsNC AUROC
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
dst <- "TCGA-LIHC_small_diff"
NC.label <- "Blood"
CRC.label <- "LIHC" #  "COAD"

### read smp table
sample.table.tcga <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/",dst,"/sample_table.txt"),check.names = F, sep = "\t", header = T)
colnames(sample.table.tcga) <- c("sample","group")
sample.table.tcga$sample <- gsub(".bam","",sample.table.tcga$sample)
rownames(sample.table.tcga) <- sample.table.tcga$sample
sample.table.tcga$group <- gsub("TCGA-","",sample.table.tcga$group)
table(sample.table.tcga$group)
sample.table.tcga[1:3,]
sample.table.tcga$group <- gsub("LAML","Blood",sample.table.tcga$group)
sample.table.tcga <- sample.table.tcga[sample.table.tcga$group %in% c(NC.label, CRC.label),] # this peak index can not classify other cancer types !
sample.table.tcga$group <- factor(sample.table.tcga$group,levels = c(NC.label,CRC.label)) # , "PAAD", "PRAD"
sample.table.tcga <- sample.table.tcga[order(sample.table.tcga$group),]

logcpm <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_all/count_matrix/cfpeakCNN_b5_d50_p1_GSE110381_CPMrowsum.txt"),check.names = F,header = T)  # log2cpm+1
rownames(logcpm) <- logcpm$gene_id # feature
logcpm <- logcpm[,2:ncol(logcpm)] ## cpm, not log2cpm !!!
max(logcpm[,1])
cpm <- logcpm

# logcpm <- logcpm[feature,]
sample.table.tcga <- sample.table.tcga[sample.table.tcga$sample %in% colnames(logcpm),]
table(sample.table.tcga$group)
positive_samples <- sample.table.tcga[sample.table.tcga$group==CRC.label,"sample"]
negative_samples <- sample.table.tcga[sample.table.tcga$group==NC.label,"sample"]
samples <- c(positive_samples, negative_samples)
# sample.table.tcga <- sample.table.tcga[match(samples,sample.table.tcga$sample),]
# group <- c(rep("positive",length(positive_samples)),rep("negative",length(negative_samples)))
# logcpm <- logcpm[,sample.table.tcga$sample]



## build ref db
nrow(df)==length(feature) # TRUE !!!
rownames(df) <- feature
feature.lab.sort <- c(rna,dna)
feature.lab.sort <- feature.lab.sort[feature.lab.sort %in% as.character(unique(df$RNA))]
feature.lab.sort <- c(feature.lab.sort, "RNA", "DNA","all") # 

feature.list <- list()
mean.list <- list()
sd.list <- list()
for(i in feature.lab.sort){
  print(i)
  if(i=="all"){
    feature.list[[i]] <- rownames(df)
  }else if(i=="RNA"){
    feature.list[[i]] <- rownames(df)[df$RNA %in% rna]
  }else if(i=="DNA"){
    feature.list[[i]] <- rownames(df)[df$RNA %in% dna]
  }else if(i=="highROC"){
    feature.list[[i]] <- rownames(df)[df$RNA %in% highROC]
  }else{
    feature.list[[i]] <- rownames(df)[df$RNA==i]
  }
  
  #extract fixed coad&blood group
  blood.mat <- cpm[feature.list[[i]], negative_samples] # cpm, not log2cpm !!!
  coad.mat <- cpm[feature.list[[i]], positive_samples]
  mean.list[[i]] <- apply(blood.mat,1,mean) # length(mean.vec)
  sd.list[[i]] <- apply(cbind(coad.mat,blood.mat),1,sd) # blood mat has many 0 : Zhu - 2021 - NC - Tissue-specific cell-free DNA degradation...
  # sd.list <- (rep(sd(as.matrix(blood.mat)),length(mean.vec)))
  #print(table(sd.list==0))
  print( table(sd.list[[i]]==0) ) # need all FALSE !!! or else enlarge cohort when calculating sd (eg: add other groups)
}




### read count mat
#x <- paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",dst,"/cfPeakCNN_GSE110381_smallDomain_diff_COADvsLAML.cpm")
#logcpm <- read.table(x,check.names = F,header = T)  # log2cpm+1
#dim(logcpm)
logcpm <- log2(logcpm+1)
#max(logcpm[,1])


logcpm.sum.list <- list()
for(i in feature.lab.sort ){
  print(i)
  if(i=="all"){
    feature.list[[i]] <- rownames(df)
  }else if(i=="RNA"){
    feature.list[[i]] <- rownames(df)[df$RNA %in% rna]
  }else if(i=="DNA"){
    feature.list[[i]] <- rownames(df)[df$RNA %in% dna]
  }else if(i=="highROC"){
    feature.list[[i]] <- rownames(df)[df$RNA %in% highROC]
  }else{
    feature.list[[i]] <- rownames(df)[df$RNA==i]
  }
  logcpm.sum.list[[i]] <- getPeakIndex(feature.lab = i,logcpm = logcpm, sample.table = sample.table.tcga, prescale = T, mean.list = mean.list, sd.list = sd.list, zscore.cutoff = 100, onlyCountZscoreOutlier = F)
}


### plot boxplot
for(y in names(logcpm.sum.list)) {
  logcpm.sum.list[[y]]$feature.lab <- y
  logcpm.sum.list[[y]]$value <- logcpm.sum.list[[y]]$mean.log2cpm.scale # sum.log2cpm.scale, mean.log2cpm
  logcpm.sum.list[[y]]$group <- factor(logcpm.sum.list[[y]]$group,levels = c(NC.label,CRC.label))
} 
logcpm.sum.df <- as.data.frame(do.call(rbind,logcpm.sum.list))
logcpm.sum.df$feature.lab <- factor(logcpm.sum.df$feature.lab,levels = feature.lab.sort)

plotViolin(logcpm.sum.df,sig.size = 4) + # 2e-16
  facet_grid(~feature.lab)
ggsave(filename = paste0("./",dst,"_CRC-idx_box_all19.pdf"),width = 48,height = 8)
plotViolin(logcpm.sum.df[logcpm.sum.df$feature.lab=="all",],sig.size = 4) # 2e-16
ggsave(filename = paste0("./",dst,"_CRC-idx_box_all.pdf"),width = 8,height = 8)

plotMultiROC(logcpm.sum.list, direction = "<", RNA_colors = RNA_colors)
ggsave(filename = paste0("./",dst,"_CRC-idx_multiroc_all19.pdf"),width = 9,height = 9)








## GSE123972_diff CRCvsNC
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
dst <- "GSE123972_diff"
NC.label <- "NC"
CRC.label <- "HCC"

# sample.table <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/",dst,"/sample_classes.txt"),check.names = F,header = T)
# sample.table$group <- "LiverDisease"
# sample.table$group[sample.table$label=="Normal"] <- "NC"
# sample.table$group[grepl("stage",sample.table$label)] <- "LIHC"
# colnames(sample.table)[1] <- "sample"
# sample.table <- sample.table[,c("sample","group")]
# write.table(sample.table,paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/",dst,"/sample_table.txt"),quote = F,row.names = F,col.names = T,sep = "\t" )

### read smp table
sample.table <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/",dst,"/sample_table.txt"),check.names = F,header = T)
colnames(sample.table) <- c("sample","group")
sample.table$group[sample.table$group=="LiverDisease"] <- "NC"
table(sample.table$group)

### read cpm
logcpm <- read.table(paste0(pre,"/output/",dst,"/call_peak_all/count_matrix/cfpeakCNN_b5_d50_p1_GSE110381_CPMrowsum.txt"),check.names = F,header = T)
# x <- paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",dst,"/filterSmp_exPeakCNN_smallDomain_diff_CRCvsNC.cpm")
# logcpm <- read.table(x,check.names = F,header = T)  # log2cpm+1
rownames(logcpm) <- logcpm$gene_id
logcpm$gene_id <- NULL
max(as.matrix(logcpm))
logcpm <- log2(logcpm+1)

logcpm.sum.list <- list()
for(i in feature.lab.sort ){
  print(i)
  if(i=="all"){
    feature.list[[i]] <- rownames(df)
  }else if(i=="RNA"){
    feature.list[[i]] <- rownames(df)[df$RNA %in% rna]
  }else if(i=="DNA"){
    feature.list[[i]] <- rownames(df)[df$RNA %in% dna]
  }else if(i=="highROC"){
    feature.list[[i]] <- rownames(df)[df$RNA %in% highROC]
  }else{
    feature.list[[i]] <- rownames(df)[df$RNA==i]
  }
  logcpm.sum.list[[i]] <- getPeakIndex(feature.lab = i,logcpm = logcpm, sample.table = sample.table, prescale = T, mean.list = mean.list, sd.list = sd.list, zscore.cutoff = 100, onlyCountZscoreOutlier = F)
}

### plot boxplot
for(y in names(logcpm.sum.list)) {
  logcpm.sum.list[[y]]$feature.lab <- y
  logcpm.sum.list[[y]]$value <- logcpm.sum.list[[y]]$mean.log2cpm.scale # sum.log2cpm.scale, mean.log2cpm
  logcpm.sum.list[[y]]$group <- factor(logcpm.sum.list[[y]]$group,levels = c(NC.label,"LIHC")) # ,"LiverDisease"
} 
logcpm.sum.df <- as.data.frame(do.call(rbind,logcpm.sum.list))
logcpm.sum.df$feature.lab <- factor(logcpm.sum.df$feature.lab,levels = feature.lab.sort)

plotViolin(logcpm.sum.df,sig.size = 4,colors = c("grey50","firebrick")) + # 2e-16 # "grey20",
  facet_grid(~feature.lab)
ggsave(filename = paste0("./",dst,"_CRC-idx_box_all19.pdf"),width = 48,height = 8)
plotViolin(logcpm.sum.df[logcpm.sum.df$feature.lab=="all",],sig.size = 4,colors = c("grey50","firebrick")) # 2e-16
ggsave(filename = paste0("./",dst,"_CRC-idx_box_all.pdf"),width = 8,height = 8)
plotViolin(logcpm.sum.df[logcpm.sum.df$feature.lab %in% c("all","RNA","DNA"),],sig.size = 6) + # 2e-16
  facet_grid(~feature.lab)
ggsave(filename = paste0("./",dst,"_CRC-idx_box_all3.pdf"),width = 8,height = 8)

plotMultiROC(logcpm.sum.list, direction = "<", RNA_colors = RNA_colors) # poor if use onlyCountZscoreOutlier=3 !!!!
ggsave(filename = paste0("./",dst,"_CRC-idx_multiroc_all19.pdf"),width = 9,height = 9)
plotMultiROC(list("all"=logcpm.sum.list[["all"]],"RNA"=logcpm.sum.list[["RNA"]],"DNA"=logcpm.sum.list[["DNA"]]), direction = "<", RNA_colors = RNA_colors) # poor if use onlyCountZscoreOutlier=3 !!!!
ggsave(filename = paste0("./",dst,"_CRC-idx_multiroc_all3.pdf"),width = 9,height = 9)
# plotMultiROC(list("all"=logcpm.sum.list[["all"]],"pri_miRNA"=logcpm.sum.list[["pri_miRNA"]])) 
# ggsave(filename = paste0("./",dst,"_CRC-idx_multi2roc_2.pdf"),width = 9,height = 9)

#all in-house data from PKU has opposite trend ...



# 3tumor-peak-index TCGA_small --> GSE71008 (GSE71008 peak) (cfpeak suppl Fig? fail)  ------------------------------------

# TCGA_small_diff CRCvsNC AUROC
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
dst <- "TCGA_small_diff3"

### read smp table
sample.table.tcga <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/archive/",dst,"/sample_table.txt"),check.names = F, sep = "\t", header = T)
#sample.table.tcga <- sample.table.tcga[,1:2]
sample.table.tcga <- sample.table.tcga[,c(2,5)]
colnames(sample.table.tcga) <- c("sample","group")
sample.table.tcga$sample <- gsub(".bam","",sample.table.tcga$sample)
rownames(sample.table.tcga) <- sample.table.tcga$sample
sample.table.tcga$group <- gsub("TCGA-","",sample.table.tcga$group)
table(sample.table.tcga$group)
sample.table.tcga <- sample.table.tcga[order(sample.table.tcga$group),]

#logcpm <- read.table("count_matrix",check.names = F,header = T) # cpm
#logcpm <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",dst,"/cfPeakCNN_GSE110381_smallDomain_diff_COADvsLAML.cpm"),check.names = F,header = T)  # log2cpm+1
#logcpm <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_all/count_matrix/cfpeakCNN_b5_d50_p1_GSE110381_CPMrowsum.txt"),check.names = F,header = T)  # log2cpm+1
#logcpm <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_all/count_matrix/cfpeakCNN_b5_d50_p1_GSE71008_CPMrowsum.txt"),check.names = F,header = T)  # log2cpm+1
logcpm <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",dst,"/OvR_cfPeakCNN_GSE71008_smallDomain_diff_COADvsR.cpm"),row.names = 1,check.names = F,header = T)
# rownames(logcpm) <- logcpm$gene_id # feature
max(logcpm[,3])
cpm <- 2^logcpm-1 ## zscore are calculated in cpm, not log2cpm !!!
#max(logcpm[,1])
# logcpm <- logcpm[feature,]
sample.table.tcga <- sample.table.tcga[sample.table.tcga$sample %in% colnames(logcpm),]
table(sample.table.tcga$group)
# Blood  COAD  PAAD  PRAD 
# 60    60    60    60



tumors <- c("COAD","LAML","PAAD","PRAD")
#feature.lab.sort <- tumors
feature.list <- list()
mean.list <- list()
sd.list <- list()
for(disease in tumors){
  print(disease)
## build ref db
  positive_samples <- sample.table.tcga[sample.table.tcga$group==disease,"sample"]
  negative_samples <- sample.table.tcga[sample.table.tcga$group!=disease,"sample"]

  feature.list[[disease]] <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_all/count_matrix/cfpeakCNN_b5_d50_p1_GSE71008.txt.diff.dcb.",disease,"vsR"))$V1

  #extract fixed coad&blood group
  blood.mat <- cpm[feature.list[[disease]], negative_samples] # cpm, not log2cpm !!!
  coad.mat <- cpm[feature.list[[disease]], positive_samples]
  mean.list[[disease]] <- apply(blood.mat,1,mean) # zscore are calculated in cpm, not log2cpm
  sd.list[[disease]] <- apply(cbind(coad.mat,blood.mat),1,sd) # blood mat has many 0 : Zhu - 2021 - NC - Tissue-specific cell-free DNA degradation...
  # sd.list <- (rep(sd(as.matrix(blood.mat)),length(mean.vec)))
  #print(table(sd.list==0))
  # print( table(sd.list[[i]]==0) ) # need all FALSE !!! or else enlarge cohort when calculating sd (eg: add other groups)
}
length(mean.list[["COAD"]])
length(sd.list[["LAML"]])
#not need keep same num of TDCPs in diff cancer
#


logcpm.sum.list <- list()
for(disease in tumors){
  # positive_samples <- sample.table.tcga[sample.table.tcga$group==disease,"sample"]
  # negative_samples <- sample.table.tcga[sample.table.tcga$group!=disease,"sample"]
  # samples <- c(positive_samples, negative_samples)
  logcpm.sum.list[[disease]] <- getPeakIndex(feature.lab = disease, logcpm = logcpm, sample.table = sample.table.tcga, prescale = T, mean.list = mean.list, sd.list = sd.list, zscore.cutoff = 100, onlyCountZscoreOutlier = 3)
}
#dim(logcpm.sum.list[["COAD"]])
#names(mean.list)
tmp <- logcpm.sum.list[["COAD"]]
#


##plot peak index score 
logcpm.sum.df <- as.data.frame(do.call(rbind,logcpm.sum.list))
table(logcpm.sum.df$group==logcpm.sum.df$peak.precursor)

NC.label <- "LAML"
logcpm.sum.df$group <- factor(logcpm.sum.df$group, levels = c("LAML","COAD","PAAD","PRAD"))
logcpm.sum.df$peak.precursor <- factor(logcpm.sum.df$peak.precursor, levels = c("LAML","COAD","PAAD","PRAD"))
p1 <- ggplot(logcpm.sum.df, aes(x=group,y=mean.log2cpm,fill=group))+ #sum.log2cpm.scale
  geom_violin(alpha=0.7)+
  geom_boxplot(width=0.1)+
  scale_fill_manual(name="Group",values = c("grey50","salmon","steelblue","purple"))+#colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(unique(logcpm.sum$group))))+)
  facet_grid(.~peak.precursor,scales = "free")+
  ylab("Scaled peak-index") +
  ggpubr::stat_compare_means(
    label.x.npc = "middle", label.y.npc = "top",
    #size =12,step.increase = 0.08,#paired = TRUE,
    aes(x=group,y=mean.log2cpm,
        label = ..p.format..), # p.signif, p.format, group = Group,
    symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01,1),
                     symbols = c( "***", "**", "*", "ns")),
    # method.args = list(alternative="greater"),
    # ref.group = peak.precursor,
    hide.ns=F, size=5,
    method = "wilcox.test"
  ) +
  theme_bw() +
  my_theme_box
p1
ggsave(paste0(dst,"_multiclass_box.pdf"),width = 24,height = 8)
# p2 <- ggplot(logcpm.sum.df, aes(x=peak.precursor,y=mean.log2cpm,fill=peak.precursor))+ #sum.log2cpm.scale
#   geom_violin(alpha=0.7)+
#   geom_boxplot(width=0.1)+
#   scale_fill_manual(name="Group",values = c("grey50","salmon","steelblue","purple"))+#colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(unique(logcpm.sum$group))))+)
#   facet_grid(.~group,scales = "free")+
#   ylab("Scaled peak-index") +
#   ggpubr::stat_compare_means(
#     label.x.npc = "middle", label.y.npc = "top",
#     #size =12,step.increase = 0.08,#paired = TRUE,
#     aes(x=group,y=mean.log2cpm,
#         label = ..p.format..), # p.signif, p.format, group = Group,
#     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01,1),
#                      symbols = c( "***", "**", "*", "ns")),
#     # method.args = list(alternative="greater"),
#     # ref.group = peak.precursor,
#     hide.ns=F, size=5,
#     method = "wilcox.test"
#   ) +
#   theme_bw() +
#   my_theme_box
# p2
# ggsave(paste0(dst,"_multiclass_byGroup_box.pdf"),width = 24,height = 8)



##plot confusion mat
logcpm.sum.tbl <- dplyr::as_tibble(logcpm.sum.df) %>%
  dplyr::group_by(sample) %>%
  dplyr::top_n(n=1,wt=mean.log2cpm) # sum.log2cpm, mean.log2cpm
logcpm.sum.tbl <- as.data.frame(logcpm.sum.tbl)

p <- plotConfusionMat(Actual = (logcpm.sum.tbl$group), Predict = (logcpm.sum.tbl$peak.precursor))
p
table(logcpm.sum.tbl$group,logcpm.sum.tbl$peak.precursor)
ggsave(plot = p, filename = paste0(dst,"_multiclass_confusion_mat.pdf"),width = 8,height = 8)

#




## GSE71008 CRCvsNC (poor, even rm batch, only use this dst for multi-classifier )
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
dst <- "GSE71008_diff"
# NC.label <- "NC"
# CRC.label <- "CRC"
### read smp table
sample.table <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/GSE71008_diff/sample_table.txt",check.names = F,header = T)
table(sample.table$disease_status)
# sample.table <- sample.table[sample.table$disease_status %in% c("Colorectal","Healthy"),]
sample.table$sample <- rownames(sample.table)
#sample.table <- sample.table[,c("source","sample")]
sample.table <- sample.table[,c("disease_status","sample")]
colnames(sample.table) <- c("group","sample")
sample.table$group <- gsub("Colorectal","COAD",sample.table$group)
sample.table$group <- gsub("Healthy","LAML",sample.table$group)
sample.table$group <- gsub("Pancreatic","PAAD",sample.table$group)
sample.table$group <- gsub("Prostate","PRAD",sample.table$group)
table(sample.table$group)
# Colorectal    Healthy Pancreatic   Prostate 
# 100         50          6         36
# COAD LAML PAAD PRAD 
# 100   50    6   36 

### read cpm (all fig6 use rowsum cpm)
# logcpm <- read.table(paste0(pre,"/output/",dst,"/call_peak_all/count_matrix/expeakCNN_b5_d50_p1_CPMrowsum.txt"),check.names = F,header = T)
# rownames(logcpm) <- logcpm$gene_id
# logcpm$gene_id <- NULL
# max(as.matrix(logcpm))
# logcpm <- log2(logcpm+1)
logcpm <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",dst,"/allSmp_rmBatch_OvR_exPeakCNN_smallDomain_diff_CRCvsR.cpm"),check.names = F,header = T)




#logcpm.sum <- getPeakIndex(logcpm = logcpm, sample.table = sample.table, prescale = T, zscore.cutoff = 100) # lower zscore.cutoff lead to lower auc in this dataset !!!  
logcpm.sum.list <- list()
tumors <- unique(sample.table$group)
for(disease in tumors){
  # positive_samples <- sample.table.tcga[sample.table.tcga$group==disease,"sample"]
  # negative_samples <- sample.table.tcga[sample.table.tcga$group!=disease,"sample"]
  # samples <- c(positive_samples, negative_samples)
  logcpm.sum.list[[disease]] <- getPeakIndex(feature.lab = disease, logcpm = logcpm, sample.table = sample.table, prescale = T, mean.list = mean.list, sd.list = sd.list, zscore.cutoff = 100, onlyCountZscoreOutlier = 3)
}
#dim(logcpm.sum.list[["COAD"]])
#names(mean.list)


logcpm.sum.df <- as.data.frame(do.call(rbind,logcpm.sum.list))
table(logcpm.sum.df$group==logcpm.sum.df$peak.precursor)

NC.label <- "LAML"
logcpm.sum.df$group <- factor(logcpm.sum.df$group, levels = c("LAML","COAD","PAAD","PRAD"))
logcpm.sum.df$peak.precursor <- factor(logcpm.sum.df$peak.precursor, levels = c("LAML","COAD","PAAD","PRAD"))
p1 <- ggplot(logcpm.sum.df, aes(x=group,y=mean.log2cpm,fill=group))+ #sum.log2cpm.scale
  geom_violin(alpha=0.7)+
  geom_boxplot(width=0.1)+
  scale_fill_manual(name="Group",values = c("grey50","salmon","steelblue","purple"))+#colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(unique(logcpm.sum$group))))+)
  facet_grid(.~peak.precursor,scales = "free")+
  ylab("Scaled peak-index") +
  ggpubr::stat_compare_means(
    label.x.npc = "middle", label.y.npc = "top",
    #size =12,step.increase = 0.08,#paired = TRUE,
    aes(x=group,y=mean.log2cpm,
        label = ..p.format..), # p.signif, p.format, group = Group,
    symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01,1),
                     symbols = c( "***", "**", "*", "ns")),
    # method.args = list(alternative="greater"),
    # ref.group = peak.precursor,
    hide.ns=F, size=5,
    method = "wilcox.test"
  ) +
  theme_bw() +
  my_theme_box
p1
ggsave(paste0(dst,"_multiclass_box.pdf"),width = 24,height = 8)
# p2 <- ggplot(logcpm.sum.df, aes(x=peak.precursor,y=sum.log2cpm,fill=peak.precursor))+ #sum.log2cpm.scale
#   geom_violin(alpha=0.7)+
#   geom_boxplot(width=0.1)+
#   scale_fill_manual(name="Group",values = c("grey50","salmon","steelblue","purple"))+#colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(unique(logcpm.sum$group))))+)
#   facet_grid(.~group,scales = "free")+
#   ylab("Scaled peak-index") +
#   ggpubr::stat_compare_means(
#     label.x.npc = "middle", label.y.npc = "top",
#     #size =12,step.increase = 0.08,#paired = TRUE,
#     aes(x=group,y=sum.log2cpm,
#         label = ..p.format..), # p.signif, p.format, group = Group,
#     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01,1),
#                      symbols = c( "***", "**", "*", "ns")),
#     # method.args = list(alternative="greater"),
#     # ref.group = peak.precursor,
#     hide.ns=F, size=5,
#     method = "wilcox.test"
#   ) +
#   theme_bw() +
#   my_theme_box
# p2
# ggsave(paste0(dst,"_multiclass_byGroup_box.pdf"),width = 24,height = 8)




logcpm.sum.tbl <- dplyr::as_tibble(logcpm.sum.df) %>%
  dplyr::group_by(sample) %>%
  dplyr::top_n(n=1,wt=mean.log2cpm) # sum.log2cpm
logcpm.sum.tbl <- as.data.frame(logcpm.sum.tbl)
# logcpm.sum.tbl <- as.data.frame(table(logcpm.sum.tbl$group,logcpm.sum.tbl$peak.precursor) )
# logcpm.sum.tbl <- reshape2::dcast(data = logcpm.sum.tbl, formula = Var1~Var2)
# rownames(logcpm.sum.tbl) <- logcpm.sum.tbl$Var1
# logcpm.sum.tbl$Var1 <- NULL
# rows: true
# cols: pred
# 
# #
p <- plotConfusionMat(Actual = (logcpm.sum.tbl$group), Predict = (logcpm.sum.tbl$peak.precursor))
p
table(logcpm.sum.tbl$group,logcpm.sum.tbl$peak.precursor)
ggsave(plot = p, filename = paste0(dst,"_multiclass_confusion_mat.pdf"),width = 8,height = 8)

# fail when TDCPs defined from TCGA
# TODO: define TDCP from plasma









# 3tumor-peak-index GSE71008 --> TCGA_small (GSE71008 peak) (cfpeak suppl Fig? fail)  ------------------------------------

## GSE71008 CRCvsNC (poor, even rm batch, only use this dst for multi-classifier )
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
dst <- "GSE71008_diff"
# NC.label <- "NC"
# CRC.label <- "CRC"
### read smp table
sample.table <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/GSE71008_diff/sample_table.txt",check.names = F,header = T)
table(sample.table$disease_status)
# sample.table <- sample.table[sample.table$disease_status %in% c("Colorectal","Healthy"),]
sample.table$sample <- rownames(sample.table)
#sample.table <- sample.table[,c("source","sample")]
sample.table <- sample.table[,c("disease_status","sample")]
colnames(sample.table) <- c("group","sample")
sample.table$group <- gsub("Colorectal","COAD",sample.table$group)
sample.table$group <- gsub("Healthy","LAML",sample.table$group)
sample.table$group <- gsub("Pancreatic","PAAD",sample.table$group)
sample.table$group <- gsub("Prostate","PRAD",sample.table$group)
table(sample.table$group)
# Colorectal    Healthy Pancreatic   Prostate 
# 100         50          6         36
# COAD LAML PAAD PRAD 
# 100   50    6   36 

### read cpm (all fig6 use rowsum cpm)
# logcpm <- read.table(paste0(pre,"/output/",dst,"/call_peak_all/count_matrix/expeakCNN_b5_d50_p1_CPMrowsum.txt"),check.names = F,header = T)
# rownames(logcpm) <- logcpm$gene_id
# logcpm$gene_id <- NULL
# max(as.matrix(logcpm))
# logcpm <- log2(logcpm+1)
logcpm <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",dst,"/allSmp_rmBatch_OvR_exPeakCNN_smallDomain_diff_CRCvsR.cpm"),check.names = F,header = T)
cpm <- 2^logcpm-1 ## zscore are calculated in cpm, not log2cpm !!!


# build ref
tumors <- c("COAD","LAML","PAAD","PRAD")
#feature.lab.sort <- tumors
feature.list <- list()
mean.list <- list()
sd.list <- list()
top <- 1000
for(disease in tumors){
  print(disease)
  ## build ref db
  positive_samples <- sample.table[sample.table$group==disease,"sample"]
  negative_samples <- sample.table[sample.table$group!=disease,"sample"]
  if(disease=="LAML"){
    disease.old <- "NC"
  }else if(disease=="COAD"){
    disease.old <- "CRC"
  }else if(disease=="PRAD"){
    disease.old <- "PRCA"
  }else if(disease=="PAAD"){
    disease.old <- "PACA"
  }
  
  feature.list[[disease]] <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_all/count_matrix/cfpeakCNN_b5_d50_p1.txt.diff.",disease.old,"vsR.top",top))$V1 # diff, dcb
  
  #extract fixed coad&blood group
  blood.mat <- cpm[feature.list[[disease]], negative_samples] # cpm, not log2cpm !!!
  coad.mat <- cpm[feature.list[[disease]], positive_samples]
  mean.list[[disease]] <- apply(blood.mat,1,mean) # zscore are calculated in cpm, not log2cpm
  sd.list[[disease]] <- apply(cbind(coad.mat,blood.mat),1,sd) # blood mat has many 0 : Zhu - 2021 - NC - Tissue-specific cell-free DNA degradation...
  # sd.list <- (rep(sd(as.matrix(blood.mat)),length(mean.vec)))
  #print(table(sd.list==0))
  # print( table(sd.list[[i]]==0) ) # need all FALSE !!! or else enlarge cohort when calculating sd (eg: add other groups)
}
length(mean.list[["COAD"]])
length(sd.list[["LAML"]])
#not need keep same num of TDCPs in diff cancer
#



#logcpm.sum <- getPeakIndex(logcpm = logcpm, sample.table = sample.table, prescale = T, zscore.cutoff = 100) # lower zscore.cutoff lead to lower auc in this dataset !!!  
logcpm.sum.list <- list()
tumors <- unique(sample.table$group)
for(disease in tumors){
  # positive_samples <- sample.table.tcga[sample.table.tcga$group==disease,"sample"]
  # negative_samples <- sample.table.tcga[sample.table.tcga$group!=disease,"sample"]
  # samples <- c(positive_samples, negative_samples)
  logcpm.sum.list[[disease]] <- getPeakIndex(feature.lab = disease, logcpm = logcpm, sample.table = sample.table, prescale = T, mean.list = mean.list, sd.list = sd.list, zscore.cutoff = 100, onlyCountZscoreOutlier = 3)
}
#dim(logcpm.sum.list[["COAD"]])
#names(mean.list)


logcpm.sum.df <- as.data.frame(do.call(rbind,logcpm.sum.list))
table(logcpm.sum.df$group==logcpm.sum.df$peak.precursor)

NC.label <- "LAML"
logcpm.sum.df$group <- factor(logcpm.sum.df$group, levels = c("LAML","COAD","PAAD","PRAD"))
logcpm.sum.df$peak.precursor <- factor(logcpm.sum.df$peak.precursor, levels = c("LAML","COAD","PAAD","PRAD"))
p1 <- ggplot(logcpm.sum.df, aes(x=group,y=mean.log2cpm,fill=group))+ #sum.log2cpm.scale
  geom_violin(alpha=0.7)+
  geom_boxplot(width=0.1)+
  scale_fill_manual(name="Group",values = c("grey50","salmon","steelblue","purple"))+#colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(unique(logcpm.sum$group))))+)
  facet_grid(.~peak.precursor,scales = "free")+
  ylab("Scaled peak-index") +
  ggpubr::stat_compare_means(
    label.x.npc = "middle", label.y.npc = "top",
    #size =12,step.increase = 0.08,#paired = TRUE,
    aes(x=group,y=mean.log2cpm,
        label = ..p.format..), # p.signif, p.format, group = Group,
    symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01,1),
                     symbols = c( "***", "**", "*", "ns")),
    # method.args = list(alternative="greater"),
    # ref.group = peak.precursor,
    hide.ns=F, size=5,
    method = "wilcox.test"
  ) +
  theme_bw() +
  my_theme_box
p1
ggsave(paste0(dst,"_multiclass_box.pdf"),width = 24,height = 8)
# p2 <- ggplot(logcpm.sum.df, aes(x=peak.precursor,y=sum.log2cpm,fill=peak.precursor))+ #sum.log2cpm.scale
#   geom_violin(alpha=0.7)+
#   geom_boxplot(width=0.1)+
#   scale_fill_manual(name="Group",values = c("grey50","salmon","steelblue","purple"))+#colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(unique(logcpm.sum$group))))+)
#   facet_grid(.~group,scales = "free")+
#   ylab("Scaled peak-index") +
#   ggpubr::stat_compare_means(
#     label.x.npc = "middle", label.y.npc = "top",
#     #size =12,step.increase = 0.08,#paired = TRUE,
#     aes(x=group,y=sum.log2cpm,
#         label = ..p.format..), # p.signif, p.format, group = Group,
#     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01,1),
#                      symbols = c( "***", "**", "*", "ns")),
#     # method.args = list(alternative="greater"),
#     # ref.group = peak.precursor,
#     hide.ns=F, size=5,
#     method = "wilcox.test"
#   ) +
#   theme_bw() +
#   my_theme_box
# p2
# ggsave(paste0(dst,"_multiclass_byGroup_box.pdf"),width = 24,height = 8)

logcpm.sum.tbl <- dplyr::as_tibble(logcpm.sum.df) %>%
  dplyr::group_by(sample) %>%
  dplyr::top_n(n=1,wt=mean.log2cpm) # sum.log2cpm
logcpm.sum.tbl <- as.data.frame(logcpm.sum.tbl)

p <- plotConfusionMat(Actual = (logcpm.sum.tbl$group), Predict = (logcpm.sum.tbl$peak.precursor))
p
table(logcpm.sum.tbl$group,logcpm.sum.tbl$peak.precursor)
ggsave(plot = p, filename = paste0(dst,"_multiclass_confusion_mat.pdf"),width = 8,height = 8)

# success when TDCPs defined from plasma






# TCGA_small_diff CRCvsNC AUROC
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
dst <- "TCGA_small_diff3"

### read smp table
sample.table.tcga <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/",dst,"/sample_table.txt"),check.names = F, sep = "\t", header = T)
#sample.table.tcga <- sample.table.tcga[,1:2]
sample.table.tcga <- sample.table.tcga[,c(2,5)]
colnames(sample.table.tcga) <- c("sample","group")
sample.table.tcga$sample <- gsub(".bam","",sample.table.tcga$sample)
rownames(sample.table.tcga) <- sample.table.tcga$sample
sample.table.tcga$group <- gsub("TCGA-","",sample.table.tcga$group)
table(sample.table.tcga$group)
sample.table.tcga <- sample.table.tcga[order(sample.table.tcga$group),]

#logcpm <- read.table("count_matrix",check.names = F,header = T) # cpm
#logcpm <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",dst,"/cfPeakCNN_GSE110381_smallDomain_diff_COADvsLAML.cpm"),check.names = F,header = T)  # log2cpm+1
#logcpm <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_all/count_matrix/cfpeakCNN_b5_d50_p1_GSE110381_CPMrowsum.txt"),check.names = F,header = T)  # log2cpm+1
#logcpm <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_all/count_matrix/cfpeakCNN_b5_d50_p1_GSE71008_CPMrowsum.txt"),check.names = F,header = T)  # log2cpm+1
logcpm <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",dst,"/OvR_cfPeakCNN_GSE71008_smallDomain_diff_COADvsR.cpm"),row.names = 1,check.names = F,header = T)
# rownames(logcpm) <- logcpm$gene_id # feature
cpm <- 2^logcpm-1 ## zscore are calculated in cpm, not log2cpm !!!
#max(logcpm[,1])
# logcpm <- logcpm[feature,]
sample.table.tcga <- sample.table.tcga[sample.table.tcga$sample %in% colnames(logcpm),]
table(sample.table.tcga$group)
# Blood  COAD  PAAD  PRAD 
# 60    60    60    60


# # build ref
# tumors <- c("COAD","LAML","PAAD","PRAD")
# #feature.lab.sort <- tumors
# feature.list <- list()
# mean.list <- list()
# sd.list <- list()
# for(disease in tumors){
#   print(disease)
#   ## build ref db
#   positive_samples <- sample.table.tcga[sample.table.tcga$group==disease,"sample"]
#   negative_samples <- sample.table.tcga[sample.table.tcga$group!=disease,"sample"]
#   
#   feature.list[[disease]] <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_all/count_matrix/cfpeakCNN_b5_d50_p1_GSE71008.txt.diff.dcb.",disease,"vsR"))$V1
#   
#   #extract fixed coad&blood group
#   blood.mat <- cpm[feature.list[[disease]], negative_samples] # cpm, not log2cpm !!!
#   coad.mat <- cpm[feature.list[[disease]], positive_samples]
#   mean.list[[disease]] <- apply(blood.mat,1,mean) # zscore are calculated in cpm, not log2cpm
#   sd.list[[disease]] <- apply(cbind(coad.mat,blood.mat),1,sd) # blood mat has many 0 : Zhu - 2021 - NC - Tissue-specific cell-free DNA degradation...
#   # sd.list <- (rep(sd(as.matrix(blood.mat)),length(mean.vec)))
#   #print(table(sd.list==0))
#   # print( table(sd.list[[i]]==0) ) # need all FALSE !!! or else enlarge cohort when calculating sd (eg: add other groups)
# }
# length(mean.list[["COAD"]])
# length(sd.list[["LAML"]])
# #not need keep same num of TDCPs in diff cancer
# #


logcpm.sum.list <- list()
for(disease in tumors){
  # positive_samples <- sample.table.tcga[sample.table.tcga$group==disease,"sample"]
  # negative_samples <- sample.table.tcga[sample.table.tcga$group!=disease,"sample"]
  # samples <- c(positive_samples, negative_samples)
  logcpm.sum.list[[disease]] <- getPeakIndex(feature.lab = disease, logcpm = logcpm, sample.table = sample.table.tcga, prescale = T, mean.list = mean.list, sd.list = sd.list, zscore.cutoff = 100, onlyCountZscoreOutlier = 3)
}
#dim(logcpm.sum.list[["COAD"]])
#names(mean.list)



##plot peak index score 
logcpm.sum.df <- as.data.frame(do.call(rbind,logcpm.sum.list))
table(logcpm.sum.df$group==logcpm.sum.df$peak.precursor)

NC.label <- "LAML"
logcpm.sum.df$group <- factor(logcpm.sum.df$group, levels = c("LAML","COAD","PAAD","PRAD"))
logcpm.sum.df$peak.precursor <- factor(logcpm.sum.df$peak.precursor, levels = c("LAML","COAD","PAAD","PRAD"))
p1 <- ggplot(logcpm.sum.df, aes(x=group,y=mean.log2cpm,fill=group))+ #sum.log2cpm.scale
  geom_violin(alpha=0.7)+
  geom_boxplot(width=0.1)+
  scale_fill_manual(name="Group",values = c("grey50","salmon","steelblue","purple"))+#colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(unique(logcpm.sum$group))))+)
  facet_grid(.~peak.precursor,scales = "free")+
  ylab("Scaled peak-index") +
  ggpubr::stat_compare_means(
    label.x.npc = "middle", label.y.npc = "top",
    #size =12,step.increase = 0.08,#paired = TRUE,
    aes(x=group,y=mean.log2cpm,
        label = ..p.format..), # p.signif, p.format, group = Group,
    symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01,1),
                     symbols = c( "***", "**", "*", "ns")),
    # method.args = list(alternative="greater"),
    # ref.group = peak.precursor,
    hide.ns=F, size=5,
    method = "wilcox.test"
  ) +
  theme_bw() +
  my_theme_box
p1
ggsave(paste0(dst,"_multiclass_box.pdf"),width = 24,height = 8)
# p2 <- ggplot(logcpm.sum.df, aes(x=peak.precursor,y=mean.log2cpm,fill=peak.precursor))+ #sum.log2cpm.scale
#   geom_violin(alpha=0.7)+
#   geom_boxplot(width=0.1)+
#   scale_fill_manual(name="Group",values = c("grey50","salmon","steelblue","purple"))+#colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(unique(logcpm.sum$group))))+)
#   facet_grid(.~group,scales = "free")+
#   ylab("Scaled peak-index") +
#   ggpubr::stat_compare_means(
#     label.x.npc = "middle", label.y.npc = "top",
#     #size =12,step.increase = 0.08,#paired = TRUE,
#     aes(x=group,y=mean.log2cpm,
#         label = ..p.format..), # p.signif, p.format, group = Group,
#     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01,1),
#                      symbols = c( "***", "**", "*", "ns")),
#     # method.args = list(alternative="greater"),
#     # ref.group = peak.precursor,
#     hide.ns=F, size=5,
#     method = "wilcox.test"
#   ) +
#   theme_bw() +
#   my_theme_box
# p2
# ggsave(paste0(dst,"_multiclass_byGroup_box.pdf"),width = 24,height = 8)



##plot confusion mat
logcpm.sum.tbl <- dplyr::as_tibble(logcpm.sum.df) %>%
  dplyr::group_by(sample) %>%
  dplyr::top_n(n=1,wt=mean.log2cpm) # sum.log2cpm, mean.log2cpm
logcpm.sum.tbl <- as.data.frame(logcpm.sum.tbl)

p <- plotConfusionMat(Actual = (logcpm.sum.tbl$group), Predict = (logcpm.sum.tbl$peak.precursor))
p
table(logcpm.sum.tbl$group,logcpm.sum.tbl$peak.precursor)
ggsave(plot = p, filename = paste0(dst,"_multiclass_confusion_mat.pdf"),width = 8,height = 8)

#also fail...
#







# overlapped diff peak number with the same trend between plasma and tissue, reanalysis 2018cfmedip fig2d/e (cfpeak fig6) ------------------------
#only consider CRCvsBlood/LAML, for TCGA has little matched CRC normal, bias exist when use GTEx colon normal, and 2018cfmedip has proved there are very little difference when using blood as ctrl
#also we consider both shuffle DMR plasma and in tissue
setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/")
#setwd("/Users/baopengfei/Desktop/lulab/tmp/projects/WCHSU-FTC/exSeek-dev")
options(stringsAsFactors = F)

pre <- "/BioII/lulab_b/baopengfei/projects"
diff <- list()
plasma.dst="PRJNA540919_diff" # GSE110381_diff, PRJNA540919_diff, lulab_oscc_plasma_diff
diff[[plasma.dst]] <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",plasma.dst,"/cfPeakCNN_smallDomain_diff_CRCvsNC.diff"))
# diff[[plasma.dst]] <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",plasma.dst,"/cfPeakCNN_smallDomain_diff_MvsL.diff"))

tissue.dst="TCGA_small_diff3" # TCGA_small_diff3,TCGA-HNSC_small_diff
diff[[tissue.dst]] <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",tissue.dst,"/cfPeakCNN_GSE110381_smallDomain_diff_COADvsLAML.diff"))
# diff[[tissue.dst]] <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",tissue.dst,"/cfPeakCNN_lulab_oscc_plasma_diff_smallDomain_MvsL.diff"))


diff[[plasma.dst]]$peak.id <- unlist(sapply(strsplit(rownames(diff[[plasma.dst]]),"|",fixed = T),"[",4))
diff[[tissue.dst]]$peak.id <- unlist(sapply(strsplit(rownames(diff[[tissue.dst]]),"|",fixed = T),"[",4))

# #only keep peaks overlapped with original TCGA consensus (rm logFC==0 might be enough)
# options (bedtools.path = "/BioII/lulab_b/baopengfei/biosoft")
# keep <- bedtoolsr::bt.intersect(a = paste0(pre,"/WCHSU-FTC/output/",tissue.dst,"/call_peak_all/cfpeakCNN/b5_d50_p1_GSE110381.bed"),
#                                 b = paste0(pre,"/WCHSU-FTC/output/",tissue.dst,"/call_peak_all/cfpeakCNN/b5_d50_p1.bed"),
#                                 u=T, s = T)
# #only keep peaks overlapped with original PRJNA540919 consensus
# keep <- bedtoolsr::bt.intersect(a = keep, # paste0(pre,"/WCHSU-FTC/output/",plasma.dst,"/call_peak_all/cfpeakCNN/b5_d50_p1_GSE110381.bed"), #
#                                 b = paste0(pre,"/WCHSU-FTC/output/","PRJNA540919_diff","/call_peak_all/cfpeakCNN/b5_d50_p1.bed"),
#                                 u=T, s = T)

# diff[[plasma.dst]] <- diff[[plasma.dst]][diff[[plasma.dst]]$peak.id %in% keep$V4,]
# diff[[tissue.dst]] <- diff[[tissue.dst]][diff[[tissue.dst]]$peak.id %in% keep$V4,]
# dim(diff[[plasma.dst]] )


# #rm highblood lowtissue peaks (optional 1)
# cutoff.crc <- 0.1
# cutoff.blood <- 0.1
# table(diff.tcga$posMeanCPM<=cutoff.crc)
# table(diff.tcga$negMeanCPM>=cutoff.blood)
# table(diff.tcga$posMeanCPM<=cutoff.crc & diff.tcga$negMeanCPM>=cutoff.blood)
# lowTissue <- diff.tcga$feature[diff.tcga$posMeanCPM<=cutoff.crc & diff.tcga$negMeanCPM>=cutoff.blood]
# length(lowTissue) # 12642
# length(diff.tcga$feature)
# diff.merge.keep <- diff.merge.keep[!(diff.merge.keep$feature %in% lowTissue),]
# # 0.13055

#select dcb peaks (optional 2)
dcb <- read.table(paste0(pre,"/WCHSU-FTC/output/",tissue.dst,"/call_peak_all/count_matrix/cfpeakCNN_b5_d50_p1_GSE110381.txt.diff.dcb"))
table(rownames(diff[[plasma.dst]]) %in% dcb$V1)
# FALSE  TRUE
# 82512  1745
# 
# keep <- (rownames(diff[[plasma.dst]]) %in% dcb$V1) # (diff[[plasma.dst]]$log2FoldChange != 0) & (diff[[tissue.dst]]$log2FoldChange != 0) #  & (rownames(diff[[plasma.dst]]) %in% dcb$V1) # dcb: all tissue log2FC>0
# table(keep)
keep <- (diff[[plasma.dst]]$log2FoldChange != 0) 
#table(keep)
diff[[plasma.dst]] <- diff[[plasma.dst]][keep,]
diff[[tissue.dst]] <- diff[[tissue.dst]][keep,]
dim(diff[[tissue.dst]])

#top.cutoff <- 50000
getTopPeakId <- function(p,top.cutoff){
  topPeak <- list()
  # p <- diff[[plasma.dst]]
  p <- p[order(p$pvalue,decreasing = F),]
  p.up <- p[p$log2FoldChange>0,]
  p.down <- p[p$log2FoldChange<0,]
  topPeak[["up"]] <- rownames(p.up)[1:min(top.cutoff,nrow(p.up))]
  topPeak[["dw"]] <- rownames(p.down)[1:min(top.cutoff,nrow(p.down))]
  return(topPeak)
}


nums <- c(100,200,300,500,1000,10000) # ,50000
for(top.cutoff in nums){
  # top.cutoff <- 50000
  print(top.cutoff)
  
plasma.top <- getTopPeakId(diff[[plasma.dst]], top.cutoff)
tissue.top <- getTopPeakId(diff[[tissue.dst]], top.cutoff)

#length(plasma.top[["up"]])
obser.upup.num <- length(intersect(plasma.top[["up"]],tissue.top[["up"]]))
obser.dwdw.num <- length(intersect(plasma.top[["dw"]],tissue.top[["dw"]]))
obser.updw.num <- length(intersect(plasma.top[["up"]],tissue.top[["dw"]]))
obser.dwup.num <- length(intersect(plasma.top[["dw"]],tissue.top[["up"]]))


## shuf 500 times
shuf.times <- 500
shuf.num <- list()
shuf.num[['upup']] <- list()
shuf.num[['dwdw']] <- list()
shuf.num[['updw']] <- list()
shuf.num[['dwup']] <- list()
set.seed(1234)
for(i in 1:shuf.times){
  # i <- 1
  set.seed(i)
  tmp.diff.plasma <- diff[[plasma.dst]]
  dim(tmp.diff.plasma)
  tmp.up.plasma <- sum(tmp.diff.plasma$log2FoldChange>0)
  tmp.dw.plasma <- sum(tmp.diff.plasma$log2FoldChange<0)
  shuf.tmp.plasma <- rownames(tmp.diff.plasma)[sample(1:nrow(tmp.diff.plasma),(tmp.up.plasma + tmp.dw.plasma),replace = F)]
  length(shuf.tmp.plasma)
  
  tmp.diff.tissue <- diff[[tissue.dst]]
  dim(tmp.diff.tissue)
  tmp.up.tissue <- sum(tmp.diff.tissue$log2FoldChange>0)
  tmp.dw.tissue <- sum(tmp.diff.tissue$log2FoldChange<0)
  shuf.tmp.tissue <- rownames(tmp.diff.tissue)[sample(1:nrow(tmp.diff.tissue),(tmp.up.tissue + tmp.dw.tissue),replace = F)]
  
  shuf.num[['upup']][[i]] <- length(intersect( shuf.tmp.plasma[1:length(plasma.top[["up"]])],                                 tissue.top[["up"]] )) # shuf.tmp.tissue[1:length(tissue.top[["up"]])]
  shuf.num[['dwdw']][[i]] <- length(intersect( shuf.tmp.plasma[(tmp.up.plasma+1):(tmp.up.plasma+length(plasma.top[["dw"]]))], tissue.top[["dw"]] )) # shuf.tmp.tissue[(tmp.up.tissue+1):(tmp.up.tissue+length(tissue.top[["dw"]]))]
  shuf.num[['updw']][[i]] <- length(intersect( shuf.tmp.plasma[1:length(plasma.top[["up"]])],                                 tissue.top[["dw"]] )) # shuf.tmp.tissue[(tmp.up.tissue+1):(tmp.up.tissue+length(tissue.top[["dw"]]))]
  shuf.num[['dwup']][[i]] <- length(intersect( shuf.tmp.plasma[(tmp.up.plasma+1):(tmp.up.plasma+length(plasma.top[["dw"]]))], tissue.top[["up"]] )) # shuf.tmp.tissue[1:length(tissue.top[["up"]])]
  # shuf.num[['upup']][[i]] <- length(intersect(shuf.tmp.plasma[1:tmp.up.plasma],shuf.tmp.tissue[1:tmp.up.tissue]))
  # shuf.num[['dwdw']][[i]] <- length(intersect(shuf.tmp.plasma[(tmp.up.plasma+1):(tmp.up.plasma+tmp.dw.plasma)],shuf.tmp.tissue[(tmp.up.tissue+1):(tmp.up.tissue+tmp.dw.tissue)]))
  # shuf.num[['updw']][[i]] <- length(intersect(shuf.tmp.plasma[1:tmp.up.plasma],shuf.tmp.tissue[(tmp.up.tissue+1):(tmp.up.tissue+tmp.dw.tissue)]))
  # shuf.num[['dwup']][[i]] <- length(intersect(shuf.tmp.plasma[(tmp.up.plasma+1):(tmp.up.plasma+tmp.dw.plasma)],shuf.tmp.tissue[1:tmp.up.tissue]))
}
shuf.upup <- do.call("c",shuf.num[['upup']])
shuf.dwdw <- do.call("c",shuf.num[['dwdw']])
shuf.updw <- do.call("c",shuf.num[['updw']])
shuf.dwup <- do.call("c",shuf.num[['dwup']])


## get z-score of overlapped number
# shuf.up.scale <- (shuf.up-mean(shuf.up))/sd(shuf.up)
mat <- cbind(shuf.upup,shuf.updw,shuf.dwdw,shuf.dwup)
mat.scale <- scale(mat)

obser.upup.num.scale <- (obser.upup.num-mean(mat[,1]))/sd(mat[,1])
obser.updw.num.scale <- (obser.updw.num-mean(mat[,2]))/sd(mat[,2])
obser.dwdw.num.scale <- (obser.dwdw.num-mean(mat[,3]))/sd(mat[,3])
obser.dwup.num.scale <- (obser.dwup.num-mean(mat[,4]))/sd(mat[,4])
obser.num <- c(obser.upup.num,obser.updw.num,obser.dwdw.num,obser.dwup.num)
#boxplot(mat)

## get pvalue
pval.vec=vector(mode="numeric", length=ncol(mat.scale))
obser.num.scale=vector(mode="numeric", length=ncol(mat.scale))
for(j in 1:ncol(mat.scale)) {
  # i <- 2
  obser.num.scale[j] <- (obser.num[j]-mean(mat[,j]))/sd(mat[,j])
  
  this.mean <- mean(mat.scale[,j])
  this.sd <- sd(mat.scale[,j])
  LOGpvalue = signif(as.double(pnorm(obser.num.scale[j], mean=this.mean, sd=this.sd, lower.tail=TRUE, log=TRUE)), digits=3)
  if ( obser.num.scale[j] > this.mean) 
  {
    LOGpvalue = signif(as.double(pnorm(obser.num.scale[j], mean=this.mean, sd=this.sd, lower.tail=FALSE, log=TRUE)), digits=3)
  }
  # message(LOGpvalue < -745)
  pvalue = signif(exp(LOGpvalue), digits=3) 
  # if (LOGpvalue < -745) pvalue=paste0("e^", LOGpvalue)
  pval.vec[j] <- pvalue
}
# yMin=10*floor(min(BiolVals_Vector)/10)
# if (yMin>=0) yMin=-10
# yMax=10*ceiling(max(BiolVals_Vector)/10)
# if (yMax<=0) yMax=10
yMin <- min(mat.scale,obser.num.scale)
yMax <- max(mat.scale,obser.num.scale)

pdf(paste0("BoxPlots_Overlaps_CancerClass_vs_Controls_TCGA_vsBlood_",top.cutoff,"_2.pdf"),width = 4.5,height = 6)
# pdf("BoxPlots_Overlaps_M_vs_L_TCGA_vsBlood.pdf",width = 4.5,height = 6) 
{
boxplot(mat.scale,
        ylim=c(yMin-4, yMax+4), ylab="Z-Scores", #las=2,
        names=c("A", "B", "C", "D"),
        cex.axis=1.5#, main=paste0("With RRBS ", Biol_Comparands[k]), cex.main=0.9
)
col.vec <- c("firebrick","grey50","firebrick","grey50")
print(obser.num.scale)
points(seq(1,4,by=1), obser.num.scale, pch=23, col=col.vec, 
       bg=col.vec, cex=1.5)
text(x=seq(1,4,by=1), y=rep(c(yMax+2, yMin-2), 2), paste0("p=",pval.vec), cex=1)
}
dev.off()
#abline(h=NegCutoffZscore, col="darkgray", lty=5)
#abline(h=PosCutoffZscore, col="darkgray", lty=5)
}

#




# DCB
#nums <- c(100,1000,10000,50000)
for(top.cutoff in nums){
  # top.cutoff <- 2000
  print(top.cutoff)
  
  # plasma.top <- getTopPeakId(diff[[plasma.dst]], top.cutoff)
  # tissue.top <- getTopPeakId(diff[[tissue.dst]], top.cutoff)
  
  plasma.top <- list()
  plasma.top[["up"]] <- rownames(diff[[plasma.dst]])[(rownames(diff[[plasma.dst]]) %in% dcb$V1) & diff[[plasma.dst]]$log2FoldChange>0]
  plasma.top[["dw"]] <- rownames(diff[[plasma.dst]])[(rownames(diff[[plasma.dst]]) %in% dcb$V1) & diff[[plasma.dst]]$log2FoldChange<0]
  # length(plasma.top[["dw"]])
  
  tissue.top <- list()
  tissue.top[["up"]] <- rownames(diff[[tissue.dst]])[(rownames(diff[[tissue.dst]]) %in% dcb$V1) & diff[[tissue.dst]]$log2FoldChange>0]
  # tissue.top[["dw"]] <- rownames(diff[[tissue.dst]])[(rownames(diff[[tissue.dst]]) %in% dcb$V1) & diff[[tissue.dst]]$log2FoldChange<0]
  # length(tissue.top[["dw"]])
  
  # table(diff[[plasma.dst]]$log2FoldChange>0 & diff[[tissue.dst]]$log2FoldChange>0) # 829/9257=0.08955
  # table(diff[[plasma.dst]]$log2FoldChange<0 & diff[[tissue.dst]]$log2FoldChange>0) # 566/5590=0.10125
  
  
  obser.upup.num <- length(intersect(plasma.top[["up"]],tissue.top[["up"]]))
  obser.dwup.num <- length(intersect(plasma.top[["dw"]],tissue.top[["up"]]))
  
  
  ## shuf 500 times
  shuf.times <- 500
  shuf.num <- list()
  shuf.num[['upup']] <- list()
  shuf.num[['dwup']] <- list()
  set.seed(1234)
  for(i in 1:shuf.times){
    # i <- 1
    set.seed(i)
    tmp.diff.plasma <- diff[[plasma.dst]]
    dim(tmp.diff.plasma)
    tmp.up.plasma <- sum(tmp.diff.plasma$log2FoldChange>0)
    tmp.dw.plasma <- sum(tmp.diff.plasma$log2FoldChange<0)
    shuf.tmp.plasma <- rownames(tmp.diff.plasma)[sample(1:nrow(tmp.diff.plasma),(tmp.up.plasma + tmp.dw.plasma),replace = F)]
    
    tmp.diff.tissue <- diff[[tissue.dst]]
    dim(tmp.diff.tissue)
    tmp.up.tissue <- sum(tmp.diff.tissue$log2FoldChange>0)
    tmp.dw.tissue <- sum(tmp.diff.tissue$log2FoldChange<0)
    shuf.tmp.tissue <- rownames(tmp.diff.tissue) # [sample(1:nrow(tmp.diff.tissue),(tmp.up.tissue + tmp.dw.tissue),replace = F)]
    
    shuf.num[['upup']][[i]] <- length(intersect( shuf.tmp.plasma[1:length(plasma.top[["up"]])],                                 tissue.top[["up"]] )) # shuf.tmp.tissue[1:length(tissue.top[["up"]])]
    # shuf.num[['dwdw']][[i]] <- length(intersect( shuf.tmp.plasma[(tmp.up.plasma+1):(tmp.up.plasma+length(plasma.top[["dw"]]))], shuf.tmp.tissue[(tmp.up.tissue+1):(tmp.up.tissue+length(tissue.top[["dw"]]))] ))
    # shuf.num[['updw']][[i]] <- length(intersect( shuf.tmp.plasma[1:length(plasma.top[["up"]])],                                 shuf.tmp.tissue[(tmp.up.tissue+1):(tmp.up.tissue+length(tissue.top[["dw"]]))] ))
    shuf.num[['dwup']][[i]] <- length(intersect( shuf.tmp.plasma[(tmp.up.plasma+1):(tmp.up.plasma+length(plasma.top[["dw"]]))], tissue.top[["up"]] )) # shuf.tmp.tissue[1:length(tissue.top[["up"]])]
    
  }
  shuf.upup <- do.call("c",shuf.num[['upup']])
  # shuf.dwdw <- do.call("c",shuf.num[['dwdw']])
  # shuf.updw <- do.call("c",shuf.num[['updw']])
  shuf.dwup <- do.call("c",shuf.num[['dwup']])
  
  
  ## get z-score of overlapped number
  # shuf.up.scale <- (shuf.up-mean(shuf.up))/sd(shuf.up)
  mat <- cbind(shuf.upup,shuf.dwup)
  mat.scale <- scale(mat)
  
  obser.upup.num.scale <- (obser.upup.num-mean(mat[,1]))/sd(mat[,1])
  # obser.updw.num.scale <- (obser.updw.num-mean(mat[,2]))/sd(mat[,2])
  # obser.dwdw.num.scale <- (obser.dwdw.num-mean(mat[,3]))/sd(mat[,3])
  obser.dwup.num.scale <- (obser.dwup.num-mean(mat[,2]))/sd(mat[,2])
  obser.num <- c(obser.upup.num,obser.dwup.num)
  #boxplot(mat)
  
  ## get pvalue
  pval.vec=vector(mode="numeric", length=ncol(mat.scale))
  obser.num.scale=vector(mode="numeric", length=ncol(mat.scale))
  for(j in 1:ncol(mat.scale)) {
    # j <- 1
    obser.num.scale[j] <- (obser.num[j]-mean(mat[,j]))/sd(mat[,j])
    
    this.mean <- mean(mat.scale[,j])
    this.sd <- sd(mat.scale[,j])
    LOGpvalue = signif(as.double(pnorm(obser.num.scale[j], mean=this.mean, sd=this.sd, lower.tail=TRUE, log=TRUE)), digits=3)
    if ( obser.num.scale[j] > this.mean) 
    {
      LOGpvalue = signif(as.double(pnorm(obser.num.scale[j], mean=this.mean, sd=this.sd, lower.tail=FALSE, log=TRUE)), digits=3)
    }
    # message(LOGpvalue < -745)
    pvalue = signif(exp(LOGpvalue), digits=3) 
    # if (LOGpvalue < -745) pvalue=paste0("e^", LOGpvalue)
    pval.vec[j] <- pvalue
  }
  # yMin=10*floor(min(BiolVals_Vector)/10)
  # if (yMin>=0) yMin=-10
  # yMax=10*ceiling(max(BiolVals_Vector)/10)
  # if (yMax<=0) yMax=10
  yMin <- min(mat.scale,obser.num.scale)
  yMax <- max(mat.scale,obser.num.scale)
  
  pdf(paste0("DCB_BoxPlots_Overlaps_CancerClass_vs_Controls_TCGA_vsBlood_",top.cutoff,"_2.pdf"),width = 2.5,height = 6)
  # pdf("DCB_BoxPlots_Overlaps_M_vs_L_TCGA_vsBlood.pdf",width = 4.5,height = 6) 
  {
    boxplot(mat.scale,
            ylim=c(yMin-4, yMax+4), ylab="Z-Scores", #las=2,
            # names=c("A", "B", "C", "D"),
            cex.axis=1.5#, main=paste0("With RRBS ", Biol_Comparands[k]), cex.main=0.9
    )
    col.vec <- c("firebrick","grey50","firebrick","grey50")
    print(obser.num.scale)
    points(seq(1,2,by=1), obser.num.scale, pch=23, col=col.vec, 
           bg=col.vec, cex=1.5)
    text(x=seq(1,2,by=1), y=rep(c(yMax+2, yMin-2), 2), paste0("p=",pval.vec), cex=1)
  }
  dev.off()
  #abline(h=NegCutoffZscore, col="darkgray", lty=5)
  #abline(h=PosCutoffZscore, col="darkgray", lty=5)
}

#




# OSCC plasma + HNSC TCGA cancer (plot each heatmap seperate) -------------------
cpm.list <- list()
diff.list <- list()
rna <- c("rRNA","pri_miRNA","lncRNA", "mRNA", "piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA')
dna <- c("intron","promoter", "enhancer", "repeats") 

## get GSE110381 diff/CPM/RPKM stats
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
dst <- "lulab_oscc_plasma_diff"

### read smp table
sample.table <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/lulab_oscc_plasma_diff/sample_table.txt",sep = "\t", check.names = F,header = T)
sample.table <- sample.table[,c("Patient","Metastasis")]
colnames(sample.table) <- c("sample","group")

### read count
count <- read.table(paste0(pre,"/output/",dst,"/call_peak_dedup/count_matrix/cfpeakCNN_b5_d50_p1.txt"),check.names = F,header = T)
rownames(count) <- count$feature
count <- count[,2:ncol(count)]
sample.table <- sample.table[sample.table$sample %in% colnames(count),]
sample.table$group <- gsub("Yes","Y",sample.table$group)
sample.table$group <- gsub("No","W",sample.table$group)
rownames(sample.table) <- sample.table$sample
# #sample.table <- sample.table[order(sample.table$group),]
# positive_samples <- sample.table[sample.table$group=="Metastasis","sample"]
# negative_samples <- sample.table[sample.table$group=="Localized","sample"]
# samples <- c(positive_samples, negative_samples)
# sample.table <- sample.table[match(samples,sample.table$sample),]
# group <- c(rep("positive",length(positive_samples)),rep("negative",length(negative_samples)))
# method <- "edger_glmlrt"
# norm_method <- "TMM"

cpm <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",dst,"/filterSmp_exPeakCNN_smallDomain_diff_YvsW.cpm"),check.names = F,header = T)
cpm$feature <- rownames(cpm)
cpm <- cpm[,c(ncol(cpm),1:(ncol(cpm)-1))]
cpm.list[[dst]] <- cpm
diff.list[[dst]] <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",dst,"/filterSmp_exPeakCNN_smallDomain_diff_YvsW.diff"),check.names = F,header = T)
diff.list[[dst]]$feature <- rownames(diff.list[[dst]])


cpm.cutoff <- 50 #no filter
#table(df.mean$log.abs.mean.coef>=coef.cutoff)
table(diff.list[[dst]]$cpm.mean<=cpm.cutoff )

df.mean.filter <- as.data.frame(diff.list[[dst]])
rownames(df.mean.filter) <- df.mean.filter$feature


### filter logcpm mat
cpm.filter <- cpm.list[[dst]]
rownames(cpm.filter) <- cpm.filter$feature
cpm.filter$feature <- NULL
cpm.filter <- cpm.filter[df.mean.filter$feature, sample.table$sample] #cpm.filter[rownames(cpm.filter) %in% df.mean.filter$feature, colnames(cpm.filter) %in% sample.table$sample]

df.mean.filter[1:3,]
df.mean.filter$tx <- unlist(sapply(strsplit(df.mean.filter$feature,"|",fixed = T),"[",4))
df.mean.filter$type <- unlist(sapply(strsplit(df.mean.filter$feature,"|",fixed = T),"[",2))
df.mean.filter$type <- gsub("_rev|_for","",df.mean.filter$type,perl = T)
table(df.mean.filter$type)

## order anno tbls
annotation_col_gse <- sample.table[colnames(cpm.filter),]

annotation_col_gse$sample <- NULL
annotation_row_gse <- data.frame( RNA = df.mean.filter[rownames(cpm.filter),"type"],
                                  row.names = rownames(cpm.filter)) # , diff=genes$d


### get DCB candidates peaks
dst <- "lulab_oscc_plasma_diff"
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/"
tmp.diff <- diff.list[[dst]][match(rownames(cpm.filter),diff.list[[dst]]$feature),]
dcb.gse <- rownames(cpm.filter)[(tmp.diff$pvalue <= 0.1)]
# table((tmp.diff$pvalue <= 0.1))
# table(tmp.diff$log2FoldChange > 0)
# table((tmp.diff$posGT1RatioCount >= max(0.05,2/7)))
# table(tmp.diff$negGT1RatioCount <= 0.1)
# table(tmp.diff$negGT1RatioCount <= max(0.1,2/7))
# dcb.gse <- rownames(cpm.filter)[(tmp.diff$pvalue <= 0.2) & (tmp.diff$log2FoldChange > 0) & (tmp.diff$posGT1RatioCount >= max(0.05,2/7)) & (tmp.diff$negGT1RatioCount <= max(0.1,2/7)) ]  # & (tmp.diff$negCent90CPM <= 0.1)
#7 Meta vs. 7 local, require minimum 2 samples in postive samples: 2/7
#table(tmp.diff$negGT1RatioCount <= 0.2)  # 44844
#table(tmp.diff$negMeanCPM <= 1)
length(dcb.gse) # p0.01:389, p0.05:1472, p0.1:2937
#87




# ## select peak by top 
# peak.list <- list()
# for (top in c(50,100,200)){
#   peak.list[[paste0("top",top)]] <- list()
# }
# disease <- "metastasis"
# disease.label <- "Y"
# normal.label <- "W"
# 
# tmp <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",dst,"/filterSmp_exPeakCNN_smallDomain_diff_",disease.label,"vs",normal.label,".diff"),sep = "\t",header = T,check.names = F,stringsAsFactors = F)
# tmp <- tmp[order(as.numeric(tmp$pvalue),decreasing = F),]
# rownames(tmp) <- unlist(sapply(strsplit(rownames(tmp),"|",fixed=T),"[",4))
# 
# # op1: fixed cutoff
# #tmp.sig <- tmp[tmp$pvalue<0.00001,] 
# 
# # op2: fixed top
# for (i in c(50,100,200)){
#   peak.id <- rownames(tmp)[1:i]
#   # print(length(peak.id))
#   peak.list[[paste0("top",i)]][[disease]] <- peak.id
# }
# peak.ids <- list()
# for (top in c(50,100,200)){
#   peak.ids[[paste0("top",top)]] <- unique(do.call("c",peak.list[[paste0("top",top)]])) 
#   print(length( peak.ids[[paste0("top",top)]] ))
# }
# 
# 
# 
# 
# 
# ## heatmap
# ref <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt") #,header=T
# ref <- as.data.frame(ref)
# rownames(ref) <- ref$transcript_id
# ref[1:3,1:3]
# 
# disease.label <- "Y"
# normal.label <- "W"
# 
# logcpm <- cpm.filter
# logcpm <- logcpm[unlist(sapply(strsplit(rownames(logcpm),"|",fixed=T),"[",4)) %in% peak.ids[["top200"]],] # rownames(diff)[diff$pvalue<=0.01]
# 
# hemoTree2 <- stats::hclust(dist(logcpm, method = "euclidean"), method = "complete")  # method: "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"
# logcpm <- logcpm[hemoTree2$labels,]
# all(rownames(logcpm)==hemoTree2$labels)
# 
# annotation_col <- sample.table[colnames(logcpm),]
# annotation_col$group[annotation_col$group=="Y"] <- "Yes"
# annotation_col$group[annotation_col$group=="W"] <- "No"
# annotation_col <- annotation_col[,c("group"),drop=F]
# table(annotation_col$group)
# 
# txid <- unlist(sapply(strsplit(rownames(logcpm),"|",fixed=T),"[",3))
# annotation_row <- data.frame( RNA = ref[txid,"transcript_type"],
#                               row.names = rownames(logcpm))
# annotation_row$RNA <- gsub("_for|_rev","",annotation_row$RNA,perl = T)
# tmp <- as.data.frame(table(annotation_row$RNA))
# tmp$lab <- paste0(tmp$Var1," (n=",tmp$Freq,")")
# annotation_row$RNA <-tmp$lab[match(annotation_row$RNA,tmp$Var1)]
# 
# RNA_colors <- list()
# for(i in 1:length(c(rna,dna))){
#   j <- c(rna,dna)[i]
#   RNA_colors[[j]] <- c(pal_nejm_adaptive()(15)[1:14],"#11838D")[i] # pal_d3_adaptive()(15)
# }
# RNA_colors <- do.call("c",RNA_colors)
# names(RNA_colors) <- tmp$lab[match(names(RNA_colors),tmp$Var1)]
# RNA_colors <- RNA_colors[!is.na(names(RNA_colors))]
# ann_colors <- list(group=c("Yes"="firebrick2","No"="orange2"),RNA=RNA_colors) # ,"salmon",seagreen,orange,chocolate
# dim(logcpm)
# pheatmap::pheatmap(mat = logcpm,
#                    annotation_col = annotation_col, #data frame that specifies the annotations shown on left side of the heatmap.
#                    annotation_row = annotation_row,
#                    annotation_colors = ann_colors, # RColorBrewer::brewer.pal(6, "Set3")
#                    border_color="NA",
#                    scale = "row", # row,none
#                    #labels_col = 3, labels_row = 6,
#                    cluster_cols = T,cluster_rows = T,
#                    # gaps_col = 24,
#                    #cutree_cols = 2,cutree_rows = 3,
#                    show_colnames=T, show_rownames=F,
#                    fontsize = 12,
#                    height = 8, width = 7,
#                    color = colorRampPalette(RColorBrewer::brewer.pal(n = 10, name = "RdBu"))(100),  # steelblue/dodgerblue4-firebrick2, viridis::viridis_pal()(100), #
#                    #fontsize_row = 5,
#                    filename = paste0("./oscc_exPeakCNN_smallDomain_diff_",disease.label,"vs",normal.label,".cpm.heatmap3.pdf")
#                    # filename = paste0("./output/",dst,"/allSmpP01_exPeak_smallDomain_diff_",disease.label,"vs",normal.label,".cpm.heatmap.pdf")
# )

#







## tcga heatmap 
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
dst <- "TCGA-HNSC_small_diff"
### read smp table
sample.table.tcga <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/TCGA-HNSC_small_diff/sample_table.txt",check.names = F, sep = "\t", header = T)
sample.table.tcga[1:3,]
#"id","metastasis","anatomic_neoplasm_subdivision"
sample.table.tcga <- sample.table.tcga[,1:2]

colnames(sample.table.tcga) <- c("sample","group")
sample.table.tcga$sample <- gsub(".bam","",sample.table.tcga$sample)
rownames(sample.table.tcga) <- sample.table.tcga$sample
table(sample.table.tcga$group)
sample.table.tcga[1:3,]
smpid <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/TCGA-HNSC_small_diff/sample_ids.txt")$V1
#write.table(sample.table.tcga[sample.table.tcga$sample %in% smpid,],"/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/TCGA_small_diff/call_peak_all/too/smps.csv",sep = ",",quote = F,row.names = F,col.names = T)
sample.table.tcga$group <- factor(sample.table.tcga$group,levels = c("Y", "N")) # , "PAAD", "PRAD"
sample.table.tcga <- sample.table.tcga[order(sample.table.tcga$group),]


### read count mat
count <- read.table(paste0(pre,"/output/",dst,"/call_peak_all/count_matrix/cfpeakCNN_b5_d50_p1_lulab_oscc_plasma.txt"),check.names = F,header = T)
rownames(count) <- count$feature
count <- count[,2:ncol(count)]
sample.table.tcga <- sample.table.tcga[sample.table.tcga$sample %in% colnames(count),]


#diff.list <- list()
cpm.tcga <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",dst,"/filterSmp_exPeakCNN_smallDomain_lulab_oscc_plasma_diff_YvsW.cpm"),row.names = 1,check.names = F,header = T)
cpm.list[[dst]] <- cpm.tcga
diff.list[[dst]] <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",dst,"/filterSmp_exPeakCNN_smallDomain_lulab_oscc_plasma_diff_YvsW.diff"),check.names = F,header = T)
diff.list[[dst]]$feature <- rownames(diff.list[[dst]])
dim(cpm.tcga)

### filter logcpm mat
cpm.tcga.filter <- cpm.tcga[df.mean.filter$feature, sample.table.tcga$sample] 
#cpm.tcga.filter <- cpm.tcga
cpm.tcga.filter <- cpm.tcga.filter[apply(cpm.tcga.filter,1,sum)>0,] # rm 0 sum
cpm.tcga.filter <- cpm.tcga.filter[apply(cpm.tcga.filter,1,sd)>0,] # rm 0 std
#table(apply(cpm.tcga.filter,1,sum)==0)

### order anno tbls
annotation_col_tcga <- sample.table.tcga[colnames(cpm.tcga.filter),]
annotation_col_tcga$sample <- NULL
annotation_row_tcga <- data.frame( RNA = df.mean.filter[rownames(cpm.tcga.filter),"type"],
                                   row.names = rownames(cpm.tcga.filter)) # , diff=genes$diff


### get DCB candidates peaks
dst <- "TCGA-HNSC_small_diff"
pre <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/"
tmp.diff <- diff.list[[dst]][match(rownames(cpm.tcga.filter),diff.list[[dst]]$feature),]
#table(tmp.diff$negCent90CPM<=0.5)
dcb.tcga <- rownames(cpm.tcga.filter)[(tmp.diff$pvalue <= 0.1)] # just diff peak is ok
# table(tmp.diff$pvalue <= 0.1)
# table(tmp.diff$negGT1RatioCount <= max(0.1,2/20)) 
# dcb.tcga <- rownames(cpm.tcga.filter)[(tmp.diff$pvalue <= 0.2) & (tmp.diff$log2FoldChange > 0) & (tmp.diff$posGT1RatioCount >= max(0.1,2/20))  & (tmp.diff$negGT1RatioCount <= max(0.1,2/20)) ]  # this version not use ex libsize
#20 vs. 20, required min 2 sampels in positive sample: tmp.diff$posGT1RatioCount >= max(0.1,2/20)
#table(tmp.diff$negGT1RatioCount <= 0.2) 
# data.table::fwrite(x = as.data.frame(dcb.tcga),file = paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_all/count_matrix/expeakCNN_b5_d50_p1_GSE110381.txt.diff.dcb"),quote = F,sep = "\t",row.names = F,col.names = F)
length(dcb.tcga) # pval0.2:3976; pval0.1:2185; pval0.05:1212; pval0.01:495; 
#490




## merge  DCB candidates peaks, plot venn (2 dst)
df <- data.frame(peak=unique(c(dcb.gse,dcb.tcga)),gse=0,tcga=0)
rownames(df) <- df$peak
df[dcb.gse,"gse"] <- 1
df[dcb.tcga,"tcga"] <- 1
df$peak <- NULL
dcb <- rownames(df)[rowSums(df)==2]
#dcb <- dcb[1:20]
dcb.peak <- unlist(sapply(strsplit(dcb,"|",fixed=T),"[",3))
dcb.RNA <- unlist(sapply(strsplit(dcb,"|",fixed=T),"[",2))
dcb.df <- data.frame(dcb=dcb,dcb.RNA=dcb.RNA,dcb.peak=dcb.peak)
dcb.df$dcb.RNA <- gsub("_rev|_for","",dcb.df$dcb.RNA,perl = T)
dcb.df$dcb.RNA <- factor(dcb.df$dcb.RNA, levels = c(rna,dna))
dcb.df <- dcb.df[order(dcb.df$dcb.RNA,decreasing = F),]
dcb <- dcb.df$dcb
dcb.peak <- dcb.df$dcb.peak
#dcb.peak[4:5]
length(dcb.peak)
dcb[duplicated(dcb.peak)]

library(ggvenn)
x <- list(dcb.gse,dcb.tcga)
names(x) <- c("GSE110381 plasma", "TCGA tissue")
ggvenn(data = x, stroke_color = "white",text_color = "black", text_size = 10,
       fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF"), # , "#868686FF"
       stroke_size = 0.5, set_name_size = 8
)
ggsave("OSCC_dcb_venn.pdf",width = 9,height = 9)




# ## heatmap of lulab_oscc_plasma
# ref <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt") #,header=T
# ref <- as.data.frame(ref)
# rownames(ref) <- ref$transcript_id
# ref[1:3,1:3]
# 
# dst <- "lulab_oscc_plasma_diff"
# disease.label <- "Y"
# normal.label <- "W"
# 
# diff <- read.table(paste0("./output/",dst,"/filterSmp_exPeakCNN_smallDomain_diff_",disease.label,"vs",normal.label,".diff"), header = T, sep="\t",check.names = F)
# logcpm <- read.table(paste0("./output/",dst,"/filterSmp_exPeakCNN_smallDomain_diff_",disease.label,"vs",normal.label,".cpm"), header = T, sep="\t",check.names = F)
# logcpm <- logcpm[matrixStats::rowSds(as.matrix(logcpm))!=0,] # rm 0 std, or meet error for pca
# 
# logcpm <- logcpm[unlist(sapply(strsplit(rownames(logcpm),"|",fixed=T),"[",4)) %in% peak.ids[["top100"]],] # rownames(diff)[diff$pvalue<=0.01]
# #table(diff$padj<=0.01)
# 
# 
# annotation_col <- sample.table.pair[colnames(logcpm),]
# annotation_col <- annotation_col[,c("group"),drop=F]
# 
# txid <- unlist(sapply(strsplit(rownames(logcpm),"|",fixed=T),"[",3))
# annotation_row <- data.frame( RNA = ref[txid,"transcript_type"],
#                               row.names = rownames(logcpm))
# annotation_row$RNA <- gsub("_for|_rev","",annotation_row$RNA,perl = T)
# tmp <- as.data.frame(table(annotation_row$RNA))
# tmp$lab <- paste0(tmp$Var1," (n=",tmp$Freq,")")
# annotation_row$RNA <-tmp$lab[match(annotation_row$RNA,tmp$Var1)]
# table(annotation_row$RNA)
# annotation_row[1:3,]
# 
# rna <- c("rRNA","pri_miRNA","lncRNA", "mRNA", "piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA')
# dna <- c("intron","promoter", "enhancer", "repeats") 
# RNA_colors <- list()
# for(i in 1:length(c(rna,dna))){
#   j <- c(rna,dna)[i]
#   RNA_colors[[j]] <- c(pal_nejm_adaptive()(15)[1:14],"#11838D")[i] # pal_d3_adaptive()(15)
# }
# RNA_colors <- do.call("c",RNA_colors)
# names(RNA_colors) <- tmp$lab[match(names(RNA_colors),tmp$Var1)]
# RNA_colors <- RNA_colors[!is.na(names(RNA_colors))]
# ann_colors <- list(group=c("Yes"="firebrick2","No"="orange2"),RNA=RNA_colors) # ,"salmon",seagreen,orange,chocolate
# 
# pheatmap::pheatmap(mat = logcpm,
#                    annotation_col = annotation_col, #data frame that specifies the annotations shown on left side of the heatmap.
#                    annotation_row = annotation_row,
#                    annotation_colors = ann_colors, # RColorBrewer::brewer.pal(6, "Set3")
#                    border_color="NA",
#                    scale = "row", # row,none
#                    #labels_col = 3, labels_row = 6,
#                    cluster_cols = T,cluster_rows = T,
#                    # gaps_col = 24,
#                    #cutree_cols = 2,cutree_rows = 3,
#                    show_colnames=T, show_rownames=F,
#                    fontsize = 12,
#                    height = 8, width = 7,
#                    color = colorRampPalette(RColorBrewer::brewer.pal(n = 10, name = "RdBu"))(100),  # steelblue/dodgerblue4-firebrick2, viridis::viridis_pal()(100), #
#                    #fontsize_row = 5,
#                    filename = paste0("./output/",dst,"/FilterSmp2P01_exPeakCNN_smallDomain_diff_",disease.label,"vs",normal.label,".cpm.heatmap3.pdf")
#                    # filename = paste0("./output/",dst,"/allSmpP01_exPeak_smallDomain_diff_",disease.label,"vs",normal.label,".cpm.heatmap.pdf")
# )
# #

## plot merge heatmap for overlapped feature
### pre-scale each matrix to prevent 3-dst-scaling together)
cpm.filter.norm <- t(scale(t(cpm.filter))) # default: zscore each column
# cpm.inhouse.filter.norm <- t(scale(t(cpm.inhouse.filter)))
cpm.tcga.filter.norm <- t(scale(t(cpm.tcga.filter)))
#colMeans(cpm.filter.norm)
#colSds(cpm.filter.norm)
#rowMeans(cpm.filter.norm)
#rowSds(cpm.filter.norm)


# txid <- unlist(sapply(strsplit(rownames(logcpm),"|",fixed=T),"[",3))
# annotation_row <- data.frame( RNA = ref[txid,"transcript_type"],
#                               row.names = rownames(logcpm))
# annotation_row$RNA <- gsub("_for|_rev","",annotation_row$RNA,perl = T)
# tmp <- as.data.frame(table(annotation_row$RNA))
# tmp$lab <- paste0(tmp$Var1," (n=",tmp$Freq,")")
# annotation_row$RNA <-tmp$lab[match(annotation_row$RNA,tmp$Var1)]
# table(annotation_row$RNA)
# annotation_row[1:3,]


RNA_colors <- list()
for(i in 1:length(c(rna,dna))){
  j <- c(rna,dna)[i]
  RNA_colors[[j]] <- c(pal_nejm_adaptive()(15)[1:14],"#11838D")[i] #pal_d3_adaptive()(15)[i]
}
RNA_colors <- do.call("c",RNA_colors)


#2dst (inhouse OSCC + TCGA)
cpm.merge <- as.data.frame(cbind(cpm.filter.norm[dcb,],cpm.tcga.filter.norm[dcb,])) # cpm.inhouse.filter.norm, cpm.filter.norm
annotation_col_merge <- rbind(annotation_col_gse[,"group",drop=F],annotation_col_tcga[,"group",drop=F])
annotation_col_merge$dataset <- c(rep("OSCC plasma",nrow(annotation_col_gse)),rep("TCGA-HNSC tissue",nrow(annotation_col_tcga)))
table(annotation_col_merge$dataset)
annotation_col_merge$group <- gsub("Y","Metastatic",annotation_col_merge$group )
annotation_col_merge$group <- gsub("W","Localized",annotation_col_merge$group )
annotation_col_merge$group <- gsub("N","Localized",annotation_col_merge$group )

annotation_row_merge <- annotation_row_gse[dcb,,drop=F]
ann_colors <- list(group=c("Metastatic"="firebrick","Localized"="grey70"),dataset=c(`OSCC plasma`="chocolate",`TCGA-HNSC tissue`="steelblue"),RNA=RNA_colors) # "Metastatic HNSC tissue"="firebrick4","Localized HNSC tissue"="grey50","salmon",seagreen,orange,chocolate
table(annotation_col_merge$dataset )
annotation_col_merge$dataset <- factor(annotation_col_merge$dataset,levels = c("OSCC plasma","TCGA-HNSC tissue"))
annotation_col_merge$group <- factor(annotation_col_merge$group,levels = c("Metastatic","Localized"))
cpm.merge[1:3,1:3]
#table(annotation_col_merge$dataset)
str(annotation_col_merge)
annotation_col_merge <- annotation_col_merge[order(annotation_col_merge$dataset,annotation_col_merge$group),]
cpm.merge <- cpm.merge[,rownames(annotation_col_merge)]
pheatmap::pheatmap(mat = cpm.merge,
                   annotation_col = annotation_col_merge, #data frame that specifies the annotations shown on left side of the heatmap.
                   annotation_row = annotation_row_merge,
                   annotation_colors = ann_colors, # RColorBrewer::brewer.pal(6, "Set3")
                   border_color="NA", # "grey20",
                   scale = "none", # row
                   #labels_col = 3, labels_row = 6,
                   cluster_cols = F,cluster_rows = T,
                   gaps_col = 14,
                   #cutree_cols = 2,cutree_rows = 3,
                   show_colnames=F, show_rownames=F,
                   fontsize = 12,
                   # height = 7,width =10,
                   height = 7,width =15,
                   color = colorRampPalette(RColorBrewer::brewer.pal(n = 10, name = "RdBu"))(100),  # colorRampPalette(RColorBrewer::brewer.pal(n = 10, name = "RdBu"))(100),  #viridis::viridis_pal()(100), #colorRampPalette(c("blue","white","red"))(100),  # steelblue/dodgerblue4-firebrick2
                   #fontsize_row = 5,
                   filename ="./merge_oscc-tcga_heatmap.pdf") # merge_gse-tcga_heatmap_2pos.pdf


#1dst: lulab
cpm.merge1 <- as.data.frame((cpm.filter.norm[dcb,])) # cpm.inhouse.filter.norm, cpm.filter.norm
annotation_col_merge1 <- annotation_col_merge[annotation_col_merge$dataset=="OSCC plasma",]
dim(cpm.merge1)
pheatmap::pheatmap(mat = cpm.merge1,
                   annotation_col = annotation_col_merge1, #data frame that specifies the annotations shown on left side of the heatmap.
                   annotation_row = annotation_row_merge,
                   annotation_colors = ann_colors, # RColorBrewer::brewer.pal(6, "Set3")
                   border_color="NA", # "grey20",
                   scale = "none", # row
                   #labels_col = 3, labels_row = 6,
                   cluster_cols = T,cluster_rows = T,
                   gaps_col = 14,
                   #cutree_cols = 2,cutree_rows = 3,
                   show_colnames=F, show_rownames=F,
                   fontsize = 12,
                   # height = 7,width =10,
                   height = 7,width =15,
                   color = colorRampPalette(RColorBrewer::brewer.pal(n = 10, name = "RdBu"))(100),  # colorRampPalette(RColorBrewer::brewer.pal(n = 10, name = "RdBu"))(100),  #viridis::viridis_pal()(100), #colorRampPalette(c("blue","white","red"))(100),  # steelblue/dodgerblue4-firebrick2
                   #fontsize_row = 5,
                   filename ="./oscc-inhouse_heatmap.pdf") # merge_gse-tcga_heatmap_2pos.pdf



## select eg. tx
for (peak in rownames(annotation_row_merge)){
  #peak <- "ENST00000385214_____1_30_49_+|pri_miRNA|ENST00000385214_____1|peak_6181|ENST00000385214_____1|30|49"  # "T295851_13316_13332_+|tucpRNA|T295851|peak_57431|T295851|13316|13332"
  # peak <- "L2a__chr16___3509264____3509549_neg_84_99_+|repeats_rev|L2a__chr16___3509264____3509549_neg|peak_48841|L2a__chr16___3509264____3509549_neg|84|99"
  for (dataset in unique(annotation_col_merge$dataset)){
    #dataset <- "GSE110381 plasma"
    tmp1.annotation_col_merge <- annotation_col_merge[annotation_col_merge$dataset==dataset,]
    tmp1 <- as.data.frame(cpm.merge[peak,annotation_col_merge$dataset==dataset,drop=F])
    for (group in unique(tmp1.annotation_col_merge$group)){
      #group <- "CRC plasma"
      # annotation_row_merge[peak,paste0(dataset,"_",group)] <- ""
      
      tmp2.annotation_col_merge <- tmp1.annotation_col_merge[tmp1.annotation_col_merge$group==group,]
      tmp2 <- as.data.frame(t(tmp1[,tmp1.annotation_col_merge$group==group,drop=F]))
      if(grepl("CRC",group)){
        tmp2 <- tmp2[order(tmp2[,1],decreasing = T),,drop=F]
      } else {
        tmp2 <- tmp2[order(tmp2[,1],decreasing = F),,drop=F]
      }
      top3 <- rownames(tmp2)[1:3]
      annotation_row_merge[peak,paste0(dataset,"_",group,"_",1:3)] <- top3
    }
  }
}
annotation_row_merge$peak <- unlist(sapply(strsplit(rownames(annotation_row_merge),"|",fixed=T),"[",4))
annotation_row_merge$chr <- unlist(sapply(strsplit(rownames(annotation_row_merge),"|",fixed=T),"[",5))
annotation_row_merge$start <- unlist(sapply(strsplit(rownames(annotation_row_merge),"|",fixed=T),"[",6))
annotation_row_merge$end <- unlist(sapply(strsplit(rownames(annotation_row_merge),"|",fixed=T),"[",7))
annotation_row_merge$width <- as.numeric(annotation_row_merge$end) - as.numeric(annotation_row_merge$start)
annotation_row_merge$idx <- 1:nrow(annotation_row_merge)
annotation_row_merge <- annotation_row_merge[,c((ncol(annotation_row_merge)-5):ncol(annotation_row_merge),1:(ncol(annotation_row_merge)-6))]
#use annotation_row_merge for IGV eg.


write.table(x = annotation_row_merge, file = "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/tmp/OSCC_annotation_row_merge.txt",quote = F,sep = "\t",row.names = F,col.names = T)

#


# test --------------------------------------------------------------------
tx.tb <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008_NCpool/call_peak_all/expeak_by_sample/b5_d50_p1/intersect/NCpool.bed",data.table = F,header = F,sep = "\t",stringsAsFactors = F)
hist(tx.tb$V10,breaks = 100)
table(tx.tb$V10>=0.05)
table(tx.tb$V9<=0.05)
max(na.omit(tx.tb$V9))

tx.tb[tx.tb$,][1:20,]
tail(tx.tb)

hist(tx.tb$V20)
hist(tx.tb$V20/(tx.tb$V3-tx.tb$V2))
table(tx.tb$V20/(tx.tb$V3-tx.tb$V2) >= 0.5)

chr.size <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/transcriptome_genome",data.table = F,header = F,sep = "\t",stringsAsFactors = F)
peak <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax/domains.bed",data.table = F,header = F,sep = "\t",stringsAsFactors = F)
table(peak$V1 %in% chr.size$V1)



bed <- data.table::fread("/Users/baopengfei/Desktop/lulab/tmp/projects/WCHSU-FTC/output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax_by_sample/csFTA-17_1.bed",data.table = F,header = F,sep = "\t",stringsAsFactors = F)
bed[1:3,]
table(duplicated(paste0(bed$V1,bed$V2,bed$V3)))



library(IRanges)
library(GenomicRanges)
bed <- data.table::fread("/Users/baopengfei/Desktop/lulab/tmp/projects/WCHSU-FTC/output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax_by_sample/csFTA-17_1.bed",data.table = F,header = F,sep = "\t",stringsAsFactors = F)




library(tidyverse)
library(dplyr)
setwd("/Users/baopengfei/Desktop/lulab/tmp/projects/WCHSU-FTC/")

## read ref
# ref <- data.table::fread("exSeek-dev/genome/hg38/tbed/10RNA_map.txt",sep = "\t",header = F,stringsAsFactors = F)
# rownames(ref) <- ref$V1
# ref[1:3,]
chr.size <- data.table::fread("exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length.txt",sep = "\t",header = T,stringsAsFactors = F)
chr.size[1:3,]


bw0 <- rtracklayer::import ("output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/tbigwig_long_RNA/csFTA-17_1.transcriptome.bigWig",format="BigWig")


bw.tbl <- as_tibble(bw0)
bw.tbl$start <- bw.tbl$start-1
bw.tbl$tx.length = chr.size$tx.length[match(bw.tbl$seqnames,chr.size$transcript_id)]
bw.tbl$type = chr.size$transcript_type[match(bw.tbl$seqnames,chr.size$transcript_id)]

bw.tbl <- bw.tbl %>% 
  group_by(seqnames) %>% 
  mutate(tx.cov = sum(width * score) / tx.length, tx.cov.zerotrunc = sum(width * score) / sum(width) ,
         ratio = tx.cov/tx.cov.zerotrunc )


library(ggplot2)
p1 <- bw.tbl %>% 
  ggplot(breaks=100,aes(x=tx.cov))+
  geom_histogram(,binwidth = 0.1)+
  xlim(c(0,10))+
  facet_grid(type~.,scales = "free_y")
ggsave(p1,'hist.pdf',width = 7, height = 21)

p2 <- bw.tbl %>% 
  ggplot(breaks=100,aes(x=tx.cov.zerotrunc))+
  geom_histogram(,binwidth = 0.1)+
  xlim(c(0,50))+
  facet_grid(type~.,scales = "free_y")
ggsave(p2, 'hist_tx.cov.zerotrunc.pdf',width = 7, height = 21)

p3 <- bw.tbl %>% 
  ggplot(breaks=100,aes(x=ratio))+
  geom_histogram(,binwidth = 0.1)+
  xlim(c(0,1))+
  facet_grid(type~.,scales = "free_y")
ggsave(p3, 'hist_tx.cov.ratio.pdf',dth = 7, height = 21)






domain <- rtracklayer::import ("output/FTC_small/call_domain_noRepeats_MQ0_localmaxMin5_noQvalue/domains_localmax_by_sample/csFTA-17_1.bed",format="bed")
domain.tbl <- as_tibble(domain)
domain.tbl$start <- domain.tbl$start-1
domain.tbl$tx.length = chr.size$tx.length[match(domain.tbl$seqnames,chr.size$transcript_id)]
domain.tbl$type = chr.size$transcript_type[match(domain.tbl$seqnames,chr.size$transcript_id)]
table(domain.tbl$type)
table(domain.tbl$seqnames %in% bw.tbl$seqnames )
domain.tbl$tx.cov <- bw.tbl$tx.cov[match(domain.tbl$seqnames,bw.tbl$seqnames)]
domain.tbl$tx.cov.zerotrunc <- bw.tbl$tx.cov.zerotrunc[match(domain.tbl$seqnames,bw.tbl$seqnames)]
domain.tbl$bg.max.ratio <- domain.tbl$tx.cov/domain.tbl$score

p4 <- domain.tbl %>% 
  ggplot(breaks=100,aes(x=bg.max.ratio))+
  geom_histogram(,binwidth = 0.1)+
  xlim(c(0,1))+
  facet_grid(type~.,scales = "free_y")
p4
#ggsave(p4, 'hist.pdf',dth = 7, height = 21)




tx.tb <- data.table::fread("genome/hg38/transcript_table/all.txt",data.table = F,header = F,sep = "\t",stringsAsFactors = F)
library(dplyr)
tx.tb <- as_tibble(tx.tb)
table(duplicated(tx.tb$V8))

tx.tb[duplicated(tx.tb$V8),]
#ENST00000607781.1
tx.tb[tx.tb$V8=="ENST00000607781.1",]



setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/")
tmp <- read.table("./test.txt")
head(tmp)
hist(tmp$V1)
tmp$V2 <- 1-log2(tmp$V1+1)
tmp$V3 <- 3+40*(1-log2(tmp$V1+1))
hist(tmp$V2)
hist(tmp$V3)

left <- read.table("tmp/left.txt")
right <- read.table("tmp/right.txt")
summary(left$V7 + right$V7)
hist((left$V7 + right$V7), xlim = c(0,50), breaks = 1000)



tmp <- read.table("/lulabdata/baopengfei/shared_reference/tRNA/metadata_tsRFun.txt",header = T)
table(duplicated(tmp$seq))




# test run-time --------------------------------
#cat ../output/GSE71008_NCpool/call_peak_all/archive/expeak_by_sample/b5_d05_p05/log/NCpool.log | grep "^2023" | cut -d"." -f 1 | sed s/"^ "/""/g > ../output/GSE71008_NCpool/call_peak_all/archive/expeak_by_sample/b5_d05_p05/log/NCpool.log2
#cat ../output/GSE71008_NCpool/call_peak_all/archive/expeak_by_sample/b5_d05_p05/log/NCpool.log | grep "start chrom" | sed s/"^ "/""/g > ../output/GSE71008_NCpool/call_peak_all/archive/expeak_by_sample/b5_d05_p05/log/NCpool.log3
tmp <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008_NCpool/call_peak_all/archive/expeak_by_sample/b5_d05_p05/log/NCpool.log2")
tmp$V3 <- paste0(tmp$V1, " ", tmp$V2, " ", "CST")
tmp[1:3,1:3]
# tmp$V1 <- as.Date(tmp$V1)

l <- list()
for(i in 1:(nrow(tmp)-1)){
  if(i%%10000==0){
    print(i)
  }
  # i <- 1
  l[[i]] <- as.numeric(difftime(tmp$V3[i+1],tmp$V3[i],units = "secs"), units = "secs")
}
df <- do.call(rbind,l)
df[100]
dim(df)
hist(df[,1],xlim = c(-50,100), breaks = 100)
summary(df[,1])
