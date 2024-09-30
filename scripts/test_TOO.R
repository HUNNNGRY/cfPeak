#sncRNA TOO


setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/")
#sample_ids <- read.table("./data/GSE163534/sample_ids.txt",header = F)
#sample_ids <- read.table("./data/GSE71008/sample_ids.txt",header = F)
sample_ids <- read.table("./data/FTC_small/sample_ids.txt",header = F)

#sample_tabl <- read.table("./data/GSE163534/sample_table_allhsa.txt",header = T,sep = '\t')
#sample_tabl <- read.table("./data/GSE71008/sample_table.txt",header = T,sep = '\t')
sample_tabl <- read.table("./data/FTC_small/sample_table.txt",header = T,sep = '\t')
sample_tabl <- sample_tabl[sample_tabl$sample %in% sample_ids$V1,]
colnames(sample_tabl)
sample_tabl <- sample_tabl[,c("sample","group")]
# colnames(sample_tabl) <- c("sample","group")
# sample_tabl <- sample_tabl[sample_tabl$group  %in% c("FTA_cf","FTC_cf"),] #  %in% c("FTA_EV","FTC_EV")
table(sample_tabl$group)
#sample_tabl <- sample_tabl[as.character(sample_tabl$group) %in% c("bone","bowel-colon","esophagus","pancreas","prostate","bladder","kidney","thyroid","liver","lung","lymph node","stomach","spleen","artery","vein"),]
#c("adipocyte","artery","vein","bladder","bone","bowel-colon","bowel-small intestine","esophagus","kidney","liver","lung","lymph node","pancreas","stomach","spleen","thyroid")

# /usr/bin/Rscript ~/projects/multi-omics-explore/scripts/run-NormCountMat.R \
# -i ../output/GSE163534/call_domain_withRepeats_all/count_matrix/domains_long_localmax2_b5_d05_p05.txt \
# -o ../output/GSE163534/call_domain_withRepeats_all/count_matrix/domains_long_localmax2_b5_d05_p05_CPM-TMM.txt

#outdir <- "../output/GSE71008/call_domain_withRepeats_all/too"
outdir <- "../output/FTC_small/call_domain_withRepeats_dedup/too_cf"
dir.create(outdir,recursive = T)
#mat <- read.table("../output/GSE71008/call_domain_withRepeats_all/count_matrix/domains_long_localmax2_b5_d05_p01_CPM-TMM.txt",header = T,row.names = 1)
mat <- read.table("../output/FTC_small/call_domain_withRepeats_dedup/count_matrix/domains_long_localmax2_b5_d05_p01_CPM-TMM.txt",header = T,row.names = 1,check.names = F)
sig <- read.table("../output/GSE163534/call_domain_withRepeats_all/too/signature_genes_1_0.5_50.txt",header = T,row.names = 1)
# mat <- mat[rownames(mat) %in% rownames(sig),]
mat <- mat[,sample_tabl$sample]
#mat[1:3,1:3]

all(sample_tabl$accession==sample_ids$V1)
#all(colnames(mat)==sample_ids$V1)

#write.table(sample_tabl,"./data/GSE163534/sample_table_filter.csv",sep = ",",quote = F,row.names = F,col.names = T)
#write.table(mat,"../output/GSE163534/call_domain_withRepeats_all/count_matrix/domains_long_localmax2_b5_d05_p05_CPM-TMM_filter.txt",sep = "\t",quote = F,row.names = T,col.names = T)


# /usr/bin/Rscript scripts/Signature_gene_matrix_v2.R \
# -m ../output/GSE163534/call_domain_withRepeats_all/count_matrix/domains_long_localmax2_b5_d05_p05_CPM-TMM_filter.txt \
# -s ./data/GSE163534/sample_table_filter.csv \
# -o ../output/GSE163534/call_domain_withRepeats_all/too \
# -n 50 -ts 0.5 -e 1 


# /usr/bin/Rscript ~/projects/multi-omics-explore/scripts/run-NormCountMat.R \
# -i ../output/GSE71008/call_domain_withRepeats_all/count_matrix/domains_long_localmax2_b5_d05_p01.txt \
# -o ../output/GSE71008/call_domain_withRepeats_all/count_matrix/domains_long_localmax2_b5_d05_p01_CPM-TMM.txt




# UMAP in R (source: wechat 降维你还在用多种R包吗？快来试试这个宝藏R包吧！)
#remotes::install_github("jlmelville/uwot")
#install.packages("tidydr")
library(tidydr)#降维
library(ggplot2)#可视化
#install.packages("mlr3")
#library(mlr3)

#tidydr::available_methods()
#x <- tidydr::dr(data = t(mat), fun = Rtsne::Rtsne,perplexity = 10)# uwot::tumap, uwot::umap, uwot::lvish, Rtsne::Rtsne, stats::prcomp
x <- tidydr::dr(data = t(mat), fun = Rtsne::Rtsne,perplexity = 5)# uwot::tumap, uwot::umap, uwot::lvish, Rtsne::Rtsne, stats::prcomp
#stats::prcomp()
autoplot(x, aes(fill=group),alpha=0.9, metadata = sample_tabl[, "group", drop=FALSE] ) +
  geom_point(size=5,alpha=0.99,shape=21,color="white")+
  theme_dr()+
  scale_fill_manual(name="Group",values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(unique(sample_tabl$group))))+ 
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
ggsave(filename = paste0(outdir,"/","tSNE_sig.pdf"), width =  9, height = 7)
#dev.off()
#x <- dr(data = iris[,1:4], fun = prcomp)
#x$drdata
# dim(mat)




# DNase -  CPM cor (holding) --------------------------------------------------------
setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/")
options(stringsAsFactors = F)

for (i in c("colon","liver","stomach")){
  #i <- "colon"
  fs <- Sys.glob(paste0("/BioII/lulab_b/baopengfei/shared_data/ENCODE/DNaseI/",i,"/*.bigWig"))
  for (j in 1:length(fs)){
    #  j <- 1
    x <- fs[j]
    #x <- "/BioII/lulab_b/baopengfei/shared_data/ENCODE/DNaseI/colon/ENCFF613BDA.bigWig"
    print(x)
    bw <- rtracklayer::import.bw(x)
    #summary(bw$score)
    bw.df <- as.data.frame(bw)
    bw.df$start <- bw.df$start-1
    bw.df$width <- "X"
    bw.df <- bw.df[,c("seqnames","start","end","width","score","strand")]
    #bw.df[1:3,]
    bw.df <- bw.df[bw.df$score>0,]
    data.table::fwrite(bw.df,paste0(x,".bed"),row.names = F,col.names = F,quote = F,sep = '\t') 
  }
}


#convert tx domain to gn domain
i="output/GSE163534/call_domain_withRepeats_all/domains_localmax_EM/b5_d05_p05.bed"
{{
  grep -v '^chr' $i | exSeek-dev/bin/tbed2gbed <(cat exSeek-dev/genome/hg38/bed/{long_DNA,long_RNA,tRNA,pri_miRNA,piRNA}_newTxID.bed ) /dev/stdin /dev/stdout
  awk 'BEGIN{{OFS="\t";FS="\t"}}/^chr/{{print $1,$2,$3,$4,$5,$6,0,0,0,1,$3-$2,0}}' $i
}} | bedtools sort \
> output/GSE163534/call_domain_withRepeats_all/domains_localmax_EM/b5_d05_p05_gn.bed

bedtools slop -s -b 10000 \
  -g /BioII/lulab_b/baopengfei/shared_reference/hg38/genome.chr.sizes \
  -i output/GSE163534/call_domain_withRepeats_all/domains_localmax_EM/b5_d05_p05_gn.bed \
  > output/GSE163534/call_domain_withRepeats_all/domains_localmax_EM/b5_d05_p05_gn_slop10k.bed

#library(rhdf5)
for i in colon liver stomach
do
cat /BioII/lulab_b/baopengfei/shared_data/ENCODE/DNaseI/colon/ENCFF305IPH.bigWig.bed | bedtools sort -i stdin > /BioII/lulab_b/baopengfei/shared_data/ENCODE/DNaseI/colon/ENCFF305IPH.bigWig.bed.sort
bedtools merge -c 5 -o mean -i /BioII/lulab_b/baopengfei/shared_data/ENCODE/DNaseI/$i/tmp.sort > /BioII/lulab_b/baopengfei/shared_data/ENCODE/DNaseI/$i/bigWig.mergebed

bedtools map \
  -a output/GSE163534/call_domain_withRepeats_all/domains_localmax_EM/b5_d05_p05_gn_slop10k.bed
  -b /BioII/lulab_b/baopengfei/shared_data/ENCODE/DNaseI/$i/mergebed \
  -c 5,5 -o sum,mean > /BioII/lulab_b/baopengfei/shared_data/ENCODE/DNaseI/$i/mergebed
done
