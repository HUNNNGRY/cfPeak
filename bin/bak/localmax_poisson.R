#! /usr/bin/env Rscript
# add poisson p stats for localmax results
# last 2204 by bpf 
# b.p.f@qq.com
# transmit form local to bioii at 220604

suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser(description='Filter domains with significant coverage')
parser$add_argument('-i', '--inputFile', type='character', default='-',
                    help='input BED file with covearge in the 5th column, default: -')
parser$add_argument('-b','--rawReadBed', type='character', 
                    help='rawread bed file for calculating count, canbe output from bedtools bamtobed, required')
parser$add_argument('-o', '--outputFile', type='character', default='-',
                    help='output BED file with adjusted p-value in the last column, default: -')
# parser$add_argument('-z', '--zeroTruncDistribution', type='character', default='negbin',
#                     help='model distribution, dist = c("poisson", "negbin", "geometric"), default: negbin')
parser$add_argument('-p', '--pvalue', type='character', default="0.01",
                    help='adjusted p-value threshold for defining significant domains, default: 0.01')
parser$add_argument('-m', '--method', type='character', default="max",
                    help='combine strategy of global tx and local lambda. options: max, min. default: max')
# parser$add_argument('-b', '--background', type='double', default=0.99,
#                     help='set coverage quantile below this quantile as background regions')
parser$add_argument('--bedtoolsPath', type='character', default="/BioII/lulab_b/baopengfei/biosoft",
                    help='bedtools path, default: /BioII/lulab_b/baopengfei/biosoft')
parser$add_argument('--chrSize', type='character', default="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/transcriptome_genome_sort_uniq_newTxID",
                    help='chrSize file path, used for shuffle random regions, default: /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/transcriptome_genome_sort_uniq_newTxID')
# parser$add_argument('--readLength', type="integer", default=20,
#                     help='average read length, default: mean(readlength)')
# parser$add_argument('--minFitRegionNum', type="integer", default=100,
#                     help='number of regions with count>0 used for fit model, default: 100')
parser$add_argument('--startSeed', type="integer", default=1234,
                    help='start seed for shuffle random regions with the same length, default: 1234')
args <- parser$parse_args()


for(i in 1:length(args)){
  v.name <- names(args)[i]
  assign(v.name,args[[i]])
  message(paste0(v.name,": ",get(v.name)))
}

pvalue <- as.numeric(pvalue)
# myfunc <- function(v1) {
#   deparse(substitute(v1))
# }
# myfunc(bedtools_path)
# get(bedtools_path)


# # test
# # #/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/FTC_long/call_domain/tbed_long_RNA/elFTA-17_L3.bed.gz
# #setwd("/Users/baopengfei/Desktop/lulab/tmp/")
# pvalue <- 0.01 # 0.000001
# inputFile <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008/call_domain_withRepeats_all/domains_localmax_by_sample_EM2/b5_d05_p01/SRR2105127.bed"
# outputFile <- "./test.bed"
# bedtoolsPath <- "/BioII/lulab_b/baopengfei/biosoft"
# # bedtoolsPath <- "/Users/baopengfei/anaconda3/bin"
# chrSize <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/transcriptome_genome_sort_uniq_newTxID"
# startSeed <- 1234
# rawReadBed <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008/call_domain_withRepeats_all/tbed_long_RNA_EM/SRR2105127.bed.gz"



# bin/localmax_poisson.R             
# --pvalue 0.1 
# --startSeed 1234 
# --bedtoolsPath /BioII/lulab_b/baopengfei/biosoft 
# --chrSize genome/hg38/chrom_sizes/all_transcript_id             
# --inputFile /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/FTC_small/call_domain/domains_localmax_by_sample/FTC-18_1.bed             
# --rawReadBed /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/FTC_small/call_domain/tbed_long_RNA/FTC-18_1.bed.gz             
# --outputFile /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/FTC_small/call_domain/domains_localmax_significant/FTC-18_1.bed


#setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/")
options(stringsAsFactors = F)
#suppressPackageStartupMessages(library(countreg))
options(bedtools.path = bedtoolsPath) # /BioII/lulab_b/baopengfei/anaconda3/envs/exvariance/bin
suppressPackageStartupMessages(library(bedtoolsr))
suppressPackageStartupMessages(library(stats))
#a,b need bed format df for write.table to temp file, for GRanges, need conver to bed.df first
#bedtoolsr was built with bedtools version 2.30.0



#ppois(20, lambda = 15, lower.tail = FALSE) # 0.08297091 or 8.3%
#1 - ppois(20, lambda = 15)        # Equivalent
#1 - sum(dpois(0:20, lambda = 15)) # Equivalent
# ppois(1:3,lambda = 5:7)
# ppois(3,lambda = 7)


GRanges2bed <- function(gr=ds){
  #1 based to 0 based
  out <- data.frame(chr=gr@seqnames, start=gr@ranges@start-1, end=gr@ranges@start+gr@ranges@width-1, name=gr$name, score=gr$score, strand=gr@strand )
  out$chr <- as.character(out$chr)
  out$name <- as.character(out$name)
  out$strand <- as.character(out$strand)
  out <- as.data.frame(out)
  return(out)
}


# run --------------------------------------------------------------------
## read peak files
if(inputFile == '-'){
  # peak <- rtracklayer::import.bed(file('stdin'))
  #bed <- read.table(file('stdin'), sep='\t', header=FALSE)
  peak.tmp <- data.table::fread(file('stdin'), sep='\t', header=FALSE, stringsAsFactors = F, data.table = F)
} else {
  #peak <- rtracklayer::import.bed(as.character(inputFile))
  peak.tmp <- data.table::fread(as.character(inputFile), sep='\t', header=FALSE, stringsAsFactors = F, data.table = F)
}
peak.tmp <- peak.tmp[,1:7]
colnames(peak.tmp) <- c("chr","start","end","name","score","strand","maxpos") # score: max cov, maxpos: max pos
if("-" %in% as.character(peak.tmp$strand)){
  message("error: - strand detected in bed input file")
}
#score is localmax coverage/count, not the same with bedtools coverage count below calculated
message(paste0('peak number: ',nrow(peak.tmp)))

# peak.df <- GRanges2bed(peak)
# peak.tmp <- peak.df#[peak$bin.idx==i,]
# plot(x=peak.tmp$maxpos-peak.tmp$start, y=peak.tmp$end-peak.tmp$maxpos,alpha=0.1)
# cor.test(peak.tmp$maxpos-peak.tmp$start, peak.tmp$end-peak.tmp$maxpos)
# hist((peak.tmp$maxpos-peak.tmp$start)/(peak.tmp$end-peak.tmp$start))
# table(peak.tmp$maxpos>=peak.tmp$start)
# table(peak.tmp$maxpos<=peak.tmp$end)



#b <- rtracklayer::import.bed(as.character(rawReadBed))
b <- data.table::fread(as.character(rawReadBed), sep='\t', header=FALSE, stringsAsFactors = F, data.table = F)
colnames(b) <- c("chr","start","end","name","score","strand")
#b <- GRanges2bed(b)
message(paste0('rawread number: ',nrow(b)))
b.wd <- mean(b$end-b$start, na.rm = T)
b.wd <- as.integer(b.wd)
message(paste0('rawread median length: ',b.wd))


#message(paste0('NA score number: ',sum(is.na(peak.tmp$score))))
#peak.tmp <- peak.tmp[!is.na(peak.tmp$score),]
#table(tmp$bin.idx)
#rand.num <- nrow(peak.tmp) # max(100000,nrow(peak.tmp)) # 5000000
#peak.tmp[1:3,]
## cal peak cov, and rm potential 0 cov peak records
peak.cov <- bedtoolsr::bt.coverage(s = T, counts = T, a = peak.tmp, b = b) # sorted = T, g = chr.size seems not working
l <- peak.cov[[ncol(peak.cov)]]>=1
message(paste0('peak cov >=1: ',sum(l),"/",nrow(peak.cov)))
#res.stat.df[i,"peak"] <- nrow(peak.cov)
#res.stat.df[i,"peak.cov1"] <- sum(l)

if (sum(l) < nrow(peak.cov)){
  peak.cov <- peak.cov[l,]
  peak.tmp <- peak.tmp[l,]
  message(paste0('peak cov >=1: ',sum(l),"/",nrow(peak.cov)))
}


#chrSize.df[1:3,]
chrSize.df <- data.table::fread(chrSize,data.table = F,sep = "\t",stringsAsFactors = F,header = F)
chrSize.df <- chrSize.df[chrSize.df$V1 %in% peak.tmp$chr,]
tx.length <- chrSize.df$V2[match(peak.tmp$chr,chrSize.df$V1)]
domain.length <- peak.tmp$end-peak.tmp$start


## get tx full length intervals & count by bedtools  
seed <- startSeed
tx.tmp <- peak.tmp
tx.tmp$start <- 0
tx.tmp$end <- tx.length 
tx.tmp <- tx.tmp[!duplicated(tx.tmp$chr),] # remove dup chr/tx
tx.tmp[1:3,]
#tx.df <- bedtoolsr::bt.shuffle(seed = seed, noOverlapping=T, i = tx.tmp, g = chrSize)  # noOverlapping = F might be more reasonable for peaks can be everywhere
tx.cov <- bedtoolsr::bt.coverage(s = T, counts = T, a = tx.tmp, b = b)


r <- domain.length/(tx.length-b.wd)
r[r>1] <- 1
r[r<0] <- 0
#hist(r, breaks = 1000)
#summary(r)
#wide <- r>=0.8 # & tx.length*r <= 10
#peak.tmp[1:3,]


## get local length intervals & count by bedtools  
#bedtools slop -i {params.tmp1} -g {params.chrsize} -b 5 -pct
peak.tmp.small <- peak.tmp  # r>0.8
tx.length.small <- tx.length # [r>0.8]
left <- peak.tmp.small$maxpos-peak.tmp.small$start
right <- peak.tmp.small$end-peak.tmp.small$maxpos
left.ratio <- left/(left+right)
#hist(r1)
# ceiling(0.4)

peak.tmp.small$start <- peak.tmp.small$start + as.integer(left.ratio*0.5*r*tx.length.small) # only consider + strand
peak.tmp.small$end <- peak.tmp.small$start + as.integer(0.5*r*tx.length.small)+1 # only consider + strand, +1 to aviod zero ratio
#local1.tmp <- bedtoolsr::bt.slop(g = chrSize.df, i = peak.tmp.small, b = -0.25, s = T, pct = T)
#colnames(local1.tmp) <- c("chr","start","end","name","score","strand")
local.cov.small <- bedtoolsr::bt.coverage(s = T, counts = T, a = peak.tmp.small, b = b)
#local.cov.small[1:3,]

# table(sum(r<0.2) >= 1)
# fc <- 2 # 2 fold extend both flank
# if (sum(r<0.2) >= 1){
# peak.tmp.long <- peak.tmp[r<0.2,]  # 
# tx.length.long <- tx.length[r<0.2]
# r.long <- r[r<0.2]
# left <- peak.tmp.long$maxpos-peak.tmp.long$start
# right <- peak.tmp.long$end-peak.tmp.long$maxpos
# left.ratio <- left/(left+right)
# peak.tmp.long$start <- peak.tmp.long$start - as.integer(left.ratio * fc * r.long * tx.length.long) # only consider + strand, expand 2 folds domain len left direction
# peak.tmp.long$start[peak.tmp.long$start<0] <- 0
# peak.tmp.long$end <- peak.tmp.long$start + as.integer((1+fc*2) * r.long * tx.length.long) # only consider + strand: right boundary is determined by left boundary and width
# local.cov.long <- bedtoolsr::bt.coverage(s = T, counts = T, a = peak.tmp.long, b = b)
# # local.cov.long[1:3,]
# }

fc <- 0.5 # 2 fold extend both flank (center pos, not maxpos)
fc.ratio <- (1/(fc*2+1))
if (sum(r<fc.ratio) >= 1){
  peak.tmp.long <- peak.tmp[r<fc.ratio,]  # 
  tx.length.long <- tx.length[r<fc.ratio]
  r.long <- r[r<fc.ratio]
  # left <- peak.tmp.long$maxpos-peak.tmp.long$start
  # right <- peak.tmp.long$end-peak.tmp.long$maxpos
  # left.ratio <- left/(left+right)
  peak.tmp.long$start <- peak.tmp.long$start - as.integer(fc * r.long * tx.length.long) # only consider + strand, expand 2 folds domain len left direction
  peak.tmp.long$start[peak.tmp.long$start<0] <- 0
  peak.tmp.long$end <- peak.tmp.long$start + as.integer((1+fc*2) * r.long * tx.length.long) # only consider + strand: right boundary is determined by left boundary and width
  local.cov.long <- bedtoolsr::bt.coverage(s = T, counts = T, a = peak.tmp.long, b = b)
  # local.cov.long[1:3,]
}


## fit poisson model
x <- as.integer(peak.cov[[ncol(peak.cov)]])
bed <- peak.tmp
domain.count <- x[match(bed$name,peak.cov$V4)]
bed$count <- domain.count # add domain count info

x2 <- as.integer((tx.cov[[ncol(tx.cov)]]))
tx.count <- x2[match(bed$chr,tx.cov$V1)]
lambda <- tx.count*r
#summary(local.lambda)

x3 <- as.integer((local.cov.small[[ncol(local.cov.small)]]))
domain.count.small <- x3[match(bed$name,local.cov.small$V4)]
#table(duplicated(local.cov.long$V1)) # no NA exist

if (sum(r<fc.ratio) >= 1){
  x4 <- as.integer((local.cov.long[[ncol(local.cov.long)]]))
  domain.count.long <- x4[match(bed$name,local.cov.long$V4)] 
  #sum(is.na(domain.count.long))==sum(r>=0.2) #NA exist for those peak with ratio >= 0.2
  local.lambda <- domain.count.long*(1/(1+2*fc))
  #table(is.na(lambda))
  #plot(x=log10(local.lambda),y=log10(lambda))
    if (method=="max"){
      lambda[(lambda<local.lambda) & !is.na(local.lambda)] <- local.lambda[lambda<local.lambda & !is.na(local.lambda)] # get max lambda
    } else if (method=="min"){
      lambda[(lambda>local.lambda) & !is.na(local.lambda)] <- local.lambda[lambda>local.lambda & !is.na(local.lambda)] # get min lambda
    } else {
      message("wrong method param, need to be max or min !")
    }
  #table(is.na(lambda))
}

# print(fit)
# fitted(fit)
# residuals(fit)
# AIC(fit1,fit2,fit3)
# BIC(fit1,fit2,fit3)
# plot(residuals(fm, type = "deviance") ~ fitted(fm))
#100, 10, 1000, 10
message('fit poisson distribution')
qvalues1 <- ppois(q = domain.count,lambda = 1 + as.integer(lambda), lower.tail = F) # read will not map to chromosome boundary, thus need subtract median reads length to get effective tx length: tx.length-b.wd
qvalues2 <- ppois(q = domain.count.small,lambda = 1 + as.integer(0.5*lambda), lower.tail = F) # read will not map to chromosome boundary, thus need subtract median reads length to get effective tx length: tx.length-b.wd
qvalues <- apply(cbind(qvalues1=qvalues1,qvalues2=qvalues2),1,min) # get min lambda (equals min pvalue) for long full-length ratio domains
#hist(lambda,xlim = c(0,10),breaks = 10000000, col="grey50")

qvalues <- p.adjust(qvalues, method='BH') # adjust pvalue by BH

#hist(qvalues2)
#table(qvalues>=0.01)
# hist(domain.count,breaks = 10000,xlim = c(0,500))
# hist(tx.length,breaks = 1000000,xlim = c(0,50000))
# hist(domain.length,breaks = 1000,xlim = c(0,100))
# hist(qvalues,breaks = 10000,xlim = c(0,0.1))
# table(qvalues<=pvalue)
#pvalue <- 0.0000001
#table(qvalues<=0.01)
bed_sig <- bed[qvalues < pvalue,]
bed_sig <- cbind(bed_sig, qvalues[qvalues < pvalue])
colnames(bed_sig)[ncol(bed_sig)] <- "qvalue"
message(paste0("sig number: ",nrow(bed_sig)))
res.df <- bed_sig # do.call(rbind,res)
#hist(bed$end-bed$start,breaks = 1000,xlim = c(0,100))
#hist(bed_sig$end-bed_sig$start,breaks = 1000,xlim = c(0,100))
#hist(tx.length,breaks = 1000000,xlim = c(0,50000))
#hist(tx.length[qvalues < pvalue],breaks = 1000000,xlim = c(0,50000))
#hist((bed$end-bed$start)/tx.length,breaks = 1000,xlim = c(0,1),ylim = c(0,1000))
#hist((bed_sig$end-bed_sig$start)/tx.length[qvalues < pvalue],breaks = 1000,xlim = c(0,1),ylim = c(0,1000)) # most high adjp are nearly full-length, cause hard to test sig with similar value

## write outfile
message('number of significant bins:', nrow(res.df))
if(outputFile == '-'){
  data.table::fwrite(res.df, row.names=FALSE, col.names=FALSE, sep='\t', quote = F)
} else {
  data.table::fwrite(res.df, outputFile, row.names=FALSE, col.names=FALSE, sep='\t', quote = F)
}
res.stat.df <- bed
res.stat.df$pvalue <- qvalues
res.stat.df$lambda <- lambda+1
colnames(res.stat.df)[5] <- "maxcount"
data.table::fwrite(res.stat.df, paste0(outputFile,".stat"), row.names=T, col.names=T, sep='\t', quote = F)



# evaluate fitness of model distribution ---------------------------------------------
#https://www.r-bloggers.com/2015/01/goodness-of-fit-test-in-r/
# the process of finding the right distribution for a set of data can be broken down into four steps:
# 1.Visualization. plot the histogram of data 
# 2.Guess what distribution would fit to the data the best 
# 3.Use some statistical test for goodness of fit 
# 4.Repeat 2 and 3 if measure of goodness is not satisfactory

# three well-known and widely use goodness of fit tests that also have nice package in R.
# 1.Chi Square test
# 2.Kolmogorov–Smirnov test
# 3.Cramér–von Mises criterion
#H0 = The data is consistent with a specified reference distribution.

# num_of_samples = 1000
# x <- rgamma(num_of_samples, shape = 10, scale = 3)
# x <- x + rnorm(length(x), mean=0, sd = .1)
# p1 <- hist(x,breaks=50, include.lowest=FALSE, right=FALSE)

## chi square test
# set.seed(1234)
# library('zoo')
# breaks_cdf <- pgamma(p1$breaks, shape=10, scale=3)
# null.probs <- rollapply(breaks_cdf, 2, function(x) x[2]-x[1])
# a <- chisq.test(p1$counts, p=null.probs, rescale.p=TRUE, simulate.p.value=TRUE)
# a$p.value

## Cramér–von Mises criterion
# set.seed(1234)
# num_of_samples = 100000
# y <- rgamma(num_of_samples, shape = 10, scale = 3)
# res <- CDFt::CramerVonMisesTwoSamples(x,y)
# p.value = 1/6*exp(-res)
# p.value

## Kolmogorov–Smirnov test
# set.seed(1234)
# num_of_samples = 1000
# x <- rgamma(num_of_samples, shape = 10, scale = 3)
# x <- x + rnorm(length(x), mean=0, sd = .1)
# plot(x)
# 
# num_of_samples = 100000
# y <- rgamma(num_of_samples, shape = 10, scale = 3)
# plot(y)

# result = ks.test(x, "dunif")
# result$p.value
# hist(x,breaks = 10000, xlim = c(0,100))

## countreg: wald test (can not compare between diff models ?)
# summary(fit1)
# lmtest::waldtest(fit, fit2)


# appendix ----------------------------------------------------------------
# chisq.test(): chi-squared test (stats)
# cut: divides the range of data vector into intervals
# cvm.test(): Cramer-von Mises test for normality (nortest)
# ecdf(): computes an empirical cumulative distribution function (stats)
# fitdistr(): Maximum-likelihood fitting of univariate distributions (MASS)
# goodfit(): fits a discrete (count data) distribution for goodness-of-fit tests (vcd)
# hist(): computes a histogram of the given data values (stats)
# ks.test(): Kolmogorov-Sminorv test (stats)
# kurtosis(): returns value of kurtosis (fBasics)
# lillie.test(): Lilliefors test for normality (nortest)
# mle(): estimate parameters by the method of maximum likelihood (stats4)
# pearson.test(): Pearson chi-square test for normality (nortest)
# plot(): generic function for plotting of R objects (stats)
# qqnorm(): produces a normal QQ plot (stats)
# qqline(), qqplot(): produce a QQ plot of two datasets (stats)
# sf.test(): test di Shapiro-Francia per la normalità (nortest)
# shapiro.test():Shapiro-Francia test for normalità (stats)
# skewness(): returns value of skewness (fBasics)
# table(): builds a contingency table (stats)



# ## plot poisson labmda stat
# #dat <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/WSQ_SMARTer_NEB/call_domain_withRepeats_all/domains_localmax_significant_EM/b5_d05_p01/NC_PKU-2392860_1.bed.stat")
# dat <- read.table("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE133684/call_domain_withRepeats_dedupByPos/domains_localmax_significant_EM/b5_d05_p01/NC_PKU-2392860_1.bed.stat")
# dat[1:3,]
# 
# ref <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt",sep = '\t',header = T,check.names = F,stringsAsFactors = F)
# table(ref$transcript_id %in% unique(dat$chr))
# ref <- ref[ref$transcript_id %in% unique(dat$chr),]
# 
# dat$RNA <- ref$transcript_type[match(dat$chr,ref$transcript_id)]
# dat$RNA <- gsub("_for|_rev","",dat$RNA,perl = T)
# table(dat$RNA)
# 
# dat.tbl <- as_tibble(dat) %>%
#   dplyr::group_by(RNA) %>%
#   dplyr::summarise( lambda.median=median(lambda), lambda.mean=mean(lambda) )
# dat.tbl$x <- dat.tbl$lambda.median
# dat.tbl$y <- apply(as.matrix(dat.tbl[,2:3]), 1, FUN = function(x) dpois(x= as.integer(x[1]), lambda = as.integer(x[1] )) ) # dpois(x=dat.tbl$lambda.median, lambda = dat.tbl$lambda.median)
# #dpois(x=2,lambda = 2)
# # plot(dat.tbl$x,dat.tbl$y)
# 
# res <- list()
# #seqs <- 1:10
# count.seqs <- 1:500
# for (l in c(rna,dna) ){
#   #l <- "pri_miRNA"
#   res[[l]] <- dpois(x=count.seqs,lambda=dat.tbl$lambda.median[dat.tbl$RNA==l])
# }
# res.df <- as.data.frame(res)
# res.df[1:3,]
# rna <- c("rRNA","pri_miRNA","lncRNA", "mRNA", "piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA')
# dna <- c("intron","promoter", "enhancer","repeats") #
# colnames(res.df) <- c(rna,dna)
# res.df$count <- count.seqs
# #plot(res.df$`1`)
# res.df2 <- as_tibble(res.df) %>%
#   tidyr::pivot_longer(cols = "rRNA":"repeats", names_to = "class", values_to = "lambda")
# 
# 
# res.df2$class <- factor(res.df2$class,levels = c(rna,dna))
# res.df2$lab <- paste0(res.df2$class,": ",res.df2$lambda)
# 
# res.df2$log10.count <- log10(res.df2$count+1)
# ggplot(res.df2, aes(x=log10.count, y=lambda, color=class, fill=class)) +
#   # geom_bar(stat="identity",position = position_dodge(width=0.7), alpha=0.8, width=0.5)+
#  # geom_histogram(aes(y=..density..), position="identity", alpha=0.5)
#   geom_density(aes(x=log10.count, y=lambda, color=class, fill=class), alpha = 0.3, stat = "identity", position = "identity") +
#   # geom_line()+
#   # ggsci::scale_fill_nejm()+
#   # ggsci::scale_color_nejm()+
#   scale_x_continuous(name = "log10(k)")+ # ,trans = "log10"
#   scale_fill_manual(name="Species",values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(c(rna,dna))))+ # c(RColorBrewer::brewer.pal(8, "Accent"),RColorBrewer::brewer.pal(9, "Set1"))
#   scale_color_manual(name="Species",values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(c(rna,dna))))+ # c(RColorBrewer::brewer.pal(8, "Accent"),RColorBrewer::brewer.pal(9, "Set1"))
#   # scale_x_discrete(label=c("Domain Ratio","MIR Ratio"))+
#   # ggrepel::geom_label_repel(data = dat.tbl,inherit.aes = F, size = 6, max.overlaps = 20,
#   #                           aes(x=x,y=y,label = paste0(RNA,": ",as.integer(lambda.median)) ), # +(0.0001)
#   #                            show.legend = FALSE #  nudge_x=(0.0001),
#   #                           ) +
#   labs(title="Poisson Prob.",x="k", y = "Poisson Prob.")+ # Overlap Ratio
#   # facet_grid(method~.)+
#   # ylim(c(0,0.5))+
#   # xlim(c(0,1))+
#   theme_bw(base_size=12) +
#   theme(#axis.ticks.x=element_blank(),
#     #strip.text.y = element_blank(),
#     #strip.text.x = element_text(face="bold",family="arial",size=20),
#     axis.title.x = element_text(size=20),
#     axis.title.y = element_text(size=20),
#     axis.text.x = element_text(size = 20,hjust = 0.5,vjust = 0.5), # ,angle = 90
#     axis.text.y = element_text(size = 20),
#     plot.title = element_text(size=20),
#     strip.text = element_text(size = 20),
#     legend.position =  "right",# c(0.9,0.6), # 
#     legend.text = element_text(size= 16),
#     legend.title= element_text(size= 16))
# ggsave("possion.pdf",width = 9,height = 6)

