#!/usr/bin/env Rscript
# run PCC between mRNA_log2TPM and sncRNA_log2CPM

# Load necessary libraries
#library(optparse)
library(argparse)
library(parallel)
library(data.table)

# Define command line arguments
parser <- ArgumentParser(description = "Run PCC between mRNA_log2TPM/clinical_Meta_num and sncRNA_log2CPM (sncRNA_log2CPM rows>5000 is not advised), nrow(smpTbl)>=ncol(log2CPM)==ncol(log2TPM). Usage: Rscript run-covarPCC.R --log2CPM mat1.txt --log2TPM mat2.txt --smpTbl sampletable.txt --covar purity -c 10 -m pearson -o Pcc_purity.txt")

# Define command line arguments
parser$add_argument("--log2CPM", type = "character", required = TRUE, help = "Input log2CPM file name, 1st column as gene index")
parser$add_argument("--log2TPM", type = "character", required = TRUE, help = "Input log2TPM/clinical_Meta_num file name, 1st column as gene/clinicalNum index")
parser$add_argument("--smpTbl", type = "character", required = TRUE, help = "Input sample.table file name, 1st column as sample index")
parser$add_argument("--covar", type = "character", default = "Null", help = "Column name in smpTbl as covariate in partial cor. run default/original cor if not provided ")
parser$add_argument("-c", "--cores", type = "integer", default = 10, help = "Number of cores to use [default: 10]")
parser$add_argument("-m", "--method", type = "character", default = "pearson", help = "Correlation method [default: pearson]")
parser$add_argument("-o", "--output", type = "character", required = TRUE, help = "Output file name")

# Parse command line arguments
args <- parser$parse_args()

# Extract arguments
logCPM <- args$log2CPM
logTPM <- args$log2TPM
sample.table <- args$smpTbl
covar <- args$covar
CORES <- args$cores
method <- args$method
output_file <- args$output


# # test code
# logCPM <- "output/TCGA_small7/HNSC_sncRNA_log2CPM_topVar.txt"
# logTPM <- "/BioII/lulab_b/baopengfei/projects/tcga/data/TCGA-HNSC_long/TCGAbiolinks_log2TPM_filter.txt"
# sample.table <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/TCGA_small7/sample_table2.txt"
# covar <- "purity"
# CORES <- 10
# method <- "pearson"
# output_file <- "test.txt"

# logCPM <- "output/TCGA_small7/STAD_sncRNA_log2CPM_topVar_rmPurity.txt"
# logTPM <- "/BioII/lulab_b/baopengfei/projects/tcga/data/STAD_long/log2TPM_filter.txt"
# sample.table <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/TCGA_small7/sample_table2.txt"
# covar <- "purity"
# CORES <- 10
# method <- "pearson"
# output_file <- "test.txt"


# Function to calculate correlation (assuming getCorVec is defined elsewhere)
## covar PCC
fun_mtx_pcr <- function(x,y,z,method="pearson"){ # z: covar
  # x <- x0
  # y <- y0
  # r12=cor((x),(y),method = method)
  # r13=cor((x),z,method = method)
  # r23=cor(z,(y),method = method)
  r12=cor.test((x),(y),method = method)$estimate
  r13=cor.test((x),z,method = method)$estimate
  r23=cor.test(z,(y),method = method)$estimate
  r123=r13%*%r23
  rup=r12-r123
  rd1=sqrt(1-r13*r13)
  rd2=sqrt(1-r23*r23)
  rd=rd1%*%rd2
  rrr=rup/rd
  return(rrr)
}

getCorVec <- function(x0,y0, pcc_type="original", covar="Null",  method="pearson",conf.level=0.95){
  # x0 <- logTPM.l[[1]]
  # y0 <- as.numeric(logCPM[1, ])
  # covar <- as.numeric(sample.table[["purity"]]) # 
  # pcc_type <- "covar"
  # conf.level=0.95
  # method="pearson"
  if(covar=="Null"){
    tmp <- cor.test(x0,y0, method = method, conf.level = conf.level ) # paired; conf.level only used for the Pearson product
    if(method=="pearson" & ("conf.int" %in% names(tmp)) ){
      res <- data.frame(confInt=paste0(round(tmp$conf.int[1],digits = 6),"|",round(tmp$conf.int[2],digits = 6)),pvalue=tmp$p.value,cor=tmp$estimate)
    }else{ # kendall or spearman
      res <- data.frame(confInt=paste0("|"),pvalue=tmp$p.value,cor=tmp$estimate)
      # res1 <- res
    }
  }else if(covar!="Null"){
    RBP_cancer_out0_log2 <- x0 # not matter for assign of x0 and y0
    gene_cancer_other_out0_log2 <- y0
    tumor_inter <- covar # sample.table[[covar]]
    # message(tumor_inter)
    pcor<-fun_mtx_pcr(RBP_cancer_out0_log2,gene_cancer_other_out0_log2,tumor_inter) # [-1,+1] partial correlation coefficient computed by fun_mtx_pcr.
    n <- length(RBP_cancer_out0_log2) # 1st var in fun_mtx_pcr input
    gn <- 1
    statistic<- pcor*sqrt((n-2-gn)/(1-pcor^2)) # used for pval cal, the test statistic calculated to assess the significance of the partial correlation. 
    p.value<- as.numeric(2*pnorm(-abs(statistic)))
    # padj <- p.adjust(p.value, method = 'BH') # vec-wise not feasiable
    RS<- -log10(p.value+0.001)*sign(pcor) # use this as rank: 
    res <- data.frame(RS=RS,pvalue=p.value,p.stat=statistic,cor=pcor)
    # res2 <- res
  }
  return(res)
}




# Main script
logTPM <- data.table::fread(logTPM, sep = "\t", header = T, stringsAsFactors = F, check.names = F, data.table = F)
rownames(logTPM) <- logTPM[,1]
logTPM[,1] <- NULL
logTPM <- logTPM[rowSums(as.matrix(logTPM))>0,]
logTPM <- logTPM[apply(as.matrix(logTPM),1,sd)>0,]

logCPM <- data.table::fread(logCPM, sep = "\t", header = T, stringsAsFactors = F, check.names = F, data.table = F)
rownames(logCPM) <- logCPM[,1]
logCPM[,1] <- NULL
logCPM <- logCPM[rowSums(as.matrix(logCPM))>0,]
logCPM <- logCPM[apply(as.matrix(logCPM),1,sd)>0,]

if(nrow(logCPM) > 5000){
  message("nrow(logCPM) > 5000, consider reduce rownum by filter topDiff or topVar")
}

sample.table <- data.table::fread(sample.table, sep = "\t", header = T, stringsAsFactors = F, check.names = F, data.table = F)
sample.table <- sample.table[!duplicated(sample.table$sample),]
rownames(sample.table) <- sample.table$sample
if(covar=="Null"){
  message("run default cor without covar")
  sample.table[[covar]] <- 1 # no covar
}else if( covar %in% colnames(sample.table) ){
  message("run partial cor, covar: ",covar)
}else{
  message(covar," not exist")
  stop()
}
sample.table <- sample.table[!is.na(sample.table[[covar]]) & sample.table[[covar]]!="" | sample.table[[covar]]!=" " & !is.nan(sample.table[[covar]]), ]
sample.table <- sample.table[intersect(colnames(logCPM), colnames(logTPM)),]
logCPM <- logCPM[,sample.table$sample]
logTPM <- logTPM[,sample.table$sample]
if(sum(grepl("peak|cluster|block",rownames(logCPM)[1],perl = T))==1){
  message("use 4th | as "," peak id")
  tmp <- unlist(sapply(strsplit(rownames(logCPM),"|",fixed = T),"[",4))
  logCPM <- logCPM[!duplicated(tmp),]
  rownames(logCPM) <- tmp[!duplicated(tmp)]
}
if(sum(grepl("ENSG",rownames(logTPM)[1],perl = T))==1){
  message("use 1st | as "," gene id")
  tmp <- unlist(sapply(strsplit(rownames(logTPM),".",fixed = T),"[",1))
  logTPM <- logTPM[!duplicated(tmp),]
  rownames(logTPM) <- tmp[!duplicated(tmp)]
}

logTPM.l <- as.list(as.data.frame(t(logTPM)))
for (g in 1:nrow(logCPM) ) {
  # g <- 1
  idx <- rownames(logCPM)[g]
  if(covar!="Null"){
    tmp <- mclapply(mc.cores = CORES, X = logTPM.l, FUN = getCorVec, 
                    y0 = as.numeric(logCPM[g, ]), method = method, 
                    conf.level = 0.95, #pcc_type = "covar",
                    covar = as.numeric(sample.table[[covar]]) )
  }else{
    tmp <- mclapply(mc.cores = CORES, X = logTPM.l, FUN = getCorVec,  
                    y0 = as.numeric(logCPM[g, ]), method = method,
                    conf.level=0.95, #pcc_type="original",
                    covar="Null")
  }
  #lapply(logTPM.l, FUN = getCorVec, y0 = as.numeric(logCPM[g, ]), method = method, conf.level = 0.95, pcc_type = "covar", covar = as.numeric(sample.table[[covar]]) )
  tmp2 <- as.data.frame(do.call(rbind, tmp))
  #tmp1 <- as.data.frame(do.call(rbind, tmp))
  #cor.test(tmp1$RS,tmp1$p.stat)
  tmp2$padj <- p.adjust(as.numeric(tmp2$pvalue), method = 'BH')
  tmp2$mRNA <- rownames(tmp2)
  tmp2$sncRNA <- idx
  output_file_tmp <- paste0(output_file, "_", idx, ".txt")
  fwrite(tmp2, output_file_tmp, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}


res.list <- list()
for (g in 1:nrow(logCPM) ) {
  # g <- 1
  idx <- rownames(logCPM)[g]
  output_file_tmp <- paste0(output_file, "_", idx, ".txt")
  res.list[[idx]] <- data.table::fread(output_file_tmp, sep = "\t", header = T, stringsAsFactors = F, check.names = F, data.table = F)
}
res.df <- as.data.frame(do.call(rbind, res.list))
idx.all <- names(res.list)
rm(list="res.list")

# Write the output to a file
fwrite(res.df, output_file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

unlink(x = paste0(output_file, "_", idx.all, ".txt") ) #

