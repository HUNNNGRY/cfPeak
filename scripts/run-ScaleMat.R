
# last 220218 by bpf
# b.p.f@qq.com

options(stringsAsFactors = F)
suppressPackageStartupMessages(library(argparse))
parser <- ArgumentParser(description='scale/transform normalized matrix')
parser$add_argument('-i', '--matrix', type='character', required=TRUE,
                    help='input CPM/TPM matrix. Rows are genes. Columns are samples.')
parser$add_argument('-s', '--sample', type='character', default="all",
                    help='tsv/txt file: pass filter colnames/samples, (default: keep all colnames/samples)')
#parser$add_argument('--byCtrl', type='logical', default=TRUE,
#                    help='whether to scale by ctrl samples (default=TRUE)')
parser$add_argument('--ctrlSample', type='character', default="normal",
                    help='control samples mean and sd value used for scaling. (default: normal, i.e. colnames or samples with HC|HD|NC|Normal|Healthy|normal|healthy character), set to other char if not needed')
parser$add_argument('--omitNAzeroSD', type='character', default="Y",
                    help='whether or not to remove rows with NA and zero std. (default: Y)')
parser$add_argument('-m', '--method', type='character', default="standard",
                    help='scaling method: standard, pareto, quantile|rank. (default: standard)')
parser$add_argument('-o', '--outfile', type='character', required=TRUE,
                    help='output file path')
args <- parser$parse_args()

mat <- args$matrix
sample <- args$sample
ctrl.sample <- args$ctrlSample
omitNA.zeroSD <- args$omitNAzeroSD
method <- args$method
outfile <- args$outfile

# # test
# mat <- "/BioII/lulab_b/baopengfei/projects/exOmics/DIP-seq/output/lulab/matrix/CPM_matrix_gene.txt"
# sample <- "/BioII/lulab_b/baopengfei/projects/multi-omics-explore/meta/lulab/paired-medip-passQCloose.txt"
# ctrl.sample <- "normal"
# omitNA.zeroSD <- "Y"
# method <- "standard"
# outfile <- "test_scaling.txt"

mat <- read.table(mat,check.names = F,sep = "\t",header = T,stringsAsFactors = F,row.names = 1)

if(sample=="all"){
  message(paste0("keep all samples from input mat"))
  sample <- colnames(mat)
} else {
  message(paste0("only keep samples from ",sample))
  sample <- read.delim(sample,header = F,stringsAsFactors = F)$V1
}


## def scaling func.
mat.scaling <- function(mat, samples, ctrl.sample, omitNA.zeroSD, method){
  #mat <- "/BioII/lulab_b/baopengfei/projects/exOmics/DIP-seq/output/lulab/matrix/CPM_matrix_gene.txt"
  #sample <- "/BioII/lulab_b/baopengfei/projects/multi-omics-explore/meta/lulab/paired-medip-passQCloose.txt"
  #sample.nc <- "/BioII/lulab_b/baopengfei/projects/multi-omics-explore/meta/lulab/multiomic-medip-NC-passQCloose.txt"
  #message(mat)
  #samples <- read.table(sample,check.names = F,sep = "\t",header = F,stringsAsFactors = F)$V1
  norm.idx <- grep("HC|HD|NC|Normal|Healthy|normal|healthy",samples)
  #samples[norm.idx]
  #samples.nc <- read.table(sample.nc,check.names = F,sep = "\t",header = F,stringsAsFactors = F)$V1
  
  tmp <- mat
  message(paste0("rows: ",dim(tmp)[1],", cols: ",dim(tmp)))

  
  ## filter samples passQC
  tmp <- tmp[,samples]
  
  ## omit NA and zeroSD rows/features
  message(paste0("omit NA and zeroSD rows/features: ",omitNA.zeroSD))
  if(omitNA.zeroSD=="Y" | omitNA.zeroSD=="y"){
    tmp <- na.omit(tmp)
    l <- apply(tmp,1,sd) #
    tmp <- tmp[l>0,] #
  }
  message(paste0("rows: ",dim(tmp)[1],", cols: ",dim(tmp)))
  
  ## unify sample ID (only for cf mulomics data)
  colnames(tmp) <- gsub("-wgs|-pico|-me","",colnames(tmp))
  
  message(paste0("normalize by method: ",method))
  message(paste0("reference samples: ",ctrl.sample))
  if(ctrl.sample=="normal" & length(norm.idx)>1){
    ## normalize by NC samples
    if(method=="pareto"){
      tmp1 <- apply(tmp, 1, function(y) (y - mean(y[norm.idx])) / (sd(y[norm.idx])^0.5) ^ as.logical(sd(y[norm.idx])))  # 主要需要括号 ！！！
    }else if(method=="standard"){
      tmp1 <- apply(tmp, 1, function(y) (y - mean(y[norm.idx])) / sd(y[norm.idx]) ^ as.logical(sd(y[norm.idx])))
    }else if(method=="quantile" | method=="rank"){
      ## seems not need norm.idx
      tmp1 <- apply(tmp, 1, function(y) dplyr::percent_rank(y)) # min_rank, dense_rank, percent_rank, cume_dist  
    }
  }else{
    ## normalize by all samples
    if(method=="pareto"){
      tmp1 <- apply(tmp, 1, function(y) (y - mean(y)) / (sd(y)^0.5) ^ as.logical(sd(y)))  # 主要需要括号 ！！！
    }else if(method=="standard"){
      tmp1 <- apply(tmp, 1, function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y)))
    }else if(method=="quantile" | method=="rank"){
      ## seems not need norm.idx
      tmp1 <- apply(tmp, 1, function(y) dplyr::percent_rank(y)) # min_rank, dense_rank, percent_rank, cume_dist  
    }
  }
  
  tmp1 <- as.data.frame(t(tmp1))
  return(tmp1)
}

## run scaling
out <- mat.scaling(mat = mat, samples = sample, ctrl.sample = ctrl.sample, omitNA.zeroSD = omitNA.zeroSD, method = method)
out <- as.data.frame(cbind(feature=rownames(out),out))
write.table(x = out, file = outfile, sep = "\t", quote = F, row.names = F, col.names = T)


