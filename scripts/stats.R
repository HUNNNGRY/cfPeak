
# stats function in R




# correlation/MutualInfo functions ------
#library(minerva)
#install.packages("minerva")
#library(energy)
#install.packages("energy")

doubleCenter <- function(x){
  centered <- x
  for(i in 1:dim(x)[1]){
    for(j in 1:dim(x)[2]){
      centered[i,j] <- x[i,j] - mean(x[i,]) - mean(x[,j]) + mean(x)
    }
  }
  return(centered)
}
distanceCovariance <- function(x,y){
  N <- length(x)
  distX <- as.matrix(dist(x))
  distY <- as.matrix(dist(y))
  centeredX <- doubleCenter(distX)
  centeredY <- doubleCenter(distY)
  calc <- sum(centeredX * centeredY)
  return(sqrt(calc/(N^2)))
}
distanceVariance <- function(x){
  return(distanceCovariance(x,x))
}
distanceCorrelation <- function(x,y){
  cov <- distanceCovariance(x,y)
  sd <- sqrt(distanceVariance(x)*distanceVariance(y))
  return(cov/sd)
}
# # Compare with Pearson's r
# x <- -10:10
# y <- x^2 + rnorm(21,0,10)
# cor(x,y) # --> 0.057
# distanceCorrelation(x,y) # --> 0.509

bootstrap <- function(x,y,reps,alpha){
  estimates <- c()
  original <- data.frame(x,y)
  N <- dim(original)[1]
  for(i in 1:reps){
    S <- original[sample(1:N, N, replace = TRUE),]
    estimates <- append(estimates, distanceCorrelation(S$x, S$y))
  }
  u <- alpha/2 ; l <- 1-u
  interval <- quantile(estimates, c(l, u))
  return(2*(dcor(x,y)) - as.numeric(interval[1:2]))
}
# Use with 1000 reps and threshold alpha = 0.05
# x <- -10:10
# y <- x^2 + rnorm(21,0,10)
# bootstrap(x,y,1000,0.05) # --> 0.237 to 0.546

permutationTest <- function(x,y,reps){
  estimates <- c()
  observed <- distanceCorrelation(x,y)
  N <- length(x)
  for(i in 1:reps){
    y_i <- sample(y, length(y), replace = T)
    estimates <- append(estimates, distanceCorrelation(x, y_i))
  }
  p_value <- mean(estimates >= observed)
  return(p_value)
}

entropy <- function(x){
  pr <- prop.table(table(x))
  H <- sum(pr * log(pr,2))
  return(-H)
}

crossEntropy <- function(x,y){
  prX <- prop.table(table(x))
  prY <- prop.table(table(y))
  H <- sum(prX * log(prY,2))
  return(-H)
}

KL_divergence <- function(x,y){
  kl <- crossEntropy(x,y) - entropy(x)
  return(kl)
}

jointDist <- function(x,y){
  N <- length(x)
  u <- unique(append(x,y))
  joint <- c()
  for(i in u){
    for(j in u){
      f <- x[paste0(x,y) == paste0(i,j)]
      joint <- append(joint, length(f)/N)
    }
  }
  return(joint)
}

marginalProduct <- function(x,y){
  N <- length(x)
  u <- unique(append(x,y))
  marginal <- c()
  for(i in u){
    for(j in u){
      fX <- length(x[x == i]) / N
      fY <- length(y[y == j]) / N
      marginal <- append(marginal, fX * fY)
    }
  }
  return(marginal)
}

mutualInfo <- function(x,y){
  joint <- jointDist(x,y)
  marginal <- marginalProduct(x,y)
  Hjm <- - sum(joint[marginal > 0] * log(marginal[marginal > 0],2))
  Hj <- - sum(joint[joint > 0] * log(joint[joint > 0],2))
  return(Hjm - Hj)
}

MIC <- function(x,y){
  N <- length(x)
  maxBins <- ceiling(N ** 0.6)
  MI <- c()
  for(i in 2:maxBins) {
    for (j in 2:maxBins){
      if(i * j > maxBins){
        next
      }
      Xbins <- i; Ybins <- j
      binnedX <-cut(x, breaks=Xbins, labels = 1:Xbins)
      binnedY <-cut(y, breaks=Ybins, labels = 1:Ybins)
      MI_estimate <- mutualInfo(binnedX,binnedY) 
      MI_normalized <- MI_estimate / log(min(Xbins,Ybins),2)
      MI <- append(MI, MI_normalized)
    }
  }
  return(max(MI))
}


#cor() computes the correlation coefficient. cor.test() test for association/correlation between paired samples.
#difference also exist for NA records

# test
# set.seed(1234)
# # Noise
# x0 <- rnorm(100,0,1)
# y0 <- rnorm(100,0,1)
# plot(y0~x0, pch = 18)
# cor(x0,y0,method = c("pearson"))
# cor(x0,y0,method = c("spearman"))
# distanceCorrelation(x0,y0)
# MIC(x0,y0)


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

# getCor <- function(i,sample.table,logcpm, pcc_var, pcc_type="original", covar="Null",  method="spearman",conf.level=0.95){
#   # i <- 1
#   # message(i)
#   ## need sample.table and logcpm
#   # pcc_type <- "original"
#   # method <- "pearson" #"spearman"
#   # conf.level=0.95
#   # covar <- "DateOperator"
#   x0 <- sample.table[[pcc_var]]
#   y0 <- as.numeric(logcpm[i,])
# # length(y0)
#     #as.data.frame(tmp)
#   if(pcc_type=="original"){
#     tmp <- cor.test(x0,y0, method = method, conf.level = conf.level ) # paired; conf.level only used for the Pearson product
#     if(method=="pearson"){
#       res <- data.frame(id=i,confInt=paste0(round(tmp$conf.int[1],digits = 6),"|",round(tmp$conf.int[2],digits = 6)),pvalue=tmp$p.value,cor=tmp$estimate)
#     }else{
#       res <- data.frame(id=i,pval=tmp$p.value,cor=tmp$estimate)
#     }
#   }else if(pcc_type=="covar"){
#     RBP_cancer_out0_log2 <- x0
#     gene_cancer_other_out0_log2 <- y0
#     tumor_inter <- sample.table[[covar]]
#     # message(tumor_inter)
#     pcor<-fun_mtx_pcr(RBP_cancer_out0_log2,gene_cancer_other_out0_log2,tumor_inter)
#     n<-length(RBP_cancer_out0_log2)
#     gn<-1
#     statistic<- pcor*sqrt((n-2-gn)/(1-pcor^2))
#     p.value<- 2*pnorm(-abs(statistic))
#     RS<- -log10(p.value+0.001)*sign(pcor) # use this as rank
#     res <- data.frame(id=i,RS=RS,pvalue=p.value,cor=statistic)
#   }
#   return(res)
# }

#



## DIY
flat2CorrMat <- function(x){
  # x should be uniq combination
  # return upper mat
  # x <- jac[,c("from","to","jaccard")]
  rowN <- unique(x[,2])
  colN <- unique(x[,1])
  comb <- c(rowN,unique(colN[!(colN %in% rowN)])) #unique(c(rowN,colN))
  tmp <- matrix(nrow = length(comb), ncol = length(comb), dimnames = list(comb, comb)  )
  for (i in 1:nrow(x)){
    # message(i)
    tmp[x[i,2],x[i,1]] <- x[i,3]
    tmp[x[i,1],x[i,2]] <- x[i,3]
  }
  return(tmp)
}
#






# significance -------
pvalue2plabel <- function(x){
  #x <- seq(0.0001,1,0.01)
  x2 <- as.numeric(x)
  x2[x>=0.1] <- "ns"
  x2[x<0.1] <- "*"
  x2[x<0.01] <- "**"
  x2[x<0.001] <- "***"
  x2[x<0.0001] <- "****"
  x2[is.na(x) | is.nan(x)] <- "na"
  # table(x2)
  return(x2)
}


metaPval <- function(p_values){
#install.packages("metap")
library(metap)
# Create a vector of p-values from multiple tests
#p_values <- c(0.03, 0.05, 0.08, 0.01)
# Use Fisher's method to combine the p-values
return(metap::sumlog(p_values)$p)
}
#* When to use Fisher's combined probability test or p.adjust?
#* Use Fisher's combined probability test when you want to assess the overall significance of a set of independent tests, without necessarily controlling for multiple comparisons.
#* Use p.adjust when you want to control the family-wise error rate or the false discovery rate when performing multiple independent tests on the same data.


# multiple fisher exact test
fisherExactFromCountTable <- function(contingency_table){
  #contingency_table contains >=2 columns as grps for comparison (e.g.: Tumor subtypes), >=2 rows each as a category label, value is count of samples 
  #iterate each pair of columns pairs, then each rows for a OvR pattern to build 2x2 matrix
  cmb <- combn(x = colnames(contingency_table), m = 2)
  results_df2 <- list()
  for(idx in 1:ncol(cmb)){
    # idx <- 5
    data <- as.data.frame(contingency_table[,c(cmb[1,idx],cmb[2,idx])])
    
    #colnames(data) <- c("grp1","grp2")
    data$color <- rownames(data)
    data <- data[, c(3,1,2)]
    # Calculate totals
    total_grp1 <- sum(data[,2])
    total_grp2 <- sum(data[,3])
    total_all <- total_grp1 + total_grp2
    
    # Initialize a list to store results
    results_list <- list()
    for (i in 1:nrow(data)) {
      # i <- 4
      color_count_grp1 <- data[,2][i]
      color_count_grp2 <- data[,3][i]
      other_count_grp1 <- total_grp1 - color_count_grp1
      other_count_grp2 <- total_grp2 - color_count_grp2
      
      # Create a 2x2 matrix for Fisher's test
      contingency_matrix <- matrix(c(color_count_grp1, other_count_grp1,
                                     color_count_grp2, other_count_grp2),
                                   nrow = 2, byrow = TRUE)
      test_result <- fisher.test(contingency_matrix)
      # print(test_result$p.value)
      # print(fisher.test(t(contingency_matrix))$p.value) # t() not change
      results_list[[data$color[i]]] <- test_result$p.value
    }
    results_df <- data.frame(
      category = names(results_list),
      grp1=cmb[1,idx],grp2=cmb[2,idx],
      p_value = unlist(results_list)
    )
    results_df$adjusted_p_value <- p.adjust(results_df$p_value, method = "BH") # not in filnal step
    results_df2[[idx]] <- results_df
  }
  results_df2 <- as.data.frame(do.call(rbind,results_df2))
  return(results_df2)
}
#


# permutation ----------
permVec <- function(x){
  x <- x[sample(1:length(x),length(x),replace = F)]
  return(x)
}
permDfByCol <- function(df){
  # df <- diff.sig
  
  tmp <- rownames(df)
  # rownames(df) <- rownames(df)[sample(1:nrow(df),nrow(df),replace = F)] # shuf rowname not enough, need shuf each col seperately
  # df <- df[tmp,]
  
  df <- lapply(df,FUN = permVec )
  df <- as.data.frame(do.call(cbind,df))
  rownames(df) <- tmp
  return(df)
}


