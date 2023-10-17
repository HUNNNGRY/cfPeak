# func. for cfPeak
# last 2305 by bpf 
# b.p.f@qq.com

# def func. ---------------------------------------------------------------
old2newTxID <- function(x){
  x <- gsub("(-)","_neg",x,fixed = T)
  x <- gsub("(+)","_pos",x,fixed = T)
  x <- gsub("::","__",x,fixed = T)
  x <- gsub(":","___",x,fixed = T)
  x <- gsub("-","____",x,fixed = T)
  x <- gsub(".","_____",x,fixed = T)
  return(x)
}
new2oldTxID <- function(x){
  # x <- txID
  x <- gsub("_neg","(-)",x,fixed = T)
  x <- gsub("_pos","(+)",x,fixed = T)
  x <- gsub("_____",".",x,fixed = T)
  x <- gsub("____","-",x,fixed = T)
  x <- gsub("___",":",x,fixed = T)
  x <- gsub("__","::",x,fixed = T)
  return(x)
}

GRange2bed <- function(gr){
  # gr <- hg38[[1]] # 1-base to 0-base
  bed <- data.frame("chr"=gr@seqnames,"start"=gr@ranges@start-1,"end"=gr@ranges@start,"name"=gr$name,"score"=gr$score,"strand"=gr@strand)
  return(bed)
}

get.mat <- function(signal,region,signal.label,region.label,up=50,down=50,bin=10,ratio=0.3,mean_mode = "coverage"){
  # signal <- "output/GSE71008_NCpool/call_peak_all/tbed_RNA_EM/NCpool.bed.gz"
  # region <- "output/GSE71008_NCpool/call_peak_all/domains_clipper_by_sample/b5_p05/NCpool.bed"
  # up=50,down=50,bin=10,ratio=0.3
  print(signal)
  print(region)
  # signal.gr <- rtracklayer::import(con=as.character(signal)) # import(.bed) got error if input bed not proper formated
  # region.gr <- rtracklayer::import(con=as.character(region))
  signal.gr <- data.table::fread(as.character(signal))
  signal.gr <- GRanges(signal.gr$V1, IRanges(signal.gr$V2, signal.gr$V3), name=paste(signal.gr$V4,signal.gr$V1,signal.gr$V2,signal.gr$V3,signal.gr$V6,sep="_"), score=signal.gr$V5, strand=signal.gr$V6)
  region.gr <- data.table::fread(as.character(region))
  region.gr <- GRanges(region.gr$V1, IRanges(region.gr$V2, region.gr$V3), name=paste(region.gr$V4,region.gr$V1,region.gr$V2,region.gr$V3,region.gr$V6,sep="_"), score=region.gr$V5, strand=region.gr$V6)
  
  chrs <- intersect(unique(seqnames(signal.gr)), unique(seqnames(region.gr)))
  # length(chrs)
  signal.gr <- signal.gr[seqnames(signal.gr) %in% chrs,]
  region.gr <- region.gr[seqnames(region.gr) %in% chrs,]
  
  
  #strand info not preserved, not flipped
  mat1 = EnrichedHeatmap::normalizeToMatrix(signal.gr, region.gr, 
                                            #value_column = "cov",  
                                            # extend = c(250, 250), # gn:250,250
                                            extend = c(up, down), # tx:50,50
                                            target_ratio = ratio, #0.3 
                                            # k = 10, # bins of target regions
                                            mean_mode = mean_mode, # "coverage",  # c("absolute", "weighted", "w0", "coverage")
                                            w = bin, # bins of extended windows
                                            keep = c(0, 0.99), # Percentiles in the normalized matrix to keep.
                                            background = 0, 
                                            smooth = TRUE # set TRUE may get negative cov 
  )
  mat1 <- as.data.frame(mat1)
  # will normalize whole matrix !
  # signal <- basename(as.character(signal))
  # signal <- unlist(sapply(strsplit(signal,".",fixed = T),"[",1))
  # region <- basename(as.character(region))
  # region <- unlist(sapply(strsplit(region,".",fixed = T),"[",1))
  
  res.tmp <- as.data.frame(cbind(bed=region.gr$name,region=rep(region.label,nrow(mat1)),signal=rep(signal.label,nrow(mat1)),mat1))
  print(nrow(res.tmp))
  print(ncol(res.tmp))
  return(res.tmp)
}

maxmin.normalize <- function(x) {  
  #x <- sweep(x, 2, apply(x, 2, min)) 
  #x <- sweep(x, 2, apply(x, 2, max), "/")  # cannot handle zero sd
  x <- apply(x, 2, function(y) (y - min(y))/(max(y)-min(y))^as.logical(sd(y)) ) # can handle zero sd, changed as 230130
  x  # (0,1)
  #2*x - 1  # (-1,1)
}
maxmin.normalize.vec <- function(y) {  
  #x <- sweep(x, 2, apply(x, 2, min)) 
  #x <- sweep(x, 2, apply(x, 2, max), "/")  # cannot handle zero sd
  x <- (y - min(y))/(max(y)-min(y))^as.logical(sd(y)) # can handle zero sd, changed as 230130
  x  # (0,1)
  #2*x - 1  # (-1,1)
}

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


## tpm
tpm <- function(count,gene.len){
  mat <- as.matrix(count)
  matrix_tpm <- 1000*mat / gene.len
  matrix_tpm <- t(t(matrix_tpm) * 1e6 / colSums(matrix_tpm))
  return(matrix_tpm)
}


fpkmToTpm <- function(fpkm)
{
  # fpkm[1:3,1:3]
  # TCGA-DM-A288-01A TCGA-QL-A97D-01A TCGA-CM-6164-01A
  # 1:        0.0000000           0.0000           0.0000
  # 2:        0.0051099           0.0000           0.0000
  # 3:        2.7694750           1.9955           3.0041
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

## def limma diff func.
#with weight
limma.trend <- function(logcpm,group){
  #logcpm = mat.norm.tmp
  #group = sample.table.tmp$source 
  suppressPackageStartupMessages(library(limma))
  suppressPackageStartupMessages(library(edgeR))
  #y <- DGEList(counts=mat, samples=samples, group=group)
  #y <- calcNormFactors(y, method=norm_method)
  model <- model.matrix(~group)
  #y <- voom(y, model, plot=FALSE)
  #fit <- lmFit(y, model)
  
  ## add weight
  we <- limma::arrayWeights(logcpm, design = model, method = "genebygene")
  fit <- lmFit(logcpm, model, weights = we)
  fit <- eBayes(fit, robust=TRUE, trend=TRUE) # , trend=TRUE (limma.trend), voom比limma-trend更适用于样本库大小不一的情况
  #fit2 <- contrasts.ft(fit)
  #fit2 <- eBayes(fit2, robust=TRUE, trend=TRUE)
  #top_table <- topTable(fit2, sort.by='none', n=Inf)
  top_table <- topTable(fit, coef=2, sort.by='none', n=Inf)
  # rename columns
  mapped_names <- colnames(top_table)
  for(i in 1:ncol(top_table)){
    if(colnames(top_table)[i] == 'logFC'){
      mapped_names[i] <- 'log2FoldChange'
    }else if(colnames(top_table)[i] == 'P.Value'){
      mapped_names[i] <- 'pvalue'
    }else if(colnames(top_table)[i] == 'adj.P.Val') {
      mapped_names[i] <- 'padj'
    }else{
      mapped_names[i] <- colnames(top_table)[i]
    }
  }
  colnames(top_table) <- mapped_names
  res <- top_table
  return(res)
}

## limma.trend.paired deprecated !!!
limma.trend.paired <- function(logcpm,group,patient){
  #logcpm <- mat.norm.tmp
  #group <- sample.table.tmp$source
  #patient <- sample.table.tmp$cell_id
  suppressPackageStartupMessages(library(limma))
  suppressPackageStartupMessages(library(edgeR))
  #y <- DGEList(counts=mat, samples=samples, group=group)
  #y <- calcNormFactors(y, method=norm_method)
  model <- model.matrix(~patient+group)
  #y <- voom(y, model, plot=FALSE)
  #fit <- lmFit(y, model)
  fit <- lmFit(logcpm, model)
  fit <- eBayes(fit, robust=TRUE, trend=TRUE) # , trend=TRUE (limma.trend), voom比limma-trend更适用于样本库大小不一的情况
  #fit2 <- contrasts.ft(fit)
  #fit2 <- eBayes(fit2, robust=TRUE, trend=TRUE)
  #top_table <- topTable(fit2, sort.by='none', n=Inf)
  top_table <- topTable(fit, coef=ncol(model), sort.by='none', n=Inf)
  # rename columns
  mapped_names <- colnames(top_table)
  for(i in 1:ncol(top_table)){
    if(colnames(top_table)[i] == 'logFC'){
      mapped_names[i] <- 'log2FoldChange'
    }else if(colnames(top_table)[i] == 'P.Value'){
      mapped_names[i] <- 'pvalue'
    }else if(colnames(top_table)[i] == 'adj.P.Val') {
      mapped_names[i] <- 'padj'
    }else{
      mapped_names[i] <- colnames(top_table)[i]
    }
  }
  colnames(top_table) <- mapped_names
  res <- top_table
  return(res)
}


## define deseq2 diff func.
#paired patient
diff <- function(mat,samples,group,patient,method,norm_method, filterType="small", featureType){
  suppressPackageStartupMessages(library(edgeR))
  mat <- mat[,samples]
  print(dim(mat))
  y <- DGEList(counts=mat, samples=samples, group=group)
  
  ## filter low expr
  message(filterType)
  counts <- edgeR::getCounts(y)
  
  if (filterType=="small"){
    keep <- rowSums(counts>=1) >= 0.5*length(samples)
  }else if(filterType=="long"){
    if(featureType=="gene"){
      gene.len = as.numeric(unlist(lapply(strsplit(rownames(counts), "|", fixed = T), function(x) x[2])))
    }else if (featureType=="domain"){
      gene.len = as.numeric(unlist(lapply(strsplit(rownames(counts), "|", fixed = T), function(x) x[7]))) - as.numeric(unlist(lapply(strsplit(rownames(counts), "|", fixed = T), function(x) x[6])))
    }
    tpm_mat <- tpm(count = counts, gene.len = gene.len)
    keep <- rowSums(tpm_mat >= 1) >= 0.5*length(samples)
    
    tpm_mat <- tpm_mat[keep,]
    logtpm <- log2(tpm_mat+1)
    print(dim(logtpm))
    print(logtpm[1:2,1:2])
  }
  #keep <- filterByExpr(y,min.count = 1, min.total.count = 5, min.prop = 0.1) 
  #min.count: Minimum count required for at least some samples.
  #min.total.count: Minimum total count required.
  #min.prop: Minimum proportion of samples in the smallest group that express the gene.
  y <- y[keep,,keep.lib.sizes=FALSE]
  print(dim(y))
  
  y <- calcNormFactors(y, method=norm_method)
  
  logcpm <- edgeR::cpm(y, log=F, normalized.lib.sizes = T,prior.count = 0)
  logcpm <- log2(logcpm+1)
  
  design <- model.matrix(~ factor(patient) + factor(group))  # partial paired mode
  #design <- model.matrix(~  group + patient)  # partial paired mode, design order has great impact
  #design <- model.matrix(~ group) 
  y <- estimateDisp(y, design) 
  if(method == 'edger_glmqlf'){
    fit <- glmQLFit(y, design)
    test <- glmQLFTest(fit, coef=ncol(design))   
    res <- topTags(test, n=nrow(mat), sort.by='none')
  } else if(method == 'edger_glmlrt'){
    fit <- glmFit(y, design)
    test <- glmLRT(fit, coef=ncol(design))  # coef should be the pos col
    res <- topTags(test, n=nrow(mat), sort.by='none')
  } else if(method == 'edger_exact'){
    test <- exactTest(y)
    res <- topTags(test, n=nrow(mat), sort.by='none')
  }
  res <- cbind(res$table, baseMean=2^(res$table$logCPM))
  # rename columns
  mapped_names <- colnames(res)
  for(i in 1:ncol(res)){
    if(colnames(res)[i] == 'logFC'){
      mapped_names[i] <- 'log2FoldChange'
    }else if(colnames(res)[i] == 'PValue'){
      mapped_names[i] <- 'pvalue'
    }else if(colnames(res)[i] == 'FDR') {
      mapped_names[i] <- 'padj'
    }else{
      mapped_names[i] <- colnames(res)[i]
    }
  }
  colnames(res) <- mapped_names
  
  if (filterType=="small"){
    normMat <- logcpm
  }else if(filterType=="long"){
    normMat <- logtpm
  }
  
  outfile <- list()
  outfile[["normMat"]] <- normMat
  outfile[["diffTable"]] <- res
  return(outfile)
}

#not paired
diff.v2 <- function(mat,samples,group,method,norm_method, filterType="small",featureType){
  suppressPackageStartupMessages(library(edgeR))
  mat <- mat[,samples]
  print(dim(mat))
  y <- DGEList(counts=mat, samples=samples, group=group)
  
  ## filter low expr
  #keep <- filterByExpr(y,min.count = 1, min.total.count = 5, min.prop = 0.1)
  message(filterType)
  counts <- edgeR::getCounts(y)
  if (filterType=="small"){
    keep <- rowSums(counts>=1) >= 0.5*length(samples)
  }else if(filterType=="long"){
    if(featureType=="gene"){
      gene.len = as.numeric(unlist(lapply(strsplit(rownames(counts), "|", fixed = T), function(x) x[2])))
    }else if (featureType=="domain"){
      gene.len = as.numeric(unlist(lapply(strsplit(rownames(counts), "|", fixed = T), function(x) x[7]))) - as.numeric(unlist(lapply(strsplit(rownames(counts), "|", fixed = T), function(x) x[6])))
    }
    tpm_mat <- tpm(count = counts, gene.len = gene.len)
    keep <- rowSums(tpm_mat >= 1) >= 0.5*length(samples)
    
    tpm_mat <- tpm_mat[keep,]
    logtpm <- log2(tpm_mat+1)
    print(dim(logtpm))
    print(logtpm[1:2,1:2])
  }
  #min.count: Minimum count required for at least some samples.
  #min.total.count: Minimum total count required.
  #min.prop: Minimum proportion of samples in the smallest group that express the gene.
  y <- y[keep,,keep.lib.sizes=FALSE]
  print(dim(y))
  
  y <- calcNormFactors(y, method=norm_method)
  
  logcpm <- edgeR::cpm(y, log=F, normalized.lib.sizes = T,prior.count = 0)
  logcpm <- log2(logcpm+1)
  
  #design <- model.matrix(~ patient + group)  # partial paired mode
  #design <- model.matrix(~  group + patient)  # partial paired mode, design order has great impact
  design <- model.matrix(~ factor(group))
  y <- estimateDisp(y, design)
  if(method == 'edger_glmqlf'){
    fit <- glmQLFit(y, design)
    test <- glmQLFTest(fit, coef=)
    res <- topTags(test, n=nrow(mat), sort.by='none')
  } else if(method == 'edger_glmlrt'){
    fit <- glmFit(y, design)
    test <- glmLRT(fit, coef=2) # coef should be the pos col, ncol(design)=2
    res <- topTags(test, n=nrow(mat), sort.by='none')
  } else if(method == 'edger_exact'){
    test <- exactTest(y)
    res <- topTags(test, n=nrow(mat), sort.by='none')
  }
  res <- cbind(res$table, baseMean=2^(res$table$logCPM))
  # rename columns
  mapped_names <- colnames(res)
  for(i in 1:ncol(res)){
    if(colnames(res)[i] == 'logFC'){
      mapped_names[i] <- 'log2FoldChange'
    }else if(colnames(res)[i] == 'PValue'){
      mapped_names[i] <- 'pvalue'
    }else if(colnames(res)[i] == 'FDR') {
      mapped_names[i] <- 'padj'
    }else{
      mapped_names[i] <- colnames(res)[i]
    }
  }
  colnames(res) <- mapped_names
  
  if (filterType=="small"){
    normMat <- logcpm
  }else if(filterType=="long"){
    normMat <- logtpm
  }
  
  outfile <- list()
  outfile[["normMat"]] <- normMat
  outfile[["diffTable"]] <- res
  return(outfile)
}



#not paired (optional with batch)
diff.v2.dcb <- function(mat,samples,group, batch="NULL", method,norm_method, filterType="NULL",featureType, getCountGt1SmpFrac="Y", getPercExpByGrp="Y"){
  filterType <- "NULL"
  suppressPackageStartupMessages(library(edgeR))
  mat <- mat[,samples]
  print(dim(mat))
  y <- DGEList(counts=mat, samples=samples, group=group)
  
  ## filter low expr
  #keep <- filterByExpr(y,min.count = 1, min.total.count = 5, min.prop = 0.1)
  message(filterType)
  counts <- edgeR::getCounts(y)
  
  if(getCountGt1SmpFrac!="NULL"){
    ## get count>1 sample freq
    counts.neg <- counts[,group=="negative"]
    gt1.neg <- rowSums(counts.neg>=1)/length(samples[group=="negative"])
    counts.neg <- NULL
    counts.pos <- counts[,group=="positive"]
    gt1.pos <- rowSums(counts.pos>=1)/length(samples[group=="positive"])
    counts.pos <- NULL
  }
  
  if (filterType=="small"){
    keep <- rowSums(counts>=1) >= 0.5*length(samples)
  }else if(filterType=="long"){
    if(featureType=="gene"){
      gene.len = as.numeric(unlist(lapply(strsplit(rownames(counts), "|", fixed = T), function(x) x[2])))
    }else if (featureType=="domain"){
      gene.len = as.numeric(unlist(lapply(strsplit(rownames(counts), "|", fixed = T), function(x) x[7]))) - as.numeric(unlist(lapply(strsplit(rownames(counts), "|", fixed = T), function(x) x[6])))
    }
    tpm_mat <- tpm(count = counts, gene.len = gene.len)
    keep <- rowSums(tpm_mat >= 1) >= 0.5*length(samples)
    
    tpm_mat <- tpm_mat[keep,]
    logtpm <- log2(tpm_mat+1)
    # print(dim(logtpm))
    # print(logtpm[1:2,1:2])
  }else{
    message("keep all, not prefilter")
    keep <- rowSums(counts) >= 0
  }
  # if (filterType=="small"){
  #   keep <- rowSums(counts>=1) >= 0.5*length(samples)
  # }else if(filterType=="long"){
  #   tpm_mat <- tpm(count = counts, 
  #                  gene.len = as.numeric(unlist(lapply(strsplit(rownames(counts), "|", fixed = T), function(x) x[2]))))
  #   keep <- rowSums(tpm_mat >= 1) >= 0.5*length(samples)
  # }else{
  #   message("keep all, not prefilter")
  #   keep <- rowSums(counts) >= 0
  # }
  #min.count: Minimum count required for at least some samples.
  #min.total.count: Minimum total count required.
  #min.prop: Minimum proportion of samples in the smallest group that express the gene.
  y <- y[keep,,keep.lib.sizes=FALSE]
  print(dim(y))
  
  y <- calcNormFactors(y, method=norm_method)
  
  cpm <- edgeR::cpm(y, log=F, normalized.lib.sizes = T,prior.count = 0)
  logcpm <- log2(cpm+1)
  
  if(getPercExpByGrp!="NULL"){
    ## get percentile level by group
    cpm.neg <- cpm[,group=="negative"]
    #dim(cpm)
    p90.neg <- apply(cpm.neg,1, function(x) quantile(x, probs = max(2/ncol(cpm.neg),0.9) )) # 0.9; min 2 samples
    p50.neg <- apply(cpm.neg,1, function(x) quantile(x, probs = max(1/ncol(cpm.neg),0.5) )) # 0.5; min 1 sample
    p10.neg <- apply(cpm.neg,1, function(x) quantile(x, probs = max(1/ncol(cpm.neg),0.1) )) # 0.1; min 1 sample
    mean.neg <- apply(cpm.neg,1, mean ) 
    cpm.neg <- NULL
    cpm.pos <- cpm[,group=="positive"]
    #dim(cpm)
    p90.pos <- apply(cpm.pos,1, function(x) quantile(x, probs = max(2/ncol(cpm.pos),0.9) )) # 0.9; min 2 samples
    p50.pos <- apply(cpm.pos,1, function(x) quantile(x, probs = max(1/ncol(cpm.pos),0.5) )) # 0.5; min 1 sample
    p10.pos <- apply(cpm.pos,1, function(x) quantile(x, probs = max(1/ncol(cpm.pos),0.1) )) # 0.1; min 1 sample
    mean.pos <- apply(cpm.pos,1, mean ) 
    cpm.pos <- NULL
  }
  
  #design <- model.matrix(~ patient + group)  # partial paired mode
  #design <- model.matrix(~  group + patient)  # partial paired mode, design order has great impact
  if(batch=="NULL"){
    design <- model.matrix(~ factor(group))
  }else if(length(batch)==length(group)){
    design <- model.matrix(~ factor(group) + batch)
  }else {
    message("not same length of batch with group")
  }
  y <- estimateDisp(y, design)
  if(method == 'edger_glmqlf'){
    fit <- glmQLFit(y, design)
    test <- glmQLFTest(fit, coef=2)
    res <- topTags(test, n=nrow(mat), sort.by='none')
  } else if(method == 'edger_glmlrt'){
    fit <- glmFit(y, design)
    test <- glmLRT(fit, coef=2) # coef should be the pos col, ncol(design)=2
    res <- topTags(test, n=nrow(mat), sort.by='none')
  } else if(method == 'edger_exact'){
    test <- exactTest(y) # no design !
    res <- topTags(test, n=nrow(mat), sort.by='none')
  }
  
  res <- cbind(res$table, baseMean=2^(res$table$logCPM))
  if(getCountGt1SmpFrac!="NULL" & getPercExpByGrp!="NULL"){
    res <- cbind(res,negGT1RatioCount=gt1.neg,posGT1RatioCount=gt1.pos,negMeanCPM=mean.neg,negCent90CPM=p90.neg,negCent50CPM=p50.neg,negCent10CPM=p10.neg, posMeanCPM=mean.pos, posCent90CPM=p90.pos,posCent50CPM=p50.pos,posCent10CPM=p10.pos)
  }else if(getCountGt1SmpFrac=="NULL" & getPercExpByGrp!="NULL"){
    res <- cbind(res,negMeanCPM=mean.neg,negCent90CPM=p90.neg,negCent50CPM=p50.neg,negCent10CPM=p10.neg, posMeanCPM=mean.pos, posCent90CPM=p90.pos,posCent50CPM=p50.pos,posCent10CPM=p10.pos)
  }else if(getCountGt1SmpFrac!="NULL" & getPercExpByGrp=="NULL"){
    res <- cbind(res,negGT1RatioCount=gt1.neg,posGT1RatioCount=gt1.pos)
  }
    
  # rename columns
  mapped_names <- colnames(res)
  for(i in 1:ncol(res)){
    if(colnames(res)[i] == 'logFC'){
      mapped_names[i] <- 'log2FoldChange'
    }else if(colnames(res)[i] == 'PValue'){
      mapped_names[i] <- 'pvalue'
    }else if(colnames(res)[i] == 'FDR') {
      mapped_names[i] <- 'padj'
    }else{
      mapped_names[i] <- colnames(res)[i]
    }
  }
  colnames(res) <- mapped_names
  #return(res)
  
  if (filterType!="long"){
    normMat <- logcpm
  }else if(filterType=="long"){
    normMat <- logtpm
  }
  
  outfile <- list()
  outfile[["normMat"]] <- normMat
  outfile[["diffTable"]] <- res
  return(outfile)
}


## diff wilcox paired 
wilcox.paired <- function(matrix_cpm,group,paired=TRUE){
  # matrix_cpm = mat.norm.tmp
  # group = sample.table.tmp$source
  # paired = FALSE
  suppressPackageStartupMessages(library(edgeR))
  test_func <- function(x){
    wilcox.test(as.numeric(x[group == levels(group)[1]]), as.numeric(x[group == levels(group)[2]]), alternative='two.sided', paired = paired)$p.value
  }
  #matrix_cpm <- cpm(mat)
  #table(group == levels(group)[2])
  #wilcox.test(as.numeric(matrix_cpm[1,group == levels(group)[1]]), as.numeric(matrix_cpm[1,group == levels(group)[2]]), alternative='two.sided', paired = paired)$p.value
  #matrix_cpm
  pvalues <- apply(matrix_cpm, 1, test_func)
  treatMeans <- apply(matrix_cpm[,which(group == levels(group)[2])], 1, mean)
  ctrlMeans <- apply(matrix_cpm[,which(group == levels(group)[1])], 1, mean)
  logFC <- treatMeans - ctrlMeans
  res <- data.frame(log2FoldChange=logFC,
                    pvalue=pvalues,
                    padj=p.adjust(pvalues, method='BH'),
                    baseMean=apply(matrix_cpm, 1, mean),
                    treatMean=treatMeans,
                    ctrlMean=ctrlMeans)
  return(res)
}




## def cancer Peak index score func.
#tail gene not correct outlier gene number, just sum up all |zsocre|>3, thus not suited for multi-cancer classify
getRawPeakIndex <- function(feature.lab,logcpm,sample.table,prescale=F,mean.list,sd.list,zscore.cutoff=100,onlyCountZscoreOutlier=F){
  #x <- paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",dst,"/cfPeakCNN_GSE110381_smallDomain_diff_COADvsLAML.cpm")
  #logcpm <- read.table(paste0(pre,"/output/",dst,"/call_peak_all/count_matrix/expeakCNN_b5_d50_p1_GSE110381_CPMrowsum.txt"),check.names = F,header = T) # cpm
  #logcpm <- read.table(x,check.names = F,header = T)  # log2cpm+1
  #rownames(logcpm) <- logcpm$gene_id # feature
  # feature.lab <- "COAD"  #"all"
  cpm <- 2^logcpm-1 # sum before log is more sensitive
  # message(table(rownames(cpm) %in% feature.list[[feature.lab]]))
  cpm <- cpm[rownames(cpm) %in% feature.list[[feature.lab]],] # some may lost some features
  #sample.table <- sample.table.tcga
  #cpm <- cpm[,2:ncol(cpm)]
  sample.table <- sample.table[sample.table$sample %in% colnames(cpm),]
  # table(sample.table$group)
  # positive_samples <- sample.table[sample.table$group==disease,"sample"]
  # negative_samples <- sample.table[sample.table$group!=disease,"sample"]
  # grpNum <- length(unique(sample.table$group))
  # if(grpNum>2){
  #   other_samples <- list()
  #   for(grp in as.character(unique(sample.table$group[!(sample.table$group %in% c(disease,NC.label))] ))   ){
  #     #grp <- "PAAD" 
  #     #grp <- "PRAD"
  #     print(grp)
  #     other_samples[[grp]] <- sample.table[sample.table$group==grp,"sample"]
  #     samples <- c(samples,other_samples[[grp]])
  #     group <- c(group,rep(grp,length(other_samples[[grp]])))
  #   }
  #   sample.table <- sample.table[match(samples,sample.table$sample),]
  # }
  cpm <- cpm[,sample.table$sample]
  message("dim(cpm): ",dim(cpm))
  
  if(prescale){
    message("prescale")
    #standardize by coad&blood group
    cpm <- (cpm-mean.list[[feature.lab]])/sd.list[[feature.lab]]
    #trim to avoid too much outlier, may not need  
    # print(max(cpm))
    # print(min(cpm))
    # zscore.cutoff <- 10
    # tmp <- cpm > zscore.cutoff
    cpm[cpm > zscore.cutoff] <- zscore.cutoff
    cpm[cpm < -zscore.cutoff] <- -zscore.cutoff
    # onlyCountZscoreOutlier <- 3
    if(onlyCountZscoreOutlier){
      message("onlyCountZscoreOutlier: ",onlyCountZscoreOutlier)
      table(cpm>=abs(onlyCountZscoreOutlier))
      cpm[cpm<abs(onlyCountZscoreOutlier)] <- 0 # abs
      cpm[cpm>=abs(onlyCountZscoreOutlier)] <- 1 # abs
      message(cpm[1:2,1:2])
    }
  } else {
    message("not prescale")
  }
  return(cpm)
}

getPeakIndex <- function(feature.lab,logcpm,sample.table,prescale=F,mean.list,sd.list,zscore.cutoff=100,onlyCountZscoreOutlier=F){
  # feature.lab <- "COAD"
  # sample.table <- sample.table.tcga
  # prescale <- T
  # # mean.list <- mean.list[[]]
  # zscore.cutoff <- 100
  cpm <- getRawPeakIndex(feature.lab=feature.lab,logcpm = logcpm,sample.table = sample.table,prescale = prescale,mean.list = mean.list,sd.list = sd.list,zscore.cutoff = zscore.cutoff,onlyCountZscoreOutlier=onlyCountZscoreOutlier)
  library(tidyr)
  library(dplyr)
  conflict_prefer("summarise", "dplyr")
  cpm$feature <- rownames(cpm)
  cpm <- as_tibble(cpm) %>% 
    pivot_longer(cols = 1:(ncol(cpm)-1), names_to = "sample", values_to = "log2cpm")   # cpm acutually !!!
  cpm$group <- sample.table$group[match(cpm$sample,sample.table$sample)]
  #table((cpm$group))
  cpm.sum <- cpm %>% 
    group_by(sample,group) %>% 
    summarise(sum.log2cpm=sum(log2cpm),mean.log2cpm=mean(log2cpm)) 
  cpm.sum$sum.log2cpm.scale <- scale(cpm.sum$sum.log2cpm)[,1] # standize scale
  cpm.sum$sum.log2cpm.scale <- maxmin.normalize.vec(cpm.sum$sum.log2cpm.scale) # min-max scale (better for visualize)
  cpm.sum$mean.log2cpm.scale <- scale(cpm.sum$mean.log2cpm)[,1] # standize scale
  cpm.sum$mean.log2cpm.scale <- maxmin.normalize.vec(cpm.sum$mean.log2cpm.scale) # min-max scale (better for visualize)
  
  # grpNum <- length(unique(sample.table$group))
  # if(grpNum==2){
  #   cpm.sum$group <- factor(cpm.sum$group,levels = c(NC.label,CRC.label))
  # }else if(grpNum>2){
  #   cpm.sum$group <- factor(cpm.sum$group,levels = c(NC.label,CRC.label, unique(sample.table$group[!(sample.table$group %in% c(CRC.label,NC.label))] )))
  # }
  cpm.sum$peak.precursor <- feature.lab
  return(cpm.sum)
}


library(ggplot2)

## theme 
bar_theme <- theme_minimal() +  # base_size=12
  theme(#axis.ticks.x=element_blank(),
    #strip.text.y = element_blank(),
    aspect.ratio = 1,
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
my_theme_point <-   theme(plot.title = element_text(size = 28,color="black",hjust = 0.5),
                    axis.title = element_text(size = 28,color ="black"),
                    axis.text = element_text(size= 24,color = "black"), #,face="bold
                    panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),
                    axis.text.x = element_text( hjust = 0.5 ), # angle = 45,
                    axis.text.y = element_text( hjust = 0 ), # angle = 45,
                    panel.grid=element_blank(),
                    legend.position = "right",#c(0.5,0.3),
                    legend.text = element_text(size= 26,color = "black"),
                    legend.title= element_text(size= 26,color = "black"))
my_theme_box <- theme(aspect.ratio = 1.5,
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
my_theme_roc <- theme(aspect.ratio = 1,
                      plot.title = element_text(size = 24,color="black",hjust = 0.5),
                      axis.title = element_text(size = 24,color ="black"), 
                      axis.text = element_text(size= 24,color = "black"),
                      axis.text.x = element_text(size= 24,color = "black"),
                      #panel.grid=element_blank(),
                      # panel.grid.major.x=element_blank(),
                      # panel.grid.minor.x = element_blank(),
                      # panel.grid.major.y = element_blank(), #element_line(color = "grey50",linetype = "dashed"), #size= 1,
                      # panel.grid.minor.y = element_blank(),
                      #panel.grid.minor.y = element_blank(),
                      panel.border = element_blank(),
                      legend.position = "right",#c(.25,.6),
                      legend.text = element_text(size= 16),
                      legend.title= element_text(size= 20),
                      strip.text.y = element_blank(),
                      strip.text.x = element_text(size=24)
)
my_theme_pie <- theme(
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
  legend.title= element_text(size= 24)
)
my_theme_volcano <-   theme(aspect.ratio=1,
                            plot.title = element_text(size = 28,color="black",hjust = 0.5),
                            axis.title = element_text(size = 28,color ="black"), 
                            axis.text = element_text(size= 24,color = "black"), #,face="bold
                            panel.grid.minor.y = element_blank(),
                            panel.grid.minor.x = element_blank(),
                            axis.text.x = element_text( hjust = 0.5 ), # angle = 45,
                            axis.text.y = element_text( hjust = 0 ), # angle = 45,
                            panel.grid=element_blank(),
                            legend.position = "right",#c(0.5,0.3),
                            legend.text = element_text(size= 24,color = "black"),
                            legend.title= element_text(size= 26,color = "black")
)

plotViolin <- function(logcpm.sum,sig.size=10,colors=c("grey50","firebrick")){
  p1 <- ggplot(logcpm.sum, aes(x=group,y=value,fill=group))+ # 
    geom_violin(alpha=0.7)+
    geom_boxplot(width=0.1)+
    scale_fill_manual(name="Group",values = colors)+#colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(unique(logcpm.sum$group))))+)
    # facet_grid(.~name,scales = "free")+
    ylab("Scaled peak-index") +
    ggpubr::stat_compare_means(
      label.x.npc = "middle", label.y.npc = "top",
      #size =12,step.increase = 0.08,#paired = TRUE,
      aes(x=group,y=value,
          label = ..p.format..), # p.signif, p.format, group = Group,
      symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01,1),
                       symbols = c( "***", "**", "*", "ns")),
      method.args = list(alternative="greater"),
      ref.group = NC.label,hide.ns=F, size=sig.size,
      method = "wilcox.test"
    ) +
    theme_bw() + 
    my_theme_box
  return(p1)
}
plotROC <- function(logcpm.sum,direction="auto"){
  library(ggplot2)
  library(pROC)
  logcpm.sum <- as.data.frame(logcpm.sum)
  rocobj <- roc(group ~ value, data=logcpm.sum, direction=direction)
  ?roc
  auc <- round(rocobj$auc,4)
  p1 <- ggroc(rocobj, colour= 'steelblue', size = 2) + #  = 'steelblue',
    ggtitle(paste0('', '(AUC = ', auc, ')')) +
    theme_minimal() + 
    my_theme_roc
  return(p1)
}
plotMultiROC <- function(logcpm.sum.list,direction="auto",RNA_colors=RNA_colors){
  library(ggplot2)
  library(pROC)
  library(purrr)
  #library
  library(pROC)
  library(ggplot2)
  library(tidyverse)
  
  logcpm.sum <- as.data.frame(do.call(rbind,logcpm.sum.list))
  #table(logcpm.sum$peak.precursor)
  logcpm.sum$sum.log2cpm <- NULL
  logcpm.sum <- reshape2::dcast(data = logcpm.sum, formula = sample+group~peak.precursor, id.var=c("sample","group"), value.var = "value")
  # logcpm.sum[1:3,]
  # logcpm.sum <- as.data.frame(logcpm.sum)
  
  rocobj.list <- list()
  # auc.list <- list()
  # rocobj <- roc(group ~ sum.log2cpm.scale + mir.sum.log2cpm.scale, data=logcpm.sum)
  for(i in names(logcpm.sum.list)){
    print(i)
    rocobj.list[[i]] <- roc(logcpm.sum[["group"]], logcpm.sum[[i]],direction=direction ) # , data=logcpm.sum
    print(rocobj.list[[i]]$auc)
    # auc.list[[i]] <- round(auc(logcpm.sum[["group"]], logcpm.sum[[i]]),4) # , data=logcpm.sum
  }
  # # extract auc
  # names(rocobj.list)
  # rocobj.list %>% 
  #   map(~tibble(AUC = .x$auc)) %>% 
  #   bind_rows(.id = "name") -> data.auc
  
  
  # generate labels labels
  tibble(name=names(rocobj.list), AUC=as.numeric(do.call(c,lapply(rocobj.list,FUN = function(y) {return(as.numeric(y$auc)) } )) ) ) %>% 
    mutate(label_long=paste0(name,", AUC=",paste(round(AUC,2))),
           label_AUC=paste0("AUC=",paste(round(AUC,2)))) -> data.labels
  
  # # plot a facet plot with AUC within plots
  # ggroc(roc.list) +
  #   facet_wrap(~name) +
  #   
  #   geom_text(data = data.labels,
  #             aes(0.5, 1, 
  #                 label = paste(label_AUC)),
  #             hjust = 1) 
  
  # plot on a single plot with AUC in labels
  p1 <- ggroc(rocobj.list, aes="colour", size = 2, alpha=0.9) + #  = 'steelblue',
    # ggtitle(paste0('', '(AUC = ', auc, ')')) +
    scale_color_manual(values = RNA_colors, labels=data.labels$label_long) +
    guides(fill = guide_legend(title = "Precursor" )) +
    theme_minimal() + 
    my_theme_roc
  
  return(p1)
}
#plotROC(logcpm.sum.list[["all"]])
plotConfusionMat<-function(Actual,Predict,colors=c("white","red4","dodgerblue3"),text.scl=5){
  actual = as.data.frame(table(Actual))
  names(actual) = c("Actual","ActualFreq")
  
  #build confusion matrix
  confusion = as.data.frame(table(Actual, Predict))
  names(confusion) = c("Actual","Predicted","Freq")
  
  #calculate percentage of test cases based on actual frequency
  
  confusion = merge(confusion, actual, by=c('Actual','Actual'))
  confusion$Percent = confusion$Freq/confusion$ActualFreq*100
  confusion$ColorScale<-confusion$Percent*-1
  confusion[which(confusion$Actual==confusion$Predicted),]$ColorScale<-confusion[which(confusion$Actual==confusion$Predicted),]$ColorScale*-1
  # confusion$Label<-paste(round(confusion$Percent,0),"%, n=",confusion$Freq,sep="")
  confusion$Label<-paste(confusion$Freq,sep="")
  confusion$Label[confusion$Label==0] <- ""
  tile <- ggplot() +
    geom_tile(aes(x=Actual, y=Predicted,fill=ColorScale),data=confusion, color="black",size=0.1) +
    labs(x="Actual",y="Predicted")
  tile = tile +
    geom_text(aes(x=Actual,y=Predicted, label=Label),data=confusion, size=text.scl, colour="black") +
    scale_fill_gradient2(low=colors[2],high=colors[3],mid=colors[1],midpoint = 0,guide='none') +
    theme_bw() +
    theme(
      aspect.ratio = 1,
      plot.title = element_text(size = 24,color="black",hjust = 0.5),
      axis.title = element_text(size = 24,color ="black"), 
      axis.text = element_text(size= 24,color = "black"),
      axis.text.x = element_text(size= 24,color = "black"),
      #panel.grid=element_blank(),
      panel.grid.major.x = element_line(color = "grey30",linetype = "dashed"),,
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey30",linetype = "dashed"), #size= 1,
      panel.grid.minor.y = element_blank(),
      #panel.grid.minor.y = element_blank(),
      panel.border = element_blank(),
      legend.position = "right",#c(.25,.6),
      legend.text = element_text(size= 24),
      legend.title= element_text(size= 24),
      strip.text.y = element_blank(),
      strip.text.x = element_text(size=24)
    )
}


plotPie <- function(x){
  #x <- c(1,2,2,3,3,4) 
  # x <- df$RNA
  df2 <- as.data.frame(table(x))
  df2 <- df2[df2$Freq>0,]
  df2 <- dplyr::as_tibble(df2) %>% 
    dplyr::group_by(x) %>% 
    dplyr::summarize(Freq=mean(Freq))
  df2$lab <- round(df2$Freq/sum(df2$Freq),digits = 3)
  df2 <- df2 %>%
    dplyr::mutate(csum = rev(cumsum(rev(Freq))),
                  pos = Freq/2 + lead(csum, 1),
                  pos = if_else(is.na(pos), Freq/2, pos))
  
  RNA.col <- data.frame(RNA=c(rna,dna),col=c(pal_nejm_adaptive()(15)[1:14],"#11838D"))
  RNA.col$RNA <- factor(RNA.col$RNA,levels = c(rna,dna))
  
  # str(df2)
  p <- ggplot(df2, aes(x="", y=Freq, fill=x)) +
    geom_bar(stat="identity", width=1, color="black") +
    coord_polar("y", start=0) +
    # geom_text(aes(label = lab), position = position_stack(vjust=0.5), size = 5, colour="white") +
    ggrepel::geom_label_repel(data = df2, col="white",
                              aes(y = pos, label = paste0(100*lab, "%")),
                              size = 8, nudge_x = 0.75, show.legend = FALSE) +
    labs(x = NULL, y = NULL, title ="") + # paste0("Total repeats peak: ",nrow(peak))
    scale_color_manual("black") +
    scale_fill_manual(name="precursor",values = RNA.col$col[RNA.col$RNA %in% df2$x] ) + 
    # scale_fill_nejm_adaptive(alpha = 0.8) +
    theme_bw() + 
    my_theme_pie
  return(p)
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


# auto color num ggsci
library("ggsci")

#' Adaptive palette (discrete).
#'
#' Create a discrete palette which will use the first n colors from
#' the supplied color values, and interpolate after n.
adaptive_pal <- function(values) {
  force(values)
  function(n = 10) {
    if (n <= length(values)) {
      values[seq_len(n)]
    } else {
      colorRampPalette(values, alpha = TRUE)(n)
    }
  }
}

## npg
pal_npg_adaptive <- function(palette = c("nrc"), alpha = 1) {
  palette <- match.arg(palette)
  
  if (alpha > 1L | alpha <= 0L) stop("alpha must be in (0, 1]")
  
  raw_cols <- ggsci:::ggsci_db$"npg"[[palette]]
  raw_cols_rgb <- col2rgb(raw_cols)
  alpha_cols <- rgb(
    raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], raw_cols_rgb[3L, ],
    alpha = alpha * 255L, names = names(raw_cols),
    maxColorValue = 255L
  )
  
  adaptive_pal(unname(alpha_cols))
}

scale_color_npg_adaptive <- function(palette = c("nrc"), alpha = 1, ...) {
  palette <- match.arg(palette)
  discrete_scale("colour", "npg", pal_npg_adaptive(palette, alpha), ...)
}

scale_fill_npg_adaptive <- function(palette = c("nrc"), alpha = 1, ...) {
  palette <- match.arg(palette)
  discrete_scale("fill", "npg", pal_npg_adaptive(palette, alpha), ...)
}

##nejm
pal_nejm_adaptive <- function(palette = c("default"), alpha = 1) {
  palette <- match.arg(palette)
  
  if (alpha > 1L | alpha <= 0L) stop("alpha must be in (0, 1]")
  
  raw_cols <- ggsci:::ggsci_db$"nejm"[[palette]]
  raw_cols_rgb <- col2rgb(raw_cols)
  alpha_cols <- rgb(
    raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], raw_cols_rgb[3L, ],
    alpha = alpha * 255L, names = names(raw_cols),
    maxColorValue = 255L
  )
  
  adaptive_pal(unname(alpha_cols))
}
scale_color_nejm_adaptive <- function(palette = c("default"), alpha = 1, ...) {
  palette <- match.arg(palette)
  discrete_scale("colour", "nejm", pal_nejm_adaptive(palette, alpha), ...)
}
scale_fill_nejm_adaptive <- function(palette = c("default"), alpha = 1, ...) {
  palette <- match.arg(palette)
  discrete_scale("fill", "nejm", pal_nejm_adaptive(palette, alpha), ...)
}

##d3
pal_d3_adaptive <- function(palette = c("category10"), alpha = 1) {
  palette <- match.arg(palette)
  
  if (alpha > 1L | alpha <= 0L) stop("alpha must be in (0, 1]")
  
  raw_cols <- ggsci:::ggsci_db$"d3"[[palette]]
  raw_cols_rgb <- col2rgb(raw_cols)
  alpha_cols <- rgb(
    raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], raw_cols_rgb[3L, ],
    alpha = alpha * 255L, names = names(raw_cols),
    maxColorValue = 255L
  )
  
  adaptive_pal(unname(alpha_cols))
}
scale_color_d3_adaptive <- function(palette = c("category10"), alpha = 1, ...) {
  palette <- match.arg(palette)
  discrete_scale("colour", "d3", pal_d3_adaptive(palette, alpha), ...)
}
scale_fill_d3_adaptive <- function(palette = c("category10"), alpha = 1, ...) {
  palette <- match.arg(palette)
  discrete_scale("fill", "d3", pal_d3_adaptive(palette, alpha), ...)
}



rna <- c("rRNA","pri_miRNA","lncRNA", "mRNA", "piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA')
dna <- c("intron", "promoter", "enhancer", "repeats")

RNA_colors <- list()
for(i in 1:length(c(rna,dna))){
  j <- c(rna,dna)[i]
  RNA_colors[[j]] <- c(pal_nejm_adaptive()(15)[1:14],"#11838D")[i] #pal_d3_adaptive()(15)[i]
}
RNA_colors[['all']] <- "grey50" #"grey50"
RNA_colors[['RNA']] <- "red2"
RNA_colors[['DNA']] <- "purple" #"#117C8E": similar with repeats
RNA_colors <- do.call("c",RNA_colors)

