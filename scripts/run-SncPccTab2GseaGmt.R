#!/usr/bin/env Rscript

# Load necessary libraries
options(stringsAsFactors = F)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(parallel)) # parallelly
#suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

# Create parser object
parser <- ArgumentParser(description = "Perform GSEA on PCC results")

# Define command line arguments
parser$add_argument("--input_file", type = "character", required = TRUE, help = "Input PCC result table file path")
parser$add_argument("--outPrefix", type = "character", default="./", help = "Output prefix for results")
parser$add_argument("--gene_set", type = "character", default="/BioII/lulab_b/baopengfei/shared_reference/geneset/PATH_ID_NAME_KEGGplusHallmarkReactomeImmReg.txt", 
                    help = "Gene set file path, longer form of gmt file")
parser$add_argument("--dataset", type = "character", default="dst", help = "dataset label")
parser$add_argument("--cell", type = "character", default="cell", help = "cell label")
parser$add_argument("--method", type = "character", default="method", help = "method label")
parser$add_argument("-c","--cores", type = "integer", default = 3, help = "Number of cores to use, lower this if killed by sys [default: 3]")
parser$add_argument('--PathSource', type='character', default="CHG,ImmReg,KEGG,MSigDB,Reactome", help='Gmt Long format tsv source selected, default: CHG,ImmReg,KEGG,MSigDB,Reactome')
parser$add_argument('--sncRES_cutoff', type='numeric', default=0.995, help='abs(sncRNA EnrichScore) cutoff, [0,1), default: 0.995 (>0.995 kept)')
parser$add_argument('--p_cutoff', type='numeric', default=0.1, help='P cutoff, (0,1], default: 0.1 (<0.1 kept)')
parser$add_argument('--p_label', type='character', default="p.adjust", help='pvalue label of column, pvalue or p.adjust  default: p.adjust')


# Parse command line arguments
args <- parser$parse_args()
for(i in 1:length(args)){
  v.name <- names(args)[i]
  assign(v.name,args[[i]])
  message(paste0(v.name,": ",get(v.name)))
}

# Extract arguments
message("start at: ",date())


# # test TCGA (2408)
# input_file <- "output/TCGA_small16/pcc_mRNA/HNSC_pearson_rmPurity.txt"
# outPrefix <- "output/TCGA_small16/pcc_mRNA/rmPurity"
# gene_set <- "/BioII/lulab_b/baopengfei/shared_reference/geneset/PATH_ID_NAME_KEGGplusHallmarkReactomeImmReg.txt"
# dataset <- "TCGA_small16"
# cell <- "HNSC"
# cores <- 10
# method <- "pearson"
# PathSource <- "CHG,ImmReg,KEGG,MSigDB"
# sncRES_cutoff <- 0.9 # 0.995
# p_label <- "p.adjust"
# p_cutoff <- 0.05 # 0.1

# # test SLE
# input_file <- "output/SLE/pcc_mRNA/Int.Monocytes_pearson_rmMonthOperator.txt.gz"
# outPrefix <- "output/SLE/pcc_mRNA/rmMonthOperator"
# gene_set <- "/BioII/lulab_b/baopengfei/shared_reference/geneset/PATH_ID_NAME_KEGGplusHallmarkReactomeImmReg.txt"
# cores <- 10
# dst <- "SLE"
# cellID <- "Int.Monocytes"
# method <- "pearson"
# PathSource <- "CHG,ImmReg,KEGG,MSigDB"
# sncRES_cutoff <- 0.9 #  0.995 # [0,1)
# p_cutoff <- 0.1 #0.05 # (0,1]
# p_label <- "pvalue" #"pvalue" # p.adjust

# Read input files
#pcc.diff <- data.table::fread(cmd=paste0("gzip -dc ",input_file), header = TRUE, stringsAsFactors = FALSE, data.table = FALSE, sep = "\t")
pcc.diff <- data.table::fread(input_file, header = TRUE, stringsAsFactors = FALSE, data.table = FALSE, sep = "\t")
print(head(pcc.diff))
#dim(pcc.diff)
# summary(pcc.diff$pvalue)
#dim(pcc.diff)
#length(unique(pcc.diff$sncRNA))
if(sum(grepl("peak|cluster|block",pcc.diff$sncRNA[1],perl = T))==1 & sum(grepl("|",pcc.diff$sncRNA[1],fixed = T))==1){
  message("use 4th | as "," peak id")
  pcc.diff$sncRNA <- unlist(sapply(strsplit(pcc.diff$sncRNA, "|", fixed = TRUE), "[", 4))
}
if(sum(grepl("ENSG",pcc.diff$mRNA[1],perl = T))==1  & sum(grepl("|",pcc.diff$mRNA[1],fixed = T))==1 ){
  message("use 1st | as "," gene id")
  pcc.diff$mRNA <- unlist(sapply(strsplit(pcc.diff$mRNA, "|", fixed = TRUE), "[", 1))
}
if("RS" %in% colnames(pcc.diff)){
  message("input is partial cor with covar")
  pcc.diff[["rank"]] <- pcc.diff$RS
}else if("cor" %in% colnames(pcc.diff)){
  message("input is cor no covar")
  pcc.diff[["rank"]] <- pcc.diff$cor
}else{
  stop("no cor or RS column !")
}
pcc.diff <- pcc.diff[order(pcc.diff$rank, decreasing = TRUE), ]
cor_tab <- split(pcc.diff, f = pcc.diff$sncRNA)
head(names(cor_tab))
#message(system('free -m'))
rm(list="pcc.diff")
gc() # put the garbage collection to free up memory before the next round of the loop
#message(system('free -m'))
#system(paste0("cat /proc/",Sys.getpid(),"/status | grep VmSize"))


m_t <- data.table::fread(gene_set, sep = "\t", header = TRUE, stringsAsFactors = FALSE, data.table = F)
PathSource <- strsplit(PathSource,",")[[1]]
m_t <- m_t[m_t$source %in% PathSource,]
m_t2g <- na.omit(m_t[, c("short_DESCRPTION", "ensembl_gene_id", "source")])
#m_t2e <- na.omit(m_t[, c("short_DESCRPTION", "ENTREZID")])
rm(list="m_t")
#message(system('free -m'))

# Function to perform GSEA for a given sncRNA
gsea.kegg.offline <- function(gene.list,m_t,qValue=0.2, pValue=0.1, adjustPval="BH"){
  # message("start gsea...")
  # gene.list <- rank.list
  # m_t <- m_t2g
  gene.list <- gene.list[!duplicated(names(gene.list))]
  tryCatch(
    {
      GSEA.KEGG <- GSEA(gene.list, 
                       TERM2GENE = m_t, 
                       nPerm  = 500,
                       minGSSize = 10,
                       maxGSSize = 500,
                       pvalueCutoff = pValue,
                       pAdjustMethod = adjustPval
      )
      GSEA.KEGG <- GSEA.KEGG@result
      # if(nrow(GSEA.KEGG)>0){
      #   GSEA.KEGG$trend <- "all"
      #   GSEA.KEGG$minus.log10.pvalue <- -log10(GSEA.KEGG$pvalue)
      #   GSEA.KEGG$minus.log10.padj <- -log10(GSEA.KEGG$p.adjust)
      #   GSEA.KEGG$minus.log10.qvalue <- -log10(GSEA.KEGG$qvalue)
      # }
      return(as.data.frame(GSEA.KEGG))
    },error = function(e) {
      message("not enough valid all.reg gene !")
      #dim[1,.]
      return(data.frame(ID="",Description="",setSize="",enrichmentScore="",NES="",pvalue="",p.adjust="",qvalues="",rank="",leading_edge="",core_enrichment="")) # ,trend="",minus.log10.pvalue="",minus.log10.padj="",minus.log10.qvalue=""
    }
  )
}
getSncGSEA <- function(sncRNA) { # , outdir = "."
  # sncRNA <- "peak_10026"
  gsea_outfile_snc <- paste0(outPrefix,"_",dataset,"_",cell,"_", sncRNA, ".txt")
  select.col <- c("ID","enrichmentScore","NES","pvalue","p.adjust","qvalues","source","sncRNA","sncRES","sig")
  
  if(file.exists(gsea_outfile_snc)){
    message("read existing kegg GSEA: ",sncRNA)
    res.tmp <- data.table::fread(gsea_outfile_snc, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE, data.table = FALSE)
    res.tmp <- res.tmp[,select.col]
  }else{
    tmp <- cor_tab[[sncRNA]]
    gene.list.ncbiid <- tmp$rank
    names(gene.list.ncbiid) <- tmp$mRNA
    gene.list.ncbiid <- na.omit(gene.list.ncbiid)
    gene.list.ncbiid <- gene.list.ncbiid[!is.na(names(gene.list.ncbiid))]
    
    set.seed(1234) # important!
    # GSEA.KEGG <- GSEA(gene.list.ncbiid,
    #                  TERM2GENE = m_t2g[, 1:2],
    #                  nPerm = 500,
    #                  minGSSize = 5,
    #                  maxGSSize = 500,
    #                  pvalueCutoff = 0.1,
    #                  pAdjustMethod = "BH",
    #                  by = "fgsea"
    # )
    res.tmp <- gsea.kegg.offline(gene.list = gene.list.ncbiid, m_t = m_t2g[, 1:2], pValue = 0.1, adjustPval = "BH")
    res.tmp <- res.tmp[res.tmp$ID!="",]
    # types <- list(GSEA.KEGG)
    # names(types) <- c("GSEA.KEGG")
    # # for (i in seq_along(types)) {
    # res.tmp <- as.data.frame(GSEA.KEGG) #types[[i]])
    
    if(nrow(res.tmp)==0){
      # message("no GSEA record for ",sncRNA)
      res.tmp <- data.frame(ID="",Description="",setSize=0,enrichmentScore=0,NES=0,pvalue=1,p.adjust=1,qvalues=1,
                            source="",sncRNA=sncRNA,sncRES=0,sig=0)
    }else{
      res.tmp$source <- m_t2g$source[match(res.tmp$ID, m_t2g$short_DESCRPTION)]
      res.tmp$sncRNA <- sncRNA
      res.tmp$sncRES <- 1 - 2 * res.tmp[[p_label]]
      res.tmp$sncRES[res.tmp$enrichmentScore < 0] <- 2 * res.tmp[[p_label]][res.tmp$enrichmentScore < 0] - 1
      res.tmp$sig <- ifelse(abs(res.tmp$sncRES) > sncRES_cutoff & res.tmp[[p_label]] < p_cutoff, 1, 0)
    }
    res.tmp <- res.tmp[,select.col]
    res.tmp <- res.tmp[res.tmp$sig == 1, ]
    
    data.table::fwrite(res.tmp, gsea_outfile_snc, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t") # prevent too many fail-rerun
  }
  
  return(NULL)
  # }
}

readSncGSEA <- function(sncRNA) {
  # sncRNA <- "peak_10026"
  gsea_outfile_snc <- paste0(outPrefix, "_",dataset,"_",cell, "_", sncRNA, ".txt")
  select.col <- c("ID","enrichmentScore","NES","pvalue","p.adjust","qvalues","source","sncRNA","sncRES","sig")
  res.tmp <- data.table::fread(gsea_outfile_snc, header = T, sep = "\t", check.names = F, quote = "", stringsAsFactors = F, data.table = F) 
  res.tmp <- res.tmp[,select.col]
  this.source <- m_t2g$source[match(res.tmp$ID, m_t2g$short_DESCRPTION)]
  res.tmp <- res.tmp[this.source %in% PathSource,]
  res.tmp <- res.tmp[res.tmp$sig == 1, ]
  return(res.tmp)
}

# Perform GSEA for each sncRNA
message("perform gsea enirch: ",date())
tmpL <- mclapply(names(cor_tab), getSncGSEA, mc.cores = cores) # outdir = output_dir,
# tmpL <- lapply(names(cor_tab)[1:3], getSncGSEA)
#saveRDS(object = cor_tab, file = paste0(outPrefix, "_",cellID, "_",method,"_", "GSEA.KEGG", ".rds.bz2"), compress = "bzip2")
#test <- readRDS("output/TCGA_small7/pcc_mRNA/rmPurity_BRCA_pearson_GSEA.KEGG.rds.bz2") #"output/TCGA_small7/pcc_mRNA/rmPurity_BRCA_pearson_GSEA.KEGG.txt")
tmp.name <- names(cor_tab)
rm(list="cor_tab")
gc()
message("done gsea enirch: ",date())

tmp.res <- list()
tmp.res <- mclapply(tmp.name, readSncGSEA, mc.cores = cores) 
tmp.res <- as.data.frame(do.call(rbind,tmp.res))
#table(tmp.res$ID!="")
tmp.res <- tmp.res[tmp.res$ID!="",]
message(paste0(dim(tmp.res)," "))
#gsea_outfile <- paste0(outPrefix,"_",dataset,"_",cell, ".txt.gz")
gsea_outfile <- paste0(outPrefix,"_",dataset,"_",cell, ".txt.gz")
message("write gsea table: ",date())
if(!(nrow(tmp.res)>0)){
  tmp.res <- data.frame("ID"="","enrichmentScore"="","NES"="","pvalue"="","p.adjust"="","qvalues"="","source"="","sncRNA"="","sncRES"="","sig"="")
}
data.table::fwrite(tmp.res, gsea_outfile, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t", compress = "gzip")

# Collect and process GSEA results
# gseaTbl <- list()
# for (sncRNA in names(cor_tab)) {
#   x <- paste0(outPrefix, "_",cellID, "_",method,"_", "GSEA.KEGG", "_", sncRNA, ".txt")
#     #paste0(output_dir, "/GSEA.KEGG_", sncRNA, ".txt")
#   path1 <- data.table::fread(x, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE, data.table = FALSE)
#   path1 <- path1[path1$sig == 1, ]
#   gseaTbl[[sncRNA]] <- path1  
# }
#gseaTbl <- as.data.frame(do.call(rbind, gseaTbl))

gseaTbl <- data.table::fread(cmd=paste0("gzip -dc ",gsea_outfile), sep = "\t", header = TRUE, check.names = FALSE, quote="", stringsAsFactors = FALSE, data.table = FALSE)
gseaTbl <- gseaTbl[gseaTbl$sig == 1, ]

SncGSEA2gmt <- function(gseaTbl) {
  tmp <- gseaTbl[, c("ID", "source", "sncRNA")]
  colnames(tmp) <- c("pathway", "source", "gene")
  
  tmp2 <- tmp %>%
    pivot_wider(id_cols = c("pathway", "source"), names_from = gene, values_from = gene) %>%
    as.data.frame()
  
  tmp2[is.na(tmp2)] <- ""
  return(tmp2)
}

if((nrow(tmp.res)>0)){
  # Save results in GMT format
  gmt_outfile <- paste0(outPrefix,"_",dataset,"_",cell,".gmt")
  gmtLong_outfile <- paste0(outPrefix,"_",dataset,"_",cell,".gmtLong")
  tmp2 <- SncGSEA2gmt(gseaTbl)
  message("write gmt table: ",date())
  data.table::fwrite(tmp2, file = gmt_outfile, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
  data.table::fwrite(gseaTbl[, c("ID", "sncRNA", "source")], file = gmtLong_outfile, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
}
unlink(x = paste0(outPrefix, "_",dataset,"_",cell,"_", tmp.name, ".txt") ) #
#unlink(x = paste0(outPrefix, "_",cellID, "_",method,"_", "GSEA.KEGG", "_", tmp.name, ".txt") ) #
