#!/usr/bin/env Rscript

# run clusterProfiler with sncRNA Gmt

#Note:
#* each cell/tissue use one *_sncRNA.gmtLong
#* topN = topN up + top N down
#* fig get basic top 20 pathway enriched as eg 
#* only offline supported for self-built gmt

options(stringsAsFactors = F)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(DOSE))
suppressPackageStartupMessages(library(enrichplot))
# suppressPackageStartupMessages(library(org.Hs.eg.db))

parser <- ArgumentParser(description='GO/KEGG/DO geneset enrichment by ORA/GSEA using diff table (network independent version)')
parser$add_argument('-i', '--matrix', type='character', required=TRUE,
                    help='input diff matrix. Rows are genes (^ENSG*, *_ENSG*). Columns are deseq2 tags')
parser$add_argument('--cutoffP', type='double', default=0.05,
                    help='Over-Representation Analysis cutoff pvalue (default=0.05, genes with pvalue <= 0.05 will be treated as top rank genes)')
parser$add_argument('--cutoffLog2FC', type='double', default=0.6,
                    help='Over-Representation Analysis cutoff abs(log2FC) (default=0.6, genes with abs(log2FC) >= 0.6 will be treated as top rank genes)')
parser$add_argument('--cutoffPadj', type='double', default=1,
                    help='Over-Representation Analysis cutoff padj (default=1, genes with padj <= 1 will be treated as top rank genes)')
parser$add_argument('--cutoffTop', type='logical', default=FALSE,
                    help='whether to select Over-Representation Analysis top gene by abs(log2FC) or pvalue rank (default=FALSE)')
parser$add_argument('--cutoffTopP', type='double', default=100,
                    help='Over-Representation Analysis top pvalue (default=100, genes with pvalue rank <= 100 will be treated as top rank genes)')
parser$add_argument('--cutoffTopFC', type='double', default=100,
                    help='Over-Representation Analysis top abs(log2FC) (default=100, genes with abs(log2FC) rank <= 100 will be treated as top rank genes)')
parser$add_argument('-p', '--pValue', type='double', default=1,
                    help='pValue on all enrichment tests to be saved in output file (default=1)')
parser$add_argument('-a', '--adjustPval', type='character', default="BH",
                    help='methods to adjust pValue, one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none". (default="BH")')
parser$add_argument('-q', '--qValue', type='double', default=1,
                    help='qValue (another FDR control method) on all enrichment tests to be saved in output file (default=1)')
parser$add_argument('--GmtLong', type='character', help='Gmt Long format tsv (path gene source)')
parser$add_argument('--PathSource', type='character', default="CHG,ImmReg,KEGG,MSigDB,Reactome", help='Gmt Long format tsv source selected, default: CHG,ImmReg,KEGG,MSigDB,Reactome')
parser$add_argument('-o', '--outdir', type='character', default="./",
                    help='outdir (default=./)')
parser$add_argument('-f', '--fig', type='logical', default=TRUE,
                    help='whether to plot figures (default=TRUE)')
parser$add_argument('--online', type='logical', default=TRUE,
                    help='whether to use online real-time updated KEGG, offline will use pre-saved PATH_ID_NAME_KEGGplusHallmark.txt file (default=TRUE)')
args <- parser$parse_args()
matrix <- args$matrix
cutoff.p <- args$cutoffP
cutoff.padj <- args$cutoffPadj
cutoff.fc <- args$cutoffLog2FC
cutoff.Top <- args$cutoffTop
cutoff.TopP <- args$cutoffTopP
cutoff.TopFC <- args$cutoffTopFC
pValue <- args$pValue
adjustPval <- args$adjustPval
qValue <- args$qValue
outdir <- args$outdir
fig <- args$fig
online <- args$online
GmtLong <- args$GmtLong
PathSource <- args$PathSource

#set.seed(1234)

# # # test TCGA
# setwd("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/")
# matrix <- "output/TCGA_small7/diff_BRCA_YesvsNo_rmPurity.diff"
# cutoff.p <- 0.05
# cutoff.padj <- 1
# cutoff.fc <- 0.5
# cutoff.Top <- T
# cutoff.TopP <- 500
# cutoff.TopFC <- 100
# pValue <- 1
# adjustPval <- "BH"
# qValue <- 1
# outdir <- "./test"
# fig <- TRUE
# GmtLong <- "output/TCGA_small7/pcc_mRNA/rmPurity_BRCA_pearson_sncRNA.gmtLong"
# PathSource <- "CHG,ImmReg,KEGG,MSigDB"

## prepare mess.
message(paste0("matrix:",matrix))
message(paste0("cutoff.p:",cutoff.p))
message(paste0("cutoff.padj:",cutoff.padj))
message(paste0("cutoff.fc:",cutoff.fc))
message(paste0("cutoff.top:",cutoff.Top))
message(paste0("cutoff.topP:",cutoff.TopP))
message(paste0("cutoff.topFC:",cutoff.TopFC))
message(paste0("pValue:",pValue))
message(paste0("adjustPval:",adjustPval))
message(paste0("qValue:",qValue))
message(paste0("outdir:",outdir))
message(paste0("fig:",fig))

if(adjustPval=="none")
{
  message("not adjust pvalue")
} else {
  message(paste0("adjust pvalue by ",adjustPval))
}



## print gene num pass cutoff based on pvalue
dir.create(outdir,recursive = T, showWarnings = F)
res <- read.table(matrix,stringsAsFactors = F,header = T,row.names = 1)


message(paste0('ORA enrich selects diff records:'), sum(res$pvalue<=cutoff.p & res$padj<=cutoff.padj & abs(res$log2FoldChange)>=cutoff.fc, na.rm = T),
        " up:",sum(res$pvalue<=cutoff.p & res$padj<=cutoff.padj & res$log2FoldChange>=cutoff.fc, na.rm = T),
        " down:",sum(res$pvalue<=cutoff.p & res$padj<=cutoff.padj & res$log2FoldChange<=-cutoff.fc, na.rm = T))

## tidy input mat
resOrder <- res[order(res$log2FoldChange,-res$pvalue,decreasing = T),]  # order rows by logFC and then pVal
resOrder$ENSG <- rownames(resOrder)
#resOrder$ENSG <- as.character(lapply(strsplit(resOrder$ENSG,".",fixed = TRUE),function(x) x[1]))

if (length(grep("^ENSG",rownames(resOrder)))>=1) {
	message("input gene id is ^ENSG*, use 1st as id, 3rd as class")
  snc <- F
	#resOrder$ENSG <- resOrder$ENSG
  resOrder$class <- as.character(lapply(strsplit(resOrder$ENSG,"|",fixed = TRUE),function(x) x[3]))
	resOrder$ENSG <- as.character(lapply(strsplit(resOrder$ENSG,"|",fixed = TRUE),function(x) x[1]))
} else if(length(grep("peak|block|cluster",rownames(resOrder),perl = T))>=1){
  message("input gene id is peak_*, use 4th as id, 7th as class")
  snc <- T
  resOrder$class <- as.character(lapply(strsplit(resOrder$ENSG,"|",fixed = TRUE),function(x) x[7]))
  resOrder$ENSG <- as.character(lapply(strsplit(resOrder$ENSG,"|",fixed = TRUE),function(x) x[4]))
}else {
  message("input gene id is not ^ENSG*")
  snc <- F
  resOrder$class <- "all"
	resOrder$ENSG <- as.character(lapply(strsplit(resOrder$ENSG,"_",fixed = TRUE),function(x) x[2]))
}
resOrder <- resOrder[!duplicated(resOrder$ENSG),]
#gene.list.ncbiid <- bitr(resOrder$ENSG, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# message("gene class: ",paste0(table(resOrder$class)," "))

#ref <- read.table("/BioII/lulab_b/baopengfei/shared_reference/hg38/gene.gtf",sep = "\t",header = T)
# library(gprofiler2)
# tmp <- gprofiler2::gconvert(query = ref$ensg, organism = "hsapiens", target = "ENTREZGENE_ACC" )
# tmp <- tmp[,c("input","target")]
# colnames(tmp) <- c("ensg","ENTREZID")
# ref$ENTREZID <- tmp$ENTREZID[match(ref$ensg,tmp$ensg)]
# ref <- ref[!duplicated(ref$ensg),]
# table(duplicated(ref$ENTREZID))
# write.table(ref,"/BioII/lulab_b/baopengfei/shared_reference/geneset/encodev27_gene.txt",quote = F,sep = "\t",col.names = T,row.names = F)
if(!snc){
  gene.list.ncbiid <- read.table("/BioII/lulab_b/baopengfei/shared_reference/geneset/encodev27_gene_v2.txt",sep = "\t",header = T)
  gene.list.ncbiid <- gene.list.ncbiid[!is.na(gene.list.ncbiid$ENTREZID),]
  gene.list.ncbiid <- gene.list.ncbiid[,c("ensg","ENTREZID")]
  gene.list.ncbiid <- gene.list.ncbiid[!duplicated(gene.list.ncbiid$ensg) & !duplicated(gene.list.ncbiid$ENTREZID),]
  colnames(gene.list.ncbiid)[1] <- "ENSG"
  
  resOrder <- dplyr::left_join(resOrder, gene.list.ncbiid)
}else{
  resOrder$ENTREZID <- resOrder$ENSG
}

resOrder <- resOrder[!duplicated(resOrder$ENTREZID) & !duplicated(resOrder$ENSG),]
resOrder <- resOrder[!is.na(resOrder$ENTREZID) & !is.na(resOrder$ENSG) & !is.na(resOrder$log2FoldChange),]
resOrder <- resOrder[order(resOrder$log2FoldChange,-resOrder$pvalue,decreasing = T),]  # order rows by logFC (no abs), then pVal 

# get all geneList (for GSEA)
gene.list.ncbiid <- resOrder$log2FoldChange
names(gene.list.ncbiid) <- resOrder$ENTREZID
gene.list.ensgid <- resOrder$log2FoldChange
names(gene.list.ensgid) <- resOrder$ENSG
message(paste0("remaining geneList lenth: ",length(gene.list.ensgid)))

## read gmt_long reference
if(!snc){
  m_t <- read.table(GmtLong,sep = "\t",header = T, stringsAsFactors = F,row.names = 1)
  #colnames(m_t)
  #m_t2e <- m_t2e[!duplicated(m_t2e$short_DESCRPTION),]
}else{
  m_t <- read.table(GmtLong,sep = "\t",header = F, stringsAsFactors = F)
  colnames(m_t) <- c("short_DESCRPTION","ensembl_gene_id","source")
  m_t$ENTREZID <- m_t$ensembl_gene_id
}
PathSource <- strsplit(PathSource,",")[[1]]
m_t <- m_t[m_t$source %in% PathSource,]
#"CHG","ImmReg","KEGG","Manually curated","MSigDB","Reactome"
m_t2g <- m_t[,c("short_DESCRPTION","ensembl_gene_id")]
#m_t2g <- m_t2g[!duplicated(m_t2g$short_DESCRPTION),]
m_t2e <- m_t[,c("short_DESCRPTION","ENTREZID")]


# get top gene (for ORA)
if (cutoff.Top){
  if(cutoff.TopP){
    cutoff.TopP <- min(nrow(resOrder),cutoff.TopP)
    message(paste0("get top gene by pvalue rank: ",cutoff.TopP))
    resOrderP <- resOrder[order(-resOrder$pvalue,abs(resOrder$log2FoldChange),decreasing = T),]  # order rows by pVal, then abs(logFC)
    resOrderP.up <- resOrderP[resOrderP$log2FoldChange>0,]
    resOrderP.down <- resOrderP[resOrderP$log2FoldChange<0,]
    ## get top gene (all)
    gene.top.ncbiid <- resOrderP$ENTREZID[1:cutoff.TopP]   # pvalue
    ## get top gene (up-regulated)
    gene.top.ncbiid.up <- resOrderP.up$ENTREZID[1:cutoff.TopP]
    ## get top gene (down-regulated)
    gene.top.ncbiid.down <- resOrderP.down$ENTREZID[1:cutoff.TopP]
    message(paste0("remaining top up gene lenth: ",length(gene.top.ncbiid.up)))
    message(paste0("remaining top down gene lenth: ",length(gene.top.ncbiid.down)))
  } else {
    message(paste0("get top gene by abs(FoldChange) rank: ",cutoff.TopFC))
    resOrderFC <- resOrder[order(abs(resOrder$log2FoldChange),-resOrder$pvalue,decreasing = T),]  # order rows by abs(logFC), then pVal
    resOrderFC.up <- resOrderFC[resOrderFC$log2FoldChange>0,]
    resOrderFC.down <- resOrderFC[resOrderFC$log2FoldChange<0,]
    ## get top gene (all)
    gene.top.ncbiid <- resOrderP$ENTREZID[1:cutoff.TopFC]   # pvalue
    ## get top gene (up-regulated)
    gene.top.ncbiid.up <- resOrderFC.up$ENTREZID[1:cutoff.TopFC]
    ## get top gene (down-regulated)
    gene.top.ncbiid.down <- resOrderFC.down$ENTREZID[1:cutoff.TopFC]
    message(paste0("remaining top up gene lenth: ",length(gene.top.ncbiid.up)))
    message(paste0("remaining top down gene lenth: ",length(gene.top.ncbiid.down)))
  }
}else{
  message(paste0("get top gene by p, padj and abs(log2FC) value."))
  ## get top gene (all)
  gene.top.ncbiid <- resOrder$ENTREZID[resOrder$pvalue<=cutoff.p & resOrder$padj<=cutoff.padj & abs(resOrder$log2FoldChange)>=cutoff.fc]   # pvalue
  ## get top gene (up-regulated)
  gene.top.ncbiid.up <- resOrder$ENTREZID[resOrder$pvalue<=cutoff.p & resOrder$padj<=cutoff.padj & resOrder$log2FoldChange>=cutoff.fc]
  ## get top gene (down-regulated)
  gene.top.ncbiid.down <- resOrder$ENTREZID[resOrder$pvalue<=cutoff.p & resOrder$padj<=cutoff.padj & resOrder$log2FoldChange<=-cutoff.fc]
  message(paste0("remaining top up gene lenth: ",length(gene.top.ncbiid.up)))
  message(paste0("remaining top down gene lenth: ",length(gene.top.ncbiid.down)))
}



# 1.Over-Representation Analysis
message("start ORA enrich...")
# 1.1 ORA all gene
ORA.KEGG <- enricher(gene.top.ncbiid,
               #universe = names(ccle.gene.spearman$g.l), 
               TERM2GENE=m_t2e,
               pvalueCutoff = pValue,
               pAdjustMethod = adjustPval,
               qvalueCutoff = qValue
               #, maxGSSize = 350
               )

# 1.2 ORA up gene
ORA.KEGG.up <- enricher(gene.top.ncbiid.up,
                        #universe = names(ccle.gene.spearman$g.l), 
                        TERM2GENE=m_t2e,
                        pvalueCutoff = pValue,
                        pAdjustMethod = adjustPval,
                        qvalueCutoff = qValue
                        #, maxGSSize = 350
)

# 1.3 ORA down gene
ORA.KEGG.down <- enricher(gene.top.ncbiid.down,
                        #universe = names(ccle.gene.spearman$g.l), 
                        TERM2GENE=m_t2e,
                        pvalueCutoff = pValue,
                        pAdjustMethod = adjustPval,
                        qvalueCutoff = qValue
                        #, maxGSSize = 350
)


# 2.Gene Set Enrichment Analysis
message("start gsea...")
GSEA.KEGG <- GSEA(gene.list.ncbiid, 
                 TERM2GENE = m_t2e, 
                 nPerm  = 500,
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = pValue,
                 pAdjustMethod = adjustPval
                 )

# 3.write table output
types <- list(ORA.KEGG,ORA.KEGG.up,ORA.KEGG.down,GSEA.KEGG)
names(types) <- c("ORA.KEGG","ORA.KEGG.up","ORA.KEGG.down","GSEA.KEGG")

message("start writing tables...")
    for(i in 1:length(types)){
      tryCatch(                   # keep running when plot error exists
        {
          write.table(types[[i]],paste0(outdir,"/", names(types)[i], ".txt"),row.names = F,col.names = T, quote = F,sep = "\t")
        },
        error = function(e) {
          message("fail export ",types)
        }
      )
    }


# plot figure output
if(fig)
{
# 4.plot
message("start plotting...")
dir.create(paste0(outdir,"/plot/"),showWarnings = F)
## 4.1.barplot for top gene enrich (not informative...)
#barplot(enrich.GO, showCategory = 10)  
#barplot(ORA.KEGG, showCategory = 10) 

## 4.2.dotplot for gene set enrich
  for (i in 1:length(types))
  {
    tryCatch(                   # keep running when plot error exists
      {
        pdf(paste0(outdir,"/plot/dotplot-",names(types)[i],".pdf"))
        print(enrichplot::dotplot(types[[i]], font.size = 9, title=names(types)[i],showCategory = 20))  # must print in for loop !
        dev.off()
      },
      error = function(e) {
        message("fail plot ",types)
      }
    )
  }
message("offline enrich done (KEGG only)")
}

