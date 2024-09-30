
# prep snc-mRNA bed6 --------------------
#pre <- "/data2/lulab1/bpf/projects/WCHSU-FTC"
diffpre="/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output"
pre="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
dst="TCGA_small7" # "SLE"
dedup="_all"

# prep sncRNA within reduced.enh
bedtools intersect -s -u -wa \
-a <(cat $pre/output/$dst/call_peak${dedup}/cfpeakCNN/b5_d50_p1_8DNA_gn.bed) -b /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/bed/enhancer.stranded.newTxID.reduced2.bed6 | grep -v "," | cut -f1-6 > \
$pre/output/$dst/call_peak${dedup}/cfpeakCNN/b5_d50_p1_gn_interEnh.bed6
chr_size="/BioII/lulab_b/baopengfei/shared_reference/hg38/genome.chr.sizes"
bedtools slop -s -l 3000 -r 3000 -g ${chr_size} \
-i $pre/output/$dst/call_peak${dedup}/cfpeakCNN/b5_d50_p1_gn_interEnh.bed6 > \
$pre/output/$dst/call_peak${dedup}/cfpeakCNN/b5_d50_p1_gn_interEnh.bed6.ext


# prep filtered sncRNA with target mRNA
cellID="BRCA" # "Naive.B"
surfix="" # _rmMonthOperator
sncRNA.mRNA <- data.table::fread(cmd=paste0("gzip -dc /BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/",dst,"/pcc_mRNA_enhancer/",cellID,"_pearson",surfix,".txt.gz"), data.table = F,sep = '\t',header = T,stringsAsFactors = F,check.names = F) # top10% of all sncRNA

ref <- data.table::fread("/BioII/lulab_b/baopengfei/shared_reference/geneset/encodev27_gene_v2.txt",header = T,data.table = F)
ref.snc <- data.table::fread(paste0(pre,"/output/",dst,"/call_peak",dedup,"/cfpeakCNN/b5_d50_p1_8DNA_gn.bed"),header = F,data.table = F)

sncRNA <- data.table::fread( paste0(pre,"/output/",dst,"/call_peak",dedup,"/cfpeakCNN/b5_d50_p1_gn_interEnh.bed6.ext"),data.table = F )
nrow(sncRNA)
sncRNA.mRNA.filter <- sncRNA.mRNA[sncRNA.mRNA$padj<0.1 & sncRNA.mRNA$cor>0.2,]
nrow(sncRNA.mRNA.filter)
#length(unique(sncRNA.mRNA$sncRNA)) # 30
#length(unique(sncRNA.mRNA.filter$sncRNA)) # 15
#length(unique(sncRNA.mRNA.filter$mRNA)) # 1595
sncRNA.mRNA.filter <- sncRNA.mRNA.filter[sncRNA$V4 %in% sncRNA.mRNA.filter$sncRNA,]
sncRNA.mRNA.filter <- sncRNA.mRNA.filter[order(sncRNA.mRNA.filter$pvalue,decreasing = F),]
#colnames(sncRNA.mRNA)

# change to gn coord
sncRNA.mRNA.filter$mRNA.chr <- ref$chr[match(sncRNA.mRNA.filter$mRNA,ref$ensg)]
sncRNA.mRNA.filter$mRNA.start <- ref$start[match(sncRNA.mRNA.filter$mRNA,ref$ensg)]
sncRNA.mRNA.filter$mRNA.end <- ref$end[match(sncRNA.mRNA.filter$mRNA,ref$ensg)]
sncRNA.mRNA.filter$sncRNA.chr <- ref.snc$V1[match(sncRNA.mRNA.filter$sncRNA,ref.snc$V4)]
sncRNA.mRNA.filter$sncRNA.start <- ref.snc$V2[match(sncRNA.mRNA.filter$sncRNA,ref.snc$V4)]
sncRNA.mRNA.filter$sncRNA.end <- ref.snc$V3[match(sncRNA.mRNA.filter$sncRNA,ref.snc$V4)]
#sncRNA.mRNA.filter <- sncRNA.mRNA.filter[,c("mRNA.chr","mRNA.start","mRNA.end","sncRNA.chr","sncRNA.start","sncRNA.end")]
table(duplicated(sncRNA.mRNA.filter))
#sncRNA.mRNA.filter <- sncRNA.mRNA.filter[!duplicated(paste0(sncRNA.mRNA.filter$mRNA.chr,sncRNA.mRNA.filter$mRNA.start,sncRNA.mRNA.filter$mRNA.end)),]
#sncRNA.mRNA.filter <- sncRNA.mRNA.filter[!duplicated(paste0(sncRNA.mRNA.filter$sncRNA.chr,sncRNA.mRNA.filter$sncRNA.start,sncRNA.mRNA.filter$sncRNA.end)),]
sncRNA.mRNA.filter <- sncRNA.mRNA.filter[sncRNA.mRNA.filter$mRNA.chr==sncRNA.mRNA.filter$sncRNA.chr,]
summary(abs(sncRNA.mRNA.filter$mRNA.start-sncRNA.mRNA.filter$sncRNA.start))
#table(abs(sncRNA.mRNA.filter$mRNA.start-sncRNA.mRNA.filter$sncRNA.start)<2000000)
distan <- abs(sncRNA.mRNA.filter$mRNA.start-sncRNA.mRNA.filter$sncRNA.start)
sncRNA.mRNA.filter <- sncRNA.mRNA.filter[(distan<2000000) & (distan>3000),]
sncRNA.mRNA.filter <- na.omit(sncRNA.mRNA.filter)
#str(sncRNA.mRNA.filter)
options(bedtools.path = "/BioII/lulab_b/baopengfei/biosoft")
library(bedtoolsr)
sncRNA.mRNA.filter[1:3,]
sncRNA.mRNA.filter.ext1 <- bedtoolsr::bt.slop(sncRNA.mRNA.filter[,c("mRNA.chr","mRNA.start","mRNA.end")],g='/BioII/lulab_b/baopengfei/shared_reference/hg38/hg38.chrom.sizes',l = 1000,r=1000)
sncRNA.mRNA.filter.ext2 <- bedtoolsr::bt.slop(sncRNA.mRNA.filter[,c("sncRNA.chr","sncRNA.start","sncRNA.end")],g='/BioII/lulab_b/baopengfei/shared_reference/hg38/hg38.chrom.sizes',l = 1000,r=1000)

data.table::fwrite(cbind(sncRNA.mRNA.filter.ext1,sncRNA.mRNA.filter.ext2),paste0(diffpre,"/",dst,"/pcc_mRNA_enhancer/",cellID,"_sncRNA_mRNA.bed6"),sep = '\t',quote = F,row.names = F,col.names = F)
#



# plot Hi-C --------------------
# BiocManager::install("HiCExperiment", ask = FALSE)
# BiocManager::install("HiCool", ask = FALSE)
# BiocManager::install("HiContacts", ask = FALSE)
# BiocManager::install("HiContactsData", ask = FALSE)
# BiocManager::install("fourDNData", ask = FALSE)
# BiocManager::install("DNAZooData", ask = FALSE)

library(dplyr)
library(ggplot2)
library(HiCExperiment)
library(HiContacts)
library(HiContactsData)
library(rtracklayer)
library(InteractionSet)


#pre <- "/data2/lulab1/bpf/projects/WCHSU-FTC"
diffpre="/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output"
pre="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
dst="TCGA_small7" # "SLE"
dedup="_all"
cellID="BRCA"

focus_region <- "chr16" # chr22:22000000-25000000|chr2:22000000-25000000
res <- 10000
chrSize <- "/BioII/lulab_b/baopengfei/shared_reference/hg38/genome.chr.sizes"


#/BioII/lulab_b/baopengfei/shared_reference/3D_genome/YueFengLab/hg38.loops/Rao_2014.GM12878.hg38.peakachu-merged.loops.chr22
loops <- paste0(diffpre,"/",dst,"/pcc_mRNA_enhancer/",cellID,"_sncRNA_mRNA.bed6") |>  # chr16
  HiCExperiment::import(format="bedpe") |> 
  InteractionSet::makeGInteractionsFromGRangesPairs() 
length(loops)
# table(loops@regions@seqnames)
loops <- loops[loops@regions@seqnames=="chr1"]
length(loops)

borders <- "/BioII/lulab_b/baopengfei/shared_reference/3D_genome/YueFengLab/hg38.TAD/GM12878_Rao_2014-raw_TADs.txt" |> 
  HiCExperiment::import(format="bed")

loops.df <- as.data.frame(loops)
tmp1 <- loops.df[,1:5]
tmp2 <- loops.df[,6:10]
colnames(tmp1) <- c("seqnames","start","end","width","strand")
colnames(tmp2) <- c("seqnames","start","end","width","strand")
loops.df2 <- as.data.frame(rbind(tmp1,tmp2))
loops.df2[1:3,]
#loops.df2 <- GRanges(loops.df2)
loops.df2 <- bedtoolsr::bt.slop(i = loops.df2, g = chrSize, l = res * 20, r = res * 20) # + res * 20: extend each side by res*20 ?
loops.df2 <- bedtoolsr::bt.sort(i = loops.df2)
loops.df2 <- bedtoolsr::bt.merge(i = loops.df2)
loops.df2$pos <- paste0(loops.df2$V1,":",loops.df2$V2,"-",loops.df2$V3)
focus_region <- paste0(loops.df2$pos[1:5],collapse = "|")
#focus_region <- "chr16:13507945-82492142"
#[1] "chr1:155409206-155815967|chr1:156378728-156937960|chr1:157174303-157576316|chr1:158343551-158755405|chr1:158861485-159348096"

#4DNFIIL851T7.mcool # GM***, lymp (53M)
#4DNFIL76YMY6.mcool  # GM***, lymp (13G)
#4DNFIOH4T1RI.mcool # Hela, (500M)
#cool_file <- "/BioII/lulab_b/baopengfei/shared_reference/3D_genome/4Dnucleome/Hela/4DNFIO8HVKOL.hic" # 
cool_file <- "/BioII/lulab_b/baopengfei/shared_reference/3D_genome/4Dnucleome/Hela/4DNFIOH4T1RI.mcool" # 
#cool_file <- "/BioII/lulab_b/baopengfei/mambaforge/envs/HiC/lib/R/library/HiCExperiment/extdata/S288C-borders.bed" #HiContactsData('yeast_wt', format = 'cool')
#> see ?HiContactsData and browseVignettes('HiContactsData') for documentation
#> loading from cache

#cooler dump -t pixels --header --join -r chr19 output.cool /BioII/lulab_b/baopengfei/shared_reference/3D_genome/4Dnucleome/Hela/4DNFIL76YMY6.mcool


#The GM12878 lymphoblastoid cell line is derived from B lymphocytes, which are a type of white blood cell responsible for producing antibodies. Among the cell types listed in your query, the ones most similar to B lymphocytes include:
#Naive.B
hic <- HiCExperiment::import(cool_file, format = 'cool', resolution = res, focus = focus_region) # , focus = focus_region
hic <- HiCExperiment::import("/BioII/lulab_b/baopengfei/shared_reference/3D_genome/Rao_2014_cell/GSE63525_HMEC_combined_30.hic", format = 'hic',resolution = 50000) # , focus = focus_region
hic

# m <- matrix(rpois(100, 5), 10, 10)
# HiCcompare::KRnorm(m)


## Horizontal matrix
HiContacts::plotMatrix(
  # hic,
  HiCExperiment::refocus(hic, focus_region),
  use.scores = 'balanced', 
  limits = c(-4, -1), # NULL
  # maxDistance = 2500000, # rect: NULL, horizon: 200000
  scale = "log10",
  loops = loops,
  borders = borders,
  dpi = 200,
  rasterize = TRUE,
  chrom_lines = TRUE,
  show_grid = T, # FALSE
  cmap = afmhotrColors(), # NULL
  caption = TRUE
)
ggsave("HiC_rect.pdf")
HiContacts::plotMatrix(
  # hic,
  HiCExperiment::refocus(hic, focus_region),
  use.scores = 'balanced', 
  limits = c(-4, -1), # NULL
  maxDistance = 2500000, # rect: NULL, horizon: 200000
  scale = "log10",
  loops = loops,
  borders = borders,
  dpi = 200,
  rasterize = TRUE,
  chrom_lines = TRUE,
  show_grid = T, # FALSE
  cmap = afmhotrColors(), # NULL
  caption = TRUE
)
ggsave("HiC_hori.pdf")


# #mcool_file <- HiContactsData('yeast_wt', format = 'mcool')
# #mcool_file <- cool_file
# #> see ?HiContactsData and browseVignettes('HiContactsData') for documentation
# #> loading from cache
# #> /BioII/lulab_b/baopengfei/shared_reference/3D_genome/YueFengLab/hg38.loops/Rao_2014.GM12878.hg38.peakachu-merged.loops
# loops <- system.file("extdata", 'S288C-loops.bedpe', package = 'HiCExperiment') |> 
#   import() |> 
#   makeGInteractionsFromGRangesPairs()
# p <- import(mcool_file, format = 'mcool', focus = 'IV') |> 
#   plotMatrix(loops = loops, limits = c(-4, -1), dpi = 120)
# p
# 
# #/BioII/lulab_b/baopengfei/shared_reference/3D_genome/YueFengLab/hg38.TAD/GM12878_Rao_2014-raw_TADs.txt
# borders <- system.file("extdata", 'S288C-borders.bed', package = 'HiCExperiment') |> 
#   import()
# p <- import(mcool_file, format = 'mcool', focus = 'IV') |> 
#   plotMatrix(loops = loops, borders = borders, limits = c(-4, -1), dpi = 120)
# p
# aggr_centros <- HiContacts::aggregate(
#   hic, targets = loops, BPPARAM = BiocParallel::SerialParam()
# )

#hic <- HiCExperiment::zoom(hic, 1000)
aggr_loops <- aggregate(hic, targets = loops, BPPARAM = BiocParallel::SerialParam(), flankingBins = 15) # 15*res
##  Going through preflight checklist...
##  Parsing the entire contact matrice as a sparse matrix...
##  Modeling distance decay...
##  Filtering for contacts within provided targets...
#aggr_loops

#> Going through preflight checklist...
#> Parsing the entire contact matrice as a sparse matrix...
#> Modeling distance decay...
#> Filtering for contacts within provided targets...
# slices(aggr_loops)
##  List of length 4
##  names(4): count balanced expected detrended
# dim(slices(aggr_loops, 'count'))
##  [1]  31  31 148
# topologicalFeatures(aggr_loops, 'targets')
plotMatrix(
    aggr_loops, 
    use.scores = 'detrended', 
    scale = 'linear', 
    limits = c(-1, 1), 
    cmap = bgrColors()
)
ggsave("HiC_aggre_sncRNA.pdf")

