# ref: https://jserizay.com/VplotR/articles/VplotR.html
#BiocManager::install("VplotR")
# library("VplotR")
# 
# 
# data(MNase_sacCer3_Henikoff2011)
# data(ABF1_sacCer3)
# summary(width(ABF1_sacCer3))
# p <- plotVmat(
#   x = MNase_sacCer3_Henikoff2011,
#   granges = ABF1_sacCer3
# )
# p
# summary(width(reg.gr))
# table(width(reg.gr)>=100)

# # Footprints
# #VplotR also implements a function to profile the footprint from MNase or ATAC-seq over sets of genomic loci. For instance, CTCF is known for its ~40-bp large footprint at its binding loci.
# p <- plotFootprint(
#   list_params[["AGO"]][[1]],
#   list_params[["AGO"]][[2]]
#   
# )
# p
# # Local fragment distribution (IGV eg.)
# #VplotR provides a function to plot the distribution of paired-end fragments over an individual genomic window.
# data(MNase_sacCer3_Henikoff2011_subset)
# genes_sacCer3 <- GenomicFeatures::genes(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene::
#                                           TxDb.Scerevisiae.UCSC.sacCer3.sgdGene
# )
# p <- plotProfile(
#   fragments = MNase_sacCer3_Henikoff2011_subset,
#   window = "chrXV:186,400-187,400", 
#   loci = ABF1_sacCer3, 
#   annots = genes_sacCer3,
#   min = 20, max = 200, alpha = 0.1, size = 1.5
# )
# p


# # Plotting fragment size distribution
# #The distribution of fragment sizes can be computed with the getFragmentsDistribution() function:
# df <- getFragmentsDistribution(
#   fragments_from_bed,
#   reg.gr
# )
# p <- ggplot(df, aes(x = x, y = y)) + 
#   geom_line() + 
#   theme_ggplot2() + 
#   labs(x = "Fragment size", y = "# of fragments")
# p



library(ggplot2)
library(VplotR)
#notes: 
#VplotR will recenter bed/prom/loci and extend/expand to region width ylim[2]-ylim[1]  
#VplotR only plot frag center in vplot



# test cfRNA on GSE71008_NC_pool hg38_tx_subsample_bam ------------------------------
## prepare reads/frag bam

# bamfile <- "/BioII/lulab_b/baopengfei/tmp/GSE71008_NCpool_11RNA_dwnsmp.bam"
# fragments <- importPEBamFiles(
#   bamfile
#   # shift_ATAC_fragments = TRUE
# )
# fragments
bam_bed <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008_NCpool/call_peak_all/tbed_11RNA_primary/NCpool.bed.gz"
fragments_from_bed <- rtracklayer::import (bam_bed)

bed <- list()
bed[["AGO"]] <- "/lulabdata/baopengfei/shared_reference/RBP/splitByRBP/AGO2_tx.bed"
bed[["RBP"]] <- "/BioII/lulab_b/baopengfei/shared_reference/RBP/splitByRBP/RBPs_tx.bed"
bed[["hotspot"]] <- "/BioII/lulab_b/baopengfei/shared_reference/RBP/POSTAR3_hotspot/all_merge_newTxID.bed"
bed[["EV"]] <- "/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/lulab/11RNAsmallDomain_EVenrich.bed"
bed[["G4"]] <- "/BioII/lulab_b/baopengfei/shared_reference/structure/quadratlas_g4grinder_tx.bed"
bed[["gold"]] <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/tmp/6dst_filter3_11RNA.bed"

reg.gr <- list()
for (i in names(bed)){
  reg.gr[[i]] <- rtracklayer::import (bed[[i]])
  print(length(reg.gr[[i]]))
}


list_params <- list()
for (i in names(bed)){
  list_params[[i]] = list(fragments_from_bed, reg.gr[[i]]) 
}

COLORSCALE_VMAT <- viridis::viridis_pal()(100) #colorRampPalette(c("blue","white","red"))(100)
for(met in c("libdepth+nloci","quantile","zscore","max","pct","none")){
  print(met)
p <- plotVmat(
  list_params, 
  cores = 6,
  nrow = 1, ncol = 6,
  normFun=met, # c("libdepth+nloci","quantile","zscore","max","pct","none")
  # breaks = NULL, 
  colors = COLORSCALE_VMAT,
  xlim = c(-150, 150),
  ylim = c(10, 65),
  return_Vmat = FALSE,
  verbose = 1,
  # main = '', 
  # xlab = 'Dist. from center',
  # ylab = 'Frag. length',
  # key = 'Score'
)
#?plotVmat
p
ggsave(plot=p, filename = paste0("./vplot_",met,".pdf"), width = 16,height = 4)
}



# Footprints
#VplotR also implements a function to profile the footprint from MNase or ATAC-seq over sets of genomic loci. For instance, CTCF is known for its ~40-bp large footprint at its binding loci.
foot.plt.list <- list()
for(i in names(bed)){
  print(i)
  foot.plt.list[[i]] <- plotFootprint(
  list_params[[i]][[1]],
  list_params[[i]][[2]],
  xlim =  c(-150, 150)
)
}
names(foot.plt.list)
foot.plt.list[["gold"]]
tmp <- cowplot::plot_grid(plotlist = foot.plt.list, nrow = 1,  common.legend = T,align = "hv", axis = "b", legend = "right")
ggplot2::ggsave(plot = cowplot::plot_grid(plotlist = tmp, ncol = 1,  common.legend = T, align = "hv", axis = "b"), #  labels = c("RBPs", "", "EV", "", "G4"), rel_widths = c(1,1,1), 
                paste0("./vplot_frag_length_line.pdf"), 
                width=10, height=16) # "_",sample




frag.plt.list <- list()
for(i in names(bed)){
  print(i)
  df <- getFragmentsDistribution(
    list_params[[i]][[1]],
    list_params[[i]][[2]],
    limits = c(0,100),cores = 6
  )
  frag.plt.list[[i]] <- ggplot(df, aes(x = x, y = y)) +
    geom_line() +
    xlim(0,100) +
    # theme_ggplot2() +
    labs(main=i, x = "Frag. size", y = "Frag. num") +
    theme_minimal() +  # base_size=12
    theme(#axis.ticks.x=element_blank(),
      #strip.text.y = element_blank(),
      aspect.ratio = 0.5,
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
}
# cowplot::plot_grid(plotlist = list(p1,p2), nrow = 1, align = "hv", axis = "b")
tmp <- cowplot::plot_grid(plotlist = frag.plt.list, nrow = 1,  common.legend = T,align = "hv", axis = "b", legend = "right")
ggplot2::ggsave(plot = cowplot::plot_grid(plotlist = frag.plt.list, nrow = 1,  common.legend = T, align = "hv", axis = "b"), #  labels = c("RBPs", "", "EV", "", "G4"), rel_widths = c(1,1,1), 
       paste0("./vplot_frag_length_line.pdf"), 
       width=12, height=6) # "_",sample

#




# test on cfDNA on ctcf bind sites  ---------------------------------
#bamfile <- "/BioII/lulab_b/baopengfei/tmp/NC_dwsmp.bam" #system.file("extdata", "ex1.bam", package = "Rsamtools")
#samtools view -bs 0.001 /BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/bam-sorted-deduped-merged/NC-passQCstringent.bam > /BioII/lulab_b/baopengfei/tmp/NC_dwsmp.bam
bamfile <- "/BioII/lulab_b/baopengfei/tmp/GSE71008_NCpool_11RNA_dwnsmp.bam"
fragments <- importPEBamFiles(
  bamfile
  # shift_ATAC_fragments = TRUE
)
fragments
#bam_bed <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008_NCpool/call_peak_all/tbed_11RNA_primary/NCpool.bed.gz"
#fragments_from_bed <- rtracklayer::import (bed)

# preapare peak bed
# # ctcf <- read.table("/lulabdata/baopengfei/shared_reference/TF/CTCF_CTCFBSDB1.0_allexp.txt",sep = "\t", row.names = F, header = F,stringsAsFactors = F, check.names = F)
# ctcf <- data.table::fread("/lulabdata/baopengfei/shared_reference/TF/CTCF_CTCFBSDB1.0_allexp.txt",sep = "\t", header = T,stringsAsFactors = F, check.names = F)
# ctcf <- ctcf[ctcf$Species=="Human",]
# #table(ctcf$Species)
# #ctcf <- ctcf[1:500,]
# ctcf$chr <- unlist(sapply(strsplit(ctcf$`chromosome Location`,":",fixed=T),"[",1))
# ctcf$start <- unlist(sapply(strsplit(ctcf$`chromosome Location`,":",fixed=T),"[",2))
# ctcf$end <- unlist(sapply(strsplit(ctcf$start,"-",fixed=T),"[",2))
# ctcf$start <- unlist(sapply(strsplit(ctcf$start,"-",fixed=T),"[",1))
# ctcf$score <- 1
# ctcf$strand <- "*"
# ctcf$width <- as.numeric(ctcf$end)-as.numeric(ctcf$start)
# reg.gr <- makeGRangesFromDataFrame (as.data.frame(ctcf[,c("chr","start","end","chromosome Location","score","strand")]))
