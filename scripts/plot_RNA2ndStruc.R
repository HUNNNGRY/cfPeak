#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(R4RNA))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(Biostrings))
#install.packages("RRNA")
suppressPackageStartupMessages(library(RRNA))
suppressPackageStartupMessages(library(msa))

parser <- ArgumentParser()

# Add arguments
parser$add_argument("-i","--inFile", help = "Input file")
parser$add_argument("-o","--outPre", help = "Output prefix")
parser$add_argument("--RNAfoldPath", default="/BioII/lulab_b/jinyunfan/anaconda3/envs/bioinfo_py36/bin/RNAfold", help = "Path to RNAfold")

# Parse command line arguments
args <- parser$parse_args()

# Use the arguments in your script
RNAfoldPath <- args$RNAfoldPath
inFile <- args$inFile
outPre <- args$outPre

# # # test
# RNAfoldPath <- "/BioII/lulab_b/jinyunfan/anaconda3/envs/bioinfo_py36/bin/RNAfold"
# #sequences <- read.fasta("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/test.fa")
# #sequences = sequences[1]
# # sequences
# inFile <- "/BioII/lulab_b/baopengfei/test.txt"
# outPre <- "test"

#



## define func.  ------------
getDotBracketEnery <- function(rna_fold_output){
  #"................(((.((((.........))))))) ( -4.10)"
  #".............. (  0.00)"
  library(stringr)
  pattern <- "-?\\d+\\.\\d+" # Define the regular expression pattern to match floating-point numbers
  matches <- str_extract_all(rna_fold_output, pattern) # Use str_extract_all to find all matches of the pattern in the string
  numeric_matches <- as.numeric(unlist(matches)) # Convert the matches to numeric values
  
  pattern2 <- "^(.*?)(?=\\s\\()" # Define the regular expression pattern to extract the desired substring
  dotbracket <- str_extract(rna_fold_output, pattern2) # Use str_extract to extract the substring
  
  return(list(dotbracket,numeric_matches))
}
#



## read in ------------
df <- data.table::fread(inFile,sep = "\t",quote = F,header = T,check.names = F,data.table = F)
#print(colnames(df))
if(!all( c("name","MFE_fasta","peak_fasta") %in% colnames(df) )){
  stop(paste("name","MFE_fasta","peak_fasta","3 least cols not found in input"))
}
df <- df[,c("name","MFE_fasta","peak_fasta")]
colnames(df) <- c("ID","extend_fasta","peak_fasta")
df$extend_fasta <- toupper(df$extend_fasta)
df$extend_fasta <- gsub("T","U",df$extend_fasta)
#df$extend_fasta <- gsub("t","u",df$extend_fasta)
df$peak_fasta <- toupper(df$peak_fasta)
df$peak_fasta <- gsub("T","U",df$peak_fasta)
df <- df[nchar(df$peak_fasta)>10,] # filter too short 
#df$peak_fasta <- gsub("t","u",df$peak_fasta)
#





## plot vienna RNA 2nd structures ------------
#test fa input
#seqFastadna2seq = function(x) {paste(getSequence(x), collapse = "")}
#sequences <- as(vapply(sequences, seqFastadna2seq, character(1)), "DNAStringSet")
#sequences <- as.data.frame(sequences)
message("Start Plot Vienna")
ct.list <- list()
coord.list <- list()
for( i in 1:nrow(df)){
  # i <- 1
  message(i)
  peakID <- df[["ID"]][i]
  ct.list[[peakID]] <- list()
  coord.list[[peakID]] <- list()
  for ( j in c("extend","peak")){
    #j <- "extend"
    message(j)

    vienna2=system2(command = RNAfoldPath, args = c("-d2", "--noLP"), input = df[[paste0(j,"_fasta")]][i], stdout = TRUE) # sequences[["x"]][2]
    vienna2=getDotBracketEnery(vienna2[2])
  
    df[[paste0(j,"_struc")]][i] <- vienna2[[1]]
    df[[paste0(j,"_struc_mfe")]][i] <- vienna2[[2]]
  
    tryCatch({
      ct.list[[peakID]][[j]] <- RRNA::makeCt(df[[paste0(j,"_struc")]][i], df[[paste0(j,"_fasta")]][i])
      coord.list[[peakID]][[j]] <- RRNA::ct2coord(ct.list[[peakID]][[j]])  
    }, error = function(e) {
      message("An error occurred: ", e$message)
      ct.list[[peakID]][[j]] <- NULL
      coord.list[[peakID]][[j]] <- NULL
    })
    # ct.list[[peakID]][[j]]=RRNA::makeCt( df[[paste0(j,"_struc")]][i], df[[paste0(j,"_fasta")]][i])
    # coord.list[[peakID]][[j]]=RRNA::ct2coord(ct.list[[peakID]][[j]])  
  }



# RNAPlot(dat,ranges,labTF=TRUE)
  message(i,df[[paste0("ID")]][i])
  ### Highlight the part within extended sequences that overlapped with peak ###
  dat <- coord.list[[i]][["extend"]]
  pdf(paste0(outPre,"",df[["ID"]][i],".vienna.pdf"))
  tryCatch({
    RNAPlot(dat,
            nt=F, # show dot
            # ranges, # highlight option1: at data.frame with specific color ...
            hl=c(df[["peak_fasta"]][i]),seqcol=c(2), # highlight option2: at specific matched sequence (case-sensitive)
            # modspec=TRUE, modp=c(1:4,43:46),mod=c(17,17,15,15,16,16,16,16),modcol=c(rep(2,2),rep(3,2),rep(4,4)), # highlight option3: show modified position/shape/color point element
            pointSize=1, lineWd = 0.75, # show point/text/line element  
            labTF=F,main=paste0(df[["ID"]][i],"\n",
                                  "peak: ",df[["peak_struc_mfe"]][i]," kcal/mol","\n","extend: ",df[["extend_struc_mfe"]][i]," kcal/mol"
                                  ) # show legend and titles
    )
    RNAPlot(dat,add=TRUE, nt=TRUE, tsize = 0.75, dp=0.75) # add nucleotides, dp: how far from the coordinates nucleotide sequence should be plotted (0-5 work best)
  }, error = function(e) {
    message("An error occurred: ", e$message)
  })
if (!is.null(dev.list())) {
  dev.off()
}


  ## plot Arch 
  message(i,df[[paste0("ID")]][i],df[[paste0("extend","_struc")]][i])
  # i <- 2
  pdf(paste0(outPre,"",df[["ID"]][i],".arch.pdf"))
  tryCatch({
    helix <- R4RNA::viennaToHelix(df[[paste0("extend","_struc")]][i])
      if(nrow(helix)>1){
        #helix.peak <- R4RNA::viennaToHelix(df[[paste0("peak","_struc")]][i])
        #helix$value <- 0.00001 # pval-like
        helix$lty <- 1
        helix$lwd <- 2
        helix$col <- "#4393C3"
        
        R4RNA::plotHelix(helix,flip = F,arrow = F,shape = "circle") 
        #plotDoubleHelix(helix,helix.peak)
        #plotOverlapHelix(helix,helix.peak)
        dev.off()
      }
  }, error = function(e) {
    message("An error occurred: ", e$message)
  })
if (!is.null(dev.list())) {
  dev.off()
}    
}




# #BiocManager::install("msa")
# library(msa)
# 
# q1 <- seqinr::as.SeqFastadna("UGAGUGGUGUUGUUGGCUGCAUUAUGAUGUUGGUUAUAUUCUG")
# q2 <- seqinr::as.SeqFastadna("AAAAUGAGUGGUGUUG")
# q3 <- seqinr::as.SeqFastadna("UUGGCUGCAU")
# q4 <- seqinr::as.SeqFastadna("GUUAUAUUCUGAAAA")
# msa::msa(c(q1,q2,q3,q4), type = "rna", order = "input", method = "ClustalOmega")
# #ClustalOmega seem to be more efficient and accurate for exact matches






# data(helix)
# # Plot helix plain
# plotHelix(known)
# # Apply global appearance options
# plotHelix(known, line = TRUE, arrow = TRUE, col = "blue", lwd = 1.5)
# # Add extra column with styling options
# known$lty <- 1:4
# known$lwd <- 1:2
# known$col <- c(rgb(1, 0, 0), "orange", "yellow", "#00FF00", 4, "purple")
# plotHelix(known)
# # Manually colour helices according to value
# helix$col <- "red"
# helix$col[which(helix$value < 1e-3)] <- "orange"
# helix$col[which(helix$value < 1e-4)] <- "green"
# helix$col[which(helix$value < 1e-5)] <- "blue"
# plotHelix(helix) 
# # Automatically creating a similar plot with legend
# coloured <- colourByValue(helix, log = TRUE, get = TRUE)
# plotHelix(coloured, line = TRUE, arrow = TRUE)
# legend("topleft", legend = attr(coloured, "legend"),
#        fill = attr(coloured, "fill"), title = "P-value", text.col = "black")	
# # Plot both helices with styles
# plotDoubleHelix(helix, known)
# # Overlap helix
# plotOverlapHelix(helix, known)





