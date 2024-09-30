
# not finished !!!!


#keep an eye on targe gene universe !!!
options(stringsAsFactors = F)

### add mir rna target label  (multiMiR)
library(multiMiR)
tmp <- multiMiR::get_multimir(mirna = unique(mirVec))@data #eg.: c("MIMAT0000072","hsa-miR-129-5p","hsa-miR-1-3p","hsa-let-7i-3p","hsa-miR-548w")
table(tmp$support_type)
table(tmp$type)
##better only keep "Functional MTI" ?: number seem suitable if filter
#tmp <- tmp[tmp$support_type=="Functional MTI",]

mirVec <- unlist(sapply(strsplit(tx$V4,"|",fixed = T),"[",1))
write.table(tmp,"/BioII/lulab_b/baopengfei/shared_reference/mirbase/hsa_target.txt",quote = F,sep = "\t",row.names = F,col.names = T)



# ### add mir rna target label  (miRanda)
# miRanda.rna <- read.table("/BioII/lulab_b/baopengfei/biosoft/miRanda-1.9-i686-linux-gnu/mir-rna.filter.txt",header = F,sep = "\t")
# miRanda.rna$cmb1 <- paste0(miRanda.rna$V1,":",miRanda.rna$V2)
# miRanda.rna$cmb2 <- paste0(miRanda.rna$V2,":",miRanda.rna$V1)
# k <- paste0(my_cor_matrix$source,":",my_cor_matrix$target)
# my_cor_matrix$relation[(k %in% unique(miRanda.rna$cmb1)) | (k %in% unique(miRanda.rna$cmb2))] <- "miRNA-mRNA.UTR3"
# table(my_cor_matrix$relation). # 173 true
# colnames(miRanda.rna) <- c("miRNA","mRNA","Tot Score","Tot Energy","Max Score","Max Energy","Strand","Len1","Len2","Positions", "cmb1","cmb2")




#############################################
# target function annotation  ---------------
#############################################

# tmp <- as.data.frame(mygene::queryMany(annotation_row_merge$enst[enst.idx],scopes="ensembl.transcript",fields=c("ensembl.gene","symbol"),species="human"))
# dim(tmp)
# table(duplicated(tmp$ensembl.gene))
tmp <- read.table("/BioII/lulab_b/baopengfei/shared_reference/mirbase/hsa_target.txt", header = T, sep = "\t", check.names = F, quote="")
ensgVec <- tmp$target_ensembl

library(dplyr)
library(biomaRt)
#mart <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl") # need online network
#saveRDS(object = mart, file = "/BioII/lulab_b/baopengfei/shared_reference/ensembl/biomaRt.ensembl.hsapiens_gene_ensembl.rds")
mart <- readRDS("/BioII/lulab_b/baopengfei/shared_reference/ensembl/biomaRt.ensembl.hsapiens_gene_ensembl.rds")
ann <- biomaRt::getBM(c("hgnc_symbol","description","chromosome_name","band","strand","start_position","end_position","ensembl_gene_id","transcript_biotype"),
                      "ensembl_gene_id", unique(ensgVec), mart, useCache = FALSE)
dim(ann)
table(duplicated(ann$hgnc_symbol))
table(duplicated(ann$ensembl_gene_id))
table(ann$transcript_biotype)
#protein_coding: 19731
#protein_coding_CDS_not_defined: 10608
#retained_intron: 10291
#nonsense_mediated_decay: 7695
#lncRNA: 294
#other: <100
unique(ann$transcript_biotype)

#rm dup by priority
# ann$transcript_biotype2 <- "other"
# ann$transcript_biotype2[ann$transcript_biotype=="protein_coding"] <- "mRNA"
ann2 <- ann[ann$transcript_biotype=="protein_coding",]

tmp$func <- ann$description[match(tmp$ensembl.gene,ann$ensembl_gene_id)]

