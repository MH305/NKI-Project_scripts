library("systemPipeR")
library(GenomicFeatures)
library(ggplot2)
library(grid)
setwd("/Users/mashrurhaidernew/Desktop/NKI/Bam/")
SS = 19
LS = 34
annotation_file <- ("/Users/mashrurhaidernew/Desktop/NKI/gencode.vM21.annotation.gff3.gz")
txdb <- makeTxDbFromGFF(file = annotation_file, format = "gff3", organism = "Mus musculus")
feat <- genFeatures(txdb, featuretype = "all", reduce_ranges = TRUE, 
                      upstream = 1000, downstream = 0, verbose = TRUE)
file_list <- list.files("/Users/mashrurhaidernew/Desktop/NKI/Bam/", pattern = '*.bam$')
#Count and plot reads of any length
listInputBam <- c()
for (i in file_list){
  srr_id <- substr(i,15,24)
  listInputBam <- c(BamFile(file = paste("/Users/mashrurhaidernew/Desktop/NKI/Bam/",i,sep = '/'),index = paste(paste("/Users/mashrurhaidernew/Desktop/NKI/Bam/",i,sep = '/'),'.bai',sep = '')))
  fc <- featuretypeCounts(bfl=BamFileList(listInputBam, yieldSize=50000), grl=feat,
                          singleEnd=TRUE, readlength=NULL, type="data.frame")
  p <- plotfeaturetypeCounts(x=fc, graphicsfile=c(paste("/Users/mashrurhaidernew/Desktop/NKI/SPR_results/",srr_id,"_FC_",as.character(SS),"-",as.character(LS),"_all.pdf",sep=""), graphicsformat="pdf", scales="free",
                               anyreadlength=TRUE, scale_length_val=NULL))
}



  

  
  
  
  

