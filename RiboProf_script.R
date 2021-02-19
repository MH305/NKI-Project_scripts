library(RiboProfiling)
library(GenomicAlignments)
require(svMisc)
setwd("~/Desktop/NKI/Bam")
file_list <- list.files(path="~/Desktop/NKI/Bam",pattern = '*.bam$')

#txdb object with annotations
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

for (i in file_list) {
  print(i) 
  name<-gsub(x = i,pattern = 'accepted_hits_',replacement = '',fixed = TRUE)
  Sys.sleep(0.01)
  aln <- readGAlignments(BamFile(file = i,index = paste(i,'.bai',sep = '')), use.names = TRUE)
  matchLenDistr <- histMatchLength(aln, 0)
  pdf(file = paste(name,'_histogram.pdf',sep = ''), width = 10, height = 10, onefile = FALSE)
  print(matchLenDistr[[2]])
  dev.off()
  #transform the GAlignments object into a GRanges object
  #(faster processing of the object)
  alnGRanges <- readsToStartOrEnd(aln, what="start")
  oneBinRanges <- aroundPromoter(txdb, alnGRanges, flankSize = 100)
  #the coverage in the TSS flanking region for the reads with match sizes 29:31
  listPromoterCov <-  readStartCov(alnGRanges, oneBinRanges, matchSize=c(26:32), fixedInterval=c(-100, 100), 
                                       charPerc = "sum", renameChr="aroundTSS")
  pdf(paste(name,"_periodicity.pdf", sep=""), width = 10, height = 10, onefile = FALSE)
  print(plotSummarizedCov(listPromoterCov))
  dev.off()
}


