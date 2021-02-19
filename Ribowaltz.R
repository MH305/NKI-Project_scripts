library(riboWaltz)
library(ggplot2)
library(gridExtra)

gencodeM21 = create_annotation(gtfpath = "/Users/mashrurhaidernew/Desktop/NKI/RiboSAK/26-32-TX.profile.counts")

gencodeM21 = create_annotation(gtfpath = "D:/NKI/gencode.vM21.annotation.gtf")


reads_list <- bamtolist(bamfolder = "D:/NKI/BAM/", annotation = gencodeM21,
                        transcript_align = TRUE, refseq_sep = '|')# needs change 

psite_data <- psite(data = reads_list,plot = TRUE,extremity = "5end")
reads_psite_list <- psite_info(reads_list, psite_data,fastapath = "D:/NKI/gencode.vM1.pc_transcripts.fa.gz",
                               fasta_genome = FALSE) # needs change 

samples = c("R10_SC_dep_3", "R10_SC_dep_4", "R2_SC_dep_1", "R2_SC_dep_2", "R7_SC_dep", "R8_SC_dep", "R10_SC_enr_3", "R10_SC_enr_4", "R7_SC_enr", "R8_SC_enr")

## Length dist ##
plots<-list()
for (sampleid in samples){
  example_length_dist_zoom <- rlength_distr(reads_list, sample = sampleid)
  plots[[sampleid]] = example_length_dist_zoom[["plot"]]
}
g<-grid.arrange(grobs=plots,nrow=3)
ggsave("all_transcripts_read_length_dist.png",g,width = 18,height = 8)

# 5' 3' end heatmaps
plots<-list()
for (sampleid in samples){
  example_ends_heatmap <- rends_heat(reads_list, gencodeM21, sample = sampleid, utr5l = 25, cdsl = 50, utr3l = 25)
  plots[[sampleid]] = example_ends_heatmap[["plot"]]
}
g<-grid.arrange(grobs=plots,nrow=4)
ggsave("5p-3p-end_heatmaps.pdf",g,width = 30,height = 30)

# Psite overlap region barplots
example_psite_region <- region_psite(reads_psite_list, gencodeM21)
ggsave("psite_region_bars.png",example_psite_region[["plot"]]+ theme(axis.text.x = element_text(angle = 90)), width = 8, height = 6)

# Region frame sepc counts
example_frames_stratified <- frame_psite_length(reads_psite_list, region = "all")
ggsave("frame_region_heatmaps.png",example_frames_stratified[["plot"]],width = 8, height = 20)
rm(example_psite_region,example_frames_stratified)
example_frames <- frame_psite(reads_psite_list, region = "all")
ggsave("frame_region_boxplots.png",example_frames[["plot"]],width = 8, height = 20)
rm(example_frames)

# Start Stop plots
plots<-list()
for (sampleid in samples){
  example_metaprofile <- metaprofile_psite(reads_psite_list, gencodeM21, sample = sampleid, utr5l = 20, cdsl = 50, utr3l = 20,  plot_title = "auto")
  plots[[sampleid]] = example_metaprofile[["plot"]]
}
g<-grid.arrange(grobs=plots,nrow=4)
ggsave("start-stop-periodicity.pdf",g,width = 30,height = 16)
rm(example_metaprofile)

