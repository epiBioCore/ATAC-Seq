#!/usr/bin/env Rscript

######################################################################################################################
#
#  Calculates the Fraction of fragments in peaks (FRiP)
#
#  usage:
#	Frip.R summary=counts.txt.summary samples=samples.txt out=encode_QC
#
# inputs:
#	summary- the summary file generated from featureCounts of the number of reads aligning to a peak
#	out- results directory
#
# output:
#	FRiP.pdf - barchart of FRiP scores
#	FRiP_scores.csv FRiP scores for all samples		
#
#
# Danielle Perley	
#
######################################################################################################################


library("tidyverse")
args <- commandArgs()
args
summary_dir <- str_split(args[grep("summary",args)],"=",simplify = T)[2]
out_dir <- str_split(args[grep("out",args)],"=",simplify = T)[2]


###calculate FrIP (Fragments in peaks)
summary_files <- list.files(path=summary_dir,pattern="summary",full.names = T)
count_summary <- map(summary_files, ~read.delim(.x,header = T,stringsAsFactors = F,check.names = F))

count_summary <- map(count_summary,function(x) {
colnames(x) <- basename(colnames(x)) %>% str_split(.,"_",simplify = T) %>% .[,1] 
rownames(x) <- x$Status
x <- x[,-1,drop = F]
return(x)})

count_summary <- do.call("cbind",count_summary)
count_summary
FRiP<-apply(count_summary,2,function(y) y["Assigned"]/sum(y))
FRiP <- round(FRiP,2)
FRiP
FRiP_df <- data.frame(coreNumber = names(FRiP),FRiP = FRiP)


#match by corenumber

ggplot(FRiP_df,aes(coreNumber,FRiP)) +
  geom_col() +
  geom_hline(yintercept = 0.2,linetype = "dashed") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust =1))


ggsave(file = file.path(out_dir,"FRiP.pdf"))
write.csv(FRiP_df,file=file.path(out_dir,"FRiP_scores.csv"),row.names = F)

