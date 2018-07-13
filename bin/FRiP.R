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
#	samples- tab-delimited file of samples with three columns sampleName, description and coreNumber.
#		sample name must be in the format of "group_rep1"  ex. Treat_rep1.
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
summary_file <- str_split(args[grep("summary",args)],"=",simplify = T)[2]
sample_file <- str_split(args[grep("samples",args)],"=",simplify = T)[2]
out_dir <- str_split(args[grep("out",args)],"=",simplify = T)[2]
sample_file

###calculate FrIP (Fragments in peaks)
count_summary <- read.delim(summary_file,header = T,stringsAsFactors = F)

colnames(count_summary) <- gsub(".*Corrected_bams\\.(\\d+)[\\.-_](\\d+)[\\.-_](\\d+).*",paste("\\1","\\2","\\3",sep = "-"),colnames(count_summary))
rownames(count_summary) <- count_summary$Status
count_summary <- count_summary[,-1]

count_summary
FRiP<-apply(count_summary,2,function(x) x["Assigned"]/sum(x))
FRiP <- round(FRiP,2)

FRiP_df <- data.frame(coreNumber = names(FRiP),FRiP = FRiP)

FRiP_df
##make a better chart, with sample Names and exp.
sample_info <- read.delim(sample_file,stringsAsFactors = F,header = F)
colnames(sample_info) <- c("sample","description","coreNumber")

sample_info
##add batch and group

sample_info$group<-gsub("_rep.*","",sample_info$sample)
#match by corenumber

FRiP_df<-inner_join(sample_info,FRiP_df)
FRiP_df
ggplot(FRiP_df,aes(sample,FRiP)) +
  geom_col(aes(fill=group)) +
  geom_hline(yintercept = 0.2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))


ggsave(file = file.path(out_dir,"FRiP.pdf"))
write.csv(FRiP_df,file=file.path(out_dir,"FRiP_scores.csv"),row.names = F)

