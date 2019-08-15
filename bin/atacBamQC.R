#!/usr/bin/env Rscript

library(ATACseqQC)
library(tidyverse)
library(scales)

###This script takes a directory of properly paired bam files
args <- commandArgs(trailingOnly=T)
args
length(args)
dir <- args[1]
outDir <- args[2]



theme_set(theme_bw())
if (!dir.exists(outDir)) {dir.create(outDir)}
args


###list of properly paired files
bams <- list.files(path = dir, pattern = "properly_paired_marked.bam$",full.names = T)
labels <- gsub("_properly_paired_sorted.bam","",basename(bams))
index <- list.files(path = dir, pattern = "properly_paired_marked.bam.bai",full.names = T)

##fragment size distribution
pdf(file = file.path(outDir,"FragmentSizeDistribution.pdf"))
fragSize <- fragSizeDist(bams,labels)
dev.off()

#library complexity
pdf(file = file.path(outDir,"LibraryComplexity.pdf"))
map2(bams,labels, function(x,y) {
	estimateLibComplexity(readsDupFreq(x))
	title(sub=y)
	})

#bamQC
cleaned="./Cleaned_bams"
if(!dir.exists(cleaned)) {dir.create(cleaned)}

cleaned_bams <- file.path(cleaned,gsub("_properly_paired_marked.bam","_cleaned.bam",basename(bams)))
l <- list(bam=bams,index=index,out=cleaned_bams)

QC <- pmap(l,function(bam,index,out){bamQC(bamfile = bam,index = index,outPath=out,doubleCheckDup=TRUE)})
names(QC) <- labels
save(QC, file=file.path(outDir,"bamQC_object.RData"))

QC_df <- data.frame(sample=names(QC),
                    total_reads=map_dbl(QC,"totalQNAMEs"),
                    proper_pairs=map_dbl(QC,"properPairRate"),
                    unmapped=map_dbl(QC,"unmappedRate"),
                    unmapped_mate=map_dbl(QC,"hasUnmappedMateRate"), 
                    duplication=map_dbl(QC,"duplicateRate"),  
                    percent_mitochondria = map_dbl(QC,"mitochondriaRate"),
                    nonRedundantFraction = map_dbl(QC,"nonRedundantFraction"),
                    PCR1=map_dbl(QC,"PCRbottleneckCoefficient_1"),
                    PCR2=map_dbl(QC,"PCRbottleneckCoefficient_2"))


write.csv(QC_df,file=file.path(outDir,"BamQCstats.csv"),row.names=F)

ggplot(QC_df,aes(sample,duplication)) +
  geom_col() +
  ylab("% Duplication") +
  scale_y_continuous(limits=c(0,1),labels = percent_format()) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

ggsave(file=file.path(outDir,"Percent_Duplication.png"))

ggplot(QC_df,aes(sample,percent_mitochondria)) +
  geom_col() +
  ylab("% mitochondrial reads") +
  scale_y_continuous(limits=c(0,1),labels = percent_format()) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))


ggsave(file=file.path(outDir,"Percent_Mitochondria.png"))


ggplot(QC_df,aes(sample,nonRedundantFraction)) +
  geom_col() +
  geom_hline(yintercept = 0.9,linetype = "dashed") +
  scale_y_continuous(limits = c(0,1)) +
  ylab("nonRedundant Fraction") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

ggsave(file=file.path(outDir,"nonRedundantFraction.png"))

ggplot(QC_df,aes(sample,PCR1)) +
  geom_col() +
  geom_hline(yintercept = 0.9,linetype = "dashed") +
  ylab("PCR bottleneck coefficient 1") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

ggsave(file=file.path(outDir,"PCR1.png"))

ggplot(QC_df,aes(sample,PCR2)) +
  geom_col() +
  geom_hline(yintercept = 3,linetype = "dashed") +
  ylab("PCR bottleneck coefficient 2") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

ggsave(file=file.path(outDir,"PCR2.png"))


