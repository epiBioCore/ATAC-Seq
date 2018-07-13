#!/usr/bin/env Rscript

#################################################################################################################
#
#   Calculates TSS enrichment scores and plots the profile around TSS
#
#   usage: 
#	TSS.R matrix=matrix.txt.gz samples=samples.txt out=encode_QC d=1000 binSize=10
#
#   inputs:
#	matrix -gziped matrix file (matrix.txt.gz) from deepTools computeMatrix function
#	samples- samples.txt file with three columns: sampleName, description, coreNumber.
#		sample names must be in the form: group_rep1, ex. Treat_rep1
#	out-results directory
#	d- the before/after value used in deepTools computetMatrix TSS
#	binSize-the binSize used in deepTools computeMatrix TSS (the default is 10)
#
#
#  outputs:
#	number_of_excluded_genes.csv - the number of genes that had 0 signal 1000 bp from TSS. 
#					Thes are exclude from the enrichment profile	
#					
#	TSS_scores.pdf - barplot of TSS enrichment score accross all sample
#	TSS_scores.csv - the TSS enrichment scores
#	enrichment_profile_by_group.pdf - enrichment profile of all samples, colored by group,or treatment
#	enrichment_profile_by_sample.pdf - enrichment profiles of each sample.
#
#
#
#
#  author: Danielle Perley
#
###################################################################################################################

library(tidyverse)

args<-commandArgs()
args

##parsing command line arguments
matrix_file <- str_split(args[grep("matrix",args)],"=",simplify = T)[2]
sample_file <- str_split(args[grep("samples",args)],"=",simplify = T)[2]
out_dir <- str_split(args[grep("out",args)],"=",simplify = T)[2]
d <- str_split(args[grep("^d",args)],"=",simplify = T)[2]
bin_size <- str_split(args[grep("binSize",args)],"=",simplify = T)[2]
annotation <- str_split(args[grep("annotation",args)],"=",simplify = T)[2]

d <- as.numeric(d)
bin_size <- as.numeric(bin_size)


mat <- read.delim(gzfile(matrix_file),skip = 2,header = F)
values<-mat[,-c(1:6)]
rownames(values) <- mat$V4
header<-read.table(gzfile(matrix_file),sep=",",nrows=1,stringsAsFactors = F)


##extract sample names from header and repeat, so that it matches the matrix
last <- ncol(header)-19
sample_ind<-8:last
sample_names<-header[,sample_ind]
sample_names$V8<-gsub("sample_labels:\\[","",sample_names$V8)
sample_names[,ncol(sample_names)] <- gsub("\\]","",sample_names[,ncol(sample_names)])

sample_names<-unlist(sample_names)

##calculate the number of times to repeat (number of bins)
##for each sample, computeMatrix counts the number of reads in each bin of x size, plus or minus y distance from TSS
total_bp <- d*2
nbin <-total_bp/bin_size
sample_name_header<-rep(sample_names,each=nbin)
colnames(values)<-sample_name_header

##sep by sample
values_by_sample <- lapply(sample_names,function(x) values[,grep(x,colnames(values))])

names(values_by_sample) <- sample_names

sapply(values_by_sample,dim)

######normalization factor, 100bp upstream and downstream  the TSS
##get the bins that correspond to the first  100 bp and last 100 bp


getNF <- function(x) {
  	
  upstream <-  1:round(100/bin_size)
  downstream <- seq(ncol(x) - (round(100/bin_size)-1),ncol(x)) 		
  ind<-c(upstream,downstream)
  rowMeans(x[,ind])
}


nfs<-lapply(values_by_sample,getNF)



###normalize values
normed<-Map(function(x,nf) {x/nf},
            x = values_by_sample,nf = nfs)

#filter out rows with NA values, that have been divided by zero
filtered_norm<-lapply(normed,function(x) {
  x[!is.na(rowSums(x)),]
})


###the number of genes with zero
zero_enrich <- lapply(normed,function(x) {
  sum(is.na(rowSums(x)))
})


zero_enrich_df <- t(as.data.frame(zero_enrich))
colnames(zero_enrich_df) <- "num_Genes"
write.csv(zero_enrich_df,file = file.path(out_dir,"number_excluded_transcripts.csv"))

##calculate average enrichment across all genes
enrichment<-lapply(filtered_norm,colMeans,na.rm = T)

enrichment_df<-as.data.frame(enrichment)

##get TSS score
TSS <- nbin/2
TSS_score<-enrichment_df[TSS,]
TSS_score <- round(TSS_score,1)


TSS_score <- as.data.frame(t(TSS_score))
colnames(TSS_score)<-"TSS_score"

TSS_score$sample <- rownames(TSS_score)


sample_info <- read.delim(sample_file, header = F, stringsAsFactors = F)
colnames(sample_info) <- c("sample","description","coreNumber")

sample_info$Group <- gsub("_rep.*","",sample_info$sample)

TSS_scores <- merge(sample_info,TSS_score)



p <- ggplot(TSS_scores,aes(sample,TSS_score)) +
     geom_col(aes(fill = Group)) +
     labs(x = "Sample Name", y = "TSS enrichment score") +
     theme_bw() +
     theme(axis.text.x = element_text(angle = 90))

if (annotation == "hg19") {
	p +
	  geom_hline(yintercept = 6) +
  	  geom_hline(yintercept = 10) 
} else if (annotation == "mm10") {
	p +
          geom_hline(yintercept = 10) +
          geom_hline(yintercept = 15) 
} else  {
	p
	}

ggsave(file = file.path(out_dir,"TSS_scores.pdf"))

write.csv(TSS_scores,file = file.path(out_dir,"TSS_scores.csv"),row.names = F)

######Plot enrichment

##first format
##making the position be the midpoint of the bin
position <- seq(-995,length.out = 200,by = 10)
enrichment_df$position <- position


enrichment_df_l <- gather(enrichment_df,sample,enrichment,-position)

##add group and batch
enrichment_df_l$Group <- sapply(enrichment_df_l$sample, function(x) sample_info[grep(x,sample_info$sample),"Group"])

ggplot(enrichment_df_l,aes(position,enrichment)) +
  geom_line(aes(group = sample,color = Group)) +
  theme_bw()

ggsave(file = file.path(out_dir,"enrichment_profile_by_group.pdf"))


###facet by sample
ggplot(enrichment_df_l,aes(position,enrichment)) +
  geom_line(aes(group = sample,color = Group)) +
  facet_wrap(~sample) +
  theme_bw()

ggsave(file = file.path(out_dir,"enrichment_profile_by_sample.pdf"))


