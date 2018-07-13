##!/bin/bash

###fastqc
###trim
##fastqc
##align
##markdups
##get insert size
##shift reads,
## call peaks with lenient p-val
## merge peaks
## make bed file with the coordinates of merged peaks
## R


###vars
hisat2Index=/local_storage/annotation_db//UCSC/mm10/Sequence/newHisat
TSSref=hg19_refseq_.bed
chromInfo=hg19_chromInfo.txt
sampleSheet=18-015_samples.txt
annotation=mm10
### in alignment section change  the sed command to clean up the file names

###### Step 1. FastQC

if [ ! -d "./FastQC" ]
then
	mkdir ./FastQC
fi 


fastqc -t 7 -o ./FastQC Raw_fastq/*gz


##### Step 2. Trim Adaptors

trimDir=./Trimmed_fastq

if [ ! -d $trimDir ]
then

	mkdir $trimDir

fi

#####
for i in `ls Raw_fastq/*.gz | cut -f2 -d "/" | cut -f1-7 -d'_' | uniq`

do

	echo "trimming adapter of sample $i" 

	trimmomatic-0.33.jar PE -threads 14 -phred33 -trimlog $outDir/${i}.log.txt \
	Raw_fastq/${i}_1.fq.gz Raw_fastq/${i}_2.fq.gz \
	$trimDir/${i}_1_paired_at.fq.gz $trimDir/${i}_1_Unpaired_at.q.gz \
	$trimDir/${i}_2_paired_at.fq.gz $trimDir/${i}_2_Unpaired_at.fq.gz \
	ILLUMINACLIP:/usr/share/Trimmomatic-0.33/adapters/NexteraPE-PE.fa:2:30:10 MINLEN:30

	echo $i trim done!

done

##now run fastqc on trimmed files


fastqc_dir="./Trimmed_fastqc"

if [ ! -d $fastqc_dir ]
then

	mkdir $fastqc_dir
fi

fastqc -o $fastqc_dir ${trimDir}/*fq.gz





##############Step 3. Align to genome 


alignDir="./Alignments"

if [ ! -d "$alignDir" ]
then

	mkdir $alignDir

fi


if [ ! -e "genome.1.ht2" ]
then

ln -s ${hisat2index}/genome* .
fi

for i in $(ls ./Trimmed_fastq/*.gz | cut -f3 -d '/' | sed 's/_[12]_*paired_at.fq.gz//g' | uniq) 
do 
	##cleanup file names
	sample=$(echo $i | sed 's/ATAC_|_USPD[0-9]+_[A-Z]+_L[0-9]//g' | sed 's/_/-/g')

	hisat2 -x genome -1 $trimDir/${i}_1_paired_at.fq.gz -2 $trimDir/${i}_2_paired_at.fq.gz \
	-p 14 --no-spliced-alignment -I 10 -X 2000 | samtools view -f2 -b - > $alignDir/${sample}.bam 2>${alignDir}/${sample}_hisat2_stderr.txt

	
	samtools sort -@ 6 -m 5G -o ${alignDir}/${sample}_sorted.bam ${alignDir}/${sample}.bam

	samtools index ${alignDir}/${sample}_sorted.bam
	
	rm ${alignDir}/${sample}.bam

	echo "$i alignment done!"

done

############## Step 4. Get Fragment Size

insertDir="./Insert_size"

if [ ! -d $insertDir ]
then

	mkdir $insertDir

fi


for f in ${alignDir}/*sorted.bam
do 

	sample=$(echo $f | cut -f3 -d "/" | sed 's/_sorted.bam//g')

	picard.jar CollectInsertSizeMetrics I=$f O=$insertDir/${sample}_insertMetrics.txt H=$insertDir/${sample}_insertHist.pdf 

done

################# Step 5. Get  Duplication metrics

dupDir="./Marked_bams"
if [ ! -d "$dupDir" ]
then

	mkdir $dupDir
fi

for f in $alignDir/*sorted.bam
do

	sample=$(echo $f | cut -f3 -d "/" | sed 's/_sorted.bam//g')

	picard.jar MarkDuplicates I=$f O=$dupDir/${sample}_marked.bam M=$dupDir/${sample}_dupInfo.txt 

done

################ Step 6. Transposition correction 


shifted="./Corrected_bams"
if [ ! -d "$shifted" ]
then

	mkdir $shifted

fi


for file in ${dupDir}/*marked.bam
do

        sample=$(echo $file | cut -f3 -d "/"| sed 's/_marked.bam//g' )
        echo correcting $sample

        date +"%m/%d/%Y %H:%M:%S"
        atacShift.py -b $file -o ${shifted}/${sample}_corrected.bam
         date +"%m/%d/%Y %H:%M:%S"      

done

################### Step 7.  Calling peaks for each sample
macs="./macs"

if [ ! -d $macs ]
then

        mkdir $macs

fi

for file in $shifted/*bam
do 

        sample=$(echo $file | cut -f3 -d "/" | sed 's/_corrected.bam//g') 
        
        echo "calling peaks for sample $sample" 
        date +"%m/%d/%Y %H:%M:%S"
        macs2 callpeak -t $file -f BAMPE --outdir $macs/${sample}_macs -n $sample -B --SPMR  -g mm -q 0.1
        date +"%m/%d/%Y %H:%M:%S"
done



################### Step 8. convert bedgraph to bigWigs

#convert bedgraphs to big wig

for f in $macs/*_macs/*pileup.bdg

do

echo $f

sample=`echo $f | cut -f3 -d "/" | sed 's/_macs//g')`

sort -k1,1 -k2,2n $f > tmp.bdg

bedGraphToBigWig tmp.bdg $chromInfo $macs/${sample}.bw

rm tmp.bdg

done



##################### Step 9. Merge peaks and counting

merged="Peak_Counts"

if [ ! -d $merged ]
then
	mkdir $merged

fi


files=$(ls $macs/*/*narrowPeak)


mergePeaks -matrix $merged/merged -venn $merged/merged_peaks_venn.txt $files  > $merged/merged_peaks.bed


echo  -e "Geneid\tChr\tStart\tEnd\tStrand" > tmp
grep -v "^#" $merged/merged_peaks.bed | awk 'BEGIN{OFS = "\t"} {print $1,$2,$3+1,$4,$5}' >> tmp

mv tmp $merged/merged_peaks.saf

featureCounts -f -F SAF -o $merged/peak_counts.txt -a $merged/merged_peaks.saf -p -d 10 -D 2000 ${shifted}/*bam 2> $merged/featureCounts.log


##################### Step 10. Encode QC

#### TSS enrichment

QC="Encode_QC"

if [ ! -d $QC ]
then

	mkdir $QC

fi

files=$(ls ${macs}/*bw)
samples=$(cut -f1 $sampleSheet | tr "\r" " ")
computeMatrix reference-point --referencePoint TSS \
                -b 1000 -a 1000 -R $TSSref \
                -S $files -out $QC/TSS_matrix.txt.gz \
                --samplesLabel $samples


TSS.R matrix=${QC}/TSS_matrix.txt.gz out=$QC d=1000 binSize=10 samples=$sampleSheet annotation=$annotion


#### FRiP


##a bin size of 10 bp is the default for computeMatrix

