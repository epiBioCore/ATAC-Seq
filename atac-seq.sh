#!/bin/bash

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
hisat2index=/local_storage/annotation_db/Mus_musculus/UCSC/mm10/Sequence/newHisat
TSSref=mm10_refseq_all_tx.bed 
chromInfo=mm10_chrominfo.txt
sampleSheet=19-064_samples.txt 
annotation=mm10


###### Step 1. FastQC

if [ ! -d "./FastQC" ]
then
	mkdir ./FastQC
fi 


fastqc -t 10 -o ./FastQC Raw_fastq/*gz


##### Step 2. Trim Adaptors

trimDir=./Trimmed_fastq

if [ ! -d $trimDir ]
then

	mkdir $trimDir

fi

#####
while read -r coreNumber name r1 r2
do

	echo "trimming adapter of sample $coreNumber" 

	java -jar /usr/share/Trimmomatic-0.33/trimmomatic-0.33.jar PE -threads 14 -phred33 -trimlog $trimDir/${coreNumber}.log.txt \
	Raw_fastq/${r1} Raw_fastq/${r2} \
	$trimDir/${coreNumber}_1_paired_at.fq.gz $trimDir/${coreNumber}_1_Unpaired_at.q.gz \
	$trimDir/${coreNumber}_2_paired_at.fq.gz $trimDir/${coreNumber}_2_Unpaired_at.fq.gz \
ILLUMINACLIP:/usr/share/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10:3:TRUE MINLEN:10 2> ${trimDir}/${coreNumber}_trimStats.txt

	echo $coreNumber trim done!

done < ${sampleSheet}


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

for i in $(ls ${trimDir}/*gz | cut -f3 -d '/' | cut -f1 -d "_" | sort | uniq) 
do 

	echo aligning $i
	hisat2 -x genome -1 $trimDir/${i}_1_paired_at.fq.gz -2 $trimDir/${i}_2_paired_at.fq.gz \
	-p 14 --no-spliced-alignment -I 10 -X 2000 | samtools view -b - > $alignDir/${i}.bam 2>${alignDir}/${i}_hisat2_stderr.txt

	samtools view -b -f2 $alignDir/${i}.bam > $alignDir/${i}_properly_paired.bam
	samtools sort -@ 6 -m 5G -o ${alignDir}/${i}_properly_paired_sorted.bam ${alignDir}/${i}.bam

	samtools index ${alignDir}/${i}_properly_paired_sorted.bam
	
	rm ${alignDir}/${i}_properly_paired.bam

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

	java -jar /usr/share/picard-tools-2.2.4/picard.jar CollectInsertSizeMetrics I=$f O=$insertDir/${sample}_insertMetrics.txt H=$insertDir/${sample}_insertHist.pdf 

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

	java -jar /usr/share/picard-tools-2.2.4/picard.jar MarkDuplicates I=$f O=$dupDir/${sample}_marked.bam M=$dupDir/${sample}_dupInfo.txt 

done

################# Step 6. BamQC


QC="./QC"
if [ ! -d "$QC" ]
then

	mkdir $QC

fi


#for file in ${dupDir}/*marked.bam
#do
#
#        sample=$(echo $file | cut -f3 -d "/"| sed 's/_marked.bam//g' )
#        echo correcting $sample
#
#        date +"%m/%d/%Y %H:%M:%S"
#        atacShift.py -b $file -o ${shifted}/${sample}_corrected.bam
#         date +"%m/%d/%Y %H:%M:%S"      
#
#done


#atacBamQC.R $dupDir $QC





################### Step 7.  Calling peaks for each sample
macs="./macs"

if [ ! -d $macs ]
then

        mkdir $macs

fi

if [ $annotation == "hg19" ]
then
	genome="hs"
elif [ $annotation == "mm10" ]
then
	genome="mm"
fi



for file in $dupDir/*bam
do 

        sample=$(basename $file | cut -f1 -d "_") 
        
        echo "calling peaks for sample $sample" 
        date +"%m/%d/%Y %H:%M:%S"
        macs2 callpeak -t $file -f BAMPE --outdir $macs/${sample}_macs -n $sample -B --SPMR  -g $genome -q 0.05
        date +"%m/%d/%Y %H:%M:%S"
done



################## Step 8. convert bedgraph to bigWigs

#convert bedgraphs to big wig

for file in $macs/*_macs/*pileup.bdg

do

$file

sample=$(basename $file | cut -f1 -d "_")

sort -k1,1 -k2,2n $file > tmp.bdg

bedGraphToBigWig tmp.bdg $chromInfo $macs/${sample}.bw

rm tmp.bdg

done



############## Step 9. Encode QC


files=$(ls ${macs}/*bw)
samples=$(for i in $files; do basename $i | cut -f1 -d "."; done)

computeMatrix reference-point --referencePoint TSS \
                -b 1000 -a 1000 -R $TSSref \
                -S $files -out $QC/TSS_matrix.txt.gz \
                --samplesLabel $samples

TSS.R matrix=${QC}/TSS_matrix.txt.gz out=$QC d=1000 binSize=10 samples=$sampleSheet annotation=$annotation


### FRiP

counts=./counts
if [ ! -d $counts ]
then
        mkdir $counts
fi

for file in ${macs}/*/*narrowPeak
do

        sample=$(basename $file | cut -f1 -d "_")
       saf=${counts}/${sample}_peaks.saf
       echo  -e "Geneid\tChr\tStart\tEnd\tStrand" > tmp
       awk 'BEGIN{OFS="\t"}{print $4,$1,$2 +1,$3,$6}' $file >> tmp
       mv tmp $saf

       featureCounts -f -F SAF -o ${counts}/${sample}_counts.txt -a $saf -p -d 10 -D 2000 Marked_bams/${sample}*bam 2> ${counts}/${sample}_featureCounts.log

done

FRiP.R summary=$counts out=$QC


#other QC
./atacBamQC.R $dupDir $QC
