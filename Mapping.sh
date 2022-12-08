#!/bin/sh
#SBATCH --job-name YMapping               
#SBATCH -e YMapping.err
#SBATCH -o YMapping.out
#SBATCH -p charlesworths
#SBATCH -c 20
#SBATCH -n 1                       
#SBATCH --time 72:00:00                
#SBATCH --mem=100gb

module purge
module load fastqc
module load trimgalore 
module load cutadapt
module load bwa
module load samtools
module load bedtools
module load java
module load picard/2.12.0

#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################

#### Specify the different Imput and other parameters 
DIR="/path/to/your/file"


REF="/path/to/your/reference/dsim_scaffold2_2019.fasta"

# Modify the fastq file name to correspond to this format : ${SAMPLE}_1.fastq.gz and ${SAMPLE}_2.fastq.gz
#Put them in a fastq/ folder. 

#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
##################################################################################################################################################################### 

cd ${DIR}

# create a directory for output file
for SAMPLE in `cat sample.txt`
do
#################################
########## Trimming ##########
#################################
mkdir fastq_Trim

if [ -f "fastq_Trim/${SAMPLE}_R1_val_1.fq.gz" ]; then
	    echo "reads already trimmed"
else	
		echo "trimming reads"

trim_galore --paired --nextera --length 75 --phred33 --no_report_file --fastqc -o fastq_Trim fastq/${SAMPLE}_R1.fastq.gz fastq/${SAMPLE}_R2.fastq.gz

fi

N=$(( $(zcat "fastq_Trim/${SAMPLE}_R1_val_1.fq.gz" | wc -l) /4 ))
M=$(( $(zcat "fastq_Trim/${SAMPLE}_R2_val_2.fq.gz" | wc -l) /4 ))
	
	echo "${SAMPLE} number of raw reads: R1: $N; R2: $M"


#################################
########## Mapping ###############
###################################
mkdir Mapping

if [ -f "Mapping/${SAMPLE}_bwa_sorted.bam.bai" ]; then
	echo "reads already mapped"
else
	echo "Mapping"
	
	
########### Index ########### 

bwa index ${REF} 

########### Mapping ########### 

mkdir Mapping

	bwa mem -t 20 ${REF} \
	fastq_Trim/${SAMPLE}_R1_val_1.fq.gz  \
	fastq_Trim/${SAMPLE}_R2_val_2.fq.gz \
	> Mapping/${SAMPLE}_bwa.sam

	samtools view -bSh Mapping/${SAMPLE}_bwa.sam  >  Mapping/${SAMPLE}_bwa.bam
	rm Mapping/${SAMPLE}_bwa.sam 

	samtools sort -@ 20 Mapping/${SAMPLE}_bwa.bam > Mapping/${SAMPLE}_bwa_sorted.bam
	rm Mapping/${SAMPLE}_bwa.bam

	samtools index Mapping/${SAMPLE}_bwa_sorted.bam
fi

N=$(( $(samtools view "Mapping/${SAMPLE}_bwa_sorted.bam" | wc -l) /2 ))
echo "${SAMPLE} got mapped using bwa total pairs mapped: $N"


###########################################################################
########################### Mark Duplicate #################################
###########################################################################

if [ -f "Mapping/${SAMPLE}_bwa_noDup.bam" ]; then
	echo "Duplicate alread removed"
else
	echo "Remove duplicate"
	
java -Xmx10g -jar /software/picard/2.12.0/picard.jar MarkDuplicates \
          REMOVE_DUPLICATES=TRUE \
          I=Mapping/${SAMPLE}_bwa_sorted.bam \
          O=Mapping/${SAMPLE}_bwa_noDup.bam \
          M=Mapping/${SAMPLE}_bwa_noDup.txt

fi 


###########################################################################
######################### Quality filtering ################################
###########################################################################


if [ -f "Mapping/${SAMPLE}_bwa_q30.bam" ]; then
	echo "reads filter"
else
	echo "Quality Filtering"


### Quality filtering - keeping multimapping reads	
samtools view -b -h -f 3 -F 4 -F 8 -F 256 -F 2048  Mapping/${SAMPLE}_bwa_noDup.bam > Mapping/${SAMPLE}_bwa_q0.bam
samtools index Mapping/${SAMPLE}_bwa_q0.bam


### Mapping quality
samtools view Mapping/${SAMPLE}_bwa_q0.bam | awk -F "\t" '{print $5}' > Mapping/${SAMPLE}_bwa_q0.txt

### Quality filtering - keeping uniquely mapping reads
samtools view -b -h -f 3 -F 4 -F 8 -F 256 -F 2048 -q 30 Mapping/${SAMPLE}_bwa_noDup.bam > Mapping/${SAMPLE}_bwa_q30.bam
samtools index Mapping/${SAMPLE}_bwa_q30.bam


fi

N=$(( $(samtools view "Mapping/${SAMPLE}_bwa_q30.bam" | wc -l) /2 ))
echo "${SAMPLE} Number of read after filtering: $N"



###########################################################################
######################### Coverage ################################
###########################################################################

mkdir Coverage

if [ -f "Coverage/${SAMPLE}_bwa_q30_coverage.txt" ]; then
	echo "Coverage : done "
else
	echo "Coverage"

bedtools genomecov -ibam Mapping/${SAMPLE}_bwa_q30.bam > Coverage/${SAMPLE}_bwa_q30_coverage.txt
fi 
done 

