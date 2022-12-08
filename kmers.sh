#!/bin/sh
#SBATCH --job-name Kmers               
#SBATCH -e Kmers.err
#SBATCH -o Kmers.out
#SBATCH -p standard
#SBATCH -c 20
#SBATCH -n 1                      
#SBATCH --time 5-00:00:00                
#SBATCH --mem=100gb


module load samtools

#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################

DIR="/path/to/your/file"


### Direction
cd ${DIR}


# create a directory for output file
mkdir Kmers
mkdir Kmers/RPM


#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
##################################################################################################################################################################### 

cd ${DIR}

#sample.txt : text file containing the list of samples

for SAMPLE in `cat sample.txt`
do


if [ -f "Kmers/RPM/${SAMPLE}.rep.rpm.total" ]; then
	    echo "TE already count"
else	


### need to concatenate R1 and R2 files before using k-seek. 	
zcat fastq_Trim/${SAMPLE}_R1_val_1.fq.gz fastq_Trim/${SAMPLE}_R1_val_2.fq.gz | perl /scratch/alarracu_lab/Cecile/k-seek-master/k_seek.pl - Kmers/${SAMPLE}


##################### Normalize the count in RPM

##### Number of read mapping for normalization
N=$(samtools view "Mapping/${SAMPLE}_bwa_q0.bam" | wc -l)
echo "total mapped reads number: $N."

##### Normalization
awk -v var="$N" -F'\t' '{OFS="\t"} {print $1, (($2*1000000/var))}' "Kmers/${SAMPLE}.rep.total" > "Kmers/RPM/${SAMPLE}.rep.rpm.total"


fi
done

#### Concatenate all sample in a one table
perl /scratch/alarracu_lab/Cecile/k-seek-master/k_compiler.pl Kmers/RPM Kmers/RPM/Y_all
