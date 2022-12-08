#!/bin/sh
#SBATCH --job-name Gene_count               
#SBATCH -e Gene_count.err
#SBATCH -o Gene_count.out
#SBATCH -p preempt
#SBATCH -c 10
#SBATCH -n 1                      
#SBATCH --time 10:00:00                
#SBATCH --mem=50gb


module load bedtools
module load samtools
module load blast

#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################



##########################################
################## INPUT #################
##########################################

#REFERENCEGFF="/gpfs/fs2/scratch/alarracu_lab/Genome/Dsimulans/dsim_gene.gtf"

REF="/path/to/your/reference"

DIR="/path/to/your/dir"


### Direction
cd ${DIR}

mkdir count_Genes

##########################################
################## Make a bed file  #################
##########################################


EXON=" "

### Make a bed file for each EXON 
grep "${EXON}" count_Genes/Exon_Coordinate.bed > count_Genes/${EXON}_Coordinate.bed



##########################################
############# Exon Coverage  #################
##########################################

### Total exon coverage for normalization 
samtools depth -b BED/Exon_Coordinate.bed -o tot_cov_Final.out -f sample.txt

cd BED


# COVERAGE FOR EACH EXON
for i in $EXON
do 
i=${i%_Coordinate.bed}

samtools depth -b count_Genes/BED/"$i"_Coordinate.bed -o /scratch/alarracu_lab/Cecile/Y_Project/count_Genes/Coverage/"$i"_cov.out -f /scratch/alarracu_lab/Cecile/Y_Project/count_Genes/sample.txt

done

