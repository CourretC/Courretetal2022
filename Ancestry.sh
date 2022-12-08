#!/bin/sh
#SBATCH --job-name Ancestry              
#SBATCH -p debug
#SBATCH -e Ancestry.err
#SBATCH -o Ancestry.out
#SBATCH -N 1
#SBATCH -c 12                    
#SBATCH --time 1:00:00               
#SBATCH --mem=10gb


module purge
module load bcftools
module load bedtools
module load htslib
module load samtools
module load blast
module load vcftools

#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################


#### Specify the different Imput and other parameters 
DIR="/path/to/your/file"

REF="/path/to/your/reference/dsim_scaffold2_2019.fasta"

VCF="var.Y.flt.hmz"

#####################################################################################################################################################################
cd $DIR



zcat Variant/${VCF}.vcf.gz | vcf-to-tab > Ancestry/${VCF}.tab


######## ####### ####### ####### ###### 
######### Creation of a Consensus
######## ####### ####### ####### ###### 

mkdir Ancestry 

# Create a bed file with the coordinate of each SNP
grep -v "#" Variant/${VCF}.vcf | awk -F'\t' '{ print $1"\t"$2-200"\t"$2+200 }' > Ancestry/${VCF}.bed

# Extract the Consensus for each SNP
bgzip -c Variant/${VCF}.vcf > Variant/${VCF}.vcf.gz
bcftools index Variant/${VCF}.vcf.gz
bcftools consensus --iupac-codes -f $REF Variant/${VCF}.vcf.gz > Ancestry/${VCF}_consensus.fasta # Make a reference genome where the SNP are coded in iupac-codes
samtools faidx Ancestry/${VCF}_consensus.fasta
bedtools getfasta -fi Ancestry/${VCF}_consensus.fasta -bed Ancestry/${VCF}.bed > Ancestry/${VCF}.SNP.consensus.fasta



### ####### ####### ####### ###### ####### ######  ####
######### Blast against the sister species ####### ###### 
######## ####### ####### ####### ###### ####### ###### 


Species="Dmauritiana"
REF_Sis="/path/to/your/reference/dmau_scaffold2_V2.fasta"

# Extract the Y chromosome from the reference.
grep "Y" -A 1 $REF_Sis > Ancestry/Y_${Species}.fasta

# Create the database for blast
makeblastdb -in Ancestry/Y_${Species}.fasta -input_type fasta -dbtype nucl -out Ancestry/Y_${Species}

# Blast
blastn -task blastn -db Ancestry/Y_${Species} -query Ancestry/${VCF}.SNP.consensus.fasta -evalue 1e-20 -qcov_hsp_perc 80 -perc_identity 85 -max_hsps 1  -max_target_seqs 1 -outfmt "6 qseqid sseqid evalue length pident qstart qseq sseq" -out Ancestry/SNP.${Species}.out


