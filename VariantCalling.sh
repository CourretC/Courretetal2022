#!/bin/sh
#SBATCH --job-name VariantCalling               
#SBATCH -p charlesworths
#SBATCH -e VariantCalling_all.err
#SBATCH -o VariantCalling_all.out
#SBATCH -N 1
#SBATCH -c 24                    
#SBATCH --time 5-00:00:00               
#SBATCH --mem=100gb


module load bcftools
module load vcftools
module load java
module load picard

#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################

#### Specify the different Imput and other parameters 
DIR="/path/to/your/file"

REF="/path/to/your/reference/dsim_scaffold2_2019.fasta"

cd $DIR

#####################################################################################################################################################################

############################################################
######################## SNP calling #######################
############################################################

mkdir Variant

## to keep only the variant
bcftools mpileup --threads 24 -a AD,DP,SP -f ${REF} -Ou -b Variant/bam_list.txt | bcftools call --threads 24 -f GQ,GP -Oz -mv -o Variant/var.vcf.gz
## keep all position including non variant site. Can be used to calculate the pi
bcftools mpileup --threads 24 -a AD,DP,SP -f ${REF} -Ou -b Variant/bam_list.txt | bcftools call --threads 24 -f GQ,GP -Oz -m -o Variant/all.vcf.gz


## Keep only the Y chromosome
vcftools --gzvcf Variant/all.vcf.gz --chr Y_Contig_1 --chr Y_Contig_100 --chr Y_Contig_11 --chr Y_Contig_113 --chr Y_Contig_121 --chr Y_Contig_126  --chr Y_Contig_13 \
--chr Y_Contig_135 --chr Y_Contig_138 --chr Y_Contig_141  --chr Y_Contig_16  --chr Y_Contig_28  --chr Y_Contig_33 \
--chr Y_Contig_36  --chr Y_Contig_37  --chr Y_Contig_47  --chr Y_Contig_69  --chr Y_Contig_92  --chr Y_Contig_94 \
--chr Y_scaffold1  --chr Y_scaffold2  --chr Y_scaffold3  --chr Y_scaffold4  --chr Y_scaffold5 --recode --out Variant/all.Y

gzip Variant/all.Y.recode.vcf

############################################################
######################## QC #######################
############################################################

## Before to apply the filtering we need to determine which treshold to apply. 
## To do so we are calculating a couple of statistique to have an idea of the quality of our data. 


## Number of unfiltered variants on the Y chromosomes. 
N=$(bcftools view -H Variant/all.Y.recode.vcf.gz | wc -l)
echo "$N unfiltered Variants on the Y chromosome"

VCF="Variant/var.Y.recode.vcf"
OUT="Variant/var.Y"

#####Calculate mean depth per individual
vcftools --vcf $VCF --depth --out $OUT


#####Calculate mean depth per site
vcftools --vcf $VCF --site-mean-depth --out $OUT


#####Calculate site quality
vcftools --vcf $VCF --site-quality --out $OUT


#####Calculate proportion of missing data per individual

vcftools --vcf $VCF --missing-indv --out $OUT


#####Calculate proportion of missing data per site

vcftools --vcf $VCF --missing-site --out $OUT


## The different outputs are then proceed in R. 


############################################################
######################## Filtering #######################
############################################################

VCF="Variant/var.Y.recode.vcf.gz"
OUT="Variant/var.Y.flt"

# set filters
MISS=0.9
QUAL=30
MIN_DEPTH=10
MAX_DEPTH=50

# perform the filtering with vcftools
vcftools --gzvcf $VCF \
--remove-indels --max-missing $MISS --minQ $QUAL \
--min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH \
--minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --out $OUT


## keep only homozygote SNP
VCF="Variant/all.Y.flt.recode.vcf"
OUT="Variant/all.Y.flt.hmz.vcf"

bcftools filter $VCF -e 'GT="het"' -Ov -o $OUT


## Number of unfiltered variants on the Y chromosomes. 
N=$(bcftools view -H $OUT | wc -l)
echo "$N filtered Variants on the Y chromosome"




################# Variant QC

# Creat a dict for the Ref
java -jar /software/picard/2.12.0/picard.jar CreateSequenceDictionary \
        R=$REF \
        O=/gpfs/fs2/scratch/alarracu_lab/Cecile/reference/Dsimulans/dsim_scaffold2_2019.dict
        
  
# Variant QC        
java -jar DISCVRSeq-1.3.10.jar VariantQC \
	--maxContigs 100 \
     -R $REF \
     -V Variant/var.Y.flt.hmz.vcf \
     -O Variant/VariantQC.html
     
      
