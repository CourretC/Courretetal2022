#!/bin/sh
#SBATCH --job-name PopStructure               
#SBATCH -p standard
#SBATCH -e PopStructure.err
#SBATCH -o PopStructure.out
#SBATCH -N 1
#SBATCH -c 12                    
#SBATCH --time 12:00:00               
#SBATCH --mem=10gb


module purge
module load plink/1.9
module load vcftools
module load perl
module load bcftools
module load python3

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

### Remove the SR and ST8 sample
vcftools --remove-indv SR_1 --remove-indv ST8_1 --gzvcf Variant/${VCF}.vcf.gz --recode --out Variant/${VCF}_v2

############################################################
######################## PCA #######################
############################################################


mkdir Plink

### Estimate de Linkage desequilibrium along the genome. 

plink --vcf Variant/${VCF}_newID_v2.recode.vcf --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--r2 --out Plink/ld



# perform linkage pruning - i.e. identify prune sites

plink --vcf Variant/${VCF}_newID_v2.recode.vcf --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.4 --out Plink/${VCF}_newID_v2.0.4
 


# prune and create pca
#plink --vcf Variant/${VCF}_newID_v2.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract Plink/${VCF}_newID_v2.0.4.prune.in \
--make-bed --pca --out Plink/${VCF}_newID_v2.0.4

#plink --vcf Variant/${VCF}_newID_v2.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# \
--make-bed --pca --out Plink/${VCF}_newID_v2.1


############################################################
######################## Phylogeny ######################## 
############################################################

## Make an aligment from the vcf file directly. 
python3 Variant/vcf2phylip.py -i Variant/${VCF}_newID_v2.recode.vcf --fasta --output-folder Variant --output-prefix ${VCF}_newID_v2




