#!/bin/sh
#SBATCH --job-name TE_count               
#SBATCH -e TE_count.err
#SBATCH -o TE_count.out
#SBATCH -p standard
#SBATCH -c 20
#SBATCH -n 1                      
#SBATCH --time 48:00:00                
#SBATCH --mem=100gb


module load htseq
module load samtools


######################################################################################################################
################################################ INPUT ###########################################################
#####################################################################################################################

REFERENCEGFF="/path/to/your/referencegff"

DIR="/path/to/your/dir"

Species="Dsim"

### Direction
cd ${DIR}

mkdir count_TE

### You need to modify the gff file to be compatible with htseq
awk '{for(i=0;++i<=NF-3;) printf $i"\t" ; print $(NF-2)}'  "$REFERENCEGFF" > ${Species}.htseq.gff


#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
##################################################################################################################################################################### 

cd ${DIR}


#sample.txt : text file containing the list of samples

for SAMPLE in `cat sample.txt`
do



if [ -f "TE_count/${SAMPLE}.out" ]; then
	    echo "TE already count"
else	
	

htseq-count -m intersection-nonempty -f bam -t similarity -i Target  Mapping/${SAMPLE}_bwa_q0.bam ${Species}.htseq.gff > count_TE/${SAMPLE}.out

##################### Remove the simple repeat
grep -v ')n' count_TE/${SAMPLE}.out | grep -v '\-rich' | grep -v 'DODECA_SAT' | grep -v '15mer_SAT' > count_TE/${SAMPLE}_filter.out
##################### Simplify TE names
sed 's/Motif://g' count_TE/${SAMPLE}_filter.out > count_TE/${SAMPLE}_filter2.out
rm count_TE/${SAMPLE}_filter.out


##################### Normalize the count in RPM

##### Number of read mapping for normalization
COUNTq0=$(samtools view -c -F 260 Mapping/${SAMPLE}_bwa_q0.bam)
NORMq0=$(echo "scale=5;(1000000/$COUNTq0)"| bc -l)

##### Normalization
awk -v a=$NORMq0 '{print $1 "\t" $2*a}' count_TE/${SAMPLE}_filter2.out | head -n -5 > count_TE/${SAMPLE}_RPM.out
rm count_TE/${SAMPLE}_filter2.out

fi
done



############## Merging all the file in one table ##############
### Make the header to add to the final table
HEADER=$(ls *_RPM.out | sed 's/_RPM.out//'| sed 's/\n/\t/')
TE="TE"
HH="$TE $HEADER"

## Merge all the TE count file in one table
paste *_RPM.out | cut -f 1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46 > TE_Merge.out  
# Add an header
echo -e $HH | cat -  TE_Merge.out > TE_Merge_Header.out
# Change the space for tabulatin in the header to be able to open it properly in R later. 
sed 's/ /\t/g' TE_Merge_Header.out > TE_Merge_Header2.out

######################################################################
