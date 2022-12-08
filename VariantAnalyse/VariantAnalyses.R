library(tidyr)
library(gdsfmt)
library(SNPRelate)
library(vcfR)
library(adegenet)
library(ape)
library(pheatmap)
library(maSigPro)
library(stringr)

# Open table with datas informations        
sample=read.csv("sample.csv", header=T, sep=";")

# Open vcr file
vcf <- read.vcfR( "var.Y.flt.hmz.vcf", verbose = FALSE )
#convert vcf in genlight format
x <- vcfR2genlight(vcf)
ploidy(x) <-1
# Convert in matrix and add samples name
aa=as.matrix(x)
rownames(aa)=sample$Name


####################################################################
################### SupFig_Haplotypes ################################
####################################################################
# Change rownames of the matrix
rownames(aa)=c("Egy12", "Egy08", "Ken75", "Tun83", "Tan96", "Zim97", "Fra93", "Aus11", "May09", "Fra09", "Mor216", "Mor12", "Sey81", "Sao01", "Cam16", "New91", "Haw07","Cal61", "Guy56", "Mad98", "Sey03")
pdf("SupFig_Haplotype.pdf")
pheatmap(aa, cluster_cols =F,legend = F, show_colnames = F, treeheight_row = 0)
dev.off()



####################################################################
######## How many SNP are common to all Y resistants
####################################################################


########## By Phenotype ##########
### Add the Phenotype information
obj <- df2genind(aa, ploidy=1, pop=sample$Phenotype)

popG=genind2genpop(obj)

# Number of share SNP betweeen Sensitive
length(which(tab(popG)["Sen",]==8)) #6985
# Number of share SNP betweeen Resistante
length(which(tab(popG)["Res",]==13)) #1123
# Fix SNP between Resistant and Sensistive allele
length(which(tab(popG)["Res",]==13&tab(popG)["Sen",]==0)) #590



########## By Groups ##########

### Add the Groupe information
obj <- df2genind(aa, ploidy=1, pop=sample$Groupe)

popG=genind2genpop(obj)

gg=tab(popG)

# Number of share SNP betweeen Groupe1
length(which(tab(popG)["G1",]==8)) # 6985
# Number of share SNP betweeen Groupe2
length(which(tab(popG)["G2",]==6)) # 6191
# Number of share SNP betweeen Groupe3
length(which(tab(popG)["G3",]==2)) #7488
# Number of share SNP betweeen Groupe4
length(which(tab(popG)["G4",]==2)) # 7644
# Number of share SNP betweeen Groupe5
length(which(tab(popG)["G5",]==3)) #6704



####################################################################
################  Locilization of the SNPs ######## ######## 
####################################################################

par(mfrow=c(1,1))
#Raw number if SNP per Contig
pie(sort(summary(x@chromosome)), main="Raw Nb SNPs")

# Normalized by the Contig length
chromSize=read.csv("dsim_scaffold2_2019.chrom.sizes", sep="\t", header=F)
rownames(chromSize)=chromSize[,1]
Contig=as.data.frame(summary(x@chromosome))
NormContig=cbind(Contig,chromSize[rownames(Contig),2])
NormContig$Norm=NormContig[,1]/NormContig[,2]
NormContigSort=NormContig[order(NormContig[,3], decreasing=T),]
pie(NormContigSort$Norm, rownames(NormContigSort), main="Normalize Nb SNPs")


############################################
########### pairwise Distant
############################################

Dist=read.csv("MEGA-pairwise-distance.csv", sep=",", row.names = 1, header=T)
# calculate the percentage of similar SNPs
Dist=100-(Dist/8190)*100

# Average distance between Sensitives
Sen=c("Al4", "C1674", "ST8", "Val2", "Waph11", "Mar19", "California", "Guyane")
mean(as.numeric(unlist(Dist[Sen, Sen])), na.rm=T) # 97.76818

# Average distance between Resistants from Cluster II
Res1=c("Al110", "T2", "Z613", "Mar216", "Oku11", "Se1")
mean(as.numeric(unlist(Dist[Res1, Res1])), na.rm=T) # 97.78022

# Average distance between Resistants from Cluster III
Res2=c("SR", "Noumea")
mean(as.numeric(unlist(Dist[Res2, Res2])), na.rm=T) # 96.12943

# Average distance between Resistants from Cluster IV
Res3=c("Rf45", "Mada196")
mean(as.numeric(unlist(Dist[Res3, Res3])), na.rm=T) # 97.58242

# Average distance between Resistants from Cluster V
Res4=c("Gth1", "BS5", "Hawai")
mean(as.numeric(unlist(Dist[Res4, Res4])), na.rm=T) # 95.77534

# Average distance between all Resistants 
mean(as.numeric(unlist(Dist[c(Res1, Res2, Res3, Res4), c(Res1, Res2, Res3, Res4)])), na.rm=T) # 76.89834



