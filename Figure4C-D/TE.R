library(devtools)
library(ggbiplot)
library(coexnet)
library(ggplot2)
library(ggpubr)
library(ggforce)
library(matrixTests)
library(RColorBrewer)
library(pheatmap)
library(dplyr)
library(tidyr)

data=read.csv("TE_Merge.out", sep="\t", header=T, dec=".", row.names = 1)


#### Isolate Batch 2
Batch2=c("Sexchro.1", "Sexchro.2", "Sexchro.6", "Sexchro.12", "Sexchro.13", "Sexchro.15", "Sexchro.7", "Sexchro.10", "Sexchro.3", "Sexchro.4","Sexchro.5", "Sexchro.8", "Sexchro.9", "Sexchro.11", "Sexchro.14", "Sexchro.16")
data.batch2=data[,Batch2]

#### Remove TE for which none of the strain have more than 10 RPM
data.batch2$Max=apply(data.batch2, 1, max)
data2=data.batch2[which(data.batch2$Max>10),1:16]


######### PCA on batch2

Groups=c(rep("I",6),rep("III",2), "IV", "V", "II", "V", "II", "V", "IV", "II")

data2.pca <- prcomp(as.matrix(t(data2)), center = TRUE,scale. = TRUE)
summary(data2.pca)

ggbiplot(data2.pca, var.axes=FALSE, obs.scale = 1, var.scale = 1, groups=Groups, scale_color_discrete=Col) + scale_colour_manual(values = c("blue", "green", "red", "orange", "purple"))


########################################################
######### Plot of most variable elements Figure 4D##############
########################################################

########## Order base on the coefficient of variation
data3=data2
Var=cofVar(as.matrix(data3))
data3_ord=data3[order(Var$cv, decreasing=F),]
aa=log10(data3_ord[145:169,])

TEnames=rownames(aa)
TEnames=replace(TEnames, 8, "new_sim_repeat")

pdf("Figure4D.pdf", paper="a4r",width = 12, height = 12)
dotchart(aa[,1],labels=TEnames, bg ="blue", pch = 17,cex=0.8, xlab="log10(RPM)")
points(aa[,7], 1:nrow(aa), col ="red",  cex=0.8, pch = 19)
points(aa[,8], 1:nrow(aa), col ="red",  cex=0.8, pch = 19)
points(aa[,9], 1:nrow(aa), col ="red",  cex=0.8, pch = 19) 
points(aa[,10], 1:nrow(aa), col ="red",  cex=0.8, pch = 19) 
points(aa[,11], 1:nrow(aa), col ="red",  cex=0.8, pch = 19) 
points(aa[,12], 1:nrow(aa), col ="red",  cex=0.8, pch = 19) 
points(aa[,13], 1:nrow(aa), col ="red",  cex=0.8, pch = 19) 
points(aa[,14], 1:nrow(aa), col ="red",  cex=0.8, pch = 19) 
points(aa[,15], 1:nrow(aa),col ="red",  cex=0.8, pch = 19)
points(aa[,16], 1:nrow(aa), col ="red",  cex=0.8, pch = 19) 
points(aa[,1], 1:nrow(aa), col ="blue",  cex=0.8, pch = 17)
points(aa[,2], 1:nrow(aa), col ="blue", cex=0.8,  pch = 17)
points(aa[,3], 1:nrow(aa), col ="blue",  cex=0.8, pch = 17)
points(aa[,4], 1:nrow(aa), col ="blue",  cex=0.8, pch = 17)
points(aa[,5], 1:nrow(aa), col ="blue",  cex=0.8, pch = 17)
points(aa[,6], 1:nrow(aa), col ="blue",  cex=0.8, pch = 17)
dev.off()



################################################
############## Wilcoxon test ###########
############################################

Res=c("Sexchro.10", "Sexchro.7","Sexchro.3", "Sexchro.4","Sexchro.5", "Sexchro.8", "Sexchro.9", "Sexchro.11", "Sexchro.14", "Sexchro.16")
Sen=c("Sexchro.1", "Sexchro.2", "Sexchro.6", "Sexchro.12", "Sexchro.13", "Sexchro.15")

WT=row_wilcoxon_twosample(data2[,Res], data2[,Sen])
WT$adjP = p.adjust(WT$pvalue,"fdr") # correct pvalue
data_RvS=data2[which(WT$adjP<0.01),]
dim(data_RvS) # None



