library(RColorBrewer)
library(matrixTests)
library(devtools)
library(ggbiplot)
library(coexnet)
library(ggplot2)
library(ggpubr)
library(ggforce)
library(rpgm)


# open table
data=read.csv("Y_all.rep.compiled", sep="\t")

# Change column and row names 
name=strsplit(colnames(data), split=".", fixed=T)
name=sapply(name, "[[", 1)
colnames(data)=name

Strain=gsub(".rep.rpm.total", "", data$lines)
row.names(data)=Strain
data=data[,-1]
data=data[,-dim(data)[2]]

############# Isolate batch2

Batch2=c("Sexchro-1", "Sexchro-2", "Sexchro-6", "Sexchro-12", "Sexchro-13", "Sexchro-15", "Sexchro-7", "Sexchro-10", "Sexchro-3", "Sexchro-4","Sexchro-5", "Sexchro-8", "Sexchro-9", "Sexchro-11", "Sexchro-14", "Sexchro-16")

Groups=c(rep("I",6),rep("III",2), "IV", "V", "II", "V", "II", "V", "IV", "II")
Col=c(rep("blue",6),rep("red",2), "orange", "purple", "green", "purple", "green", "purple", "orange", "green")

data2=data[Batch2,]

Sums=data2%>%summarise_if(is.numeric, max)
data2=data2[,-which(Sums<5)] # remove kemers with less than 5 counts
Sums[order(Sums, decreasing=T)] #order based on the count

##### PCA Figure 5A
data2.pca <- prcomp(as.matrix(data2), center = TRUE, scale. = TRUE)
ggbiplot(data2.pca, var.axes=FALSE, obs.scale = 1, var.scale = 1, groups=Groups, scale_color_discrete=Col) + scale_colour_manual(values = c("blue", "green", "red", "orange", "purple"))


################################################
############## Wilcoxon test ###########
############################################
Res=c("Sexchro-7", "Sexchro-10", "Sexchro-3", "Sexchro-4","Sexchro-5", "Sexchro-8", "Sexchro-9", "Sexchro-11", "Sexchro-14", "Sexchro-16")
Sen=c("Sexchro-1", "Sexchro-2", "Sexchro-6", "Sexchro-12", "Sexchro-13", "Sexchro-15")


WT=col_wilcoxon_twosample(data2[Res,], data2[Sen,], correct=FALSE)
data_RvS=data2[,which(WT$pvalue<0.01)] # isolate the kmers that are significantly differents betw Res and Sen


##########################################
############ Dot plot Figure 5C ############
##########################################
aa=log10(data_RvS) #log10 transformation
aa[aa==-Inf]=0 # remove the 0
aa=t(aa) #rotate the matrix

dotchart(aa[,1], bg ="blue", pch = 17,cex=0.8, xlab="log10(RPM)", xlim=c(-2, 4))
points(aa[,7], 1:nrow(aa), col ="red",  cex=0.8, pch = 19)
points(aa[,8], 1:nrow(aa), col ="red",  cex=0.8, pch = 19)
points(aa[,9], 1:nrow(aa), col ="red",  cex=0.8, pch = 19) #col ="yellow"
points(aa[,10], 1:nrow(aa), col ="red",  cex=0.8, pch = 19) # col ="magenta"
points(aa[,11], 1:nrow(aa), col ="red",  cex=0.8, pch = 19) #col ="green"
points(aa[,12], 1:nrow(aa), col ="red",  cex=0.8, pch = 19) #col ="magenta"
points(aa[,13], 1:nrow(aa), col ="red",  cex=0.8, pch = 19) #col ="green"
points(aa[,14], 1:nrow(aa), col ="red",  cex=0.8, pch = 19) #col ="magenta"
points(aa[,15], 1:nrow(aa),col ="red",  cex=0.8, pch = 19)# col ="yellow"
points(aa[,16], 1:nrow(aa), col ="red",  cex=0.8, pch = 19) #col ="green"

points(aa[,1], 1:nrow(aa), col ="blue",  cex=0.8, pch = 17)
points(aa[,2], 1:nrow(aa), col ="blue", cex=0.8,  pch = 17)
points(aa[,3], 1:nrow(aa), col ="blue",  cex=0.8, pch = 17)
points(aa[,4], 1:nrow(aa), col ="blue",  cex=0.8, pch = 17)
points(aa[,5], 1:nrow(aa), col ="blue",  cex=0.8, pch = 17)
points(aa[,6], 1:nrow(aa), col ="blue",  cex=0.8, pch = 17)


####################################################################
############## Sup Figure  : Correlation Resistance vs TE copy number ###########
#####################################################################


#Open the table
data.Res=read.csv("~/Dropbox/Y_Project/Phenotypage/Resitance.csv", sep=";")

### Calculate the total
data.Res$total=data.Res$femelle+data.Res$males

### Calculate the female ratio
data.Res$ratio=data.Res$femelle/data.Res$total

## transformation in arcsin
data.Res$Trans=asin(sqrt(data.Res$ratio))

###transform the Strain column in factor
data.Res$Strain=as.factor(data.Res$Strain)


MeanStrain = as.data.frame(tapply(data.Res$Trans, data.Res$Strain, mean))
MeanStrain=MeanStrain[c(1:10,12),]
Strain=c("Sexchro-8", "Sexchro-9", "Sexchro-4", "Sexchro-16", "Sexchro-5", "Sexchro-3", "Sexchro-14","Sexchro-11", "Sexchro-10", "Sexchro-7", "Sexchro-15")
MeanStrain=cbind(Strain,MeanStrain)
rownames(MeanStrain)=Strain

Kmers=colnames(data_RvS)
Resistance=MeanStrain[Strain,2]

Data_Cor=cbind(Resistance, data_RvS[Strain,])
Data_Cor=as.data.frame(Data_Cor)

Data_Cor<-gather(Data_Cor, key="measure", value="value", Kmers)

Data_Cor$value=as.numeric(as.character(Data_Cor$value))
Data_Cor$Resistance=as.numeric(as.character(Data_Cor$Resistance))

pdf("Sup_FigureXX_Correlation_kmers.pdf", paper="a4r",width = 12, height = 12)

#### Correlation plot
for(i in 1:2){
  g = ggscatter(Data_Cor, x="Resistance", y="value" ,
                add = "reg.line", conf.int = TRUE, size=1, 
                cor.coef = TRUE, cor.method = "spearman", cor.coef.size=3,
                xlab = "Resistance", ylab = "kmers CN")+
    facet_wrap_paginate(. ~measure, scales = "free",ncol = 4, nrow = 5, page=i)
  
  print(g)
}

dev.off()

################ ################ ################ 
################ Figure 5B ############## ###########
################ ################ ################ 

TARGET="AAACAAT"
par(mfrow=c(1,1))
x=unlist(data_RvS[Res,TARGET])
z=unlist(data_RvS[Sen,TARGET])
boxplot(x,z, col=c("red", "blue"), main=TARGET)

data_TARGET=Data_Cor[which(Data_Cor$measure==TARGET),]
data_TARGET$Resistance=as.numeric(as.character(data_TARGET$Resistance))

ggscatter(data_TARGET, x="Resistance", y="value",
          add = "reg.line", conf.int = TRUE, add.params = list(color = "blue", fill = "lightgray"),
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Resistance", ylab = "kmers CN")



