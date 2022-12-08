### Exon

# table with the number of copies of each exon 
exon=read.table("exon_number.csv", header=T, sep=";")

## Calculate the normalization factor base on the mean coverage
data=read.table('Count/tot_cov_Final.out')
Norm=colMeans(data[,-1])[-1]

################### Mean coverage and Normalization
Count_gene=NULL

for(i in c(1:82)){

Exon=as.character(exon[i,3])
Nb=as.integer(exon[i,4])

# Open the exon coverage
data2=read.table(paste0('Count/',Exon,'_cov.out'))
data2=data2[,3:25]

# Normalize
data2_norm=data2/Norm

# Mean by column and Multiply by the number of annotated exon. 
Count=colMeans(data2_norm)*Nb
Name=c('Egy12',	'Egy08',	'Ken75',	'Sey81',	'Tun83',	'Tan96',	'Zim97',	'Fra93',	'Aus11',	'May09',	'Fra09','Mor16','Mor12','Sey81b', 'Sao01','Cam16','New91','Haw07','Cal61','Guy56','Mad98','Tun83b','Sey03')
names(Count)=Name

Count2=Count[c('Tun83b','Fra93','Aus11','Mor12','Cal61','Guy56','Sey81b','May09','Fra09','Mor16','Sao01','Cam16','New91','Haw07','Mad98','Sey03')]
Count2=c(as.character(exon[i,1]), exon[i,2], Count2)

Count_gene=rbind(Count_gene, Count2)
}


####### Figure 4A

pdf("Figure4A.pdf")

dotchart(as.numeric(Count_gene[,3]), labels=Count_gene[,2], groups = as.factor(Count_gene[,1]), pch = 4, xlab="Exon number", cex=0.6, pt.cex=0.6, xlim=c(0,20))
for(i in c(3:8)){
  #wdy
  points(as.numeric(Count_gene[77:82,i]), 1:6, col ="blue", pch = 4, cex=0.8, xlab="Count")
  #ARY
  points(as.numeric(Count_gene[71:76,i]), 97:102, col ="blue", pch = 4, cex=0.8, xlab="Count")
  #Ppr-Y
  points(as.numeric(Count_gene[65:70,i]), 14:19, col ="blue", pch = 4, cex=0.8, xlab="Count")
  #CCY
  points(as.numeric(Count_gene[59:64,i]), 89:94, col ="blue", pch = 4, cex=0.8, xlab="Count")
  #Pp1-Y2
  points(as.numeric(Count_gene[58,i]), 22, col ="blue", pch = 4, cex=0.8, xlab="Count")
  #Pp1-Y1
  points(as.numeric(Count_gene[57,i]), 25, col ="blue", pch = 4, cex=0.8, xlab="Count")
  #ORY
  points(as.numeric(Count_gene[49:56,i]), 28:35, col ="blue", pch = 4, cex=0.8, xlab="Count")
  #PRY
  points(as.numeric(Count_gene[46:48,i]), 9:11, col ="blue", pch = 4, cex=0.8, xlab="Count")
  #kl5
  points(as.numeric(Count_gene[29:45,i]), 38:54, col ="blue", pch = 4, cex=0.8, xlab="Count")
  #kl3
  points(as.numeric(Count_gene[13:28,i]), 57:72, col ="blue", pch = 4, cex=0.8, xlab="Count")
  #kl2
  points(as.numeric(Count_gene[1:12,i]), 75:86, col ="blue", pch = 4, cex=0.8, xlab="Count")
  
}

for(i in c(9:18)){
  #wdy
  points(as.numeric(Count_gene[77:82,i]), 1:6, col ="red", pch = 4, cex=0.8, xlab="Count")
  #ARY
  points(as.numeric(Count_gene[71:76,i]), 97:102, col ="red", pch = 4, cex=0.8, xlab="Count")
  #Ppr-Y
  points(as.numeric(Count_gene[65:70,i]), 14:19, col ="red", pch = 4, cex=0.8, xlab="Count")
  #CCY
  points(as.numeric(Count_gene[59:64,i]), 89:94, col ="red", pch = 4, cex=0.8, xlab="Count")
  #Pp1-Y2
  points(as.numeric(Count_gene[58,i]), 22, col ="red", pch = 4, cex=0.8, xlab="Count")
  #Pp1-Y1
  points(as.numeric(Count_gene[57,i]), 25, col ="red", pch = 4, cex=0.8, xlab="Count")
  #ORY
  points(as.numeric(Count_gene[49:56,i]), 28:35, col ="red", pch = 4, cex=0.8, xlab="Count")
  #PRY
  points(as.numeric(Count_gene[46:48,i]), 9:11, col ="red", pch = 4, cex=0.8, xlab="Count")
  #kl5
  points(as.numeric(Count_gene[29:45,i]), 38:54, col ="red", pch = 4, cex=0.8, xlab="Count")
  #kl3
  points(as.numeric(Count_gene[13:28,i]), 57:72, col ="red", pch = 4, cex=0.8, xlab="Count")
  #kl2
  points(as.numeric(Count_gene[1:12,i]), 75:86, col ="red", pch = 4, cex=0.8, xlab="Count")
  
}

dev.off()



######## Figure4B - CK2B

data2=read.table(paste0('Count/PCKR_mod_cov.out'))
Nb=108

# Open the exon coverage
data2=data2[,3:25]

# Normalize
data2_norm=data2/Norm

# Multiply by the number of annotated exon. 
Count=colMeans(data2_norm)*Nb
Name=c('Egy12',	'Egy08',	'Ken75',	'Sey81',	'Tun83',	'Tan96',	'Zim97',	'Fra93',	'Aus11',	'May09',	'Fra09','Mor16','Mor12','Sey81b', 'Sao01','Cam16','New91','Haw07','Cal61','Guy56','Mad98','Tun83b','Sey03')
names(Count)=Name

Count2_CK2B=Count[c('Tun83b','Fra93','Aus11','Mor12','Cal61','Guy56','Sey81b','May09','Fra09','Mor16','Sao01','Cam16','New91','Haw07','Mad98','Sey03')]
Groups=c(rep("Sen", 6), rep("Res",10))

CK2B_Sen=Count2[c(1:6)]
#boxplot(CK2B_Sen)

vioplot(Count2~Groups, col=c("red", "blue"), xlab="Phenotype", ylab="Copy number estimate")
stripchart(Count2~Groups,              # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)  


######## Figure 4B - SRPK2

data2=read.table(paste0('Count/SRPK2_mod_cov.out'))
Nb=68

# Open the exon coverage
data2=data2[,3:25]

# Normalize
data2_norm=data2/Norm

# Multiply by the number of annotated exon. 
Count=colMeans(data2_norm)*Nb
Name=c('Egy12',	'Egy08',	'Ken75',	'Sey81',	'Tun83',	'Tan96',	'Zim97',	'Fra93',	'Aus11',	'May09',	'Fra09','Mor16','Mor12','Sey81b', 'Sao01','Cam16','New91','Haw07','Cal61','Guy56','Mad98','Tun83b','Sey03')
names(Count)=Name

Count2_SRPK2=Count[c('Tun83b','Fra93','Aus11','Mor12','Cal61','Guy56','Sey81b','May09','Fra09','Mor16','Sao01','Cam16','New91','Haw07','Mad98','Sey03')]
Groups=c(rep("Sen", 6), rep("Res",10))

boxplot(Count2~Groups, col=c("red", "blue"), xlab="Phenotype", ylab="Copy number estimate")
stripchart(Count2~Groups,              # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)  




######## Figure 4B - SRPK1
data2=read.table(paste0('Count/SRPK1_mod_cov.out'))
Nb=108

# Open the exon coverage
data2=data2[,3:25]

# Normalize
data2_norm=data2/Norm

# Multiply by the number of annotated exon. 
Count=colMeans(data2_norm)*Nb
Name=c('Egy12',	'Egy08',	'Ken75',	'Sey81',	'Tun83',	'Tan96',	'Zim97',	'Fra93',	'Aus11',	'May09',	'Fra09','Mor16','Mor12','Sey81b', 'Sao01','Cam16','New91','Haw07','Cal61','Guy56','Mad98','Tun83b','Sey03')
names(Count)=Name

Count2_SRPK1=Count[c('Tun83b','Fra93','Aus11','Mor12','Cal61','Guy56','Sey81b','May09','Fra09','Mor16','Sao01','Cam16','New91','Haw07','Mad98','Sey03')]
Groups=c(rep("Sen", 6), rep("Res",10))

boxplot(Count2~Groups, col=c("red", "blue"), xlab="Phenotype", ylab="Copy number estimate", ylim=c(50,180))
stripchart(Count2~Groups,              # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)  



Strain=rep(names(Count2_CK2B),3)
Count=c(Count2_CK2B,  Count2_SRPK1, Count2_SRPK2)
Gene=c(rep("CK2B", 16), rep("SRPK1", 16), rep("SRPK2", 16))
Groups=c(rep("Sen", 6), rep("Res",10),rep("Sen", 6), rep("Res",10), rep("Sen", 6), rep("Res",10))

data=data.frame(Strain,Groups, Gene, Count)


boxplot(Count~Groups*Gene, data=data, col=c("red", "blue"))
stripchart(Count~Groups*Gene,              # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)  


ggplot(data,aes(Gene,Count,fill=Groups))+
  geom_boxplot(alpha=0.1) + 
  geom_point(aes(col=Groups),position=position_jitterdodge())



