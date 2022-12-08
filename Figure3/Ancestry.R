library(stringr)
library(tidyr)


########## ####### ####### ####### ###### 
### Script R permettant de récupérer de manière automatique l'allele correspondant dans le Subject
######## ####### ####### ####### ###### 


data=read.table("var.Y.flt.hmz_newID_v2.recode.SNP.Dmauritiana.out")

data2=NULL
for (i in 1:dim(data)[1]){
  aa=strsplit(as.vector(data$V7[i]),split="")
  aa=unlist(aa)
  Start=data$V6[i]
  bb <- rep(NA_integer_, length(aa))
  bb[aa != "-"] <-seq(from = Start, by = 1, length.out = sum(aa != "-"))
  query=which(bb==200)
  if (length(query) > 0) {
    subject = str_sub(as.character(data$V8[i]), query, query) # On recupère le SNP
    totsubject=c(as.character(data[i,1]), subject)
    data2=rbind(data2, totsubject)
  }
}


write.table(data2, "var.Y.flt.hmz_newID_v2.recode.SNP.Dmauritiana.csv", quote=F, sep="\t", row.names=F, col.names=F)



######## ####### ####### ####### ###### 
### Pour les simulans
######## ####### ####### ####### ###### 

data = read.table("var.Y.flt.hmz_newID_v2.recode.tab")
k=0
for (i in c(4:24)) {
  k=k+1
  data=separate(data, paste0("V",i), c(paste0("geno", k), paste0("phred",k)), sep="/")
}

data[,1]=paste0(data[,1], ":", data[,2]-200,"-", data[,2]+200)
data=data[,c(1,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44)]


colnames(data)= c("SNP", "Al110","Al4","BS5", "C1674","California", "Gth1", "Guyane","Hawai","Mada196", "Mar19","Mar216","Noumea","Oku11", "Rf45","SR", "ST8","Se1", "T2", "Val2", "Waph11", "Z613")
write.table(data, "var.Y.flt.hmz_newID_v2.SNP.Dsim.csv", quote=F, sep="\t", row.names=F, col.names=F)



######## ####### ####### ####### ###### 
### Merge des 4 espèces
######## ####### ####### ####### ###### 

dataSech=read.table("var.Y.flt.hmz_newID_v2.recode.SNP.Dsechellia.csv")
dataMau=read.table("var.Y.flt.hmz_newID_v2.recode.SNP.Dmauritiana.csv")
dataMelano=read.table("var.Y.flt.hmz_newID_v2.recode.SNP.Dmelanogaster.csv")
dataSim = read.table("var.Y.flt.hmz_newID_v2.SNP.Dsim.csv")
colnames(dataSech)=c("SNP", "Sech")
colnames(dataMau)=c("SNP", "Mau")
colnames(dataMelano)=c("SNP", "Melano")
colnames(dataSim)= c("SNP", "Al110","Al4","BS5", "C1674","California", "Gth1", "Guyane","Hawai","Mada196", "Mar19","Mar216","Noumea","Oku11", "Rf45","SR", "ST8","Se1", "T2", "Val2", "Waph11", "Z613")

data=merge(dataSim, dataMau, by="SNP", all.x=T)
data=merge(data, dataSech, by="SNP", all.x=T)
data=merge(data, dataMelano, by="SNP", all.x=T)

length(which(as.character(dataMau$Mau)!="-")) #6138
length(which(as.character(dataSech[,2])!="-")) #6126
length(which(as.character(dataMelano[,2])!="-")) #3909



write.table(data, "var.Y.flt.hmz_newID_v2.SNP.SimMauSechMel.csv", quote=F, sep="\t", row.names=T, col.names=F)



########### recreat a fasta alignment file with each individual. 


data[is.na(data)] <- "-"
dd=NULL
for(i in 2:24){
  data[,i]=str_replace(data[,i],'\\.', "-")
  aa=paste(">", colnames(data)[i], sep="")
  dd=rbind(dd,aa)
  seq=paste(data[,i], sep="", collapse="")
  dd=rbind(dd,as.vector(seq))
}

write.table(dd, "var.Y.flt.hmz_newID_v2.SNP.SimMauSech.fasta", quote=F, sep="\t", row.names=F, col.names=F)




# Number on allele share witht the resistant haplotype

Res=c("Al110","BS5","Gth1","Hawai","Mada196","Mar216","Noumea","Oku11", "Rf45","SR","Se1", "T2", "Z613")
Sen=c("Al4","California" , "ST8", "Val2" , "C1674", "Guyane", "Mar19", "Waph11")
aa=NULL
for(i in 1:dim(data)[1]){
  if(length(unique(as.character(data[i,Res])))==1 && length(unique(as.character(data[i,Sen])))==1 && unique(as.character(data[i,Res]))!=unique(as.character(data[i,Sen])) ){
    aa=c(aa,data$SNP[i])

  }
}

data_Sim=data[aa,c("Al110","Al4","Mau", "Sech", "Melano")]

length(which(data_Sim$Mau==data_Sim$Al110))# 370
length(which(data_Sim$Mau==data_Sim$Al4))# 15
length(which(data_Sim$Mau=="-"))# 155

length(which(data_Sim$Sech==data_Sim$Al110))# 364
length(which(data_Sim$Sech==data_Sim$Al4))# 20
length(which(data_Sim$Sech=="-"))# 152

length(which(data_Sim$Melano==data_Sim$Al110))# 223
length(which(data_Sim$Melano==data_Sim$Al4))# 15
length(which(data_Sim$Melano=="-"))# 301



