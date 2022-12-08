library(ggplot2)
library(ggpubr)
library(agricolae)
library(datasets)
library(multcompView)
library(dplyr)


#Open the table
data=read.csv("Resitance.csv", sep=";")

### Calculate the total
data$total=data$femelle+data$males

### Calculate the femelle ratio
data$ratio=data$femelle/data$total

## transformation in arcsin
data$Trans=asin(sqrt(data$ratio))

###transform the Strain column in factor
data$Strain=as.factor(data$Strain)

#Reorder
data$Strain = with(data, reorder(Strain, Trans, mean))

################################
######## Anova test ########
################################

sink("anova.txt")
anova_test=summary(aov(data$Trans~data$Strain))
print(anova_test)
sink()

################################
######## Kruskal-Wallis test ########
################################
sink("KruskalWallis.txt")
KW=lapply(split(data, data$Strain), function(d) { kruskal.test(Trans~groupe, data=d) })
print(KW)
sink()

################################
######## Wilcoxon test ########
################################
sink("Wilcoxon.txt")
WT=lapply(split(data, data$Strain), function(d) { wilcox.test(d$Trans,split(data, data$Strain)$YST8$Trans)})
print(WT)
sink()


################################
######## Tukey test ########
################################

model=glm(Trans~Strain, data=data)
ANOVA=aov(model)
TUKEY=TukeyHSD(ANOVA)
cld <- multcompLetters4(ANOVA, TUKEY)


Tk <- group_by(data, Strain) %>%
  summarise(mean=mean(Trans), quant = quantile(Trans, probs = 0.75)) %>%
  arrange(desc(mean))

cld <- as.data.frame.list(cld$Strain)
cld$Strain=rownames(cld)
Tk$cld <- cld$Letters



################################
######## Plot ########
################################

col1=c(rep('red', 3), rep('red', 7), 'red', rep('red', 2), 'blue')
Name=c("SaoTome", "Cameroun", "France", "Seychelles 2003", "Morocco", "Mayotte", "Tanzania", "Egypt", "Madagascar", "Zimbabwe", "Hawaii", "New Caledonia", "Seychelles 1981", "Tunsia" )

pdf("Resistance_test.pdf", paper="a4r",width = 12, height = 12)
# Violin plot basic
p <- ggplot(data, aes(x=Strain, y=ratio, fill=Strain)) + 
  geom_violin(trim=TRUE) +
  scale_fill_manual(values=col1) +
  ylim(0.3, 1.1) +
  labs(y = "Proportion of female in the progeny") +
  stat_summary(fun="mean", geom="crossbar",width = 0.5,colour = "black") +
  geom_text(data=cld, aes(x=Strain, y=0.8, label = Letters), size = 3, vjust=-20) +
  theme(axis.line = element_line(colour = "black"),legend.position="none", panel.background = element_blank() )+
  scale_x_discrete(labels=Name) +
  theme(axis.text.x=element_text(color = "black", size=10, angle=30, hjust=0.9)) +
  geom_point(position = position_jitter(seed = 1, width = 0.1),  size=1)
p

dev.off()


################################
######## Summary Table ########
################################

WT_pvalue=signif(unlist(lapply(split(data, data$Strain), function(d) { wilcox.test(d$Trans,split(data, data$Strain)$YST8$Trans)$p.value })),digits=3)

KW_chi_squared=round(unlist(lapply(split(data, data$Strain), function(d) { kruskal.test(Trans~groupe, data=d)$statistic })),3)
KW_pvalue=round(unlist(lapply(split(data, data$Strain), function(d) { kruskal.test(Trans~groupe, data=d)$p.value })),3)

Ymin=round(unlist(lapply(split(data, data$Strain), function(d) min(d$ratio))),3)
Ymax=round(unlist(lapply(split(data, data$Strain), function(d) max(d$ratio))),3)
Ymean=round(unlist(lapply(split(data, data$Strain), function(d) mean(d$ratio))),3)

data_sum=cbind(Ymin, Ymax, Ymean, WT_pvalue, as.numeric(KW_chi_squared), KW_pvalue)

colnames(data_sum)=c("Ymin", "Ymax", "Ymean", "Wilcoxon pvalue", "Kruskal Wallis Chi Squared", "Kruskal Wallis p.value")

write.table(data_sum, "Statistique_Summary.csv", quote=F, sep=";", col.names = T)

