# load tidyverse package
library(tidyverse)

# read in data
pca <- read_table(paste0("./var.Y.flt.hmz.eigenvec"), col_names = FALSE)
eigenval <- scan(paste0("./var.Y.flt.hmz.eigenval"))

# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

## add phenotype color
COL=c("green", "blue", "blue", "blue", "green","green", "blue","blue","orange","purple","green","blue","red","purple","green","red", "purple","blue","blue","orange","green")
pca <- as_tibble(data.frame(pca, COL))

# first convert to percentage variance explained
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)

NAME=c("Al110", "Al4", "C1674","ST8_1", "T2", "Z613","Val2", "WAPH11", "Rf45", "Gth1", "Mar21", "Mar19", "SR_2", "SaoTome", "Oku11", "Noumea", "Hawai", "Californie", "Guyane", "Madagascare", "Se1") 

pdf("PCA.pdf")
# plot pca
b <- ggplot(pca, aes(PC1, PC2)) + geom_point(size = 3, col=COL, position = position_jitter(w = 0.01, h = 0.01))#
b <- b + scale_colour_manual(values = COL)
b <- b + coord_equal() + theme_light()
b <- b+ xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
b
dev.off()





