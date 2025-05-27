#############
# Libraries #
#############

library(tidyverse)
library(ggpubr)

###################
# Expression data #
###################

expr <- read.csv2("raw data/TMM_TCGA_HNSC_counts_NoNormal_log2_filtered.csv", sep = ";", as.is = T, check.names = F)
rownames(expr) <- expr[,1]
expr <- expr[,-1]
expr <- as.data.frame(t(expr))
expr <- rownames_to_column(expr, "id")

expr$id <- gsub("-01A", "",expr$id)
expr$id <- gsub("-01B", "",expr$id)

expr <- expr[!grepl("-06",expr$id),]

dup <- as.data.frame(duplicated(expr$id))

info <- read.csv2("raw data/Clinical_info_HNSC.csv", sep = ";", as.is = T, check.names = F)
head(info)

###############
# Subset gene #
###############

#TRIM32 = ENSG00000119401
#MYC =  ENSG00000136997

genes <- expr %>% select(c("id", "ENSG00000119401","ENSG00000136997"))

colnames(genes)[c(2,3)] <- c("TRIM32", "MYC")

merged <- merge(genes, info, by.x = "id", by.y = "Patient ID")


p1 <- ggscatter(merged, x = "TRIM32", y = "MYC",
                title = " TRIM32 vs MYC (n = 495)",
                color = "black", shape = 21, size = 2,fill="red",alpha=0.4, 
                add = "reg.line",  
                add.params = list(color = "#D73027", fill = "#D73027"), 
                conf.int = TRUE, 
                cor.coef = TRUE,
                cor.coef.size = 5,
                cor.coeff.args = list(method = "pearson", label.sep = "\n"),
                xlab = "TRIM32", ylab = "MYC")+ 
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face ="italic"),  
        axis.text.x = element_text(size=18), axis.ticks.x=element_blank(), 
        axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 18))

p1

png("trim32_vs_myc_pearson.png", res= 200, height = 1800, width = 1200)
print(p1)  
dev.off() 

pdf("trim32_vs_myc_pearson.pdf", height = 5, width = 5)
print(p1)
dev.off()