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

###############
# Subset gene #
###############

#TRIM32 = ENSG00000119401
#MAGED1 = ENSG00000179222

genes <- expr %>% select(c("id", "ENSG00000119401","ENSG00000179222"))

colnames(genes)[c(2,3)] <- c("TRIM32", "MAGED1")


p1 <- ggscatter(genes, x = "TRIM32", y = "MAGED1",
                     title = "TRIM32 vs MAGED1",
                     color = "black", shape = 21, size = 1.5,fill="#D73027",alpha=0.4, # Points color, shape and size
                     add = "reg.line",  # Add regressin line
                     add.params = list(color = "#D73027", fill = "#D73027"), # Customize reg. line
                     conf.int = TRUE, # Add confidence interval
                     cor.coef = TRUE,
                     cor.coef.size = 5,# Add correlation coefficient. see ?stat_cor
                     cor.coeff.args = list(method = "pearson", label.sep = "\n"),
                     xlab = "TRIM32", ylab = "MAGED1")+
      theme(legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face ="italic"), #Italic if it is a gene. 
      axis.text.x = element_text(size=10), axis.ticks.x=element_blank(), 
      axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 11),
      axis.text.y = element_text(size = 10))

p1


p2 <- ggscatter(genes, x = "MAGED1", y = "TRIM32",
          title = "TRIM32 vs MAGED1",
          color = "black", shape = 21, size = 1.5,fill="#D73027",alpha=0.4, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "#D73027", fill = "#D73027"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE,
          cor.coef.size = 5,# Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.sep = "\n"),
          xlab = "MAGED1", ylab = "TRIM32")+ 
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face ="italic"), #Italic if it is a gene. 
        axis.text.x = element_text(size=10), axis.ticks.x=element_blank(), 
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 11),
        axis.text.y = element_text(size = 10))

p2

png("trim32_vs_maged1_pearson.png", res= 200, height = 1800, width = 1200)
print(p1)  
dev.off() 

pdf("trim32_vs_maged1_pearson.pdf", height = 10, width = 8)
print(p1)
dev.off()

png("maged1_vs_trim32_pearson.png", res= 200, height = 1800, width = 1200)
print(p2)  
dev.off() 

pdf("maged1_vs_trim32_pearson.pdf", height = 10, width = 8)
print(p2)
dev.off()
