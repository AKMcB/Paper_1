#############
# Libraries #
#############

library(ggpubr)
library(ggplot2)
library(tidytext)
library(Hmisc)
library(data.table)
library(tidyverse)

########################
# Read expression file #
########################

pancan <- as.data.frame(fread("EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena"))

df <- subset(pancan, pancan$sample == "TRIM32")
rownames(df) <- df$sample
df$sample <- NULL
df <- tibble::rownames_to_column(as.data.frame(t(df)), "id")
df$id <- gsub('\\.',"-",df$id)


ann <- readr::read_tsv("Survival_SupplementalTable_S1_20171025_xena_sp")
ann <- ann[1:3]
colnames(ann)[1] <- "id"

merged <- merge(df, ann, by.x="id", by.y = "id")
merged <- merged %>% 
  filter(!grepl('-11', id))
merged <- merged[complete.cases(merged),]
merged$TRIM32 <- as.numeric(merged$TRIM32)
merged <- merged[order(merged$TRIM32, decreasing = TRUE),]
#head(merged)

colnames(merged)[4] <- "cancer_type"

fac <- with(merged, reorder(cancer_type, TRIM32, median, decreasing = T))
merged$cancer_type <- factor(merged$cancer_type, levels = levels(fac))


p <- ggboxplot(merged, x="cancer_type", y="TRIM32",outlier.shape = NA,
          fill ="cancer_type",
          palette = c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090", "#FFD92F", "#E0F3F8", "#ABD9E9",
                      "#74ADD1", "#4575B4", "#8DA0CB", "#E78AC3", "#762A83", "#9970AB", "#C2A5CF", "#E7D4E8",
                      "#F7F7F7", "#D9F0D3", "#A6DBA0", "#5AAE61", "#1B7837", "#4DAC26", "#80CDC1", "#35978F",
                      "#7FC97F", "#008837", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", 
                      "#A6761D", "#666666", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
                      "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A51A3", "#FFFF99", "#B15928", "#0571B0", "#CC4C02"), 
          ylab = "TRIM32 Expression", xlab = "Cancer Type") +
          guides(fill=guide_legend(title="Primary Disease"))+
          scale_x_reordered()+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face ="italic"), #Italic if it is a gene. 
        axis.text.x = element_text(vjust = 0.5,angle = 90, size = 10), axis.ticks.x=element_blank(), 
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 11),
        axis.text.y = element_text(size = 10))+
          stat_compare_means(method = "anova", 
                     label.x.npc = "center",
                     label.y.npc = "top", vjust = 1.0, hjust = 0.5)


p

pdf("PanCan_TRIM32_2.pdf", height = 4, width = 7)
print(p)
dev.off()

png("PanCan_TRIM32_2.png", res = 200, height = 1500, width = 1800)
print(p)
dev.off()

fwrite(merged, file = "TRIM32_exp_pancan_cancer_types.csv")



