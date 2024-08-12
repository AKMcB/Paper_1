#############
# Libraries #
#############
library(tibble)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(tidytext)
library(data.table)

########################
# Read expression file #
########################
t32 <- read.csv2("CCLE_TRIM32_exp_info_Cancer_Codes.csv", sep = ";", as.is = T, check.names = F)


t32 <- t32[complete.cases(t32), ] #Remove NA
dup<-t32[duplicated(t32$ID),] #Check for duplicate samples
t32 <- t32[!t32$OncotreeCode %in% c("ACC", "CHOL"),] #Only one sample with these codes, have to be removed

t32$TRIM32<- as.numeric(t32$TRIM32)
t32 <-t32[order(t32t$TRIM32, decreasing = TRUE),]
str(t32)


fac <- with(t32, reorder(OncotreeCode, TRIM32, median, order = TRUE))
t32$OncotreeCode <- factor(t32$OncotreeCode, levels = levels(fac))

t32 <- t32[!t32$Differentiation %in% c("Mixed", "Non-Epithelial"),] #We are only interested in epthelial cells

colnames(t32)[8] <- "Cancer Code"

p <- ggboxplot(t32, x="Cancer Code", y="TRIM32", outlier.shape = NA,
               palette = c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090", "#FFD92F", "#E0F3F8", "#ABD9E9",
                           "#74ADD1", "#4575B4", "#8DA0CB", "#E78AC3", "#762A83", "#9970AB", "#C2A5CF", "#E7D4E8",
                           "#F7F7F7", "#D9F0D3", "#A6DBA0", "#5AAE61", "#1B7837", "#4DAC26", "#80CDC1", "#35978F",
                           "#7FC97F", "#008837", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", 
                           "#A6761D", "#666666", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
                           "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A51A3", "#FFFF99", "#B15928", "#0571B0", "#CC4C02"), 
               ylab = "TRIM32 Expression",
               ggtheme = theme_pubr(legend = "right")) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face ="italic"), #Italic if it is a gene. 
        axis.text.x = element_text(size=10, angle = 90), axis.ticks.x=element_blank(), 
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 10))+
  geom_jitter(shape = 21,stroke=0.5,size=2, aes(fill= t32$`Cancer Code`, alpha=0.5), position = position_jitter(width = 0.3, height = 0.5))+
  scale_x_reordered()+ 
  stat_compare_means(method = "anova", 
                     label.x.npc = "center",
                     label.y.npc = "top", vjust = 1.0, hjust = 0.5)
p

pdf("TRIM32_exp_CCLE_2.pdf", height = 4, width = 7)
print(p)
dev.off()

png("TRIM32_exp_CCLE_2.png", res = 200, height = 1500, width = 1800)
print(p)
dev.off()

fwrite(t32, file= "CCLE_TRIM32_exp_info_Cancer_Codes_filtered.csv")
