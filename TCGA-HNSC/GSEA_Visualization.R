#############
# Libraries #
#############

library(data.table)
library(fgsea)
library(ggplot2)
library(dplyr)
library(ggpubr)

#####################
# Read gsea results #
#####################

high<-as.data.frame(fread("gsea_report_for_High_1714027534091.tsv"))
low <- as.data.frame(fread("gsea_report_for_Low_1714027534091.tsv"))
combined <- rbind(high, low)

x <- combined
x <- x %>% arrange(desc(NES))

top <- head(x, n=20)
bottom <- tail(x, n=20)

x <- rbind(top, bottom)


plot <- ggplot(x %>% mutate(NAME = gsub("GOBP_", "", NAME)), aes(reorder(NAME, NES), NES)) +
  geom_col(aes(fill = `FDR q-val` < 0.05)) +
  coord_flip() +
  scale_fill_manual(values = c("#00BFC4", "#F8766D")) + 
  labs(x = "Pathway", y = "Normalized Enrichment Score",
       title = "High vs Low TRIM32") +
  theme(
    text = element_text(size = 15, face = "bold"),  # Adjust the size to your desired value
    axis.text.x = element_text(angle = 90, hjust = 1)) 

plot


dev.off()
#To save
ggexport(plot, filename = "GSEA_NOROC_High_vs_Low_HEV_GOBP.png",res=200,width = 3200, height = 2000)
ggexport(plot, filename = "GSEA_NOROC_High_vs_Low_HEV_GOBP.pdf",width = 25, height = 15)
