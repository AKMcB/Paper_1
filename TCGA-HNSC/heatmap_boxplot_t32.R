#############
# Libraries #
#############

library(tidyverse)
library(ggpubr)

###################
# Expression data #
###################

expr <- read.csv2("raw data/TMM_TCGA_HNSC_counts_NoNormal_log2_filtered_hgnc.csv", sep = ";", as.is = T,check.names = F)
expr <- expr[,-c(1,3)]

dup <-expr[duplicated(expr$hgnc_symbol)|duplicated(expr$hgnc_symbol, fromLast=TRUE),]
expr <- distinct(expr, hgnc_symbol, .keep_all = TRUE)
expr <- expr[,-1]

rownames(expr) <- expr[,1]
expr <- expr[,-1]
expr <- as.data.frame(t(expr))
expr <- rownames_to_column(expr, "id")


expr <- expr[!grepl("-06A",expr$id),]
expr$id <- gsub("-01A", "",expr$id)
expr$id <- gsub("-01B", "",expr$id)


dup <- as.data.frame(duplicated(expr$id))


expr <- expr %>% select(c("id", "TRIM32", "MYC"))

#######################
# Heatmap column info #
#######################

ann2 <- read.csv2("Heatmap/cluster_column_order_Heatmap__myc_trim32_tcga_location_subtype_3_cluster_2025_02_27.csv", as.is = T, check.names = F)
ann2 <- ann2[,-1]

merged <- merge(expr, ann2, by.x="id", by.y = "Patient_ID")

dup <- as.data.frame(duplicated(merged$id))

merged$Cluster<- factor(merged$Cluster, 
                        levels= c("cluster1", "cluster2", "cluster3"), 
                        labels = c("Cluster 1", "Cluster 2", "Cluster 3"))

my_comparisons <- list( c("Cluster 1", "Cluster 2"),
                        c("Cluster 1", "Cluster 3"),
                        c("Cluster 2", "Cluster 3")) 
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                    symbols = c("****", "***", "**", "*", "ns"))

merged_long <- pivot_longer(merged, cols = c(TRIM32, MYC), 
                            names_to = "gene", values_to = "expression")

p <- ggplot(merged, aes(x = Cluster, y = TRIM32, fill = Cluster))+
  geom_point(alpha=0.5,position = position_jitter(width = 0.3, height = 0.5), shape= 21, size= 6)+
  geom_boxplot(fill = "white", alpha = 0.8, outlier.shape = NA) +
  labs(y = expression(paste(italic("TRIM32")~"(log2+1)")), 
       x = "Sample Type", 
       title = expression(paste(italic("TRIM32"), " Expression TCGA-HNSC Heatmap Clusters")))+ 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face ="bold"), 
        axis.text.x = element_text(size=10), axis.ticks.x=element_blank(), 
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 11, face ="italic"), 
        axis.text.y = element_text(size = 10), 
        panel.background = element_rect(fill = "white",
                                        colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid",
                                 colour = "black")) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#A0D636","#DF2DE0","#333ED4"))+ 
  geom_signif(comparisons = my_comparisons,map_signif_level = T, textsize=5, y_position = c(6.0, 6.3, 6.6))
p


png("heatmap/trim32_expr_3_cluster_2025_02_28.png", res= 200, height = 1800, width = 1200)
print(p)  
dev.off() 

pdf("heatmap/trim32_expr_3_cluster_2025_02_28.pdf", height = 10, width = 8)
print(p)
dev.off()