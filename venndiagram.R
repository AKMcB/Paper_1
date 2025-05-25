############
# Libaries #
############

library(tidyverse)
library(VennDetail)
library(ggvenn)
library(fgsea)
library(data.table)
library(ggpubr)

###################
# Read gene lists #
###################
#cluster row order form TCGA-HNSC heatmap 
tcga <- read.csv2("cluster_rows_order_Heatmap_myc_trim32_tcga_location_subtype_3_cluster_2025_02_27.csv")
tcga <- subset(tcga, tcga$Cluster == "cluster1")

#Cluster row order for noroc 
noroc <- read.csv2("2025_03_04_cluster_rows_trim32_expr_noroc_myc_target_genes.csv")
noroc <- subset(noroc, noroc$Cluster == "cluster1")

#top.table from KO vs WT TRIM32 HSC3
ko <- read.csv2("2025_03_05_limma_trim32_WT_1_3_KO_3_6.csv", sep=",")

genes_1 <- as.data.frame(gmtPathways("../gene_lists/HALLMARK_MYC_TARGETS_V1.v2023.2.Hs (1).gmt"))
genes_2 <- as.data.frame(gmtPathways("../gene_lists/HALLMARK_MYC_TARGETS_V2.v2023.2.Hs (1).gmt"))

colnames(genes_1) <- "gene"
colnames(genes_2) <- "gene"

comb <- rbind(genes_1, genes_2)
dup <-comb[duplicated(comb$gene),]
comb <- dplyr::distinct(comb, .keep_all = T)

ko <- subset(ko, ko$X %in% comb$gene)
ko$logFC <- as.numeric(ko$logFC)
col <- "logFC"
#subset for the negative values as they indicate an upragulation in WT 
sorted_df <- ko[order(ko[[col]]), ]
subset_df <- head(sorted_df, 50)

ko <- subset(ko, ko$logFC <0)
ko <- subset(ko, ko$adj.P.Val < 0.05) #onlt the significant genes

tcga <- tcga[,"GeneID"]
noroc <- noroc[,"GeneID"]
subset_df <- subset_df[,"X"]

##################
## Venn Diagram ##
##################

overlap_list <- list("TCGA-HNSC" = tcga, "NOROC" = noroc, "KO cs WT HSC3" = subset_df)
overlap_plot <- venndetail(overlap_list) #identify overlapping terms
overlap_df <- result(overlap_plot)

#fwrite(overlap_df, "2025_03_20_overlap_tcga_noroc_KO_v2.csv", row.names = TRUE)
x = list(tcga,noroc, subset_df)

names(x) <- c("TCGA-HNSC", "NOROC", "KO vs WT HSC3")

plot <- ggvenn(x,set_name_size = 3,
               fill_color = c("#E41A1C", "#377EB8", "#4DAF4A"),stroke_size = 0.5)

plot

dev.off()

ggexport(plot, filename = "2025_03_20_tcga_noroc_KO_venndiagram_myc_targets.png",res=200,width = 2500, height = 2000)
ggexport(plot, filename = "2025_03_20_tcga_noroc_KO_venndiagram_myc_targets.pdf",width = 15, height = 15)
