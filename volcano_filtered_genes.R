library(tidyverse)
library(fgsea)
library(EnhancedVolcano)
library(data.table)
#####################
## Expression file ##
#####################

expr <- as.data.frame(fread("2025_03_07_linear_reg_noroc_trim32_high_vs_low.csv"))

genes_1 <- as.data.frame(gmtPathways("gene_lists/HALLMARK_MYC_TARGETS_V1.v2023.2.Hs (1).gmt"))
genes_2 <- as.data.frame(gmtPathways("gene_lists/HALLMARK_MYC_TARGETS_V2.v2023.2.Hs (1).gmt"))

colnames(genes_1) <- "gene"
colnames(genes_2) <- "gene"

comb <- rbind(genes_1, genes_2)
dup <-comb[duplicated(comb$gene),]
comb <- dplyr::distinct(comb, .keep_all = T)

expr_1 <- subset(top.table, rownames(top.table) %in% comb$gene)

fwrite(expr, "limma_trim32_WT_1_2_3_KO_1_3_6_myc_2025_02_26.csv", row.names = T)

v <- EnhancedVolcano(expr_1,
                     lab = rownames(expr_1),
                     labSize = 6.0,
                     labCol = 'black',
                     labFace = 'bold',
                     title = "NOROC: High vs Low TRIM32",
                     x = 'logFC',
                     y = 'adj.P.Val',
                     pCutoff = 0.05,
                     colAlpha = 0.5,
                     FCcutoff = 0.5,
                     ylim = c(0,3), 
                     legendLabels = c('Not sig.', 'Log (base 2) FC', 'Adjusted p-value',
                                      'Adjusted p-value & Log (base 2) FC'),
                     caption = paste(
                       "Total significant genes (Padj): ",
                       sum(expr_1$adj.P.Val < 0.05), ", Up: ",
                       sum(expr_1$adj.P.Val < 0.05 & expr_1$logFC > 0), " Down: ",
                       sum(expr_1$adj.P.Val < 0.05 & expr_1$logFC < 0)))
v

png("2025_03_07_volcano_noroc_tri.png", res = 200, height = 1800, width = 1500)
print(v)
dev.off()

pdf("figures/volcano/volvanoplot_trim32_WT_1_3_KO_3_6_PI3K_AKT_MTOR.pdf", height = 8, width = 8)
print(v)
dev.off()

