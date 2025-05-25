library(tidyverse)
library(fgsea)
library(EnhancedVolcano)
library(data.table)

###################
# Expression file #
###################

expr <- as.data.frame(fread("limma_trim32_WT_1_3_KO_3_6_v2.csv"))

genes_1 <- as.data.frame(gmtPathways("gene_lists/HALLMARK_MYC_TARGETS_V1.v2023.2.Hs (1).gmt"))
genes_2 <- as.data.frame(gmtPathways("gene_lists/HALLMARK_MYC_TARGETS_V2.v2023.2.Hs (1).gmt"))

colnames(genes_1) <- "gene"
colnames(genes_2) <- "gene"

comb <- rbind(genes_1, genes_2)
dup <-comb[duplicated(comb$gene),]
comb <- dplyr::distinct(comb, .keep_all = T)

expr_1 <- subset(expr, rownames(expr) %in% comb$gene)

fwrite(expr, "limma_trim32_WT_1_2_3_KO_1_3_6_myc_2025_02_26.csv", row.names = T)

expr_1 <- fread("limma_trim32_WT_1_2_3_KO_1_3_6_myc_2025_02_26.csv")
expr_1 <- expr_1[,-1]

test <- as.data.frame(expr_1$logFC>0)

v <- EnhancedVolcano(expr_1,
                     lab = expr_1$V1,
                     labSize = 6.0,
                     labCol = 'black',
                     labFace = 'bold',
                     title = "NOROC: High vs Low TRIM32",
                     x = 'logFC',
                     y = 'adj.P.Val',
                     pCutoff = 0.05,
                     colAlpha = 0.5,
                     FCcutoff = 1,
                     ylim = c(0,3), 
                     legendLabels = c('Not sig.', 'Log (base 2) FC', 'Adjusted p-value',
                                      'Adjusted p-value & Log (base 2) FC'),
                     caption = paste(
                       "Total significant genes (Padj): ",
                       sum(expr_1$adj.P.Val < 0.05), ", Up: ",
                       sum(expr_1$adj.P.Val < 0.05 & expr_1$logFC > 0), " Down: ",
                       sum(expr_1$adj.P.Val < 0.05 & expr_1$logFC < 0)))
v

png("volcanoplot_trim32_WT_1_3_KO_3_6_myc_2025_02_26.png", res = 200, height = 1800, width = 1500)
print(v)
dev.off()

pdf("volcanoplot_trim32_WT_1_3_KO_3_6_myc_2025_02_26.pdf", height = 8, width = 8)
print(v)
dev.off()

