
#############
# Libraries #
#############

library(tidyverse)
library(readxl)
library(ComplexHeatmap)
library(dendsort)
library(fgsea)
library(circlize)
library(RColorBrewer)

###################
# Expression file #
###################

raw <- as.data.frame(read_xlsx("hsc3_t47d_counts.xlsx"))

raw <- distinct(raw, Name, .keep_all = T)
rownames(raw) <- raw$Name

raw <- dplyr::select(raw, contains("CPM"))

raw <- dplyr::select(raw, contains("HSC3"))

colnames(raw) <- gsub("- CPM", "", colnames(raw))
 

colnames(raw) <- c("HSC3_KO_1", "HSC3_KO_3", "HSC3_KO_6", "HSC3_WT_1", 
                    "HSC3_WT_2", "HSC3_WT_3")

raw$Gene <- rownames(raw)

genes <- as.data.frame(gmtPathways("gene_lists/HALLMARK_PI3K_AKT_MTOR_SIGNALING.v2023.2.Hs (1).gmt"))

housekeeping <- c("GAPDH", "B2M", "ACTB", "TBP", "RP",
                  "RRN18S", "DIMT1", "TUBA1A", "GUSB", "PGK1", "RPL13A",
                  "PUM1", "DHX9", "MZT2B", "UBXN4", "LARP1", "TAF2", "STX5",
                  "SYMPK", "TMEM11", "SDHA", "TFRC", "HMBS")

hs <- subset(raw, (Gene %in% genes$HALLMARK_PI3K_AKT_MTOR_SIGNALING))
hs$Gene <- NULL

###########
# Heatmap #
###########

# Compute row variances
row_variances <- apply(hs, 1, var)

# Identify rows with zero variance
zero_variance_rows <- which(row_variances == 0)

# Display zero variance rows
print(zero_variance_rows)


# Compute column variances
column_variances <- apply(hs, 2, var)

# Identify columns with zero variance
zero_variance_columns <- which(column_variances == 0)

# Display zero variance columns
print(zero_variance_columns)

#Remove rows with zero variance
rows_to_remove <- c("CXCR4","ADCY2", "GNA14", "FGF6", "PRKCB", "PLCB1")
hs <- hs[!rownames(hs) %in% rows_to_remove,]

#scaled
hs_scaled <-t(scale(t(hs)))
h <- hs_scaled

#Log2 transformed
h <- as.matrix(log2(1+hs))


h_breaks <- seq(min(h), max(h), length.out = 10)
h_breaks


#Reposition the breaks in quantile positions
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}


h_breaks <- quantile_breaks(h, n = 11)
h_breaks


row_dist <- as.dist(1-cor(t(h), method = "pearson"))
row_hc <- hclust(row_dist, method = "complete")
Rowv=dendsort(as.dendrogram(row_hc), isRevers=TRUE, type = "average") 


col_fun <- colorRamp2(quantile_breaks(h, n = 11), c("#1e90ff", "#0479f5", "#0467d2", "#0356af", "#02458c","#023369","#012246", "#130202", "#6b0002","#b20004","#ff0000"))



ht <- Heatmap(h,col = col_fun,
              cluster_columns = F,
              name = "Expression Values",
              show_heatmap_legend = T ,
              #top_annotation = ha,
              #left_annotation = hr,
              show_column_names = T,
              show_row_names = F,
              cluster_rows = Rowv,
              row_title_gp = gpar(fontsize=12),
              column_title_gp = gpar(fontsize=12),
              height = unit(12, "cm"),
              width = unit(10.55, "cm"),
              column_dend_reorder = F,
              row_dend_reorder = T,
              show_row_dend = T,
              border = T,
              column_dend_height = unit(2, "cm"),
              column_names_rot = 45,
              #legend_height = unit(4, "cm"),
              column_names_gp = grid::gpar(fontsize = 7,fontface="bold"),
              column_split = rep(c("KO", "KO", "KO", "WT", "WT", "WT")),
              heatmap_legend_param = list(title="Scaled Expression", legend_height=unit(3, "cm"),title_gp=gpar(fontsize=10, fontface="bold"),labels_gp = gpar(fontsize=8)))

draw(ht, merge_legend=TRUE, padding = unit(c(2, 2, 2, 2), "mm"))

png("heatmap_pi3k_akt_mtor_hallmark_wt_vs_ko_hsc3_cpm_scaled.png", res= 200, height = 1800, width = 1200)
print(ht)  
dev.off() 

pdf("heatmap_pi3k_akt_mtor_hallmark_wt_vs_ko_hsc3_cpm_scaled.pdf", height = 10, width = 8)
print(ht)
dev.off()
