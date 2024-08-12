#############
# Libraries #
#############
library(tidyverse)
library(readxl)
library(ComplexHeatmap)
library(dendsort)
library(circlize)
library(RColorBrewer)

###################
# Expression file #
###################

raw <- as.data.frame(read_xlsx("hsc3_t47d_counts.xlsx"))

raw <- distinct(raw, Name, .keep_all = T)
rownames(raw) <- raw$Name
raw <- dplyr::select(raw, contains("CPM"))

raw <- dplyr::select(raw, contains("T47D"))

raw <- dplyr::select(raw, contains("HSC3"))

colnames(raw) <- gsub("- CPM", "", colnames(raw))
#raw <- raw[,-6]
colnames(raw) <- c("T47D_KO_15", "T47D_KO_21", "T47D_KO_29", "T47D_KO_33", 
                    "T47D_WT_1", "T47D_WT_2", "T47D_WT_3") 

colnames(raw) <- c("HSC3_KO_1", "HSC3_KO_3", "HSC3_KO_6", "HSC3_WT_1", 
                    "HSC3_WT_2", "HSC3_WT_3")

raw$Gene <- rownames(raw)

housekeeping <- c("GAPDH", "B2M", "ACTB", "TBP", "RP",
                  "RRN18S", "DIMT1", "TUBA1A", "GUSB", "PGK1", "RPL13A",
                  "PUM1", "DHX9", "MZT2B", "UBXN4", "LARP1", "TAF2", "STX5",
                  "SYMPK", "TMEM11", "SDHA", "TFRC", "HMBS")

hs <- subset(raw, (Gene %in% housekeeping))
hs$Gene <- NULL

###########
# Heatmap #
###########

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

pdf("Housekeeping_genes_hsc3_wt_vs_ko.pdf",  width=13,height=16)
ht <- Heatmap(h,col = col_fun,
              cluster_columns = F,
              name = "Expression Values",
              show_heatmap_legend = T ,
              #top_annotation = ha,
              #left_annotation = hr,
              show_column_names = T,
              show_row_names = T,
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
              heatmap_legend_param = list(title="Log2(1+CPM)", legend_height=unit(3, "cm"),title_gp=gpar(fontsize=10, fontface="bold"),labels_gp = gpar(fontsize=8)))

draw(ht, merge_legend=TRUE, padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()
