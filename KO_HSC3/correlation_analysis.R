library(dplyr)
library(ggpubr)
library(limma)
library(magrittr)
library(edgeR)
library(readxl)
library(pheatmap)
library(data.table)
library(tidyverse)
################
## Count file ##
################

expr <- as.data.frame(fread("hsc_wt_ko_counts.csv"))
expr <- expr[,-c(2,3)]
expr <- distinct(expr,expr$Identifier,.keep_all = TRUE)
expr$`expr$Identifier` <- NULL

cell <- as.data.frame(fread("CCLE_RNAseq_genes_counts_20180929.gct.gz"))
cell <- dplyr::select(cell, contains(c("Name", "Description", "HSC3")))
cell$Name <- gsub("\\..*", "", cell$Name)

expr <- merge(expr, cell, by.x = "Identifier", by.y = "Name")
test <- expr[duplicated(expr$Identifier),] #Check for duplicates

expr <- expr[,-c(2,9)]

expr <- column_to_rownames(expr, "Identifier")

#For hsc3 data
colnames(expr) <- c("HSC3_KO_1", "HSC3_KO_3", "HSC3_KO_6", "HSC3_WT_1", 
                    "HSC3_WT_2", "HSC3_WT_3", "HSC3_CCLE")

###############
## Meta data ##
###############

meta <- as.data.frame(fread("metadata_hsc3.csv"))
meta <- meta[,-1]

ccle_info <- data.frame(name = "HSC3_CCLE",
                        treatment = "WT")

meta <- rbind(meta, ccle_info)

#get same patients in expr as the clinical file
filtered_colnames <- intersect(names(expr), meta$name)
expr <- expr[, filtered_colnames, drop = FALSE]


#Make sure the group info is in the same order as the expr file
meta <- column_to_rownames(meta, "name")
meta <- meta %>% arrange(treatment)
expr<- expr[,rownames(meta)]
all(rownames(meta) == colnames(expr))


treatment <- as.factor(meta$treatment)

##########################
## Correlation analysis ##
##########################

d0 <- DGEList(expr, group = treatment) #Converting to DGE list for edgeR
d0$genes <- data.frame(Symbol=rownames(d0))

#By default it filters out genes that has total counts below 10
keep.exprs <- filterByExpr(d0, group=treatment) 
d0 <- d0[keep.exprs,, keep.lib.sizes=FALSE]
dim(d0)

d0 <- calcNormFactors(d0, method = "TMM")

cpm_data <- cpm(d0, log = TRUE)

cor_matrix <- cor(cpm_data, method = "pearson")

library(ComplexHeatmap)
ht <- ComplexHeatmap::pheatmap(cor_matrix,
                               cluster_rows = F,
                               cluster_cols = F,
                               clustering_distance_rows = 'correlation',
                               clustering_distance_cols = 'correlation',
                               clustering_method = 'complete',
                               #column_names_rot = 45,
                               heatmap_legend_param = list(title="Correlation", legend_height=unit(3, "cm"),title_gp=gpar(fontsize=10, fontface="bold"),labels_gp = gpar(fontsize=8)))
ht

pdf("Cor_heatmap_hsc3_wt_vs_ko_ccle.pdf", width = 6,height = 5)
draw(ht, padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()


library(tidyverse)
dput(names(df))
dput(names(cor_matrix))

treated_samples <- colnames(cor_matrix)[meta$treatment == "KO"]
control_samples <- colnames(cor_matrix)[meta$treatment == "WT"]


average_cor_treated <- rowMeans(cor_matrix[, treated_samples])
average_cor_control <- rowMeans(cor_matrix[, control_samples])

# Create a dataframe
plot_data <- data.frame(
  Average_Correlation_Treated = average_cor_treated,
  Average_Correlation_Control = average_cor_control
)
plot_data$Sample <- rownames(plot_data)

#For hsc3 data
plot_data$Condition <- c("HSC3_KO_1", "HSC3_KO_3", "HSC3_KO_6", "HSC3_WT_1", 
                         "HSC3_WT_2", "HSC3_WT_3")
# dot plot
dotplot <- ggplot(plot_data, aes(x = Average_Correlation_Treated, y = Average_Correlation_Control, color = Condition)) +
  geom_point(size=5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
  xlab("Average Correlation with KO") +
  ylab("Average Correlation with WT") +
  theme_pubr()+
  theme(legend.text = element_text(size = 8))
dotplot

pdf("Cor_dotplot_hsc3_wt_vs_ko.pdf", height = 4, width = 6)
print(dotplot)
dev.off()
