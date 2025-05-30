library(dplyr)
library(ggpubr)
library(limma)
library(edgeR)
library(ComplexHeatmap)
library(data.table)
library(tidyverse)
################
## Count file ##
################

expr <- as.data.frame(fread("hsc_wt_ko_counts.csv"))
expr <- expr[,-c(2,3)]
expr <- distinct(expr,expr$Identifier,.keep_all = TRUE)
expr$`expr$Identifier` <- NULL


test <- expr[duplicated(expr$Identifier),] #Check for duplicates

expr <- expr[,-2]

expr <- column_to_rownames(expr, "Name")


#For hsc3 data
colnames(expr) <- c("HSC3_KO_1", "HSC3_KO_3", "HSC3_KO_6", "HSC3_WT_1", 
                    "HSC3_WT_2", "HSC3_WT_3")


###############
## Meta data ##
###############

meta <- as.data.frame(fread("metadata_hsc3.csv"))
meta <- meta[,-1]

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


ht <- pheatmap(cor_matrix,
               cluster_rows = F,
               cluster_cols = F,
               clustering_distance_rows = 'correlation',
               clustering_distance_cols = 'correlation',
               clustering_method = 'complete',
               #column_names_rot = 45,
               heatmap_legend_param = list(title="Correlation", legend_height=unit(3, "cm"),title_gp=gpar(fontsize=10, fontface="bold"),labels_gp = gpar(fontsize=8)))
ht

pdf("Cor_heatmap_hsc3_wt_vs_ko.pdf", width = 6,height = 5)
draw(ht, padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()
