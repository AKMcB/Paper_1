#############
# Libraries #
#############

library(scATOMIC)
library(plyr)
library(dplyr)
library(data.table)
library(randomForest)
library(caret)
library(parallel)
library(reticulate)
library(Rmagic)
library(Matrix)
library(Seurat)
library(agrmt)
library(cutoff.scATOMIC)
library(copykat)
library(ggplot2)
library(harmony)

load("merged_all_samples.RData")

merged_seurat <- JoinLayers(merged_seurat)

# Subset based on either "CD45p" or "CD45n"
subset_seurat <- merged_seurat[, grepl("CD45n", merged_seurat@meta.data$orig.ident)]

# Create a new column 
subset_seurat@meta.data$sample_ID <- sub("_.*", "", subset_seurat@meta.data$sample)
subset_seurat@meta.data$sample_ID

#################
# Normalization #
#################

options(future.globals.maxSize = 2 * 1024^3) 
combined <- SCTransform(subset_seurat, return.only.var.genes = F)

##########################
# Find variable features #
##########################

combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(combined), 10)

# plot variable features
plot1 <- VariableFeaturePlot(combined)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

###########
# Scaling # 
###########

all.genes <- rownames(combined)
combined <- ScaleData(combined, features = all.genes)

#save(combined,file= "all_samples_cd45p_norm_scaled.RData")
#load("all_samples_cd45n_norm_scaled.RData")
rm(merged_seurat, subset_seurat, plot1, plot2)
#Linear reduction 
combined <- RunPCA(combined, features = VariableFeatures(object = combined))

DimPlot(combined, reduction = "pca", group.by = "sample_ID")

###########################
# Batch effect correction #
###########################

combined_batch <- RunHarmony(combined, "sample_ID")

harmony.embeddings <- Embeddings(combined_batch, reduction = "harmony")

p1 <- DimPlot(object = combined_batch, reduction = "harmony", pt.size = .1, group.by = "sample_ID")
p2 <- VlnPlot(object = combined_batch, features = "harmony_1", group.by = "sample_ID",  pt.size = .1)
p1+p2
combined <- combined_batch

#################
# Cluster cells #
#################

combined <- combined %>% 
  FindNeighbors(reduction = "harmony") %>% 
  FindClusters(resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(combined), 5)

combined <- RunUMAP(combined, reduction = "harmony", dims = 1:20)

p1 <- DimPlot(combined, reduction = "umap", group.by = "sample_ID")
p1

pdf("batch_corr_dimplot_cd45n_all_samples.pdf", height = 8, width = 8)
print(p1)
dev.off()

save(combined,file= "all_samples_cd45n_norm_scaled_clustered_batch_corr.RData")
#load("all_samples_cd45n_norm_scaled_clustered.RData")

# Identify the sample
highlight_sample <- "HN17"  # Replace this with your sample ID

# Create a new column in metadata
combined$highlight <- ifelse(combined$sample_ID == highlight_sample, highlight_sample, "Other")

# Plot using DimPlot
DimPlot(combined, group.by = "highlight", cols = c("red", "grey")) + 
  ggtitle(paste("Highlighting", highlight_sample))

########################
# Cell type Annotation #
########################
#scAtomic
sparse_matrix <- combined@assays$SCT$counts

cell_predictions <- run_scATOMIC(sparse_matrix)
gc()

result <- create_summary_matrix(prediction_list = cell_predictions,
                                use_CNVs = FALSE, modify_results = TRUE,
                                mc.cores = 8, raw_counts = sparse_matrix, 
                                min_prop = 0.5)


scAtomic_result <- as.data.frame(table(result$scATOMIC_pred))
write.csv2(scAtomic_result, "scatomic_annotation_results_batch_corr.csv")
patient_id <- combined@meta.data$sample_ID

tree_results_non_interactive <- scATOMICTree(
  predictions_list = cell_predictions,
  summary_matrix = result,
  interactive_mode = F,
  save_results = T)

combined <- AddMetaData(combined, result)

d <- DimPlot(combined, group.by = "scATOMIC_pred") 
FeaturePlot(combined, features = c("TRIM32", "MYC"))
saveRDS(combined, "annotated_all_samples.RDS")
pdf("batch_corr_scatomic_ann_cd45n_all_samples.pdf", height = 8, width = 12)
print(d)
dev.off()

#SingelR: HPCA 
library(SingleR)
library(celldex)

singler_ref <- HumanPrimaryCellAtlasData()
table(singler_ref$label.main)
singler_ref.sub <- singler_ref[,singler_ref$label.main %in% c("Monocyte", "Macrophage", "DC",
                                                              "B_cell", "Fibroblasts", "T_cells",
                                                              "Epithelial_cells", 
                                                              "Endothelial_cells", "NK_cell", 
                                                              "Smooth_muscle_cells")]

singler_results_HPCA_sub <- SingleR::SingleR(
  test = GetAssayData(combined, assay = 'SCT', slot = 'data'),
  ref = singler_ref.sub,
  labels = singler_ref.sub@colData@listData$label.main
)

prediction_results <- as.data.frame(table(singler_results_HPCA_sub$labels))
write.csv2(prediction_results, "predictions_hpca_all_samples_cd45n_batch_corr.csv")

plotScoreHeatmap(singler_results_HPCA_sub)
combined@meta.data$cell_type_singler_HPCA_sub <- singler_results_HPCA_sub@listData$labels

p1 <- DimPlot(combined, group.by = "cell_type_singler_HPCA_sub")
p2<- FeaturePlot(combined, features = c("TRIM32", "MYC"))

pdf("batch_corr_HPCA_ann_cd45n_all_samples.pdf", height = 8, width = 12)
print(p1)
print(p2)
dev.off()

subset_cell <- subset(combined, sample_ID == highlight_sample)

DimPlot(subset_cell, group.by = "cell_type_singler_HPCA_sub") +
  ggtitle(paste("Cell Annotations for", highlight_sample))

celltype_counts <- table(subset_cell$cell_type_singler_HPCA_sub) 

# Convert the table to a data frame 
celltype_df <- as.data.frame(celltype_counts)
colnames(celltype_df) <- c("CellType", "Count")

# bar plot
ggplot(celltype_df, aes(x = CellType, y = Count, fill = CellType)) +
  geom_bar(stat = "identity") +
  ggtitle(paste("Cell Type Distribution for", highlight_sample)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels if needed

library("SCpubr")

barplot <- do_BarPlot(combined,
           group.by = "sample_ID",
           split.by = "cell_type_singler_HPCA_sub",
           position = "fill")

pdf("batch_corr_hpca_ann_barplot_cd45n_all_samples.pdf", height = 8, width = 8)
print(barplot)
dev.off()

#SingleR: BPE
singler_ref <- BlueprintEncodeData()
table(singler_ref$label.main)
singler_ref.sub <- singler_ref[,singler_ref$label.main %in% c("Monocytes", "Macrophages", "DC",
                                                              "B-cells", "Fibroblasts", "CD8+ T-cells",
                                                              "CD4+ T-cells", "Epithelial cells","Endothelial cells", "NK cells", 
                                                              "Skeletal muscle", "Smooth muscle")]

singler_results_BPE_sub <- SingleR::SingleR(
  test = GetAssayData(combined, assay = 'SCT', slot = 'data'),
  ref = singler_ref.sub,
  labels = singler_ref.sub@colData@listData$label.main
)

plotScoreHeatmap(singler_results_BPE_sub)
combined@meta.data$cell_type_singler_BPE_sub <- singler_results_BPE_sub@listData$labels

p1 <- DimPlot(combined, group.by = "cell_type_singler_BPE_sub")
p2 <- FeaturePlot(combined, features = c("TRIM32", "MYC"))

pdf("batch_corr_bpe_ann_cd45n_all_samples.pdf", height = 8, width = 12)
print(p1)
print(p2)
dev.off()


barplot <- do_BarPlot(combined,
                      group.by = "sample_ID",
                      split.by = "cell_type_singler_BPE_sub",
                      position = "fill")


barplot
pdf("batch_corr_bpe_ann_barplot_cd45n_all_samples.pdf", height = 8, width = 8)
print(barplot)
dev.off()

prediction_results <- as.data.frame(table(singler_results_BPE_sub$labels))
write.csv2(prediction_results, "predictions_bpe_all_samples_cd45n_batch_corr.csv")


saveRDS(combined, "annotated_all_samples_cd45n_batch_corr.RDS")


f <- FeaturePlot(combined, features = c("CD4", "CD14", "CD19","CD69", "CD80", "CD86" ))
f
v <- VlnPlot(combined, features = c("CD4", "CD14", "CD19", "CD69", "CD80", "CD86"), group.by = "sample_ID")
v
ggsave(f, filename = "featureplot_immune_markers_cd45p_batch_corr", device = "png", height = 10, 
       width = 10, units = "in")   

ggsave(v, filename = "violinplot_immune_markers_cd45p_batch_corr", device = "png", height = 10, 
       width = 10, units = "in")   
