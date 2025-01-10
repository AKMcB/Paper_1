#libraries 
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
combined <- readRDS("annotated_all_samples_cd45n_batch_corr.RDS")

#Annotate samples before extraction
#SingleR: BPE
singler_ref <- BlueprintEncodeData()
table(singler_ref$label.main)
singler_ref.sub <- singler_ref[,singler_ref$label.main %in% c("Monocytes", "Macrophages", "DC",
                                                              "B-cells", "Fibroblasts", "CD8+ T-cells",
                                                              "CD4+ T-cells", "Epithelial cells","Endothelial cells", "NK cells", 
                                                              "Skeletal muscle", "Smooth muscle")]

singler_results_BPE_sub <- SingleR::SingleR(
  test = GetAssayData(combined,assay = 'SCT', slot = 'data'),
  ref = singler_ref.sub,
  labels = singler_ref.sub@colData@listData$label.main
)

plotScoreHeatmap(singler_results_BPE_sub)
combined@meta.data$cell_type_singler_BPE_sub <- singler_results_BPE_sub@listData$labels
DimPlot(combined, group.by = "cell_type_singler_BPE_sub")

############################
# Extract epithelial cells #
############################

combined <- combined[, grepl("Epithelial cells", combined@meta.data$cell_type_singler_BPE_sub)]

DimPlot(combined, group.by = "cell_type_singler_BPE_sub")

###################
## Normalization ##
###################
options(future.globals.maxSize = 2 * 1024^3) 
combined <- SCTransform(combined, return.only.var.genes = F)

############################
## Find variable features ##
############################

combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(combined), 10)

###########
# Scaling #
###########

all.genes <- rownames(combined)
combined <- ScaleData(combined, features = all.genes)

#Linear reduction 
combined <- RunPCA(combined, features = VariableFeatures(object = combined))

DimPlot(combined, reduction = "pca", group.by = "pan_cancer_cluster")

###########################
# Batch effect correction #
###########################

combined <- RunHarmony(combined, "sample_ID")

harmony.embeddings <- Embeddings(combined_batch, reduction = "harmony")

p1 <- DimPlot(object = combined, reduction = "harmony", pt.size = .1, group.by = "sample_ID")
p2 <- VlnPlot(object = combined, features = "harmony_1", group.by = "sample_ID",  pt.size = .1)
p1+p2

###################
## Cluster cells ##
###################

combined <- combined %>% 
  FindNeighbors(reduction = "harmony") %>% 
  FindClusters(resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(combined), 5)

combined <- RunUMAP(combined, reduction = "harmony", dims = 1:20)

p1 <- DimPlot(combined, reduction = "umap", group.by = "pan_cancer_cluster")
p1

pdf("batch_corr_dimplot_cd45n_epithlial_cancer_vs_normal.pdf", height = 8, width = 8)
print(p1)
dev.off()

saveRDS(combined,"cd45n_epithelial_cells.RDS")
#####################
# Feature & vlnplot #
#####################

f <- FeaturePlot(combined, features = c("TRIM32", "MYC"))
f
v <- VlnPlot(combined, features = c("TRIM32", "MYC"), group.by = "pan_cancer_cluster")
v

pdf("featureplot_trim32_myc_cancer_vs_normal_cd45n.pdf", width = 12, height = 8)
print(f)
print(v)
dev.off()

saveRDS(combined, "annotated_all_samples_cd45n_batch_corr_epithelial_cells.RDS")
