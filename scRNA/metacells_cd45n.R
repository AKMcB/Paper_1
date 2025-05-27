#############
# Libraries #
#############

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(SeuratData)
library(hdWGCNA)
library(WGCNA)
set.seed(1234)

###################
# Read RDS object # 
###################

combined <- readRDS("annotated_all_samples_cd45n_batch_corr_epithelial_cells.RDS")

#Functions from https://divingintogeneticsandgenomics.com/post/how-to-do-gene-correlation-for-single-cell-rnaseq-data-part-1/ 
matrix_to_expression_df<- function(x, obj){
  df<- x %>%
    as.matrix() %>% 
    as.data.frame() %>%
    tibble::rownames_to_column(var= "gene") %>%
    tidyr::pivot_longer(cols = -1, names_to = "cell", values_to = "expression") %>%
    tidyr::pivot_wider(names_from = "gene", values_from = expression) %>%
    left_join(obj@meta.data %>% 
                tibble::rownames_to_column(var = "cell"))
  return(df)
}

get_expression_data<- function(obj, assay = "RNA", slot = "data", 
                               genes = NULL, cells = NULL){
  if (is.null(genes) & !is.null(cells)){
    df<- GetAssayData(obj, assay = assay, slot = slot)[, cells, drop = FALSE] %>%
      matrix_to_expression_df(obj = obj)
  } else if (!is.null(genes) & is.null(cells)){
    df <- GetAssayData(obj, assay = assay, slot = slot)[genes, , drop = FALSE] %>%
      matrix_to_expression_df(obj = obj)
  } else if (is.null(genes & is.null(cells))){
    df <- GetAssayData(obj, assay = assay, slot = slot)[, , drop = FALSE] %>%
      matrix_to_expression_df(obj = obj)
  } else {
    df<- GetAssayData(obj, assay = assay, slot = slot)[genes, cells, drop = FALSE] %>%
      matrix_to_expression_df(obj = obj)
  }
  return(df)
}

#######################
# Construct metacells #
#######################
#Create a slow for the results be saved in
combined <- SetupForWGCNA(
  combined,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "correlation" # the name of the hdWGCNA experiment
)

# construct metacells in each group.
combined <- MetacellsByGroups(
  seurat_obj = combined,
  group.by = c("pan_cancer_cluster", "sample_ID"), # specify the columns in seurat_obj@meta.data to group by
  k = 10, # nearest-neighbors parameter
  max_shared = 5, # maximum number of shared cells between two metacells
  ident.group = 'pan_cancer_cluster' # set the Idents of the metacell seurat object
)

##############################################
# Processing for the meta-cell seurat object #
##############################################

# extract the metacell object 
combined_metacell <- GetMetacellObject(combined)

#normalize
options(future.globals.maxSize = 2 * 1024^3) 
combined_metacell<- NormalizeData(combined_metacell,normalization.method = "LogNormalize", scale.factor = 10000)

#find variable features
combined_metacell<- FindVariableFeatures(combined_metacell, selection.method = "vst", nfeatures = 2000)

#scale data
all.genes <- rownames(combined_metacell)
combined_metacell <- ScaleData(combined_metacell, features = all.genes)

#Linear reduction:PCA
combined_metacell <- RunPCA(combined_metacell, features = VariableFeatures(object = combined_metacell))

DimPlot(combined_metacell, reduction = "pca", group.by = "pan_cancer_cluster", seed= 1234)

#Correct for batch effect
combined_metacell <- RunHarmony(combined_metacell, "sample_ID")

#harmony.embeddings <- Embeddings(combined_batch, reduction = "harmony")

p1 <- DimPlot(object = combined_metacell ,reduction = "harmony", pt.size = .1, 
              group.by = "sample_ID")
p2 <- VlnPlot(object = combined_metacell, features = "harmony_1", 
              group.by = "sample_ID",  pt.size = .1)
p1+p2

#################
# Cluster cells #
#################

combined_metacell <- combined_metacell %>% 
  FindNeighbors(reduction = "harmony") %>% 
  FindClusters(resolution = 0.5)

combined_metacell <- RunUMAP(combined_metacell, reduction = "harmony", dims = 1:20, seed.use = 1234)
Idents(combined_metacell)<- combined_metacell$pan_cancer_cluster

p1<- DimPlot(combined, reduction = "umap", label = FALSE, 
             repel = TRUE, group.by = "pan_cancer_cluster", seed = 1234) +
  ggtitle("single cell")


p2<- DimPlot(combined_metacell, reduction = "umap", label = FALSE, 
             repel = TRUE,group.by = "pan_cancer_cluster", seed = 1234 ) + 
  ggtitle("meta cell")


dim <- p1 + p2
dim

pdf("dimplot_metacells_vs_singlecells_pan_cluster_v2.pdf", height = 6, width = 12)
print(p1)
print(p2)
dev.off()

pdf("dimplot_metacells_vs_singlecells_pan_cluster_combined_v2.pdf", height = 6, width = 12)
print(dim)
dev.off()

#Featureplot 

f1 <- FeaturePlot(combined, features = c("TRIM32", "MYC")) &
  scale_colour_gradientn(colours = c("#cfe2f3","#cb3110"))
f1

f2 <- FeaturePlot(combined_metacell, features = c("TRIM32", "MYC"))&
  scale_colour_gradientn(colours = c("#cfe2f3","#cb3110"))
f2

pdf("featureplot_metacells_vs_singlecells_trim32_myc_v2.pdf", height = 6, width = 12)
print(f1)
print(f2)
dev.off()

#Function from https://www.biostars.org/p/458261/
vp_case1 <- function(gene_signature, file_name, test_sign){
  plot_case1 <- function(signature, y_max = NULL){
    VlnPlot(combined_metacell, features = signature,
            pt.size = 0.1, 
            group.by = "pan_cancer_cluster", 
            y.max = y_max # add the y-axis maximum value - otherwise p-value hidden
    ) + stat_compare_means(comparisons = test_sign, label = "p.signif")
  }
  plot_list <- list()
  y_max_list <- list()
  for (gene in gene_signature) {
    plot_list[[gene]] <- plot_case1(gene)
    y_max_list[[gene]] <- max(plot_list[[gene]]$data[[gene]]) # get the max no. for each gene
    plot_list[[gene]] <- plot_case1(gene, y_max = (y_max_list[[gene]] + 1) )
  }
  cowplot::plot_grid(plotlist = plot_list)
  file_name <- paste0(file_name, "_r.pdf")
  ggsave(file_name, width = 12, height = 6)
}

gene_sig <- c("TRIM32", "MYC")
comparisons <- list(c("Cancer", "Normal"))
vp_case1(gene_signature = gene_sig, file_name = "vlnplot_metacells_trim32_myc_cell_types_v3", test_sign = comparisons)

v1 <- VlnPlot(combined, features = c("TRIM32", "MYC"), group.by = "pan_cancer_cluster")
v1

pdf("vlnplot_singlecells_trim32_myc_cell_types_v2.pdf", height = 6, width = 12)
print(v1)
dev.off()

v2 <- VlnPlot(combined_metacell, features = c("TRIM32", "MYC"), group.by = "pan_cancer_cluster")
v2

pdf("vlnplot_metacells_trim32_myc_cell_types_v2.pdf", height = 6, width = 12)
print(v2)
dev.off()

pdf("vlnplot_metacells_vs_singlecells_trim32_myc_v2.pdf", height = 6, width = 12)
print(v1)
print(v2)
dev.off()

##########################
# Correlation: Metacells #
##########################
genes <- c("TRIM32", "MYC")

#single cells
expression_data<- get_expression_data(combined, genes = genes)

# https://github.com/LKremer/ggpointdensity
# ggpubr to add the correlation
library(ggpubr)
p1 <- ggscatter(expression_data, x ="TRIM32", y="MYC",
                title = "Single Cells",
                color= "black", shape=21, size=1, fill="#D73027",alpha = 0.4,
                add = "reg.line", 
                add.params = list(color = "#D73027", fill= "#D73027"), 
                conf.int = T,
                cor.coef = T,
                cor.coef.size = 5,
                cor.coef.args = list(method = "pearson", label.sep ="\n"),
                xlab = "TRIM32", ylab = "MYC") + 
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "italic"),
        axis.text.x = element_text(size = 18), axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 10))


p1 <- ggplot(expression_data, aes(x= TRIM32, y = MYC)) + 
  geom_smooth(method="lm", color ="black") +
  labs(title = "Single Cells")+
  geom_point(size = 0.8, color= "black") +
  facet_wrap(~pan_cancer_cluster) +
  ggpubr::stat_cor(method = "pearson")+ 
  theme_classic()

p2 <- ggplot(get_expression_data(combined_metacell, genes = genes), 
             aes(x= TRIM32, y = MYC)) + 
  labs(title = "Metacells")+
  geom_smooth(method="lm", color= "black") +
  geom_point(size = 0.8, color="black", aes(alpha= 0.5)) +
  facet_wrap(~pan_cancer_cluster) +
  ggpubr::stat_cor(method = "pearson")+
  theme_classic()

p2 <- ggscatter(get_expression_data(combined_metacell, genes = genes), x ="TRIM32", y="MYC",
                title = "Metacells",
                color= "black", shape=21, size=1, fill="#D73027",alpha = 0.4,
                add = "reg.line", 
                add.params = list(color = "#D73027", fill= "#D73027"), 
                conf.int = T,
                cor.coef = T,
                cor.coef.size = 5,
                cor.coef.args = list(method = "pearson", label.sep ="\n"),
                xlab = "TRIM32", ylab = "MYC") + 
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "italic"),
        axis.text.x = element_text(size = 18), axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 10))

p2

pdf("correlation_metacells_vs_singlecells_trim32_vs_MYC_v2.pdf", width = 5, height = 5)
print(p1)
print(p2)
dev.off()