#############
# Libraries #
#############

library(Seurat)
library(tidyverse)
library(limma)
library(clustree)
library(cowplot)
library(sctransform)
library(SeuratWrappers)
library(scCustomize)
library(SCpubr)
library(patchwork)
library(reticulate)
library(harmony)
library(msigdbr)
library(SingleCellExperiment)
library(scDblFinder)
library(optparse)

############
# Optparse #
############

option_list <- list(
  make_option(c("f", "--sample_file"),
              default = NULL,
              help = "the main directory of the folders"))

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

SAMPLE_FILE <- opt$sample_file
SAMPLE_FILE <- "GSE164690_RAW/"

########################
# Create Seurat Object #
########################

main_folder <- list.files(path = SAMPLE_FILE, pattern = "HN")
folders <- c("PBL", "CD45p", "CD45n")
high_feature_threshold <- 5000 # From Article
low_feature_threshold <- 200 #From Article


#create an empty list to store objects
seurat_objects <- list()

for (i in 1:length(main_folder)) {
  for (i2 in 1:length(folders)) {
    sample_name <- paste(main_folder[i], sep = "_", folders[i2])
    high_thresh <- 5000
    
    data_dir <- file.path(paste(SAMPLE_FILE, main_folder[i],folders[i2], sep = "/"))
    data <- Read10X(data.dir = data_dir)
    project_name <- sample_name
    seurat_obj <- CreateSeuratObject(counts = data, project = project_name)
    seurat_obj <- Add_Mito_Ribo(seurat_obj, species = "human")
    
    median_percent_MT <- median(seurat_obj$percent_mito)
    mad_percent_MT <- mad(seurat_obj$percent_mito)
    high_threshold_MT <- median_percent_MT + 3 * mad_percent_MT
    
    assign(
      paste0(sample_name, "_qc"),
      QC_Plot_UMIvsGene(
        seurat_object = seurat_obj,
        meta_gradient_name = "percent_mito",
        meta_gradient_low_cutoff = high_threshold_MT, 
        low_cutoff_gene = low_feature_threshold,
        high_cutoff_gene = high_thresh
      )
    )
    
    #subset object based on quality control values
    seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > low_feature_threshold & nFeature_RNA < high_thresh & percent_mito < high_threshold_MT)
    seurat_obj$sample <- sample_name
    
    #identify and remove the doublets
    sce <- as.SingleCellExperiment(DietSeurat(seurat_obj))
    sce <- scDblFinder(sce)
    seurat_obj$Multiplet <- sce$scDblFinder.class
    seurat_obj <- subset(seurat_obj, subset = Multiplet == "singlet")
    
    #store object in the list
    seurat_objects[[sample_name]] <- seurat_obj }
}

merged_seurat <- Reduce(function(x,y) merge(x,y), seurat_objects)

save(merged_seurat,file= "merged_all_samples.RData")

load("merged_all_samples.RData")

#Merge only the cd45p and cd45n object 
selected_names <- names(seurat_objects)[grepl("CD45p|CD45n", names(seurat_objects))]
selected_names