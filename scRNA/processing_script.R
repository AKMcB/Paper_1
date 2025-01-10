###############
## Libraries ##
###############

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

#Set variables from command line arguments 
SAMPLE_FILE <- opt$sample_file
SAMPLE_FILE <- "GSE164690_RAW/"
##########################
## Create Seurat Object ##
##########################

main_folder <- list.files(path = SAMPLE_FILE, pattern = "HN")
folders <- c("PBL", "CD45p", "CD45n")
high_feature_threshold <- 5000 # From Article
low_feature_threshold <- 200 #From Article


# Initialize an empty list to store Seurat objects
seurat_objects <- list()

# Loop over samples
for (i in 1:length(main_folder)) {
  for (i2 in 1:length(folders)) {
    sample_name <- paste(main_folder[i], sep = "_", folders[i2])
    high_thresh <- 5000
    
  data_dir <- file.path(paste(SAMPLE_FILE, main_folder[i],folders[i2], sep = "/"))
  data <- Read10X(data.dir = data_dir)
  project_name <- sample_name
  seurat_obj <- CreateSeuratObject(counts = data, project = project_name)
  
  # Add mitochondrial and ribosomal information
  seurat_obj <- Add_Mito_Ribo(seurat_obj, species = "human")
  
  # Calculate mitochondrial percentage thresholds
  median_percent_MT <- median(seurat_obj$percent_mito)
  mad_percent_MT <- mad(seurat_obj$percent_mito)
  high_threshold_MT <- median_percent_MT + 3 * mad_percent_MT
  
  # Generate QC plot
  assign(
    paste0(sample_name, "_qc"),
    QC_Plot_UMIvsGene(
      seurat_object = seurat_obj,
      meta_gradient_name = "percent_mito",
      meta_gradient_low_cutoff = high_threshold_MT,  # Using fixed low threshold from author
      low_cutoff_gene = low_feature_threshold,
      high_cutoff_gene = high_thresh
    )
  )
  
  # Subset Seurat object based on QC metrics, using MAD-based threshold for mitochondrial content
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > low_feature_threshold & nFeature_RNA < high_thresh & percent_mito < high_threshold_MT)
  seurat_obj$sample <- sample_name
  # Additional QC to identify and remove doublets
  sce <- as.SingleCellExperiment(DietSeurat(seurat_obj))
  sce <- scDblFinder(sce)
  seurat_obj$Multiplet <- sce$scDblFinder.class
  seurat_obj <- subset(seurat_obj, subset = Multiplet == "singlet")
  
  # Store the Seurat object in the list
  seurat_objects[[sample_name]] <- seurat_obj
  
  # Print completion message
  print(paste0("Processed: ", sample_name))
}
}

merged_seurat <- Reduce(function(x,y) merge(x,y), seurat_objects)

save(merged_seurat,file= "merged_all_samples.RData")

load("merged_all_samples.RData")

#Merge only the cd45p and cd45n object 
selected_names <- names(seurat_objects)[grepl("CD45p|CD45n", names(seurat_objects))]
selected_names


