# Script for PhD Project 2

## The folders are divided based on what data has been used: 
- **CCLE**: Scripts for the data from the cancer cell line encyclopedia 
- **PANCAN**: Scripts for the data from the Pan Cancer Atlas project 
- **TCGA-HNSC**: Scripts for the data from the The Cancer Genome Atlas program 
- **KO-HSC3**: Scripts for  the data from the KO HSC3 cells 
- **scRNA**: Scripts for the data from the scRNA-Sequencing data 

## CCLE
- **CCLE_boxplots_all_lines.R**
    - This script explore the CCLE data and the expression of a specific gene, and all cell lines are included 

## PANCAN 
- **PanCan_Exp_boxplot.R** 
    - This script explore the PanCan data and the expression of a specific gene, and include all the 33 cancer types. 

## TCGA-HNSC 
- **Survival_script_HNSC_t32.R** 
   - This script explores the different survival parameters such as overall survival, disease-specific survival, and progression free survival in TCGA-HNSC, based on the expression values of a specific gene. 
- **Pearson_t32.R** 
    - This script caluclates the pearson correlation between two genes in the TCGA-HNSC data. 

## KO-HSC3 
- **limma_hsc3.R** 
    - This script performs LIMMA analysis and GSEA of the LIMMA results. Gene lists for the GSEA can be downloaded from [here](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp). 
- **housekeeping_genes_hsc3.R**
    - This scripts creates a heatmap to visualize the expression of common housekeeping genes in the WT/KO HSC3 clones. 
- **correlation_analysis.R** 
    - This script calculates the correlation between the WT/KO HSC3 clones, and creates a heatmap showing the correlation values. 

## scRNA 
- **processing_script.R**
    - This script process all the samples from the scRNA data and performed quality control. 
- **cd45n_p_processing.R**
    - This script subset the dataset either on cd45n or cd45p cells. It performs also annotation of the cells.
- **epithelial_cell_cd45.R**
    - This subset the annotated data for epithelial cells. 
- **metacells_cd45n.R**
    - This script creates metacells of the single cells based on KNN.
