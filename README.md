# Script for Paper I

#### The folders are divided based on what data has been used: 
- **CCLE**: Scripts for the data from the cancer cell line encyclopedia 
- **PANCAN**: Scripts for the data from the Pan Cancer Atlas project 
- **TCGA-HNSC**: Scripts for the data from the The Cancer Genome Atlas program 
- **noroc**: Scripts for the data from the NOROC study
- **KO-HSC3**: Scripts for  the data from the KO HSC3 cells 
- **scRNA**: Scripts for the data from the scRNA-Sequencing data 

#### CCLE 
- **CCLE_boxplots_all_lines.R**
    - This script explore the CCLE data and the expression of a specific gene, and all cell lines are included 
- **2025_03_03_subset_cell_lines.R**
    - This scripts explores the CCLE data in HNSC and the expression of a specific gene.

#### PANCAN 
- **PanCan_Exp_boxplot.R** 
    - This script explore the PanCan data and the expression of a specific gene, and include all the 33 cancer types. 

#### TCGA-HNSC 
- **Survival_script_HNSC_t32.R** 
   - This script explores the different survival parameters such as overall survival, disease-specific survival, and progression free survival in TCGA-HNSC, based on the expression values of a specific gene. 
- **Pearson_t32.R** 
    - This script caluclates the pearson correlation between two genes in the TCGA-HNSC data. 
- **TRIM32_exp_Normal_VS_Tumor_HNSC.R**
    - This scripts explores TRIM32 expression in normal versus tumor tissue. 
- **clinical_table_tcga_hnsc.r**
    -Create a clinical table based on TRIM32 expression and clinica features. 
- **heatmap_locations_t32.R**
    - Creates a heatmap based on the TCGA-HNSC patients and MYC target genes from the Molecular Signatures Database. HNSC location, HPV-status and MYC expression is added as annotations. 
- **heatmap_boxplot_t32.R**
    - Creates boxplots showing the expression of genes in the heatmap clusters created in the script heatmap_locations_t32.R.
- **limma_hnsc.R**
    - LIMMA analysis of TCGA-HNSC based on TRIM32 expression.
- **updated_hnsc_locations.R**
    - Refinement of TCGA-HNSC cancer sites

#### noroc 
- **Correlation_one gene to others.R**
    - Performs correlation analysis between TRIM32 and all proteins in the NOROC data. 
- **heatmap.R**
    - Creates a heatmap showing the expression of proteins in the NOROC study filtered for MYC target genes from the Molecular Signtaures Database. 
- **boxplot_heatmap_clusters.R**
    - Creates boxplots showing the expression of genes in the heatmap clusters creates in the script heatmap.R. 
- **lmfit_noroc.R**
    - Performs a linear regression analysis on the NOROC data based on TRIM32 expression.
- **clinical_table.r**
    - Create a clinical table based on the TMA score of TRIM32
- **Survival_script_NOROC.R**
    - Explores disease-specific survival of NOROC patients based on TRIM32 and MYC TMA score

#### KO-HSC3 
- **limma_hsc3.R** 
    - This script performs LIMMA analysis and GSEA of the LIMMA results. Gene lists for the GSEA can be downloaded from [here](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp). 
- **housekeeping_genes_hsc3.R**
    - This scripts creates a heatmap to visualize the expression of common housekeeping genes in the WT/KO HSC3 clones. 
- **correlation_analysis.R** 
    - This script calculates the correlation between the WT/KO HSC3 clones, and creates a heatmap showing the correlation values. 

#### scRNA 
- **processing_script.R**
    - This script process all the samples from the scRNA data and performed quality control. 
- **cd45n_p_processing.R**
    - This script subset the dataset either on cd45n or cd45p cells. It performs also annotation of the cells.
- **epithelial_cell_cd45.R**
    - This subset the annotated data for epithelial cells. 
- **metacells_cd45n_new.R**
    - This script creates metacells of the single cells based on KNN. 

- **venndiagram.R**
    - Creates a venndiagram comparing genes between TCGA-HNSC, NOROC, and TRIM32 KO model. 
- **volcano_filtered_genes.R**
    - Creates a volcano plot based on the LIMMA results and filter the results based on the MYC target genes from the Molecular Signtaures Database
