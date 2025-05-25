#############
# Libraries #
#############

library(dplyr)
library(Seurat)
library(ggplot2)
library(ComplexHeatmap)
library(SeuratData)
library(RColorBrewer)
set.seed(1234)

###################
# Read RDS object # 
###################

combined <- readRDS("annotated_all_samples_cd45n_batch_corr_epithelial_cells.RDS")

combined_cancer <- combined[, grepl("Cancer", combined@meta.data$pan_cancer_cluster)]

library(edgeR)
y <- Seurat2PB(combined_cancer, sample ="sample_ID")

library(limma)
y <- calcNormFactors(y, method = "TMM")
lcpm <- cpm(y, log=T)

library(ggpubr)
y$samples$Sample_Names <- rownames(y$samples)

p1 <- plotMDS(lcpm)
p1

corrected_lcpm <- removeBatchEffect(lcpm, batch = y$samples$sample)

p2 <- plotMDS(corrected_lcpm)

raw_cor <- cor(lcpm)

pheatmap(raw_cor, main = "without correction")

corrected_cor <- cor(corrected_lcpm)
pheatmap(corrected_cor, main = "with correction")


corrected_lcpm_t32 <- corrected_lcpm[rownames(corrected_lcpm) %in% c("TRIM32", "MYC"),]
corrected_lcpm_t32 <- as.data.frame(t(corrected_lcpm_t32))

p1 <- ggscatter(corrected_lcpm_t32, x ="TRIM32", y="MYC",
                title = "Single Cells: TRIM32 vs MYC Cancer Cells",
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
p1


pdf("2025_03_19_correlation_trim32_myc_single_cells.pdf", height = 5, width = 5)
print(p1)
dev.off()
