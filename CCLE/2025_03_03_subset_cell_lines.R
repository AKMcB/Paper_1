#############
# Libraries #
#############
library(tidyverse)
library(data.table)
library(ggrepel)

########################
# Read expression data #
########################

expr <- fread("Raw Data/OmicsExpressionProteinCodingGenesTPMLogp1.csv")
colnames(expr) <- trimws(gsub("\\(\\d+\\)", "", colnames(expr)))
expr <- expr[,c("V1", "TRIM32", "MYC")]

info <- read.csv2("Raw Data/Model_edited.csv")
unique(info$OncotreeLineage)
info <- subset(info, info$OncotreeLineage == "Head and Neck")

merged <- merge(expr, info, by.x = "V1", by.y = "ModelID")

#fwrite(merged, "2025_03_03_trim32_myc_expr_in_hsnc_cell_lines.csv", row.names = TRUE)

#Create boxplots 
t <- ggplot(merged, aes(x = OncotreeLineage, y = TRIM32)) +
  geom_point(aes(fill = CellLineName), shape = 21, size = 6) +
  geom_boxplot(alpha = 0.5) +  # Adjust alpha for better visibility
  geom_text_repel(data = subset(merged, CellLineName == "HSC-3"),
                  aes(label = CellLineName),
                  size = 5, hjust = -1)+ 
  labs(title = "Expression in HNSCC Cell Lines", y = " TRIM32 Expression (log2+1)")+ 
  theme_bw()


pdf("2025_03_07_trim32_expr_hnsc.pdf", heigh = 10, width =10)
print(t)
dev.off()
