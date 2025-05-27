#############
# Libraries #
##############

library(tibble)
library(ggpubr)
library(data.table)
library(dplyr)

###################
# Expression data #
###################

expr <- as.data.frame(fread("TMM_TCGA_HNSC_counts_log2.csv"))
expr_1 <- expr 

normal_id <- read.csv2("Normal_Patient_ID.csv", as.is = T, check.names = F)
normal <- subset(expr_1, expr_1$V1 %in% normal_id$id)
normal <- select(normal, c(V1, ENSG00000119401))


expr_1 <- expr_1[-grep('-11A', expr_1$V1),, drop = FALSE]
expr_1 <- select(expr_1, c(V1, ENSG00000119401))

#Change the patient id before merging
normal$sample_type <- "Normal"
expr_1$sample_type <- "Tumor"

#Merge by rows 
merge <- rbind(expr_1, normal)

dup <- as.data.frame(duplicated(merge$V1)) #no duplicates

my_comparisons <- list( c("Normal", "Tumor") ) 
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                    symbols = c("****", "***", "**", "*", "ns"))

################
# Make boxplot #
################
p <- ggplot(merge, aes(x = sample_type, y = ENSG00000119401, fill = sample_type))+
  geom_point(alpha=0.5,position = position_jitter(width = 0.3, height = 0.5), shape= 21, size= 3)+
  geom_boxplot(fill = "white", alpha = 0.8, outlier.shape = NA) +
  labs(y = expression(paste(italic("TRIM32")~"(log2+1)")), 
       x = "Sample Type", 
       title = expression(paste(italic("TRIM32"), " Expression in Normal vs Tumor Tissue")))+ 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face ="bold"), 
        axis.text.x = element_text(size=10), axis.ticks.x=element_blank(), 
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 11, face ="italic"),
        axis.text.y = element_text(size = 10), 
        panel.background = element_rect(fill = "white",
                                        colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid",
                                 colour = "black")) +
  scale_fill_manual(values = c("green", "red"))+ 
  geom_signif(comparisons = my_comparisons,map_signif_level = T, textsize=5)
p

png("../expression/TRIM32_expr_normal_vs_tumor.png", res = 200 ,height = 1500, width = 1500)
print(p)
dev.off()

pdf("../expression/TRIM32_expr_normal_vs_tumor.pdf", height = 500/72, width = 500/72)
print(p)
dev.off()
