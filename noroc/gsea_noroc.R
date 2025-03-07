
library(tidyverse)
library(clusterProfiler)
library(msigdbr) 
library(org.Hs.eg.db) # Change based on species (e.g., org.Mm.eg.db for mouse)
library(fgsea)
library(EnhancedVolcano)
#########################
# Read expression file ##
#########################

expr <- read.csv2("NOROC_TMT_norm_outlier_removed_patient_removed.csv", sep = ";", as.is = T,check.names = F)
expr <- expr[,-c(2,3)]
expr <- distinct(expr,expr$`Gene names`,.keep_all = T) #to remove duplicates 
expr$`Gene names` <- gsub("\\;.*","",expr$`Gene names`)
#be aware that the distinct function creates a column at the end of the df 
#remove it before continuing
expr$`expr$\`Gene names\`` <- NULL 

rownames(expr) <- expr[,1]
expr$`Gene names` <- NULL

str(expr)
###################
# make group info #
###################

trim32 <- as.data.frame(t(expr))

trim32 <- dplyr::select(trim32, TRIM32)

hist(trim32$TRIM32)
shapiro.test(trim32$TRIM32)

trim32$TRIM32_level <- ifelse(trim32$TRIM32 >= median(trim32$TRIM32), 'High', "Low")

trim32 <- rownames_to_column(trim32, "ID")
#Define groups 
high <- trim32$ID[trim32$TRIM32_level == "High"]
low <- trim32$ID[trim32$TRIM32_level == "Low"]


##################
# Perform t-test #
##################

# Initialize a vector to store p-values from the differential expression tests
t_test_results <- numeric(nrow(expr))
# Perform tests for each gene
for (i in 1:nrow(expr)) {
  gene <- expr[i, ]
  gene_name <- rownames(gene)
  # Convert to numeric
  gene <- as.numeric(gene)
  
  # Extract data for each cluster and remove NAs
  data_high <- expr[,high]
  data_high <- data_high[gene_name,]
  data_high <- as.numeric(unlist(data_high))
  
  data_low <- expr[,low]
  data_low <- data_low[gene_name,]
  data_low <- as.numeric(unlist(data_low))
  # Check if there are enough observations for the test
  if (length(data_high) < 2 || length(data_low) < 2) {
    t_test_results[i] <- NA  # Not enough data to perform the test
    next
  }
  
  # Perform Shapiro-Wilk test on combined data
  combined_data <- c(data_high, data_low)
  combined_data <- unlist(combined_data)
  if (length(combined_data) >= 3 && length(combined_data) <= 5000) {
    shapiro_result <- shapiro.test(combined_data)
    is_normal <- shapiro_result$p.value > 0.05
  } else {
    not_normal <- shapiro_result$p.value < 0.05  # Assume non-normal if not enough data
  }
  
  # Perform the appropriate test based on normality
  if (is_normal) {
    # Use t-test for normally distributed genes
    test_result <-  t.test(data_high, data_low)
  } else {
    # Use Wilcoxon test for non-normally distributed genes
    test_result <- wilcox.test(data_high, data_low)
  }
  
  # Store the p-value
  t_test_results[i] <- test_result$p.value
}
# View the first few p-values
head(t_test_results)


p_adj <- p.adjust(t_test_results, method = "BH")


fold_changes <- rowMeans(expr[, low]) - rowMeans(expr[, high])


# Assuming expr is a matrix with genes as rows and samples as columns
mean_low <- rowMeans(expr[, low])
mean_high <- rowMeans(expr[, high])
# Fold Change
fold_changes <- mean_high / mean_low
# Log2 Fold Change
log2_fold_change <- log2(fold_change)


results_df <- data.frame(Gene = rownames(expr),
                         FC = fold_changes,
                         pvalue = t_test_results,
                         padj = p_adj)


results_df <- subset(results_df, rownames(results_df) != "TRIM32")

fwrite(results_df, "../GSEA/for_manuscript/gsea_noroc_high_vs_low_trim32_new_fold_change.csv", row.names = T)

results_df <- fread("GSEA/for_manuscript/gsea_noroc_high_vs_low_trim32.csv")
summary(results_df)

v <- EnhancedVolcano(results_df,
                     lab = results_df$Gene,
                     labSize = 6.0,
                     labCol = 'black',
                     labFace = 'bold',
                     title = "NOROC: High vs Low TRIM32",
                     x = 'FC',
                     y = 'padj',
                     xlim = c(-10,10),
                     ylim = c(0,22),
                     FCcutoff = 2,
                     colAlpha = 0.5)
v

pdf("../GSEA/for_manuscript/volcanoplot_noroc_trim32_high_vs_low_new_fold_change_2025_02_27.pdf", heigh = 10, width =10)
print(v)
dev.off()

png("../GSEA/for_manuscript/volcanoplot_noroc_trim32_high_vs_low_new_fold_change_2025_02_27.png",res= 200, heigh = 2000, width =1800)
print(v)
dev.off()


########
# GSEA #
########

#Rank genes based on fold change and significance 
indegree_rank <- results_df %>%
  mutate(signed_rank_stats = sign(FC) * -log10(padj)) %>%
  arrange(desc(signed_rank_stats))

rownames(indegree_rank) <- indegree_rank$Gene

gene_list<- indegree_rank$signed_rank_stats
names(gene_list)<- rownames(indegree_rank)

## Run fgsea
pathways <- gmtPathways("../gene_lists/h.all.v2023.2.Hs.symbols.gmt")
set.seed(123)
fgseaRes <- fgsea(pathways, gene_list, minSize=15, maxSize=500, nperm=1000)
head(fgseaRes)

# Subset to pathways with FDR < 0.05
sig <- fgseaRes[fgseaRes$padj < 1,]
# Top 10 pathways enriched in patients with KO 
sig$pathway[sig$NES > 0][1:10]
# Top 10 pathways enriched in patients with WT
sig$pathway[sig$NES < 0][1:10]

#fwrite(sig, "gsea_gobp_highvslow_t45_subtype_corr_arranged.csv", row.names = TRUE)

dat <- data.frame(fgseaRes)
# Settings
fdrcut <- 0.1 # FDR cut-off to use as output for significant signatures
dencol_neg <- "blue" # bubble plot color for negative ES KO
dencol_pos <- "red" # bubble plot color for positive ES WT
signnamelength <- 8 # set to remove prefix from signature names (2 for "GO", 4 for "KEGG", 8 for "REACTOME")
asp <- 3 # aspect ratio of bubble plot
charcut <- 100 # cut signature name in heatmap to this nr of characters

# Make signature names more readable
a <- as.character(dat$pathway) # 'a' is a great variable name to substitute row names with something more readable
for (j in 1:length(a)){
  a[j] <- substr(a[j], signnamelength+2, nchar(a[j]))
}
#a <- tolower(a) # convert to lower case (you may want to comment this out, it really depends on what signatures you are looking at, c6 signatures contain gene names, and converting those to lower case may be confusing)
for (j in 1:length(a)){
  if(nchar(a[j])>charcut) { a[j] <- paste(substr(a[j], 1, charcut), "...", sep=" ")}
} # cut signature names that have more characters than charcut, and add "..."
a <- gsub("_", " ", a)
dat$NAME <- a
# Determine what signatures to plot (based on FDR cut)
dat2 <- dat[dat[,"padj"]<fdrcut,]
dat2 <- dat2[order(dat2[,"padj"]),] 
dat2$signature <- factor(dat2$NAME, rev(as.character(dat2$NAME)))
# Determine what labels to color
sign_neg <- which(dat2[,"NES"]<0)
sign_pos <- which(dat2[,"NES"]>0)
# Color labels
signcol <- rep(NA, length(dat2$signature))
signcol[sign_neg] <- "blue" # text color of negative signatures
signcol[sign_pos] <- "red" # text color of positive signatures
signcol <- rev(signcol) # need to revert vector of colors, 
#because ggplot starts plotting these from below


# Plot bubble plot
g<-ggplot(dat2, aes(x=padj,y=signature,size=size))+
  geom_point(aes(fill=NES), shape=21)+
  theme_bw()+ # white background, needs to be placed before the "signcol" line
  xlim(0, fdrcut)+
  scale_size_area(max_size=10,guide="legend")+
  scale_fill_gradient2(low=dencol_neg, high=dencol_pos) + 
  scale_color_identity(labels = c("blue" = "sign_neg", "red" = "sign_pos"), guide = "legend")+
  ggtitle("Hallmark: High vs Low TRIM32")+ 
  theme(axis.text.y = element_text(colour=signcol),
        plot.title = element_text(hjust = 0.5, face ="bold"))+
  theme(aspect.ratio=asp, axis.title.y=element_blank()) # test aspect.ratio

g

png("GSEA/for_manuscript/hallmark_noroc_trim32_high_vs_low_2025_02_26.png", res = 200, height = 1800, width = 1500)
print(g)
dev.off()

pdf("GSEA/for_manuscript/hallmark_noroc_trim32_high_vs_low_2025_02_26.pdf", height = 10, width = 8)
print(g)
dev.off()

