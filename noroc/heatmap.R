
#############
# Libraries #
#############
library(ComplexHeatmap)
library(RColorBrewer)
library(gplots)
library(circlize)
library(dendextend)
library(tidyverse)
library(clValid)
library(dendsort)
library(fgsea)

#########################
# Read expression file ##
#########################

expr <- read.csv2("Raw data/NOROC_TMT_norm_outlier_removed_patient_removed.csv", sep = ";", as.is = T,check.names = F)
expr <- expr[,-c(2,3)]
expr <- distinct(expr,expr$`Gene names`,.keep_all = T) #to remove duplicates 
expr$`Gene names` <- gsub("\\;.*","",expr$`Gene names`)
#be aware that the distinct function creates a column at the end of the df 
#remove it before continuing
expr$`expr$\`Gene names\`` <- NULL 

####################
# Make annotations #
####################

ann <- subset(expr, expr$`Gene names` %in% c("TRIM32"))
expr <- expr[expr$ != "TRIM32", ]

genes_1 <- as.data.frame(gmtPathways("gene_lists/HALLMARK_MYC_TARGETS_V1.v2023.2.Hs (1).gmt"))
genes_2 <- as.data.frame(gmtPathways("gene_lists/HALLMARK_MYC_TARGETS_V2.v2023.2.Hs (1).gmt"))

colnames(genes_1) <- "gene"
colnames(genes_2) <- "gene"

comb <- rbind(genes_1, genes_2)
dup <-comb[duplicated(comb$gene),]
comb <- dplyr::distinct(comb, .keep_all = T)

expr <- subset(expr, expr$`Gene names` %in% comb$gene)
rownames(expr) <- expr[,1]
expr$`Gene names` <- NULL

rownames(ann) <- ann[,1]
ann$`Gene names` <-NULL
ann <- as.data.frame(t(ann))
ann$d <- "d"

ann<- ann[colnames(expr),]
all(rownames(ann) == colnames(expr))

################################
# Calculate cluster stability #
###############################

#Patients has to be in rownames therefore use t(merged) in the function 
#This will not change the df and the genes will still be in rownames when 
#converting the df into a matrix
internal <- clValid(t(expr), method = "complete", metric = "correlation", clMethods = "hierarchical", nClust = 2:10, validation = "internal")
plot(internal, legend = FALSE)


##################
# Create Heatmap #
##################

#The genes should be in the rownames
h <- as.matrix(expr) #Create a matrix of the main expression file 
str(h)

#Check skewness of the data
#Should be a sharp peak if the data is normalized
dat <- data.frame(values = as.numeric(h))

ggplot(dat, aes(values)) + geom_density(bw = "sj")


##plot based on quantiles of the expression values
#First find the breaks- 10 breaks from the lowest to the highest expression 
h_breaks <- seq(min(h), max(h), length.out = 10)
h_breaks

ann$TRIM32 <- as.numeric(ann$TRIM32)
z <- ann$TRIM32
z_breaks <- seq(min(z), max(z), length.out = 10)
z_breaks

#Reposition the breaks in quantile positions
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

h_breaks <- quantile_breaks(h, n = 11)
h_breaks

z_breaks <- quantile_breaks(z, n = 11)
z_breaks

hist(h)

#Make clusters based on their pearson correlation
row_dist <- as.dist(1-cor(t(h), method = "pearson"))
col_dist <- as.dist(1-cor(h, method = "pearson"))
col_hc <- hclust(col_dist, method = "complete")
row_hc <- hclust(row_dist, method = "complete")

#use the dendsort package and reorder the clustering
Rowv=dendsort(as.dendrogram(row_hc), isRevers=TRUE, type = "average") 
Colv=dendsort(as.dendrogram(col_hc), type = "average")

#color gradient
col_fun <- colorRamp2(quantile_breaks(h, n = 11), c("#1E90FFFF", "#0479F5FF", "#0467D2FF", "#0356AFFF", "#02458CFF", "#023369FF", "#012246FF", "#130202FF", "#6B0002FF", "#B20004FF", "#FF0000FF"))

#0%         10%         20%         30%         40%         50%         60%         70%         80%         90%        100% 
#"#1E90FFFF" "#0479F5FF" "#0467D2FF" "#0356AFFF" "#02458CFF" "#023369FF" "#012246FF" "#130202FF" "#6B0002FF" "#B20004FF" "#FF0000FF" 
a = colorRamp2(c(0, 5, 10), c("blue", "white", "red"))
ha = HeatmapAnnotation(foo = 1:10, col = list(foo = col_fun))


#Get the annotation files. Make sure this file is in the same order as the expression matrix file.
#For example, order of genes/patients in expr == order of genes/patients in annotation file.
#Otherwise the annotation will show misleading information.


ha <- HeatmapAnnotation ("TRIM32"= ann$TRIM32,
                         col = list("TRIM32"= col_fun),
                         simple_anno_size = unit(0.30, "cm"),
                         border = T,
                         annotation_name_gp= gpar(fontsize = 8.5,fontface="bold"),
                         show_legend = T,
                         annotation_legend_param = list(title_gp = gpar(fontsize=10, fontface="bold"), labels_gp=gpar(fontsize=8)))


ht <- Heatmap(h,col = col_fun,
              cluster_columns = Colv,
              name = "Expression Values",
              show_heatmap_legend = T ,
              top_annotation = ha,
              #left_annotation = hr,
              show_column_names = F,
              show_row_names = F,cluster_rows = Rowv,
              row_title_gp = gpar(fontsize=10),
              column_title_gp = gpar(fontsize=10),
              height = unit(12, "cm"),
              width = unit(10.55, "cm"),
              column_dend_reorder = T,
              row_dend_reorder = T,
              show_row_dend = T,
              border = T,
              column_dend_height = unit(2, "cm"),
              column_names_rot = 90,
              #legend_height = unit(4, "cm"),
              row_names_gp = grid::gpar(fontsize= 7.5),
              column_names_gp = grid::gpar(fontsize = 10,fontface="bold"),
              heatmap_legend_param = list(title="Scaled Expression", legend_height=unit(3, "cm"),title_gp=gpar(fontsize=10, fontface="bold"),labels_gp = gpar(fontsize=8)),
              row_split = 2, column_split = 3,
              column_gap = unit(c(0.7,0.7,0.7), "mm"),
              row_gap = unit(c(0.7,0.7), "mm"))

#Use padding to make it look better (in my opinion)
draw(ht, merge_legend=TRUE, padding = unit(c(2, 2, 2, 2), "mm"))

dev.off()

#To find the gene and patient orders, draw the heatmap object first so that it will not change with every run
ht <- draw(ht)

#For getting the column orders
for (i in 1:length(column_order(ht)))   if (i == 1) {
  clu <- t(t(colnames(expr[,column_order(ht)[[i]]])))
  out <- cbind(clu, paste("cluster", i, sep=""))
  colnames(out) <- c("Patient_ID", "Cluster")   } else {
    clu <- t(t(colnames(expr[,column_order(ht)[[i]]])))
    clu <- cbind(clu, paste("cluster", i, sep=""))
    out <- rbind(out, clu)   } 

out

write.csv2(out, "2025_03_04_cluster_column_order_trim32_expr_noroc_myc_target_genes.csv")


#For getting the gene orders
for (i in 1:length(row_order(ht))){   if (i == 1) {
  clu <- t(t(row.names(expr[row_order(ht)[[i]],])))
  out <- cbind(clu, paste("cluster", i, sep=""))
  colnames(out) <- c("GeneID", "Cluster")   } else {
    clu <- t(t(row.names(expr[row_order(ht)[[i]],])))
    clu <- cbind(clu, paste("cluster", i, sep=""))
    out <- rbind(out, clu)   } 
}

out

write.csv2(out, file = "2025_03_04_cluster_rows_trim32_expr_noroc_myc_target_genes.csv")

pdf("2025_03_04_noroc_heatmap_trim32_myc_target_genes.pdf", height = 8, width = 8)
print(ht)
dev.off()

png("2025_03_04_noroc_heatmap_trim32_myc_target_genes.png", res=200, height = 2000, width = 2000)
print(ht)
dev.off() 





