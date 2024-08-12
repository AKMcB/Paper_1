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

########################
# Read expression data #
########################

#Import the data
#The HNSC data is TMM normalized
#There are no gene or patient duplicates in the expression data
expr <- read.csv2("raw data/TMM_TCGA_HNSC_counts_NoNormal_log2_filtered.csv", sep = ';', header = TRUE, check.names = F)

rownames(expr) <- expr[,1]
expr <- expr[,-1]
expr <- as.data.frame(t(expr))
expr <- rownames_to_column(expr, "id")

expr$id <- gsub("-01A", "",expr$id)
expr$id <- gsub("-01B", "",expr$id)

expr <- expr[!grepl("-06",expr$id),]

dup <- as.data.frame(duplicated(expr$id))

#################
# Clinical info #
#################

ann <- read.csv2("raw data/Clinical _info_HNSC_V2_no_duplicates.csv", sep = ";",
                 as.is = T, check.names=F)

#Subset based on location 
ann <- ann[, c(2,3,71,120)]
ann <- ann[-c(528,529),]
#write.csv2(ann, "hnsc_locations_report.csv")

expr <- expr[expr$id %in% ann$case_submitter_id, ]
ann <- ann[ann$case_submitter_id %in% expr$id,]


rownames(expr) <- expr[,1]
expr <- expr[,-1]
expr <- as.data.frame(t(expr))

expr <- tibble::rownames_to_column(expr, "gene")
trim <- subset(expr, expr$gene %in% c("ENSG00000119401", "ENSG00000136997" ))
rownames(trim) <- trim[,1]
trim$gene <- NULL
trim <- as.data.frame(t(trim))
colnames(trim) <- c("TRIM32", "MYC")

expr <- expr[expr$gene != "ENSG00000119401", ]
expr <- expr[expr$gene != "ENSG00000136997", ]
rownames(expr) <- expr[,1]
expr <- expr[,-1]

# Use table() to count occurrences of each unique string in the column
string_counts <- table(ann$site_of_resection_or_biopsy)

# Print the counts
print(string_counts)

patterns_pahrynx <- c("Oropharynx, NOS","Posterior wall of oropharynx",
                      "Hypopharynx, NOS", "Tonsil, NOS", "Base of tongue, NOS")

# Define the patterns you want to replace as a regular expression

pattern_regex <- paste(patterns_pahrynx, collapse = "|")

# Define the replacement string
replacement <- "Pharynx, NOS"


ann$site_of_resection_or_biopsy <- gsub(pattern_regex, replacement, ann$site_of_resection_or_biopsy)
string_counts <- table(ann$site_of_resection_or_biopsy)
print(string_counts)

patterns_tongue <- c("Ventral surface of tongue, NOS","Border of tongue")
pattern_regex <- paste(patterns_tongue, collapse = "|")
replacement <- "Tongue, NOS"

ann$site_of_resection_or_biopsy <- gsub(pattern_regex, replacement, ann$site_of_resection_or_biopsy)
string_counts <- table(ann$site_of_resection_or_biopsy)
print(string_counts)

patterns_gum <- c("Lower gum","Upper gum", "Hard palate", "Palate, NOS")
pattern_regex <- paste(patterns_gum, collapse = "|")
replacement <- "Gum, NOS"


ann$site_of_resection_or_biopsy <- gsub(pattern_regex, replacement, ann$site_of_resection_or_biopsy)
string_counts <- table(ann$site_of_resection_or_biopsy)
print(string_counts)

patterns_floor <- c("Anterior floor of mouth")
pattern_regex <- paste(patterns_floor, collapse = "|")
replacement <- "Floor of mouth, NOS"


ann$site_of_resection_or_biopsy <- gsub(pattern_regex, replacement, ann$site_of_resection_or_biopsy)
string_counts <- table(ann$site_of_resection_or_biopsy)
print(string_counts)

patterns_lar <- c("Supraglottis")
pattern_regex <- paste(patterns_lar, collapse = "|")
replacement <- "Larynx, NOS"

ann$site_of_resection_or_biopsy <- gsub(pattern_regex, replacement, ann$site_of_resection_or_biopsy)
string_counts <- table(ann$site_of_resection_or_biopsy)
print(string_counts)

#Remove lip since that can be classified at melanoma 
#Remove retromolar areas as this site can not be defined further 

ann <- subset(ann,!ann$site_of_resection_or_biopsy %in% c("Lip, NOS", "Retromolar area"))
unique(ann$site_of_resection_or_biopsy)

write.csv2(ann, "tcga_hnsc_updated_locations_03112023.csv")

#-------------------------------------------------------------------------------

ann <- read.csv2("raw data/tcga_hnsc_updated_locations_03112023.csv")
ann <- ann[,-1]
str(ann)

expr <- expr[expr$id %in% ann$case_submitter_id, ]
ann <- ann[ann$case_submitter_id %in% expr$id,]

ann$trim32 <- expr$ENSG00000119401
ann$trim32_expression <- ifelse(expr$ENSG00000119401 >= median(expr$ENSG00000119401), 'High', "Low")

rownames(expr) <- expr[,1]
expr <- expr[,-1]
expr <- as.data.frame(t(expr))

#MAke sure the patient order is the same in expr file and ann file
rownames(ann) <- ann$case_submitter_id
ann$case_submitter_id <- NULL
ann<- ann[colnames(expr),]
all(rownames(ann) == colnames(expr))#Check if they are in same order. 


#Z-normalize
expr <- 2^(expr)  #linearize log2 data 
expr <- t(scale(t(expr)))
expr <- as.data.frame(expr)


###########
# Heatmap #
###########

#Complexheatmap requires a matrix
#The genes should be in the rownames!!

#Patients has to be in rownames
internal <- clValid::clValid(t(expr), method = "complete", metric = "correlation", clMethods = "hierarchical", nClust = 2:10, validation = "internal")
plot(internal, legend = FALSE)

h <- as.matrix(expr) #Create a matrix of the main expression file 

#Check skewness of the data
#Should be a sharp peak if the data is normalized
dat <- data.frame(values = as.numeric(h))

ggplot(dat, aes(values)) + geom_density(bw = "SJ")


##We will plot based on quantiles of the expression values
#First find the breaks- 10 breaks from the lowest to the highest expression 
h_breaks <- seq(min(h), max(h), length.out = 10)
h_breaks

ann$trim32 <- as.numeric(ann$trim32)#define ann32 column as numeric 
z <- ann$trim32
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
#Dendsort package for clustering

library(dendsort)

#Make clusters based on their pearson correlation
row_dist <- as.dist(1-cor(t(h), method = "pearson"))
col_dist <- as.dist(1-cor(h, method = "pearson"))
col_hc <- hclust(col_dist, method = "complete")
row_hc <- hclust(row_dist, method = "complete")


#We can use the dendsort package and reorder the clustering
Rowv=dendsort(as.dendrogram(row_hc), isRevers=TRUE, type = "average") 
Colv=dendsort(as.dendrogram(col_hc), type = "average")

#Choose the color gradient here. You can change the colors by googling R color codes.
#Run these code to add the option 
col_fun <- colorRamp2(quantile_breaks(h, n = 11), c("#1e90ff", "#0479f5", "#0467d2", "#0356af", "#02458c","#023369","#012246", "#130202", "#6b0002","#b20004","#ff0000"))

col_fun_2 <- colorRamp2(quantile_breaks(z_breaks, n = 11), c("#77D1B5","#18BDB0","#00A6AE","#008CB0","#006AA8", "#130202","#69000C", "#8F0B1B", "#B71729" , "#E02136", "#FF7078"))

library(circlize)
a = colorRamp2(c(0, 5, 10), c("blue", "white", "red"))
ha = HeatmapAnnotation(foo = 1:10, col = list(foo = col_fun))


#Get the annotation files. Make sure this file is in the same order as the expression matrix file.
#For example, order of genes/patients in expr == order of genes/patients in annotation file.
#Otherwise the annotation will show misleading information.
unique(ann$site_of_resection_or_biopsy)

ha <- HeatmapAnnotation("Location"= ann$site_of_resection_or_biopsy,
                        "TRIM32" = ann$trim32,
                        col = list("Location" = c(
                          "Tongue, NOS" = "#1F78C8",
                          "Overlapping lesion of lip, oral cavity and pharynx" = "#ff0000",
                          "Gum, NOS" = "#33a02c",
                          "Pharynx, NOS" = "#6A33C2",
                          "Larynx, NOS" = "#ff7f00",
                          "Floor of mouth, NOS" = "#FB6496",
                          "Mouth, NOS" = "#b2df8a",
                          "Cheek mucosa" = "#C814FA"),
                          "TRIM32" = col_fun_2),
                        simple_anno_size = unit(0.30, "cm"),
                        annotation_name_side = "left",
                        border = T,
                        annotation_name_gp= gpar(fontsize = 12,fontface="bold"),
                        show_legend = T,
                        annotation_legend_param = list(title_gp = gpar(fontsize=12, fontface="bold"), labels_gp=gpar(fontsize=12)))


#Time to make the heatmap. Put it in an object (ht) for further analysis -> finding the patient and gene orders
#Play around a bit with the different parameters of the heatmap. If it has more information, change the height and width.

ht <- Heatmap(h,col = col_fun,
              cluster_columns = Colv,
              name = "Expression Values",
              show_heatmap_legend = T ,
              top_annotation = ha,
              #left_annotation = hr,
              show_column_names = F,
              show_row_names = F,cluster_rows = Rowv,
              row_title_gp = gpar(fontsize=12),
              column_title_gp = gpar(fontsize=12),
              height = unit(12, "cm"),
              width = unit(10.55, "cm"),
              column_dend_reorder = T,
              row_dend_reorder = T,
              show_row_dend = T,
              border = T,
              column_dend_height = unit(2, "cm"),
              column_names_rot = 90,
              #legend_height = unit(4, "cm"),
              row_names_gp = grid::gpar(fontsize= 12),
              column_names_gp = grid::gpar(fontsize = 12,fontface="bold"),
              heatmap_legend_param = list(title="Scaled Expression", legend_height=unit(3, "cm"),title_gp=gpar(fontsize=12, fontface="bold"),labels_gp = gpar(fontsize=12)),
              row_split = 2, column_split = 3,
              column_gap = unit(c(0.7,0.7,0.7), "mm"),
              row_gap = unit(c(0.7,0.7), "mm"))

#Print the heatmap as a pdf in your local drive. Use padding to make it look better (in my opinion)
draw(ht, merge_legend=TRUE, padding = unit(c(2, 3, 2, 2), "mm"))

#You have to run dev.off() before checking the pdf in your local drive
dev.off()

#To find the gene and patient orders, draw the heatmap object first so that it will not change with every run
ht <- draw(ht)


ggexport(ht, filename = "heatmap_hnsc_locations_T32_exp_v2.pdf", height = 8, width = 12)
ggexport(ht, filename = "heatmap_hnsc_locations_T32_exp_v2.png", res=200, height =1500, width = 2200)


#For getting the column orders
for (i in 1:length(column_order(ht)))   if (i == 1) {
  clu <- t(t(colnames(expr[,column_order(ht)[[i]]])))
  out <- cbind(clu, paste("cluster", i, sep=""))
  colnames(out) <- c("Patient_ID", "Cluster")   } else {
    clu <- t(t(colnames(expr[,column_order(ht)[[i]]])))
    clu <- cbind(clu, paste("cluster", i, sep=""))
    out <- rbind(out, clu)   } 

  out

write.csv2(out, "cluster_column_order_Heatmap__HNSC_locations_T32_exp_v2.csv")


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

write.csv2(out, file = "cluster_rows_order_Heatmap__HNSC_locations_T32_exp_v2.csv")
