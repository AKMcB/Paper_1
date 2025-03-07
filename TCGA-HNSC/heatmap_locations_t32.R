
#Import the data
#The HNSC data is TMM normalized
#There are no gene or patient duplicates in the expression data
expr <- read.csv2(file = 'TMM_TCGA_HNSC_counts_NoNormal_log2_filtered_hgnc.csv', sep = ';', header = TRUE, check.names = F)
expr <- expr[,-c(1,3)]

dup <-expr[duplicated(expr$hgnc_symbol)|duplicated(expr$hgnc_symbol, fromLast=TRUE),]
expr <- distinct(expr, hgnc_symbol, .keep_all = TRUE)
expr <- expr[,-1]

rownames(expr) <- expr[,1]
expr <- expr[,-1]
#expr <- as.data.frame(t(expr))
#Z-normalize
expr <- 2^(expr)  #linearize log2 data 
#expr <- (scale(expr))
expr <- t(scale(t(expr))) #patients in columns

expr <- as.data.frame(expr)
write.csv2(expr, "TMM_TCGA_HNSC_counts_NoNorma_scaled_filtered_hgnc.csv")
#----------------------------------------------------------------------------------
#Import the libraries
library(ComplexHeatmap)
library(RColorBrewer)
library(gplots)
library(circlize)
library(dendextend)
library(dplyr)
library(ggplot2)
library(tibble)
library(fgsea)

#Extract the genes and TRIM proteins first otherwise the file becomes too big 
# z-normalize the TCGA before making the heatmap 
#Prepare the expr file 
#The expression only contain ensembl ID, so gene symbol must be added before use. 
setwd("C:/Users/abe186/UiT Office 365/O365-Phd Anne - General/TCGA-Head and neck/TMM conversion")
expr <- read.csv2("TMM_TCGA_HNSC_counts_NoNorma_scaled_filtered_hgnc.csv", sep = ";", as.is = T,check.names = F)
expr <- expr[complete.cases(expr),] #Remove NAs
rownames(expr) <- expr[,1]
expr <- expr[,-1]
expr <- as.data.frame(t(expr))
expr <- tibble::rownames_to_column(expr, "id")

expr <- expr[!grepl("-06A",expr$id),]

expr$id <- substr(expr$id, 1, nchar(expr$id) - 4)
expr <- distinct(expr,id,.keep_all = T)

rownames(expr) <- expr[,1]
expr <- expr[,-1]
expr <- as.data.frame(t(expr))

expr <- tibble::rownames_to_column(expr, "gene")
trim <- subset(expr, expr$gene %in% c("TRIM32", "MYC" ))
rownames(trim) <- trim[,1]
trim$gene <- NULL
trim <- as.data.frame(t(trim))
#colnames(trim) <- c("TRIM32", "MYC")

expr <- expr[expr$gene != "TRIM32", ]
expr <- expr[expr$gene != "MYC", ]
rownames(expr) <- expr[,1]
expr <- expr[,-1]

#####################
setwd("C:/Users/abe186/UiT Office 365/O365-Phd Anne - General/TCGA-Head and neck/raw data/")
ann <- read.csv2("tcga_hnsc_updated_locations_03112023.csv")
ann <- ann[,-1]
str(ann)

ann2 <- read.csv2("Clinical_info_HNSC.csv")
ann2 <- dplyr::select(ann2, c(Patient.ID, Subtype))

ann <- merge(ann, ann2, by.x = "case_submitter_id", by.y = "Patient.ID")

expr <- as.data.frame(t(expr))
expr <- rownames_to_column(expr, "id")

expr <- expr[expr$id %in% ann$case_submitter_id, ]
ann <- ann[ann$case_submitter_id %in% expr$id,]

rownames(expr) <- expr$id
expr$id <- NULL
expr <- as.data.frame(t(expr))

#MAke sure the patient order is the same in expr file and ann file
rownames(ann) <- ann$case_submitter_id
ann$case_submitter_id <- NULL
ann<- ann[colnames(expr),]
all(rownames(ann) == colnames(expr))#Check if they are in same order. 


trim<- trim[colnames(expr),]
all(rownames(trim) == colnames(expr))


#Filter gene list 

genes_1 <- as.data.frame(gmtPathways("../gene_lists/HALLMARK_MYC_TARGETS_V1.v2023.2.Hs (1).gmt"))
genes_2 <- as.data.frame(gmtPathways("../gene_lists/HALLMARK_MYC_TARGETS_V2.v2023.2.Hs (1).gmt"))

colnames(genes_1) <- "gene"
colnames(genes_2) <- "gene"

comb <- rbind(genes_1, genes_2)
dup <-comb[duplicated(comb$gene),]
comb <- dplyr::distinct(comb, .keep_all = T)

expr <- subset(expr, rownames(expr) %in% comb$gene)
expr <- as.matrix(expr)

#Complexheatmap needs a matrix
#Values should have "." not ","

#Patients has to be in rownames for cluster validation 
library(clValid)
internal <- clValid::clValid(t(expr), method = "complete", metric = "correlation", clMethods = "hierarchical", nClust = 2:10, validation = "internal")
plot(internal, legend = FALSE)

h <- as.matrix(expr) #Create a matrix of the main expression file 

#Check skewness of the data
#Should be a sharp peak if the data is normalized
dat <- data.frame(values = as.numeric(h))

ggplot(dat, aes(values)) + geom_density(bw = "sj")


##We will plot based on quantiles of the expression values
#First find the breaks- 10 breaks from the lowest to the highest expression 
h_breaks <- seq(min(h), max(h), length.out = 10)
h_breaks

trim$TRIM32 <- as.numeric(trim$TRIM32)#define TRIM32 column as numeric 
z <- trim$TRIM32
z_breaks <- seq(min(z), max(z), length.out = 10)
z_breaks

trim$MYC <- as.numeric(trim$MYC)#define TRIM32 column as numeric 
x <- trim$MYC
x_breaks <- seq(min(x), max(x), length.out = 10)
x_breaks
#Reposition the breaks in quantile positions
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

h_breaks <- quantile_breaks(h, n = 11)
h_breaks

z_breaks <- quantile_breaks(z, n = 11)
z_breaks

x_breaks <- quantile_breaks(x, n = 11)
x_breaks
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

col_fun_2 <- colorRamp2(quantile_breaks(h, n = 11), hcl.colors(11, palette = "Zissou 1", alpha = NULL, rev = FALSE, fixup = TRUE))

library(circlize)
a = colorRamp2(c(0, 5, 10), c("blue", "white", "red"))
ha = HeatmapAnnotation(foo = 1:10, col = list(foo = col_fun))


#Get the annotation files. Make sure this file is in the same order as the expression matrix file.
#For example, order of genes/patients in expr == order of genes/patients in annotation file.
#Otherwise the annotation will show misleading information.
table(ann$site_of_resection_or_biopsy)
table(ann$Subtype)

ha <- HeatmapAnnotation("Location"= ann$site_of_resection_or_biopsy,
                        "Subtype" = ann$Subtype,
                        "TRIM32" = trim$TRIM32,
                        "MYC" = trim$MYC,
                        col = list("Location" = c(
                          "Tongue, NOS" = "#1F78C8",
                          "Overlapping lesion of lip, oral cavity and pharynx" = "#ff0000",
                          "Gum, NOS" = "#33a02c",
                          "Pharynx, NOS" = "#6A33C2",
                          "Larynx, NOS" = "#ff7f00",
                          "Floor of mouth, NOS" = "#FB6496",
                          "Mouth, NOS" = "#b2df8a",
                          "Cheek mucosa" = "#C814FA"),
                          "Subtype" = c("HNSC_HPV+" = "red", 
                                        "HNSC_HPV-"= "blue"),
                          "TRIM32" = col_fun_2, 
                          "MYC" = col_fun_2),
                        simple_anno_size = unit(0.30, "cm"),
                        annotation_name_side = "left",
                        border = T,
                        annotation_name_gp= gpar(fontsize = 12,fontface="bold"),
                        show_legend = T,
                        annotation_legend_param = list(title_gp = gpar(fontsize=12, fontface="bold"), labels_gp=gpar(fontsize=12)))


hr <- rowAnnotation("Up/Down"=annotation2$`Up/Down regulated in cell lines`,
                    "Gene Cluster - Cell Lines" = annotation2$`Gene cluster cell lines`,
                    col=list("Gene Cluster - Cell Lines"=c("Cluster A"="cadetblue2",
                                                           "Cluster B"="burlywood2",
                                                           "Cluster C"="chocolate2",
                                                           "#N/A" = "#CCCCCC"),
                             "Up/Down"=c("Down" = "#1e90ff",
                                         "Up" = "#ff0000")),
                    simple_anno_size = unit(0.30, "cm"),
                    border = T,
                    annotation_name_gp= gpar(fontsize = 8.5,fontface="bold"),
                    show_legend = T,
                    annotation_legend_param = list(title_gp = gpar(fontsize=10, fontface="bold"), labels_gp=gpar(fontsize=8)))

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
              column_gap = unit(c(0.7,0.7, 0.7), "mm"),
              row_gap = unit(c(0.7,0.7), "mm"))

#Print the heatmap as a pdf in your local drive. Use padding to make it look better (in my opinion)
draw(ht, merge_legend=TRUE, padding = unit(c(2, 2, 2, 2), "mm"))

#You have to run dev.off() before checking the pdf in your local drive
dev.off()

#To find the gene and patient orders, draw the heatmap object first so that it will not change with every run
ht <- draw(ht)

setwd("C:/Users/abe186/UiT Office 365/O365-Phd Anne - General/TCGA-Head and neck/Heatmap")


library(ggpubr)
ggexport(ht, filename = "myc_trim32_tcga_location_subtype_3_cluster_2025_02_27.pdf", height = 8, width = 12)
ggexport(ht, filename = "myc_trim32_tcga_location_subtype_3_cluster_2025_02_27.png", res=200, height =1500, width = 2200)

#For getting the column orders
for (i in 1:length(column_order(ht)))   if (i == 1) {
  clu <- t(t(colnames(expr[,column_order(ht)[[i]]])))
  out <- cbind(clu, paste("cluster", i, sep=""))
  colnames(out) <- c("Patient_ID", "Cluster")   } else {
    clu <- t(t(colnames(expr[,column_order(ht)[[i]]])))
    clu <- cbind(clu, paste("cluster", i, sep=""))
    out <- rbind(out, clu)   } 

out

write.csv2(out, "cluster_column_order_Heatmap__myc_trim32_tcga_location_subtype_3_cluster_2025_02_27.csv")


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

write.csv2(out, file = "cluster_rows_order_Heatmap_myc_trim32_tcga_location_subtype_3_cluster_2025_02_27.csv")


#You can also create an interactive heatmap
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("InteractiveComplexHeatmap")
library(InteractiveComplexHeatmap)

htShiny(ht)
