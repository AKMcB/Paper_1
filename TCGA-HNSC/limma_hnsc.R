#############
# Libraries #
#############

library(edgeR)
library(data.table)
library(tidyverse)
library(ggpubr)
library(Glimma)
library(biomaRt)
library(Homo.sapiens)
library(EnhancedVolcano)

###################
# Expression file #
###################

expr <- as.data.frame(fread("TMM conversion/TCGA-HNSC.htseq_counts.tsv"))
expr$Ensembl_ID <- gsub("\\..*","", expr$Ensembl_ID)

rownames(expr) <- expr[,1]
expr <- expr[,-1]

#explore the expr file 
test <- as.data.frame(colnames(expr))
#500 primary cancer samples 
#2 metastatic -06 
#44 normal samples -11 

expr <- as.data.frame(t(expr))
expr <- rownames_to_column(expr,"id")

expr <- expr[!grepl("-06",expr$id),]
expr <- expr[!grepl("-11",expr$id),]

expr$id <- gsub("-01A", "",expr$id)
expr$id <- gsub("-01B", "",expr$id)

dup <- as.data.frame(duplicated(expr$id))

#Xena browser uploads the data after log2(x+1) trasformation. 
#So converting it to linear scale to get the raw count data
expr <- expr[,-c(60485:60489)]

rownames(expr) <- expr[,1]
expr$id <- NULL
expr <- as.data.frame(t(expr))

expr <- 2^expr

############################
# Create group information #
############################

t32_1 <- read.csv2("raw data/TMM_TCGA_HNSC_counts_NoNormal_log2_filtered.csv", sep = ";", as.is = T, check.names = F)
rownames(t32_1) <- t32_1[,1]
t32_1 <- t32_1[,-1]
t32_1 <- as.data.frame(t(t32_1))
t32_1 <- rownames_to_column(t32_1, "id")

t32_1$id <- gsub("-01A", "",t32_1$id)
t32_1$id <- gsub("-01B", "",t32_1$id)

t32_1 <- t32_1[!grepl("-06",t32_1$id),]

dup <- as.data.frame(duplicated(t32_1$id))

rownames(t32_1) <- t32_1[,1]
t32_1 <- t32_1[,-1]
t32_1 <- as.data.frame(t(t32_1))

filtered_colnames <- intersect(names(expr), names(t32_1))
df1_filtered <- expr[, filtered_colnames, drop = FALSE]
expr <- df1_filtered

#Make sure the group info is in the same order as the expr file
t32_1 <- rownames_to_column(t32_1, "id")
t32_2 <- subset(t32_1, t32_1$id == "ENSG00000119401")
rownames(t32_2) <- t32_2[,1] 
t32_2$id <- NULL
t32_2 <- as.data.frame(t(t32_2))

#Make sure the files match the clinical file 
info <- read.csv2("raw data/Clinical_info_HNSC.csv", sep = ";", as.is = T, check.names = F)
info <- info[,c(2,3)]

filtered_colnames <- intersect(names(expr), info$`Patient ID`)
df1_filtered <- expr[, filtered_colnames, drop = FALSE]
expr <- df1_filtered

filtered_colnames <- intersect(rownames(t32_2), info$`Patient ID`)
df1_filtered <- t32_2[filtered_colnames, ,drop = FALSE]
t32_2 <- df1_filtered

t32_2$TRIM32_expression <- ifelse(t32_2$ENSG00000119401 >= median(t32_2$ENSG00000119401), 'High', "Low")

top <- t32_2 %>% filter(t32_2$ENSG00000119401 > quantile(t32_2$ENSG00000119401, 0.75))
bottom <- t32_2 %>% filter(t32_2$ENSG00000119401 < quantile(t32_2$ENSG00000119401, 0.25))

combined <- rbind(top, bottom)
t32_2 <- combined

t32_2 <- t32_2 %>% arrange(TRIM32_expression)
expr<- expr[,rownames(t32_2),]
all(rownames(t32_2) == colnames(expr))

table(t32_2$TRIM32_expression)

group <- as.factor(t32_2$TRIM32_expression)

####################
# gene annotations #
####################
raw_data <- expr
raw_data <- rownames_to_column(raw_data, "gene")
ensembl <- raw_data$gene


mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

G_list <- getBM(mart= mart, filters= "ensembl_gene_id", values= ensembl, attributes= c("ensembl_gene_id",
                                                                                      "hgnc_symbol", "entrezgene_id")
)

G_list$entrezgene_id <- NULL

#order the genes in G-list_1 in the same order as in expr file 
expr <- rownames_to_column(expr, "id")
expr <- merge(expr, G_list, by.x = "id", by.y = "ensembl_gene_id")
expr <- expr %>% relocate("hgnc_symbol", .after = "id")

dup <-expr[duplicated(expr$id)|duplicated(expr$id, fromLast=TRUE),]

expr <- distinct(expr, id, .keep_all = TRUE)

#most of the hgnc symbols come up as duplicates because of missing values
#most of the duplicates have low counts, ~5/30 have counts higher than 10
dup <-expr[duplicated(expr$hgnc_symbol)|duplicated(expr$hgnc_symbol, fromLast=TRUE),]

expr <- distinct(expr, hgnc_symbol, .keep_all = TRUE)

dup <-expr[duplicated(expr$id)|duplicated(expr$id, fromLast=TRUE),]

rm(expr_1, dup, G_list_1, not_empty, raw_data)

expr_1 <- expr[,-1]
expr_1 <- column_to_rownames(expr_1, "hgnc_symbol")

###########
## LIMMA ##
###########

#Create DGEList object
d0 <- DGEList(expr_1)
d0

#Filtering
keep.exprs <- filterByExpr(d0, group= group) #By default it filters out genes that has total counts below 10 in the minimum number of samples
d0 <- d0[keep.exprs,, keep.lib.sizes=FALSE]
dim(d0)

#Normalizing the actual data now and visualizing
d0 <- calcNormFactors(d0, method = "TMM")
lcpm <- cpm(d0, log=TRUE)

myShapes <- c(0,1)
myColors <- c("red", "blue")

col_group <- myColors[group]
shape_group <- myShapes[group]

plotMDS(lcpm, col=col_group, pch=shape_group,dim.plot = c(1,2))

legend("topright", col=c("red","blue"),pch = c(0,1), legend=c("High TRIM32","Low TRIM32"))
 
design <- model.matrix(~0+group)

colnames(design) <- gsub("group","", colnames(design))
design

v <- voom(d0, design, plot=TRUE)
vfit <- lmFit(v, design)

contr <- makeContrasts(High - Low, levels = colnames(coef(vfit)))
contr

vfit <- contrasts.fit(vfit, contr)
efit <- treat(vfit)
plotSA(efit, main="Final model: Mean-variance trend")
summary(decideTests(efit))

top.table <- topTable(efit, sort.by = "P", n = Inf, coef = 1)
fwrite(top.table, "limma/limma_high_vs_low_trim32_2025_02_26.csv", row.names = T)

top.table <- subset(top.table, rownames(top.table) != "TRIM32")

v <- EnhancedVolcano(top.table,
                lab = rownames(top.table),
                labSize = 6.0,
                labCol = 'black',
                labFace = 'bold',
                title = "High TRIM32 vs Low TRIM32",
                x = 'logFC',
                y = 'adj.P.Val',
                colAlpha = 0.5)
v

pdf("limma/for manuscript/volcanoplot_limma_trim32_top_bottom_2025_02_26.pdf", heigh = 10, width =10)
print(v)
dev.off()

png("limma/for manuscript/volcanoplot_limma_trim32_top_bottom_2025_02_26.png",res= 200, heigh = 2000, width =1800)
print(v)
dev.off()

########
# GSEA #
########
library(fgsea)
top.table <- fread("limma/for manuscript/limma_top_bottom_trim32_2025_02_26.csv", sep = ",")

top.table$genes <- rownames(top.table)
str(top.table)
indegree_rank <- top.table %>%
  mutate(signed_rank_stats = sign(logFC) * -log10(P.Value)) %>%
  arrange(desc(signed_rank_stats))

rownames(indegree_rank) <- indegree_rank$genes

gene_list<- indegree_rank$signed_rank_stats
names(gene_list)<- rownames(indegree_rank)

## Run fgsea
pathways <- gmtPathways("gene_lists/h.all.v2023.2.Hs.symbols.gmt")
set.seed(123)
fgseaRes <- fgsea(pathways, gene_list, minSize=15, maxSize=500, nperm=1000)
head(fgseaRes)

# Subset to pathways with FDR < 0.05
sig <- fgseaRes[fgseaRes$padj < 0.05,]
# Top 10 pathways enriched in patients with KO 
sig$pathway[sig$NES > 0][1:10]
# Top 10 pathways enriched in patients with WT
sig$pathway[sig$NES < 0][1:10]

#fwrite(sig, "gsea_gobp_highvslow_t45_subtype_corr_arranged.csv", row.names = TRUE)

dat <- data.frame(fgseaRes)
# Settings
fdrcut <- 0.05 # FDR cut-off to use as output for significant signatures
dencol_neg <- "blue" # bubble plot color for negative ES KO
dencol_pos <- "red" # bubble plot color for positive ES WT
signnamelength <- 3 # set to remove prefix from signature names (2 for "GO", 4 for "KEGG", 8 for "REACTOME")
asp <- 3 # aspect ratio of bubble plot
charcut <- 100 # cut signature name in heatmap to this nr of characters

# Make signature names more readable
a <- as.character(dat$pathway) 
for (j in 1:length(a)){
  a[j] <- substr(a[j], signnamelength+2, nchar(a[j]))
}
for (j in 1:length(a)){
  if(nchar(a[j])>charcut) { a[j] <- paste(substr(a[j], 1, charcut), "...", sep=" ")}}
a <- gsub("_", " ", a)
dat$NAME <- a
#Determine top 10 positive and negative NES signatures
top_pos <- dat[dat$NES > 0, ]
#top_pos <- top_pos[order(-top_pos$NES), ][1:10, ]  # Order by NES descending and take top 10
top_neg <- dat[dat$NES < 0, ]
#top_neg <- top_neg[order(top_neg$NES), ][1:10, ]  # Order by NES ascending and take top 10
# Combine the top positive and negative NES signatures
dat2 <- rbind(top_pos, top_neg)
dat2$signature <- factor(dat2$NAME, levels = rev(dat2$NAME))

# Determine what signatures to plot (based on FDR cut)
dat2 <- dat2[dat2[,"padj"]<fdrcut,]
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
  theme_bw()+ 
  xlim(0, fdrcut)+
  scale_size_area(max_size=10,guide="legend")+
  scale_fill_gradient2(low=dencol_neg, high=dencol_pos) + 
  scale_color_identity(labels = c("blue" = "sign_neg", "red" = "sign_pos"), guide = "legend")+
  ggtitle("Hallmark: High vs Low TRIM32")+ 
  theme(axis.text.y = element_text(colour=signcol),
        plot.title = element_text(hjust = 0.5, face ="bold"))+
  theme(aspect.ratio=asp, axis.title.y=element_blank())

g

pdf("limma/for manuscript/hallmark_limma_trim32_top_bottom_2025_02_26.pdf", heigh = 10, width =10)
print(g)
dev.off()

png("limma/for manuscript/hallmark_limma_trim32_top_bottom_2025_02_26.png",res= 200, heigh = 2000, width =1800)
print(g)
dev.off()
