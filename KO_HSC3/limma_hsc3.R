#############
# Libraries #
#############

library(edgeR)
library(data.table)
library(tidyverse)
library(ggpubr)
library(Glimma)
library(EnhancedVolcano)
library(fgsea)

###################
# Expression file #
###################

expr <- as.data.frame(fread("hsc_wt_ko_counts.csv"))
expr <- expr[,-c(2:4)]
expr <- distinct(expr,expr$Name,.keep_all = TRUE)
expr$`expr$Name` <- NULL
expr <- column_to_rownames(expr, "Name")

colnames(expr) <- c("HSC3_KO_1", "HSC3_KO_3", "HSC3_KO_6", "HSC3_WT_1", 
                    "HSC3_WT_2", "HSC3_WT_3") 

expr <- subset(expr, select = c(HSC3_KO_3,HSC3_KO_6,HSC3_WT_1,
                                HSC3_WT_3))

############################
# Create group information #
############################

meta <- as.data.frame(fread("metadata_hsc3.csv"))
meta <- meta[,-1]

meta <- subset(meta, meta$name %in% c("HSC3_KO_3","HSC3_KO_6","HSC3_WT_1",
                                      "HSC3_WT_3"))
#get same patients in expr as the clinical file
filtered_colnames <- intersect(names(expr), meta$name)
expr <- expr[, filtered_colnames, drop = FALSE]

#Make sure the group info is in the same order as the expr file
rownames(meta) <- NULL
meta <- column_to_rownames(meta, "name")
expr<- expr[,rownames(meta),]
all(rownames(meta) == colnames(expr))

treatment <- as.factor(meta$treatment)

#########
# LIMMA #
#########

#Create DGEList object
d0 <- DGEList(expr)
d0

#Filtering
keep.exprs <- filterByExpr(d0, group= treatment) #By default it filters out genes that has total counts below 10 in the minimum number of samples
d0 <- d0[keep.exprs,, keep.lib.sizes=FALSE]
dim(d0)

#Normalizing the actual data now and visualizing
d0 <- calcNormFactors(d0, method = "TMM")
lcpm <- cpm(d0, log=TRUE)

#Expression values between the conditions 
house <- as.data.frame(lcpm)
house <- as.data.frame(t(house))
house <- rownames_to_column(house, "id")
meta <- rownames_to_column(meta, "id")
house <- merge(house, meta, by = "id")

my_comparisons <- list( c("KO", "WT")) 
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                    symbols = c("****", "***", "**", "*", "ns"))

ggplot(house,aes(x = treatment, y = MYC)) +
  geom_boxplot(aes(fill = treatment), alpha = 0.5, outlier.shape = NA)+
  geom_point() +
  labs(y = "Expression (log2 cpm)", title = "MYC") + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"), 
        axis.text.x = element_text(size = 10), 
        axis.ticks.x = element_blank(), 
        axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 11), 
        axis.text.y = element_text(size = 10), 
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black")) +
  scale_fill_manual(values = c("#333ED4", "red"))+ 
  geom_label(aes(label=id))+
  geom_signif(comparisons = my_comparisons,map_signif_level = T, y_position = c(8.1), textsize=8)

#boxplot(lcpm,col=group)
col_group = treatment
myColors <- c("red", "green")
levels(col_group) <- myColors
col_group <- as.character(col_group)
plotMDS(lcpm, col = col_group) 
legend("topleft", col=c("red","green"),pch = c(0,1), legend=c("KO","WT"))
boxplot(lcpm,col=col_group)

design <- model.matrix(~0+treatment)

colnames(design) <- gsub("treatment","", colnames(design))
design

v <- voom(d0, design, plot=TRUE)
vfit <- lmFit(v, design)

contr <- makeContrasts(KO - WT, levels = colnames(coef(vfit)))
contr

vfit <- contrasts.fit(vfit, contr)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")
summary(decideTests(efit))

top.table <- topTable(efit, sort.by = "P", n = Inf, coef = 1, adjust.method = "BH")
fwrite(top.table, "limma_results/2025_03_05_limma_trim32_WT_1_3_KO_3_6.csv", row.names = TRUE)

test <- subset(top.table, rownames(top.table) %in% "TRIM32")

v <- EnhancedVolcano(top.table,
                       lab = rownames(top.table),
                       labSize = 6.0,
                       labCol = 'black',
                       labFace = 'bold',
                       title = "KO vs WT TRIM32",
                       x = 'logFC',
                       y = 'adj.P.Val',
                       pCutoff = 0.05,
                       colAlpha = 0.5, 
                     ylim = c(0,4.5))
v

png("figures/volcano/volcanoplot_trim32_WT_1_3_KO_3_6_v2.png", res = 200, height = 1800, width = 1500)
print(v)
dev.off()

pdf("figures/volcano/volvanoplot_trim32_WT_1_3_KO_3_6_v2.pdf", height = 8, width = 8)
print(v)
dev.off()

########
# GSEA #
########
#indegree_rank <- setNames(object=top.table[,"t"], rownames(top.table))

top.table <- fread("limma_trim32_WT_1_3_KO_3_6_v2.csv")
#top.table$V1 <- rownames(top.table)

indegree_rank <- top.table %>%
  mutate(signed_rank_stats = sign(logFC) * -log10(P.Value)) %>%
  arrange(desc(signed_rank_stats))

rownames(indegree_rank) <- indegree_rank$V1

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

dat <- data.frame(fgseaRes)
# Settings
fdrcut <- 0.05 # FDR cut-off to use as output for significant signatures
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

# Determine top 10 positive and negative NES signatures
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
  theme_bw()+ # white background, needs to be placed before the "signcol" line
  xlim(0, fdrcut)+
  scale_size_area(max_size=10,guide="legend")+
  scale_fill_gradient2(low=dencol_neg, high=dencol_pos) + 
  scale_color_identity(labels = c("blue" = "sign_neg", "red" = "sign_pos"), guide = "legend")+
  ggtitle("Hallmark: KO vs WT TRIM32")+ 
  theme(axis.text.y = element_text(colour=signcol),
        plot.title = element_text(hjust = 0.5, face ="bold"))+
  theme(aspect.ratio=asp, axis.title.y=element_blank()) # test aspect.ratio

g

png("figures/hallmark_trim32_WT_1_3_KO_3_6_v2.png", res = 200, height = 1800, width = 1500)
print(g)
dev.off()

pdf("figures/hallmark_trim32_WT_1_3_KO_3_6_v2.pdf", height = 10, width = 8)
print(g)
dev.off()
