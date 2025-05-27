#############
# Libraries #
#############
library(edgeR)
library(limma)
library(tidyverse)
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

trim32 <- setorder(trim32,cols = TRIM32_level) #all high patients first

rownames(trim32) <- NULL
trim32 <- column_to_rownames(trim32, "ID")
expr<- expr[,rownames(trim32),]
all(rownames(trim32) == colnames(expr))


level <- as.factor(trim32$TRIM32_level)

#####################
# Linear regression #
#####################

design <- model.matrix(~0+level)
colnames(design) <- gsub("level","", colnames(design))
design

contr.matrix <- makeContrasts(
  res = High - Low,
  levels = colnames(design))
contr.matrix

myShapes <- c(0,1)
myColors <- c("red", "blue")

col_group <- myColors[level]
shape_group <- myShapes[level]

plotMDS(expr, col=col_group, pch=shape_group,dim.plot = c(1,2))

legend("topright", col=c("red","blue"),pch = c(0,1), legend=c("High TRIM32","Low TRIM32"))

vfit <- lmFit(expr, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")
summary(decideTests(efit))

top.table <- topTable(efit, sort.by = "P", n = Inf, coef = 1, adjust.method = "BH")

fwrite(top.table, "2025_03_07_linear_reg_noroc_trim32_high_vs_low.csv", row.names = TRUE)

top.table_filt <- subset(top.table, top.table$adj.P.Val < 0.05)

v <- EnhancedVolcano(top.table,
                     lab = rownames(top.table),
                     labSize = 6.0,
                     labCol = 'black',
                     labFace = 'bold',
                     title = "NOROC: High vs Low TRIM32",
                     x = 'logFC',
                     y = 'adj.P.Val',
                     pCutoff = 0.05,
                     colAlpha = 0.5, 
                     ylim = c(0,4.5),
                     FCcutoff = 0.5,
                     legendLabels = c('Not sig.', 'Log (base 2) FC', 'Adjusted p-value',
                                      'Adjusted p-value & Log (base 2) FC'),
                     caption = paste(
                       "Total significant genes (Padj): ",
                       sum(top.table$adj.P.Val < 0.05), ", Up: ",
                       sum(top.table$adj.P.Vall < 0.05 & top.table$logFC > 0), " Down: ",
                       sum(top.table$adj.P.Val < 0.05 & top.table$logFC < -0)))
v

pdf("2025_03_07_volcanoplot_noroc_trim32_high_vs_low.pdf", heigh = 10, width =10)
print(v)
dev.off()

png("2025_03_04_volcanoplot_noroc_trim32_high_vs_low.png",res= 200, heigh = 2000, width =1800)
print(v)
dev.off()

########
# GSEA #
########

indegree_rank <- setNames(object=top.table[,"t"], rownames(top.table))

top.table <- fread("2025_03_07_linear_reg_noroc_trim32_high_vs_low.csv")


#top.table$V1 <- rownames(top.table)

indegree_rank <- top.table %>%
  mutate(signed_rank_stats = sign(logFC) * -log10(P.Value)) %>%
  arrange(desc(signed_rank_stats))

#rownames(indegree_rank) <- indegree_rank$V1

gene_list<- indegree_rank$signed_rank_stats
names(gene_list)<- indegree_rank$V1

## Run fgsea
pathways <- gmtPathways("h.all.v2023.2.Hs.symbols.gmt")
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
signnamelength <- 8 # set to remove prefix from signature names (2 for "GO", 4 for "KEGG", 8 for "REACTOME")
asp <- 3 # aspect ratio of bubble plot
charcut <- 100 # cut signature name in heatmap to this nr of characters

# Make signature names more readable
a <- as.character(dat$pathway) 
for (j in 1:length(a)){
  a[j] <- substr(a[j], signnamelength+2, nchar(a[j]))
}
for (j in 1:length(a)){
  if(nchar(a[j])>charcut) { a[j] <- paste(substr(a[j], 1, charcut), "...", sep=" ")}
}
a <- gsub("_", " ", a)
dat$NAME <- a

#Determine top 10 positive and negative NES signatures
top_pos <- dat[dat$NES > 0, ]
top_pos <- top_pos[order(-top_pos$NES), ][1:10, ]  # Order by NES descending and take top 10
top_neg <- dat[dat$NES < 0, ]
top_neg <- top_neg[order(top_neg$NES), ][1:10, ]  # Order by NES ascending and take top 10
#Combine the top positive and negative NES signatures
dat2 <- rbind(top_pos, top_neg)
dat2$signature <- factor(dat2$NAME, levels = rev(dat2$NAME))

#Determine what signatures to plot (based on FDR cut)
dat2 <- dat2[dat2[,"padj"]<fdrcut,]
dat2 <- dat2[order(dat2[,"padj"]),] 
dat2$signature <- factor(dat2$NAME, (as.character(dat2$NAME)))
#Determine what labels to color
sign_neg <- which(dat2[,"NES"]<0)
sign_pos <- which(dat2[,"NES"]>0)


signcol <- rep(NA, length(dat2$signature))
signcol[sign_neg] <- "blue" # text color of negative signatures
signcol[sign_pos] <- "red" # text color of positive signatures
signcol <- rev(signcol) # need to revert vector of colors, 
#because ggplot starts plotting these from below


#bubble plot
g<-ggplot(dat2, aes(x=padj,y=signature,size=size))+
  geom_point(aes(fill=NES), shape=21)+
  theme_bw()+ 
  xlim(0, fdrcut)+
  scale_size_area(max_size=10,guide="legend")+
  scale_fill_gradient2(low=dencol_neg, high=dencol_pos) + 
  scale_color_identity(labels = c("blue" = "sign_neg", "red" = "sign_pos"), guide = "legend")+
  ggtitle("NOROC: High vs Low TRIM32")+ 
  theme(axis.text.y = element_text(colour=signcol),
        plot.title = element_text(hjust = 0.5, face ="bold"))+
  theme(aspect.ratio=asp, axis.title.y=element_blank()) # test aspect.ratio

g


pdf("2025_03_07_gsea_noroc_trim32_high_vs_low.pdf", height = 10, width = 8)
print(g)
dev.off()
