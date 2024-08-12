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

expr <- as.data.frame(fread("hsc_wt_ko_counts.csv"))
expr <- expr[,-c(2:4)]
expr <- distinct(expr,expr$Name,.keep_all = TRUE)
expr$`expr$Name` <- NULL
expr <- column_to_rownames(expr, "Name")

colnames(expr) <- c("HSC3_KO_1", "HSC3_KO_3", "HSC3_KO_6", "HSC3_WT_1", 
                    "HSC3_WT_2", "HSC3_WT_3") 

############################
# Create group information # 
############################

meta <- as.data.frame(fread("metadata_hsc3.csv"))
meta <- meta[,-1]

#get same patients in expr as the clinical file
filtered_colnames <- intersect(names(expr), meta$name)
expr <- expr[, filtered_colnames, drop = FALSE]


#Make sure the group info is in the same order as the expr file
meta <- column_to_rownames(meta, "name")
meta <- meta %>% arrange(treatment)
expr<- expr[,rownames(meta)]
all(rownames(meta) == colnames(expr))


meta$sample <- rownames(meta)
test_meta <- meta[-3,]
test_meta$lane <- c("3", "4", "4", "4", "4", "1")

treatment <- as.factor(meta$treatment)

test_treat <- as.factor(test_meta$treatment)
lane <- as.factor(test_meta$lane)

#########
# LIMMA #
#########

test <- expr[,-3]

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

house <- as.data.frame(lcpm)
house <- as.data.frame(t(house))
house <- rownames_to_column(house, "id")
meta <- rownames_to_column(meta, "id")
house <- merge(house, meta, by = "id")


ggplot(house,aes(x = treatment, y = ACTB)) +
  geom_boxplot(aes(fill = treatment), alpha = 0.5, outlier.shape = NA)+
  geom_point() +
  labs(y = "Expression (cpm)", title = "ACTB") + 
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
  geom_label(aes(label=id))

#boxplot(lcpm,col=group)
col_group = treatment
myColors <- c("red", "green")
levels(col_group) <- myColors
col_group <- as.character(col_group)
plotMDS(lcpm, col = col_group) 
legend("topright", col=c("red","green"),pch = c(0,1), legend=c("KO","WT"))
boxplot(lcpm,col=col_group)

design <- model.matrix(~0+treatment)



colnames(design) <- gsub("treatment","", colnames(design))
design


v <- voom(d0, design, plot=TRUE)
vfit <- lmFit(v, design)


contr <- makeContrasts(WT - KO, levels = colnames(coef(vfit)))
contr


vfit <- contrasts.fit(vfit, contr)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")
summary(decideTests(efit))

top.table <- topTable(efit, sort.by = "P", n = Inf, coef = 1)
#fwrite(top.table, "limma_trim32_wt_vs_ko.csv")


v <- EnhancedVolcano(top.table,
                     lab = rownames(top.table),
                     labSize = 3.0,
                     pointSize = 2,
                     labCol = 'black',
                     labFace = 'italic',
                     title = "WT vs KO TRIM32",
                     x = 'logFC',
                     y = 'adj.P.Val',
                     pCutoff = 0.05,
                     FCcutoff = 0.5,
                     colAlpha = 0.5)

