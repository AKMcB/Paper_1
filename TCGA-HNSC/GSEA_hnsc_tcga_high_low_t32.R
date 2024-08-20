
##############
##Libraries ##
##############

library(dplyr)
library(tibble)
library(caroline)
library(dplyr)
#####################
## Expression data ##
#####################

expr <- read.csv2("raw data/TMM_TCGA_HNSC_counts_NoNormal_log2_filtered.csv", sep = ";", as.is = T, check.names = F)
rownames(expr) <- expr[,1]
expr <- expr[,-1]
expr <- as.data.frame(t(expr))
expr <- rownames_to_column(expr, "id")

expr$id <- gsub("-01A", "",expr$id)
expr$id <- gsub("-01B", "",expr$id)

expr <- expr[!grepl("-06",expr$id),]

dup <- as.data.frame(duplicated(expr$id))

ann <- read.csv2("raw data/Clinical_info_HNSC.csv", sep = ";",
                 as.is = T, check.names=F)

expr <- expr[expr$id %in% ann$`Patient ID`, ]

#####################
## GSEA dataframes ##
#####################

expr$trim32_expression <- ifelse(expr$ENSG00000119401 >= median(expr$ENSG00000119401), "High", "Low")

#Check if the groups are correct
median(expr$ENSG00000119401)
test <- expr %>% select("ENSG00000119401", "trim32_expression")


high <- subset(expr, expr$trim32_expression == "High")
low <- subset(expr, expr$trim32_expression == "Low")

rownames(high) <- high[,1]
high$id <- NULL
high <- as.data.frame(t(high))

rownames(low) <- low[,1]
low$id <- NULL
low <- as.data.frame(t(low))

colnames(high) <- paste("High",colnames(high), sep= "_")

colnames(low) <- paste("Low",colnames(low), sep= "_") 

#Merge the files 
merged <- cbind(high, low)

#Remove possible NAs
merged<- merged[complete.cases(merged),]

#Create phenotype file
pheno_merged <- merged[11997, ]

merged <- merged[-11997,]


setwd("C:/Users/abe186/UiT Office 365/O365-Bioinformatikk TRIM27 - General/TCGA-Head and neck/GSEA/")

write.csv2(merged, "gsea_hnsc_tcga_high_vs_low_t32_495pas.csv") #Save phenotype file 
write.csv2(pheno_merged, "gsea_hnsc_tcga_high_vs_low_t32_pheno_495pas.csv")

