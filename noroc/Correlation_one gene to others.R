#############
# Libraries #
#############

library(tidyverse)

########################
# Read expression file #
########################

expr <- read.csv2("Raw data/NOROC_TMT_norm_outlier_removed_patient_removed.csv", sep = ";", as.is = T,check.names = F)
expr <- expr[,-c(2,3)]
expr <- distinct(expr,expr$`Gene names`,.keep_all = T) #to remove duplicates 
expr$`Gene names` <- gsub("\\;.*","",expr$`Gene names`)
#be aware that the distinct function creates a column at the end of the df 
#remove it before continuing
expr$`expr$\`Gene names\`` <- NULL 


rownames(expr) <- expr[,1]
expr$`Gene names` <- NULL
expr <- as.data.frame(t(expr))

#Find the correalation between trims
T32 <- t(cor(expr$TRIM32,expr))
T32 <- as.data.frame(T32)

#Change the column names
colnames(T32) <- "TRIM32"

#Write it as csv2
write.csv2(T32, "2025_03_05_trim32_corr_all_proteins_noroc.csv")


