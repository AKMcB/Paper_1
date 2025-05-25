#############
# Libraries #
#############

library(survival)
library(survminer)
library(tidyverse)
library(dplyr)
library(lubridate)
library(grid)
library(gridExtra)
library(reade)

########################
# Read expression data #
########################
expr <- read.csv2("TMM conversion//TMM_TCGA_HNSC_counts_NoNormal_log2_filtered.csv", sep = ";", as.is = T, check.names = F)

#check if there is any gene duplicates 
rownames(expr) <- expr[,1]
expr <- expr[,-1]
expr <- tibble::rownames_to_column(expr, "genes")
dup<-expr[duplicated(expr$genes),] #check if there is any gene duplicates 

#check if there is any patient duplicates 
rownames(expr) <- expr[,1]
expr <- expr[,-1]
expr <- as.data.frame(t(expr))
expr <- tibble::rownames_to_column(expr, "id")
dup<-expr[duplicated(expr$id),]

trims <- expr %>% select(c("id", "ENSG00000119401", "ENSG00000136997"))

trims$id <- gsub("-01A", "",trims$id)
trims$id <- gsub("-01B", "",trims$id)

trims <- trims[!grepl("-06",trims$id),]

dup <- as.data.frame(duplicated(trims$id))

#################
# Clinical info #
#################

info <- read.csv2("Clinical_info_HNSC.csv", sep = ";", as.is = T, check.names = F)
head(info)

info<- info %>% select(c("Patient ID", "Months of disease-specific survival", "Disease-specific Survival status", 
                         "Overall Survival (Months)", "Overall Survival Status", "Progress Free Survival (Months)",
                         "Progression Free Status"))

dss <- info %>% select(c("Patient ID", "Months of disease-specific survival", "Disease-specific Survival status"))

pfi <- info %>% select(c("Patient ID", "Progress Free Survival (Months)","Progression Free Status"))

os <- info %>% select(c("Patient ID", "Overall Survival (Months)", "Overall Survival Status"))

#Merging data
trim_dss <- merge(trims, dss, by.x = "id", by.y = "Patient ID")

trim_pfi <- merge(trims, pfi, by.x = "id", by.y = "Patient ID")

trim_os <- merge(trims, os, by.x = "id", by.y = "Patient ID")


trim_dss$TRIM32_expression <- ifelse(trim_dss$ENSG00000119401 >= median(trim_dss$ENSG00000119401), 'High', "Low")

trim_pfi$TRIM32_expression <- ifelse(trim_pfi$ENSG00000119401 >= median(trim_pfi$ENSG00000119401), 'High', "Low")

trim_os$TRIM32_expression <- ifelse(trim_os$ENSG00000119401 >= median(trim_os$ENSG00000119401), 'High', "Low")


trim_dss$MYC_expression <- ifelse(trim_dss$ENSG00000136997 >= median(trim_dss$ENSG00000136997), 'High', "Low")

trim_pfi$MYC_expression <- ifelse(trim_pfi$ENSG00000136997 >= median(trim_pfi$ENSG00000136997), 'High', "Low")

trim_os$MYC_expression <- ifelse(trim_os$ENSG00000136997 >= median(trim_os$ENSG00000136997), 'High', "Low")


trim_dss <- trim_dss[complete.cases(trim_dss),] #495 -> 470

trim_pfi <- trim_pfi[complete.cases(trim_pfi),] # 495 -> 494

trim_os <- trim_os[complete.cases(trim_os),] # 495 -> 494


trim_dss<- subset(trim_dss, trim_dss$`Months of disease-specific survival` > 0)

trim_pfi<- subset(trim_pfi, trim_pfi$`Progress Free Survival (Months)` > 0)

trim_os<- subset(trim_os, trim_os$`Overall Survival (Months)` > 0)



trim_dss$`Disease-specific Survival status` <- gsub("\\:.*","",trim_dss$`Disease-specific Survival status`)

trim_pfi$`Progression Free Status` <- gsub("\\:.*","",trim_pfi$`Progression Free Status`)

trim_os$`Overall Survival Status` <- gsub("\\:.*","",trim_os$`Overall Survival Status`)



trim_dss$`Months of disease-specific survival` <- as.numeric(trim_dss$`Months of disease-specific survival`)

trim_pfi$`Progress Free Survival (Months)`<- as.numeric(trim_pfi$`Progress Free Survival (Months)`)

trim_os$`Overall Survival (Months)` <- as.numeric(trim_os$`Overall Survival (Months)`)



trim_dss$years<- trim_dss$`Months of disease-specific survival`/12 #Make DSS.months into years

trim_pfi$years<- trim_pfi$`Progress Free Survival (Months)`/12

trim_os$years<- trim_os$`Overall Survival (Months)`/12



trim_dss$`Disease-specific Survival status`<- as.integer(trim_dss$`Disease-specific Survival status`)

trim_pfi$`Progression Free Status` <- as.integer(trim_pfi$`Progression Free Status`)

trim_os$`Overall Survival Status` <- as.integer(trim_os$`Overall Survival Status`)



survival_dss = Surv(time= trim_dss$years, event = trim_dss$`Disease-specific Survival status`)

survival_pfi = Surv(time= trim_pfi$years, event = trim_pfi$`Progression Free Status`)

survival_os = Surv(time= trim_os$years, event = trim_os$`Overall Survival Status`)


survival_fit_dss<- survfit(formula = survival_dss ~ trim_dss$TRIM32_expression, data = trim_dss)

survival_fit_pfi<- survfit(formula = survival_pfi ~ trim_pfi$TRIM32_expression, data = trim_pfi)

survival_fit_os<- survfit(formula = survival_os ~ trim_os$TRIM32_expression, data = trim_pfi)


#######
# DSS #
#######

dss <- ggsurvplot(fit= survival_fit_dss, 
                      pval = TRUE, 
                      surv.median.line = "hv", 
                      xlab = "DSS (Years)", 
                      ylab = "DSS Probability",
                      ylim=c(0.0,1), 
                      xlim=c(0,18),
                      palette = c("#E69F00","#0072B2", "#009E73", "#CC79A7"),
                      pval.coord=c(0.1,0.12), 
                      break.x.by= 5,         
                      conf.int = T, 
                      conf.int.alpha = c(1), 
                      conf.int.style="step",
                      risk.table = T,
                      risk.table.height = 0.15,
                      ncensor.plot = TRUE,
                      ncensor.plot.height = 0.15,
                      #legend = c(0.5, 0.95),
                      legend.labs=c("High","Low"),
                      legend.title= "TRIM32")

dss

pdf("Kaplan-Meier/trim32_hnsc_tcga_dss.pdf", width = 8, height = 7, onefile = F)
print(dss)
dev.off()

#######
# PFI #
#######

pfi <- ggsurvplot(fit= survival_fit_pfi, 
                  pval = TRUE, 
                  surv.median.line = "hv", 
                  xlab = "PFI (Years)", 
                  ylab = "PFI Probability",
                  ylim=c(0.0,1), 
                  xlim=c(0,18),
                  palette = c("#E69F00","#0072B2", "#009E73", "#CC79A7"),
                  pval.coord=c(0.1,0.12), 
                  break.x.by= 5,         
                  conf.int = T, 
                  conf.int.alpha = c(1), 
                  conf.int.style="step",
                  risk.table = T,
                  risk.table.height = 0.15,
                  ncensor.plot = TRUE,
                  ncensor.plot.height = 0.15,
                  #legend = c(0.5, 0.95),
                  legend.labs=c("High","Low"),
                  legend.title= "TRIM32")

pfi

pdf("Kaplan-Meier/trim32_hnsc_tcga_pfi.pdf", width = 8, height = 7, onefile = F)
print(pfi)
dev.off()

######
# OS #
######

os <- ggsurvplot(fit= survival_fit_os, 
                  pval = TRUE, 
                  surv.median.line = "hv", 
                  xlab = "OS (Years)", 
                  ylab = "OS Probability",
                  ylim=c(0.0,1), 
                  xlim=c(0,18),
                  palette = c("#E69F00","#0072B2", "#009E73", "#CC79A7"),
                  pval.coord=c(0.1,0.12), 
                  break.x.by= 5,         
                  conf.int = T, 
                  conf.int.alpha = c(1), 
                  conf.int.style="step",
                  risk.table = T,
                  risk.table.height = 0.15,
                  ncensor.plot = TRUE,
                  ncensor.plot.height = 0.15,
                  #legend = c(0.5, 0.95),
                  legend.labs=c("High","Low"),
                  legend.title= "TRIM32")

os

pdf("Kaplan-Meier/trim32_hnsc_tcga_os.pdf", width = 8, height = 7, onefile = F)
print(os)
dev.off()
                        
