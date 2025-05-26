#############
# Libraries #
#############
library(tidyverse)
library(survival)
library(survminer)
library(colorspace)
library(haven)
library(readxl)

#Combine the forth clinical file 
info <- read_sav("NOROC TMA_patient data (2).sav")

#7 patients missing-unknown reason 
unique(info$Hospital)
symbol <- c("Bergen" = "B", "Oslo Universitetssykehus, Oslo" = "O", "Troms??" = "T")
info$Hospital <- as.factor(info$Hospital)

info$tma_id <- symbol[info$Hospital]
info <- info %>% relocate(tma_id,.after = Hospital)
info$id_tma <- paste(info$tma_id,"-",info$Hospital_ID_number)
info$id_tma <- gsub(" ", "", info$id_tma)
info <- info %>% relocate(id_tma,.after = Hospital)
info <- info[,-c(3,4)]
info <- info[,c(2:4)]

info <- haven::as_factor(info)

#Add scoring information 
info2 <- read.csv2("t32_scores_level_noroc.csv")
info2 <- info2[,-c(2:4)]

difference <- setdiff(info2$ID,info$id_tma)
# Filter df1 based on the difference
result <- info[info$id_tma %in% difference, ]

#2 patients missing -B12 and O157- too few cells
info3 <- read_xlsx("../myc_scoring/myc_scoring.xlsx", sheet = 8)
info3 <- dplyr::select(info3, c(id_sample,myc_value))

# 3 patients missing- B-12, T-4, T-45- too few cells 
info_comb <- merge(info, info2, by.x = "id_tma", by.y = "ID") #131 patients B-12 was in the first info file
info_comb <- merge(info_comb, info3, by.x ="id_tma", by.y = "id_sample") #131 patients - B-12, t4, T45 was already gone from the first file

#remove patients: pallitative,in situ, missing pT, retracted
remove_pats <- c("T-8", "T-28","T-61", "T-63", 
                 "O-10","O-15", "O-25", "O-67",
                 "O-188") #9

info_comb <- info_comb[!(info_comb$id_tma %in% remove_pats),] #131-123, 8 removed
#T-28 were already removed, the rest was removed

#Remove B-99 due to missing data 
info_comb <- info_comb[-20,] #122 patients

#unique(info_comb$DAAR)
#info_comb$DAAR <- as.factor(info_comb$DAAR)
# Create a mapping vector
mapping <- c("Dead of other Causes 5 years" = 0, "Dead of other Cancer 5 years" = 0, 
             "Alive 5 years" = 0, "Dead of Disease 5 years" = 1)

info_comb <- info_comb %>%
  mutate(dss = recode(info_comb$DAAR, !!!mapping))

info_comb$Follow_up_months <- pmin(info_comb$Follow_up_months, 61)

#optional
replacement_vector <- c("high" = "High/Moderate", "moderate" = "High/Moderate",
                        "low" = "Low/Negative", "negative" ="Low/Negative")
info_comb$t32_level <- replacement_vector[info_comb$t32_level]
str(info_comb)
info_comb$t32_level<- factor(info_comb$t32_level,
                              levels = c("high", "moderate", "low", "negative"),
                              labels = c("High/Moderate", "High/Moderate", "Low/Negative", "Low/Negative"))


#remove DSS_months at 0 value 
DSS <- subset(info_comb, info_comb$Follow_up_months > 0) #119

#3 patients removed- O234, O-65, O-79

#Convert DSS/PFI.time into years 
DSS$years <- DSS$Follow_up_months/12

#Remove NA form the df
DSS <- DSS[complete.cases(DSS),]

hmt32 <- subset(DSS, DSS$t32_level == "High/Moderate")
lnt32 <- subset(DSS, DSS$t32_level == "Low/Negative")

table(hmt32$myc_value)
#Define survival
survival = Surv(time= DSS$years, event = DSS$dss)

survival_fit<- survfit(formula = survival ~ DSS$t32_level + DSS$myc_value, data = DSS)

ggsurvplot(fit = survival_fit)

p <- ggsurvplot(fit= survival_fit, 
                  pval = TRUE, 
                  surv.median.line = "hv", 
                  xlab = "DSS (Years)", 
                  ylab = "DSS Probability",
                  ylim=c(0.25,1), 
                  xlim=c(0,5.2),
                  #palette = c("#E69F00","#0072B2","#009E73","#CC79A7","#F64D41","#6D0026"),
                  pval.coord=c(0.1,0.27), 
                  break.x.by= 1,         
                  conf.int = T, 
                  conf.int.alpha = c(0.3), 
                  conf.int.style="step",
                  risk.table = T,
                  risk.table.height = 0.15,
                  ncensor.plot = TRUE,
                  ncensor.plot.height = 0.25,
                  #legend = c(0.5, 0.95),
                  legend.title= "All Patients", 
                tables.theme = theme(axis.text.x = element_blank(), axis.ticks.x=element_blank(), 
                                     axis.title.x = element_blank(), axis.title.y = element_blank(),
                                     axis.text.y = element_text(size = 10)),
                font.x=10, font.y=11, font.tickslab=10)

p


setwd("C:/Users/abe186/UiT Office 365/O365-PhD Anne - General/NOROC Project- New/trim32_scoring/figures/")

pdf("2025_03_19_km_myc_score_noroc.pdf", height = 9, width = 9, onefile = F)
print(p)
dev.off()

png("2025_03_18_km_trim32_score_noroc_2_groups.png", res=200, height = 2000, width = 1800)
print(p)
dev.off()

