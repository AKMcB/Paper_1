#####################
# Loading libraries #
#####################

library(data.table)
library(dplyr)
library(table1)
library(tidyverse)

###########################
# Reading expression file #
###########################

expr <- read.csv2("TMM conversion/TMM_TCGA_HNSC_counts_NoNormal_log2_filtered.csv", sep = ";", as.is = T, check.names = F)

#check if there is any gene duplicates 
rownames(expr) <- expr[,1]
expr <- expr[,-1]
expr <- tibble::rownames_to_column(expr, "genes")
dup<-expr[duplicated(expr$genes),] #check if there is any gene duplicates 

#check if there is any patient duplicates 
rownames(expr) <- expr[,1]
expr <- expr[,-1]
expr <- as.data.frame(t(expr))
expr <- tibble::rownames_to_column(expr, "Patient ID")
dup<-expr[duplicated(expr$id),]

trims <- expr %>% select(c("Patient ID", "ENSG00000119401", "ENSG00000136997"))

trims <- trims[!grepl("-06",trims$`Patient ID`),]
trims$`Patient ID` <- gsub("-01A", "",trims$`Patient ID`)
trims$`Patient ID` <- gsub("-01B", "",trims$`Patient ID`)


dup <- as.data.frame(duplicated(trims$id))

#######################################
# Subset clinical file with expr file #
#######################################

info <- read.csv2("raw data/Clinical_info_HNSC.csv", sep = ";", as.is = T, check.names = F)
head(info)
str(info)

info2 <- read.csv2("raw data/tcga_hnsc_updated_locations_03112023.csv")
info2 <- info2[,c(3,5)]
colnames(info2)[1] <- "Patient ID"

trims <- left_join(trims, info2, by = "Patient ID")


merged <- merge(trims, info, by ="Patient ID")
merged <- merged[,-c(5,6)]
head(merged)

merged$TRIM32_expression <- ifelse(merged$ENSG00000119401 >= median(merged$ENSG00000119401), 'High', "Low")

merged <- merged[, !names(merged) %in% c(
  "American Joint Committee on Cancer Publication Version Type",
  "Cancer Type",
  "TCGA PanCanAtlas Cancer Type Acronym",
  "Cancer Type Detailed",
  "Last Communication Contact from Initial Pathologic Diagnosis Date",
  "Birth from Initial Pathologic Diagnosis Date",
  "Date Calculated Day Value",
  "Prior Diagnosis",
  "Race Category",
  "Number of Samples Per Patient",
  "Sample Type",
  "Somatic Status",
  "Tissue Prospective Collection Indicator",
  "Tissue Retrospective Collection Indicator",
  "Tissue Source Site",
  "Tissue Source Site Code",
  "Tumor Disease Anatomic Site",
  "Last Alive Less Initial Pathologic Diagnosis Date Calculated Day Value",
  "Other Patient ID",
  "Tumor Type",
  "Patient Weight",
  "Neoadjuvant Therapy Type Administered Prior To Resection Text",
  "International Classification of Diseases for Oncology, Third Edition ICD-O-3 Histology Code",
  "International Classification of Diseases for Oncology, Third Edition ICD-O-3 Site Code",
  "Classification",
  "Informed consent verified",
  "In PanCan Pathway Analysis",
  "New Neoplasm Event Post Initial Therapy Indicator",
  "Oncotree Code",
  "Person Neoplasm Cancer Status"
)]

merged$`Neoplasm Disease Stage American Joint Committee on Cancer Code` <- factor(merged$`Neoplasm Disease Stage American Joint Committee on Cancer Code`,
                                                                                  levels = c("STAGE I","STAGE II", "STAGE III", "STAGE IVA", "STAGE IVB", "STAGE IVC" ),
                                                                                  labels = c("Stage I", "Stage II", "Stage III", "Stage IV", "Stage IV", "Stage IV"))

merged$`Disease Free Status` <- factor(merged$`Disease Free Status`,
                                       levels = c("0:DiseaseFree","1:Recurred/Progressed"),
                                       labels = c("Disease Free", "Recurred/Progressed"))

merged$`Disease-specific Survival status` <- factor(merged$`Disease-specific Survival status`,
                                                    levels = c("0:ALIVE OR DEAD TUMOR FREE","1:DEAD WITH TUMOR"),
                                                    labels = c("Alive/Tumor free", "Dead with tumor"))

merged$`Overall Survival Status`<- factor(merged$`Overall Survival Status`,
                                          levels = c("0:LIVING","1:DECEASED"),
                                          labels = c("Alive", "Deceased"))

merged$`Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code` <- factor(merged$`Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code`,
                                                                                             levels = c("N0","N1", "N2", "N2A", "N2B", "N2C", "N3", "NX"),
                                                                                             labels = c("N0", "N1", "N2", "N2", "N2", "N2", "N3", "NX"))

merged$`American Joint Committee on Cancer Tumor Stage Code`<- factor(merged$`American Joint Committee on Cancer Tumor Stage Code`,
                                                                      levels = c("T1","T2", "T3", "T4", "T4A", "T4B", "TX", "T0"),
                                                                      labels = c("T1", "T2", "T3", "T4", "T4", "T4", "TX", "T0"))

merged$`Progression Free Status`<- factor(merged$`Progression Free Status`,
                                          levels = c("0:CENSORED","1:PROGRESSION"),
                                          labels = c("Censored", "Progression"))

merged$Subtype<- factor(merged$Subtype,
                        levels = c("HNSC_HPV+","HNSC_HPV-"),
                        labels = c("HPV+", "HPV-"))

#############################
# Analyze clinical features #
#############################
# Define the check_normality function
check_normality <- function(data) {
  # Step 2: Identify numerical columns
  numerical_columns <- sapply(data, is.numeric)
  
  # Step 3: Perform Shapiro-Wilk test on each numerical column
  normality_results <- lapply(names(data)[numerical_columns], function(col_name) {
    col <- data[[col_name]]
    shapiro_test <- shapiro.test(col)
    return(list(
      Column = col_name,
      W = shapiro_test$statistic,
      p_value = shapiro_test$p.value,
      Normal = shapiro_test$p.value > 0.05
    ))
  })
  
  # Step 4: Identify normally distributed columns
  vars_normal <- names(data)[numerical_columns][sapply(normality_results, function(result) result$Normal)]
  return(vars_normal)
}

# Run the normality test to determine normally distributed columns
vars_normal <- check_normality(merged)

# Define the rndr function based on normally distributed columns
rndr <- function(x, name, ...) {
  cont <- ifelse(name %in% vars_normal, "Mean (SD)", "Median (Min, Max)")
  render.default(x, name, render.continuous = c("", cont), ...)
}

pvalue <- function(x, ...) {
  x <- x[-length(x)] # Remove "overall" group
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length))) # nolint
  if (is.numeric(y)) {
# check for normality using Shapiro-Wilk test and variance equality using var.test
    swtest <- shapiro.test(y)$p.value
    vartest <- var.test(y ~ g)$p.value
    if (all(swtest > 0.05, vartest > 0.05)) {
# If both p-values are >0.05 we assume normality and use t.test with equal variances
      p <- t.test(y ~ g, var.equal=TRUE)$p.value
      symbol <- "T-Equal"
    } else {
# Otherwise, we assume normality with unequal variances and use t.test with var.equal=F
      p <- t.test(y ~ g, var.equal=FALSE)$p.value
      symbol <- "T-Unequal"
    } 
    if (all(swtest <= 0.05)) {
      # Otherwise, reject normality and use Wilcoxon test
      p <- wilcox.test(y ~ g, exact=FALSE)$p.value
      symbol <- "Wil"
    }} 
  if (is.factor(y)) {
# For categorical variables, perform a chi-squared test of independence
    p <- tryCatch({
# Si warning chi-test, rC)aliser un Fisher
      chisq.test(table(y, g))$p.value
    }, warning = function(w) {
# If a warning appears (due to expected cell counts <5 or too few observations), use Fisher's Exact test
# If there is a warning message, print it and change the condition
      print(paste("Warning message:", w))
      return(NULL)
    })
    if (!is.null(p)) {
      p <- chisq.test(table(y, g))$p.value
      symbol <- "Chi"
    } else
    {
      p <- fisher.test(table(y, g), simulate.p.value=TRUE)$p.value
      symbol <- "Fish"
    }
  }
# Format the p-value, using an HTML entity for the less-than sign.
# The initial empty string places the output on the line below the variable label.
 c(paste0(sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)), "<sup>", symbol, "</sup>"))
}

caption <- c("Table 1: Association of clinical parameters 
between TRIM32 expression groups")

footnote <- c("ns= Not Significant &nbsp
               <=0.05 = * &nbsp;
               <0.01 = ** &nbsp;
               <0.001 = *** &nbsp;
               <0.0001 = ****  &nbsp;")


table <- table1(~ Sex+
                `Diagnosis Age`+
                 Subtype+
                `site_of_resection_or_biopsy`+
                `Neoplasm Histologic Grade`+
                `Neoplasm Disease Stage American Joint Committee on Cancer Code`+
                `American Joint Committee on Cancer Tumor Stage Code`+
                `Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code`+
                `American Joint Committee on Cancer Metastasis Stage Code`+
                `Overall Survival Status`+
                `Disease-specific Survival status`+
                `Progression Free Status`+
                 `Disease Free Status`+
                `Aneuploidy Score`+
                `Buffa Hypoxia Score`+
                `Ragnum Hypoxia Score`+ 
                `Winter Hypoxia Score`
                | TRIM32_expression, data = merged,
                extra.col = list(`p-value` = pvalue), caption = caption,
                footnote = footnote, extra.col.pos = 3,
                render = rndr, render.categorical = "FREQ (PCTnoNA%)",
                topclass = " Rtable1-times ")

print(table)

#Make lookup table for the p-value function
symbols <- c("****", "***", "**", "*", "ns")
cutoffs <- c(0.0001, 0.001, 0.01, 0.05)
lookup_table <- setNames(symbols, cutoffs)

#Change the p-value function so that the symbol are not included
pvalue <- function(x, ...) {
  x <- x[-length(x)] # Remove "overall" group
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    # check for normality using Shapiro-Wilk test and variance equality using var.test
    swtest <- shapiro.test(y)$p.value
    vartest <- var.test(y ~ g)$p.value
    if (all(swtest > 0.05, vartest > 0.05)) {
      # If both p-values are >0.05 we assume normality and use t.test with equal variances
      p <- t.test(y ~ g, var.equal=TRUE)$p.value
      symbol <- lookup_table[findInterval(p, cutoffs)+1]
    } else {
      # Otherwise, we assume normality with unequal variances and use t.test with var.equal=F
      p <- t.test(y ~ g, var.equal=FALSE)$p.value
      symbol <- lookup_table[findInterval(p, cutoffs)+1]
    } 
    if (all(swtest <= 0.05)) {
      # Otherwise, reject normality and use Wilcoxon test
      p <- wilcox.test(y ~ g, exact=FALSE)$p.value
      symbol <- lookup_table[findInterval(p, cutoffs)+1]
    }} 
  if (is.factor(y)) {
    # For categorical variables, perform a chi-squared test of independence
    p <- tryCatch({
      # Si warning chi-test, rC)aliser un Fisher
      chisq.test(table(y, g))$p.value
    }, warning = function(w) {
      # If a warning appears (due to expected cell counts <5 or too few observations), use Fisher's Exact test
      # If there is a warning message, print it and change the condition
      print(paste("Warning message:", w))
      return(NULL)
    })
    if (!is.null(p)) {
      p <- chisq.test(table(y, g))$p.value
      symbol <- lookup_table[findInterval(p, cutoffs)+1]
    } else 
    {
      p <- fisher.test(table(y, g), simulate.p.value=TRUE)$p.value
      symbol <- lookup_table[findInterval(p, cutoffs)+1]
    }
  } 
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  c(paste0(sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)), "<sup>", symbol, "</sup>"))
}

label(merged$`Diagnosis Age`)    <- "Age"
label(merged$site_of_resection_or_biopsy) <- "Tumor Site"
label(merged$`Neoplasm Histologic Grade`)             <- "Histologic Grade"
label(merged$`Neoplasm Disease Stage American Joint Committee on Cancer Code` )    <- "Tumor Stage"
label(merged$`American Joint Committee on Cancer Tumor Stage Code`)    <- "Tumor Stage Code"
label(merged$`Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code`) <- "Lymph Node Code"
label(merged$`American Joint Committee on Cancer Metastasis Stage Code`) <- "Metastasis Stage Code"



table <- table1(~ Sex+
                  `Diagnosis Age`+
                  Subtype+
                  `site_of_resection_or_biopsy`+
                  `Neoplasm Histologic Grade`+
                  `Neoplasm Disease Stage American Joint Committee on Cancer Code`+
                  `American Joint Committee on Cancer Tumor Stage Code`+
                  `Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code`+
                  `American Joint Committee on Cancer Metastasis Stage Code`+
                  `Aneuploidy Score`+
                  `Buffa Hypoxia Score`+
                  `Ragnum Hypoxia Score`+ 
                  `Winter Hypoxia Score`
                | TRIM32_expression, data = merged,
                extra.col = list(`p-value` = pvalue), caption = caption,
                footnote = footnote, extra.col.pos = 3,
                render = rndr, render.categorical = "FREQ (PCTnoNA%)",
                topclass = " Rtable1-times ")

print(table)


#########################
# Write the output file #
#########################
#https://github.com/benjaminrich/table1/issues/33                                                      
write_table1 <- function(x,                   # a table1 object
                         file,                # path to output .pdf file
                         landscape = FALSE,   # landscape print?
                         scale = 0.5,           # scaling factor
                         width = 8.5,         # width of resulting pdf (in.)
                         height = 11,         # height of resulting pdf (in.)
                         dump_html = TRUE) {  # delete intermediate html files?
  
  x <- htmltools::HTML(x)
  src <- system.file(package = "table1", "table1_defaults_1.0")
  default.style <- htmltools::htmlDependency("table1", "1.0",
                                             src = src,
                                             stylesheet = "table1_defaults.css")
  
  x <- htmltools::div(class = "Rtable1", default.style, x)
  x <- htmltools::browsable(x)
  
  htmltools::save_html(x, file = gsub(".pdf", ".html", file))
  pagedown::chrome_print(input = gsub(".pdf", ".html", file),
                         output = file,
                         options = list(landscape = landscape,
                                        scale = scale,
                                        paperWidth = width,
                                        paperHeight = height))
  if(dump_html) unlink(gsub(".pdf", ".html", file))
}

write_table1(table, "2025_03_19_clinical_table_hnsc_trim32_high_vs_low_v3.pdf")
