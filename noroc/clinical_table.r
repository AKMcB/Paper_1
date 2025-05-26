#######################
## Loading libraries ##
######################

library(data.table)
library(tidyverse)
library(table1)
library(htmltools)
library(pagedown)
library(readxl)
library(haven)
####################
##Load input files##
####################
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

info <- haven::as_factor(info)

info2 <- fread("2025_03_18_clinical_info_combined.csv")
info2 <- haven::as_factor(info2)
info2 <- info2[,c(3:5, 8,16,29, 30,40,42 ,50)]

difference <- setdiff(info$id_tma,info2$id_tma)
# Filter df1 based on the difference
result <- info[info$id_tma %in% difference, ]

info <- merge(info, info2, by= "id_tma") #131 patients- O-10 missing from info2

#Add scoring information 
info3 <- read.csv2("t32_scores_level_noroc.csv")
info3 <- info3[,-c(2:4)]

#2 patients missing -B12 and O157- too few cells

info4 <- read_xlsx("../myc_scoring/myc_scoring.xlsx", sheet = 8)
info4 <- dplyr::select(info4, c(id_sample,myc_value))

# 3 patients missing- B-12, T-4, T-45- too few cells 

info_comb <- merge(info, info3, by.x = "id_tma", by.y = "ID") #130 patients B-12 was in the first info file
info_comb <- merge(info_comb, info4, by.x ="id_tma", by.y = "id_sample") #130 patients - B-12, t4, T45 was already gone from the first file

#remove patients: pallitative,in situ, missing pT, retracted
remove_pats <- c("T-8", "T-28","T-61", "T-63", 
                 "O-10","O-15", "O-25", "O-67",
                 "O-188") #9

info_comb <- info_comb[!(info_comb$id_tma %in% remove_pats),] #131-123, 8 removed
#T-28 were already removed, the rest was removed

#Remove B-99 due to missing data 
info_comb <- info_comb[-20,] #122 patients

info_comb <- subset(info_comb, info_comb$Follow_up_months > 0)

####################
##Analyze the data##
####################
str(info_comb)

info_comb <- info_comb %>% 
  mutate(across(where(is.factor), as.character)) %>%  
  mutate(across(where(is.character), ~ na_if(., "Missing")))

info_comb <- info_comb %>% 
  mutate(across(where(is.factor), as.character)) %>%  
  mutate(across(where(is.character), ~na_if(.,"Missing/ not evaluable")))

info_comb <- info_comb %>% 
  mutate(across(where(is.factor), as.character)) %>%  
  mutate(across(where(is.character), ~na_if(.,"Missing/not evaluable")))

info_comb <- info_comb %>% 
  mutate(across(where(is.factor), as.character)) %>%  
  mutate(across(where(is.character), ~na_if(.,"Unknown/ Not evaluable")))

info_comb <- info_comb %>% 
  mutate(across(where(is.factor), as.character)) %>%  
  mutate(across(where(is.character), ~na_if(.,"Missing/unknown")))

info_comb <- info_comb %>% 
  mutate(across(where(is.factor), as.character)) %>%  
  mutate(across(where(is.character), ~na_if(.,"")))

#input_tab <- input_tab %>% mutate(across(where(is.character), ~na_if(.,"Missing")))

unique(info_comb$Smoking)
info_comb$Alcohol <- factor(info_comb$Alcohol, 
                            levels = c("Never", "Current", "Seldom", "Moderately", "Heavy", "Former alcholic abuse"),
                            labels = c("Never", "Current", "Current", "Current", "Current", "Current"))


info_comb$cT <- factor(info_comb$cT,
                       levels = c("cT1", "cT2", "cT3", "cT4", "cT4a"),
                       labels = c("T1", "T2", "T3", "T4", "T4"))

info_comb$cN <- factor(info_comb$cN,
                       levels = c("cNx", "cN0", "cN1" ,"cN2", "cN2a", "cN2b", "cN2c" ,"cN3"),
                       labels = c("Missing","N0", "N+", "N+", "N+", "N+", "N+", "N+"))


info_comb$pT_8th_edition_TNM <- factor(info_comb$pT_8th_edition_TNM,
                               levels = c("pT1", "pT2", "pT3"),
                               labels = c("T1", "T2", "T3"))
unique(info_comb$pN)
info_comb$pN <- factor(info_comb$pN,
                       levels = c("pN0", "pN1","pN2", "pN2a", "pN2b", "pN3"),
                       labels = c("N0", "N+", "N+", "N+", "N+", "N+"))

info_comb <- info_comb %>% 
  mutate(pT_cT = if_else(is.na(info_comb$pT_8th_edition_TNM), 
                         info_comb$cT, info_comb$pT_8th_edition_TNM))

info_comb <- info_comb %>% 
  mutate(pN_cN = if_else(is.na(info_comb$pN),
                         info_comb$cN, info_comb$pN))


info_comb$Keratinization <- factor(info_comb$Keratinization,
                                   levels = c("Highly keratinized (>50 % of tumor areal)",
                                              "Moderately keratinized (20-50% of tumor areal)",
                                              "No keratinisation (0-5 % of tumor areal)",
                                              "Minimal keratinisation (5-20% of tumor areal)"),
                                   labels = c("Highly keratinized (>50 % of tumor areal)",
                                              "Moderately keratinized (20-50% of tumor areal)",
                                              "No to minimal keratinisation (0-20% of tumor areal)", 
                                              "No to minimal keratinisation (0-20% of tumor areal)"))

unique(info_comb$Keratinization)
info_comb$NuclearPolymorphism <- factor(info_comb$NuclearPolymorphism,
                                        levels = c("Extreme nuclear polymorphism (in 75-100% of cells)",
                                                   "Abundant nuclear polymorphism (in 50-75% of the cells)",
                                                   "Moderately abundant nuclear polymorphism (in 25-50% of the cells)",
                                                   "Little nuclear polymorphism (in less than 25% of the cells)"),
                                        labels = c("Extreme nuclear polymorphism (in 75-100% of cells)",
                                                   "Abundant nuclear polymorphism (in 50-75% of the cells)",
                                                   "Little to moderate nuclear polymorphism (in 0-50% of the cells)", 
                                                   "Little to moderate nuclear polymorphism (in 0-50% of the cells)"))


unique(info_comb$Staging_101022)


#There is one with 93, but I made NA as there is no information about stage
info_comb$Staging_101022 <- factor(info_comb$Staging_101022,
                                   levels = c("Stage I", "Stage I eller II", "Stage II", "Stage III", "Stage IV"),
                                   labels = c("Stage I", "Stage I/II", "Stage II", "Stage III", "Stage IV"))

info_comb$t32_level <- factor(info_comb$t32_level,
                              levels = c("high", "moderate", "low", "negative"),
                              labels = c("High/Moderate", "High/Moderate", "Low/Negative", "Low/Negative"))

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
  vars.normal <- names(data)[numerical_columns][sapply(normality_results, function(result) result$Normal)]
  
  return(vars.normal)
}

# Run the normality test to determine normally distributed columns
vars.normal <- check_normality(info_comb)

rndr <- function(x, name, ...) {
  cont <- ifelse(name %in% vars.normal, "Mean (SD)", "Median (Min, Max)")
  render.default(x, name, render.continuous=c("", cont), ...)
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
      p <- t.test(y ~ g,var.equal=FALSE)$p.value
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

caption <-  c("Table 1: Association of clinical parameters and TRIM32 TMA Score")

footnote <- c("T-Equal = t.test with equal variance
               T-Unequal = t.test with unequal variance
               ns= not significant
               SD = standard deviation 
               OS = Overall Survival 
               DSS = Death Specific Survival")

colnames(info_comb)
str(info_comb$Age_at_diagnose)
table <- table1(~ Gender+
                  Age_at_diagnose+
                  Staging_101022+
                  pT_cT+
                  pN_cN+
                  DAAR+
                  Alcohol+
                  Smoking+
                  Differentiation+
                  Woolgar_Depht_Invasion+
                  Tumor_Thickness_mm+
                  Diameter_mm+
                  Keratinization+
                  NuclearPolymorphism+
                  myc_value|t32_level,
                data = info_comb,
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

label(info_comb$Staging_101022)    <- "Stage"
label(info_comb$Age_at_diagnose) <- "Age"
label(info_comb$pT_cT)             <- "Tumor Stage"
label(info_comb$pN_cN )      <- "Lymph Node"
label(info_comb$DAAR)            <- "Cause of Death"
label(info_comb$Woolgar_Depht_Invasion) <- "Woolgar Depht Invasion"
label(info_comb$Tumor_Thickness_mm) <- "Tumor Thickness (mm)"
label(info_comb$Diameter_mm) <- "Diameter (mm)"
label(info_comb$NuclearPolymorphism) <- "Nuclear Polymorphism"
label(info_comb$myc_value)         <- "MYC level"

fwrite(info_comb, "testing_rmd_file.csv", row.names = T)
table <- t1flex(table1(~ Gender+
                  Age_at_diagnose+
                  Staging_101022+
                  pT_cT+
                  pN_cN+
                  DAAR+
                  Alcohol+
                  Smoking+
                  Differentiation+
                  Woolgar_Depht_Invasion+
                  Tumor_Thickness_mm+
                  Diameter_mm+
                  Keratinization+
                  NuclearPolymorphism+
                  myc_value|t32_level, data = info_comb,
                extra.col = list(`p-value` = pvalue), caption = caption,
                footnote = footnote, extra.col.pos = 3,
                render = rndr, render.categorical = "FREQ (PCTnoNA%)",
                topclass = " Rtable1-times "))

print(table)

# Define section properties with standard page size
sect_properties <- officer::prop_section(
  page_size = officer::page_size(
    orient = "portrait",
    width = 8.5, height = 11), # Standard Letter size
  type = "continuous")
# Create a flextable and adjust its properties
table_df <- as.data.frame(table)
my_table <- flextable::flextable(table_df)
my_table <- flextable::fontsize(my_table, size = 8) 
my_table <- flextable::padding(my_table, padding.top = 0, padding.bottom = 0)
my_table <- flextable::set_table_properties(my_table, width = 1, layout = "autofit")
# Save the table as a Word document
my_table %>% 
  save_as_docx(path = "test.docx", pr_section = sect_properties)
?save_as_docx
#########################
##Write the output file##
#########################
write_table1 <- function(x,                   # a table1 object
                         file,                # path to output .pdf file
                         landscape = FALSE,   # landscape print?
                         scale = 1,           # scaling factor
                         width = 8.5,         # width of resulting pdf (in.)
                         height = 17,         # height of resulting pdf (in.)
                         dump_html = TRUE) {  # delete intermediate html files?
  
  x <- htmltools::HTML(x)
  src <- system.file(package = "table1", "table1_defaults_1.0")
  default.style <- htmltools::htmlDependency("table1", "1.0",
                                             src = src,
                                             stylesheet = "table1_defaults.css",
                                             head = '<link rel="icon" href="data:,">')
  
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

write_table1(table, "2025_03_18_clinical_table_myc_t32_scoring_all_patients.pdf")

