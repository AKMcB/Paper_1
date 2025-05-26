library(tidyverse)

setwd("C:/Users/abe186/UiT Office 365/O365-Bioinformatikk TRIM27 - General/TCGA-Head and neck/TMM conversion")
expr <- read.csv2("TMM_TCGA_HNSC_counts_NoNorma_scaled_filtered.csv", sep = ";", as.is = T,check.names = F)
expr <- expr[complete.cases(expr),] #Remove NAs
rownames(expr) <- expr[,1]
expr <- expr[,-1]
#expr <- as.data.frame(t(expr))
expr <- tibble::rownames_to_column(expr, "id")

expr$id <- substr(expr$id, 1, nchar(expr$id) - 4)
expr <- distinct(expr,id,.keep_all = T)

setwd("C:/Users/abe186/UiT Office 365/O365-Bioinformatikk TRIM27 - General/TCGA-Head and neck/Heatmap/")
ann <- read.csv2("Clinical _info_HNSC_V2_no_duplicates.csv", sep = ";",
                 as.is = T, check.names=F)
#Subset based on location 
ann <- ann[, c(2,3,71,120)]
ann <- ann[-c(528,529),]
#write.csv2(ann, "hnsc_locations_report.csv")

expr <- expr[expr$id %in% ann$case_submitter_id, ]
ann <- ann[ann$case_submitter_id %in% expr$id,]


# Use table() to count occurrences of each unique string in the column
string_counts <- table(ann$site_of_resection_or_biopsy)

# Print the counts
print(string_counts)

patterns_pahrynx <- c("Oropharynx, NOS","Posterior wall of oropharynx",
                      "Hypopharynx, NOS", "Tonsil, NOS", "Base of tongue, NOS")

# Define the patterns you want to replace as a regular expression

pattern_regex <- paste(patterns_pahrynx, collapse = "|")

# Define the replacement string
replacement <- "Pharynx, NOS"


# Use gsub() to replace multiple patterns with the replacement string
ann$site_of_resection_or_biopsy <- gsub(pattern_regex, replacement, ann$site_of_resection_or_biopsy)
string_counts <- table(ann$site_of_resection_or_biopsy)
print(string_counts)

patterns_tongue <- c("Ventral surface of tongue, NOS","Border of tongue")

# Define the patterns you want to replace as a regular expression

pattern_regex <- paste(patterns_tongue, collapse = "|")

# Define the replacement string
replacement <- "Tongue, NOS"


# Use gsub() to replace multiple patterns with the replacement string
ann$site_of_resection_or_biopsy <- gsub(pattern_regex, replacement, ann$site_of_resection_or_biopsy)
string_counts <- table(ann$site_of_resection_or_biopsy)
print(string_counts)

patterns_gum <- c("Lower gum","Upper gum", "Hard palate", "Palate, NOS")

# Define the patterns you want to replace as a regular expression

pattern_regex <- paste(patterns_gum, collapse = "|")

# Define the replacement string
replacement <- "Gum, NOS"


# Use gsub() to replace multiple patterns with the replacement string
ann$site_of_resection_or_biopsy <- gsub(pattern_regex, replacement, ann$site_of_resection_or_biopsy)
string_counts <- table(ann$site_of_resection_or_biopsy)
print(string_counts)

patterns_floor <- c("Anterior floor of mouth")

# Define the patterns you want to replace as a regular expression

pattern_regex <- paste(patterns_floor, collapse = "|")

# Define the replacement string
replacement <- "Floor of mouth, NOS"


# Use gsub() to replace multiple patterns with the replacement string
ann$site_of_resection_or_biopsy <- gsub(pattern_regex, replacement, ann$site_of_resection_or_biopsy)
string_counts <- table(ann$site_of_resection_or_biopsy)
print(string_counts)

patterns_lar <- c("Supraglottis")

# Define the patterns you want to replace as a regular expression

pattern_regex <- paste(patterns_lar, collapse = "|")

# Define the replacement string
replacement <- "Larynx, NOS"


# Use gsub() to replace multiple patterns with the replacement string
ann$site_of_resection_or_biopsy <- gsub(pattern_regex, replacement, ann$site_of_resection_or_biopsy)
string_counts <- table(ann$site_of_resection_or_biopsy)
print(string_counts)

#Remove lip since that can be classified at melanoma 
#Remove retromolar areas as this site can not be defined further 

ann <- subset(ann,!ann$site_of_resection_or_biopsy %in% c("Lip, NOS", "Retromolar area"))
unique(ann$site_of_resection_or_biopsy)

write.csv2(ann, "../Heatmap/tcga_hnsc_updated_locations_03112023.csv")
