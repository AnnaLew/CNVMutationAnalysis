# Load the necessary library
library(data.table)

#load the data
cnv_data_processed <- read.csv("data/processed/cnv_data_processed.csv")
mutation_data_merged <- read.csv("data/preprocessed/mutation_data_merged_simplified.csv")

#convert data.frames to data.tables
cnv_data_processed <- data.table(cnv_data_processed)
mutation_data_merged <- data.table(mutation_data_merged)

# Merge mut_sample and cnv_sample data.tables based on "Hugo_Symbol", "bcr_patient_barcode", and "acronym"
merged_data <- merge(cnv_data_processed, mutation_data_merged,
                      by = c("Hugo_Symbol", "bcr_patient_barcode", "acronym"),
                      all = TRUE, row.names = FALSE, suffixes = c("", ""))

# Select only the columns you want to keep
merged_data <- select(merged_data, -bcr_patient_barcode)

# Write merged data.table to a CSV file
dir.create("data/merged", recursive = TRUE)
write.csv(merged_data, file = "data/merged/merged_data.csv", row.names = FALSE)