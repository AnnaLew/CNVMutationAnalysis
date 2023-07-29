# In case you are running files on a cluster you need to load R module 
# system("module load R/4.0.2-foss-2018a-bare")

# Load the necessary library

library(data.table)

#load the data
cnv_data_processed <- read.csv("data/processed/cnv_data_processed.csv")
mutation_data_merged <- read.csv("data/preprocessed/mutation_data_merged.csv")

#convert data.frames to data.tables
cnv_data_processed <- data.table(cnv_data_processed)
mutation_data_merged <- data.table(mutation_data_merged)

# Merge mut_sample and cnv_sample data.tables based on "Hugo_Symbol", "bcr_patient_barcode", and "acronym"
merged_data <- merge(cnv_data_processed, mutation_data_merged,
                      by = c("Hugo_Symbol", "bcr_patient_barcode", "acronym"),
                      all = TRUE, row.names = FALSE, suffixes = c("", ""))

# Select only the columns you want to keep
merged_data <- merged_data[, c("Hugo_Symbol", "acronym", "CNV", "Mutation")]

# Write merged data.table to a CSV file
write.csv(merged_data, file = "data/merged/merged_data.csv", row.names = FALSE)