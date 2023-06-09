# In case you are running files on a cluster you need to load R module 
# system("module load R/4.0.2-foss-2018a-bare")

# Load the necessary library

library(tidyverse)

#Load the data

cnv_data_preprocessed = read.csv("cnv_data_preprocessed.csv", head = TRUE, row.names = 1)
acronyms <- readRDS("acronyms.rds")

# Make a table with Hugo symbol, barcode and CNV 

cnv_data_processed <- pivot_longer(cnv_data_preprocessed, cols = c(-Hugo_Symbol), names_to = "bcr_patient_barcode", values_to = "CNV")

# Convert my_df to data.table and append column 
setDT(cnv_data_processed) 
cnv_data_processed[, acronym := rep(acronyms, each = nrow(cnv_data_preprocessed))]

# Save processed cnv data
write.csv(cnv_data_processed, "cnv_data_processed.csv")
