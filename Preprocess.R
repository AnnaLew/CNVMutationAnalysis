# Load the necessary library
library(dplyr) 

# Read all the files
mutation_data = read.table("data/original/mc3.v0.2.8.PUBLIC.maf.gene_vclass_HGVSp_sample_likelyDriverLoose.tsv", head = TRUE)
patient_data = read.csv("data/original/clinical_PANCAN_patient_with_followup.tsv", header = TRUE, sep = "\t")
cnv_data = read.csv("data/original/all_thresholded.by_genes_whitelisted.tsv", header = TRUE, sep = "\t")

cnv_data_preprocessed <- cnv_data
names(cnv_data_preprocessed)[1] ="Hugo_Symbol"

# Add the tumor type information from patient_data to mutation_data
mutation_data_merged <- merge(mutation_data, patient_data, by="bcr_patient_barcode", all.x = TRUE, all.y = FALSE) %>%
  select(bcr_patient_barcode, acronym, Hugo_Symbol, Variant_Classification)

mutation_data_merged$bcr_patient_barcode <- gsub(pattern = "-", replacement= ".", mutation_data_merged$bcr_patient_barcode)

# Rename the "Variant_Classification" column to "Mutation"
mutation_data_merged_simplified <- mutation_data_merged
colnames(mutation_data_merged_simplified)[colnames(mutation_data_merged_simplified) == "Variant_Classification"] <- "Mutation"

# Modify the "Mutation" column in the mutation_data_merged_simplified table, changing any non-NA value to 1
mutation_data_merged_simplified$Mutation[!is.na(mutation_data_merged_simplified$Mutation)] <- 1

#Shorten and edit the names in cnv_data so that they can be matched with patient_data
names(cnv_data_preprocessed) <- substring(names(cnv_data_preprocessed), 1, 12)
names(cnv_data_preprocessed) <- gsub(pattern = "[.]", replacement= "-", names(cnv_data_preprocessed))

#Copy column names to first row, so that it will be easier to work on it (colnames should be unique)
cnv_data_preprocessed <- rbind(colnames(cnv_data_preprocessed), cnv_data_preprocessed)

#Find acronyms corresponding to the patient barcodes and overwrite them
patient_barcodes <- unique(patient_data$bcr_patient_barcode)

for (colnum in 4:ncol(cnv_data_preprocessed)) {
  if (cnv_data_preprocessed[1, colnum] %in% patient_barcodes) {
    cnv_data_preprocessed[1, colnum] = patient_data$acronym[which(cnv_data_preprocessed[1, colnum] == patient_data$bcr_patient_barcode)]
  } else {
    cnv_data_preprocessed[1, colnum] = NA
  }
}

# Remove columns that did not have matching acronyms
cnv_data_preprocessed <- cnv_data_preprocessed[ , colSums(is.na(cnv_data_preprocessed[1,]))==0]

# Remove two unnecesary columns
cnv_data_preprocessed <- cnv_data_preprocessed[, -c(2:3)]

# Change the colnames to match the ones in mutations data
colnames(cnv_data_preprocessed) <- gsub(pattern = "-", replacement= ".", colnames(cnv_data_preprocessed))

# Create list with all acronyms in the right order
acronyms <- as.character(as.vector(cnv_data_preprocessed[1,]))
acronyms <- acronyms[-1]

# Remove the first row
cnv_data_preprocessed <- cnv_data_preprocessed[-1, ]

#Save merged mutation data
dir.create("data/preprocessed", recursive = TRUE)
write.csv(mutation_data_merged, "data/preprocessed/mutation_data_merged.csv", row.names = FALSE)
write.csv(mutation_data_merged_simplified, "data/preprocessed/mutation_data_merged_simplified.csv", row.names = FALSE)

# Save preprocessed cnv data
write.csv(cnv_data_preprocessed, "data/preprocessed/cnv_data_preprocessed.csv")

# Save all acronyms
saveRDS(acronyms, file = "acronyms.rds")