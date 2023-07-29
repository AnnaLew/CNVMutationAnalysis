# Load the necessary libraries

library(dplyr)
library(tidyr)

# Read in the merged data from the CSV file
merged_data <- read.csv("data/merged/merged_data.csv", stringsAsFactors = FALSE)

# Create a new data frame to hold the counts for each combination
counts <- merged_data %>%
  group_by(acronym, Hugo_Symbol) %>%
  summarize(mut_only = sum(Mutation == 1 & (CNV %in% c(0, NA)), na.rm = TRUE),
            cnv_positive = sum(is.na(Mutation) & CNV %in% c(1, 2), na.rm = TRUE),
            cnv_negative = sum(is.na(Mutation) & CNV %in% c(-1, -2), na.rm = TRUE),
            nothing = sum(is.na(Mutation) & CNV == 0, na.rm = TRUE),
            cnv_positive_mut = sum(Mutation == 1 & CNV %in% c(1, 2), na.rm = TRUE),
            cnv_negative_mut = sum(Mutation == 1 & CNV %in% c(-1, -2), na.rm = TRUE))

# Merge the counts and all_combinations data frames, filling in 0 for any missing counts
counted_pairs <- counts %>%
  replace(is.na(.), 0) %>%
  filter(mut_only + cnv_positive + cnv_negative + nothing + cnv_positive_mut + cnv_negative_mut >= 5)

# Write the final table to a CSV file
dir.create("data/countedPairs", recursive = TRUE)
write.csv(counted_pairs, "data/countedPairs/counted_pairs.csv", row.names = FALSE)
