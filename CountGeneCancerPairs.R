# In case you are running files on a cluster you need to load R module 
# system("module load R/4.0.2-foss-2018a-bare")

# Load the necessary libraries

library(dplyr)
library(tidyr)

# Read in the merged data from the CSV file
merged_data <- read.csv("merged_data_nox.csv", stringsAsFactors = FALSE)

# Create a new data frame to hold the counts for each combination
counts <- merged_data %>%
  group_by(acronym, Hugo_Symbol) %>%
  summarize(mut_only = sum(Mutation == 1 & (CNV %in% c(0, NA))),
            cnv_positive = sum(is.na(Mutation) & CNV %in% c(1, 2)),
            cnv_negative = sum(is.na(Mutation) & CNV %in% c(-1, -2)),
            nothing = sum(is.na(Mutation) & CNV == 0),
            cnv_positive_mut = sum(Mutation == 1 & CNV %in% c(1, 2)),
            cnv_negative_mut = sum(Mutation == 1 & CNV %in% c(-1, -2)))

# Create a data frame with all possible combinations of "acronym" and "Hugo_Symbol"
all_combinations <- merged_data %>%
  select(acronym, Hugo_Symbol) %>%
  distinct()

# Merge the counts and all_combinations data frames, filling in 0 for any missing counts
final_table <- merge(all_combinations, counts, by = c("acronym", "Hugo_Symbol"), all.x = TRUE) %>%
  replace(is.na(.), 0)

# Write the final table to a CSV file
write.csv(final_table, "final_table.csv", row.names = FALSE)