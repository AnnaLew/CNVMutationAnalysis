# Fisher's exact test gene enrichment

# Asses whether gene/cancer pair is enriched for co-occuring mutations and CNVs 
# (done separately for positive and negative CNVs)

#                       cases with mut	           cases without mutation
#  cases with cnv	     cnv_positive_mut	                 cnv_positive    
#  cases without cnv	        mut                           nothing

# Load the necessary library

library(dplyr) 

# Load the data
combinations_counts <- read.csv("combinations_counts.csv", head = TRUE)

# Check for missing values
sum(is.na(combinations_counts))

# Remove all rows with NA values
combinations_counts <- na.omit(combinations_counts)

# Remove rows with zeros in columns with co-occuring mutations and CNVs
# Create separate dataframes for positive and negative CNVs
positive <- combinations_counts %>%
  filter(cnv_positive_mut > 0)
negative <- combinations_counts %>%
  filter(cnv_negative_mut > 0)

# Perform Fisher's exact test for each gene-cancer pairs

for (i in 1:nrow(positive)) {
  matrix_values <- matrix(c(positive[i, "cnv_positive_mut"], positive[i, "cnv_positive"], 
                            positive[i, "mut_only"], positive[i, "nothing"]),
                          nrow = 2)
  
  fisher_result <- fisher.test(matrix_values)
  positive[i, "p.value"] <- fisher_result$p.value
  positive[i, "odds"] <- fisher_result$estimate

}

for (i in 1:nrow(negative)) {
  matrix_values <- matrix(c(negative[i, "cnv_negative_mut"], positive[i, "cnv_negative"], 
                            negative[i, "mut_only"], positive[i, "nothing"]),
                          nrow = 2)
  
  fisher_result <- fisher.test(matrix_values)
  negative[i, "p.value"] <- fisher_result$p.value
  negative[i, "odds"] <- fisher_result$estimate
}

# Check number of p.values below 0.1
sum(positive$p.value < 0.1)
sum(negative$p.value < 0.1)

# Adjust p-values
positive$p.adj <- p.adjust(positive$p.value, method = "BH")
negative$p.adj <- p.adjust(negative$p.value, method = "BH")

# Check number of p.adj.values below 0.05
sum(positive$p.adj < 0.1)
sum(negative$p.adj < 0.1)

# Save files
dir.create("data/fishersTest", recursive = TRUE)
write.csv(positive, "data/fishersTest/fishers_positive.csv", row.names = FALSE)
write.csv(negative, "data/fishersTest/fishers_negative.csv", row.names = FALSE)

