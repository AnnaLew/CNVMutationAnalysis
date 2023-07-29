#Import needed library
library(ggplot2)

# Read the file
combinations_counts = read.csv("combinations_counts.csv", head = TRUE)

# Remove NA values
combinations_counts_clean <- combinations_counts[!is.na(combinations_counts$acronym), ]

# Find top 10 gene cancer pairs
top10 <- head(combinations_counts_clean[
  order(combinations_counts_clean$cnv_positive_mut + 
          combinations_counts_clean$cnv_negative_mut, 
        decreasing = TRUE),], 10)

# Transform top_10 into data frame
top_10 <- data.frame(top10)

# Merge the first two columns with a vertical bar separator
top_10$type <- paste(top_10[,1], top_10[,2], sep = " | ")

# Remove the original columns
top_10 <- top_10[, -c(1, 2)]

# Move the merged_columns column to the front
top_10 <- top_10[, c(ncol(top_10), 1:(ncol(top_10)-1))]

# Calculate the proportion of each value in a row
proportions_df <- round(as.data.frame(t(apply(top_10[,2:7], 1, function(row) row / sum(row)))), 4)

# Add type column to the proportion dataframe
rownames(proportions_df) <- top_10$type

# First, convert the row names into a separate column of the dataframe
proportions_df$gene <- rownames(proportions_df)
rownames(proportions_df) <- NULL

# Then, reshape the dataframe into a long format 
proportions_long <- tidyr::gather(proportions_df, key = "category", value = "proportion", -gene)

# Create "all data" plot 
# Create the plot using ggplot
ggp1 <- ggplot(proportions_long, aes(x = gene, y = proportion, fill = category)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "All data", x = "Gene", y = "Proportions") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")) +
  theme_classic() +
  theme(legend.position = "bottom")

ggp1

# Save "all data" plot
# Save the plot as a PNG file
ggsave("data/graphics/All_data.png", plot = ggp1, dpi = 300)

# Extract the data of interest for "Proportions by Gene" plot 
#REMOVE "NOTHING" COLUMN

top_10_only_mc <- top_10[,c(1:4,6:7)]

# Calculate the proportion of each value in a row
proportions_only_mc <- round(as.data.frame(t(apply(top_10_only_mc[,2:6], 1, function(row) row / sum(row)))), 4)

# Add type column to the proportion dataframe
rownames(proportions_only_mc) <- top_10_only_mc$type

# First, convert the row names into a separate column of the dataframe
proportions_only_mc$gene <- rownames(proportions_only_mc)
rownames(proportions_only_mc) <- NULL

# Then, reshape the dataframe into a long format
proportions_only_mc_long <- tidyr::gather(proportions_only_mc, key = "category", value = "proportion", -gene)

# Create the plot using ggplot
ggp2 <- ggplot(proportions_only_mc_long, aes(x = gene, y = proportion, fill = category)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Proportions by Gene", x = "Gene", y = "Proportions") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")) +
  theme_classic() +
  theme(legend.position = "bottom"),
    axis.text.x = element_text(angle = 45, size = 8, hjust = 1))

ggp2

# Save "Proportions by Gene" plot
# Save the plot as a PNG file
ggsave("data/graphics/Nothing_removed_plot.png", plot = ggp2, dpi = 300)

## Extract the data of interest for log transformed "Proportions by Gene" plot 
#LOG TRANSFORM 

top_10_only_mc_log <- top_10_only_mc

top_10_only_mc_log[,c(2:6)] <- log(top_10_only_mc_log[,c(2:6)])

top_10_only_mc_log <- rapply(top_10_only_mc_log, function(x) ifelse(is.infinite(x), 0, x), how = "replace")

# Calculate the proportion of each value in a row
proportions_only_mc_log <- round(as.data.frame(t(apply(top_10_only_mc_log[,2:6], 1, function(row) row / sum(row)))), 4)

# Add type column to the proportion dataframe
rownames(proportions_only_mc_log) <- top_10_only_mc_log$type

# First, convert the row names into a separate column of the dataframe
proportions_only_mc_log$gene <- rownames(proportions_only_mc_log)
rownames(proportions_only_mc_log) <- NULL

# Then, reshape the dataframe into a long format
proportions_only_mc_long_log <- tidyr::gather(proportions_only_mc_log, key = "category", value = "proportion", -gene)

# Create log transformed "Proportions by Gene" plot 
# Create the plot using ggplot
ggp3 <- ggplot(proportions_only_mc_long_log, aes(x = gene, y = proportion, fill = category)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Proportions by Gene", x = "Gene", y = "Proportions") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")) +
  theme_classic() +
  theme(legend.position = "bottom"),
    axis.text.x = element_text(angle = 45, size = 8, hjust = 1))

ggp3

# Save log transformed "Proportions by Gene" plot 
# Save the plot as a PNG file
ggsave("data/graphics/Log_transformed_original.png", plot = ggp3, dpi = 300)

# Extract the data of interest for merged "Proportions by Gene" plot 
# Merge cnv and cnv_mut
top_10_only_mc_log_merged <- top_10_only_mc_log
top_10_only_mc_log_merged$cnv <- top_10_only_mc_log_merged$cnv_positive + top_10_only_mc_log_merged$cnv_negative
top_10_only_mc_log_merged$cnv_mut <- top_10_only_mc_log_merged$cnv_positive_mut + top_10_only_mc_log_merged$cnv_negative_mut
top_10_only_mc_log_merged <- top_10_only_mc_log_merged[,c(1,2,7,8)]

# Calculate the proportion of each value in a row
proportions_only_mc_log_shorten <- round(as.data.frame(t(apply(top_10_only_mc_log_merged[,2:4], 1, function(row) row / sum(row)))), 4)

# Add type column to the proportion dataframe
rownames(proportions_only_mc_log_shorten) <- top_10_only_mc_log_merged$type

# First, convert the row names into a separate column of the dataframe
proportions_only_mc_log_shorten$gene <- rownames(proportions_only_mc_log_shorten)
rownames(proportions_only_mc_log_shorten) <- NULL

# Then, reshape the dataframe into a long format 
proportions_only_mc_long_log_shorten <- tidyr::gather(proportions_only_mc_log_shorten, key = "category", value = "proportion", -gene)

# Reorder the levels of the category variable, so it looks better on the graph
proportions_only_mc_long_log_shorten$category <- factor(
  proportions_only_mc_long_log_shorten$category,
  levels = c("cnv_mut", "cnv", "mut_only")

# Create merged "Proportions by Gene" plot 
# Create the plot using ggplot
ggp4 <- ggplot(proportions_only_mc_long_log_shorten, aes(x = gene, y = proportion, fill = category)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Proportions by Gene", x = "Gene", y = "Proportions") +
  scale_fill_manual(values = c("#E69F00", "#009E73", "#0072B2")) +
  theme_classic() +
  theme(legend.position = "bottom")

ggp4

# Save merged "Proportions by Gene" plot 
# Save the plot as a PNG file
ggsave("data/graphics/Log_transformed_merged.png", plot = ggp4, dpi = 300)

# Volcano plots 

# Load the data
fishers_negative <- read.csv("data/fishersTest/fishers_negative.csv", head = TRUE)

# Create volcano negative plot
volcano_plot_negative <- ggplot(fishers_negative, aes(x = odds, y = p.adj)) +
  geom_point(color = "blue", size = 3) +
  labs(title = "Volcano Plot",
       x = "Odds",
       y = "Adjusted p-value")

# Save the plot
ggsave("data/graphics/volcano_fishers_negative.png", plot = volcano_plot_negative, dpi = 300)

# Create volcano positive plot

# Load the data
fishers_positive <- read.csv("data/fishersTest/fishers_positive.csv", head = TRUE)

# Create volcano positive plot
volcano_plot_positive <- ggplot(fishers_positive, aes(x = odds, y = p.adj)) +
  geom_point(color = "blue", size = 3) +
  labs(title = "Volcano Plot",
       x = "Odds",
       y = "Adjusted p-value")

# Save the plot
ggsave("data/graphics/volcano_fishers_positive.png", plot = volcano_plot_positive, dpi = 300)
