#==========================================================================
#
# Unsupervised hierarchical tree clustering of samples
#
#==========================================================================

# User Configuration - Modify these paths as needed
base_path <- "YOUR_BASE_PATH_HERE"  # Set your base path here

output_path <- file.path(base_path, "Figures")
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

library_list <- c("RColorBrewer","svglite","tidyverse","rlang", "corrplot", "dendextend")
lapply(library_list, require, character.only=TRUE)

# Define ggplot theme upfront
tissue_plot <- sort(c("whole_blood","spleen","omental_at","heart","kidney","lung","adrenal","thymus","thyroid", "pituitary", "liver", "skeletal_muscle"))
extended_palette <- setNames(colorRampPalette(brewer.pal(8, "Dark2"))(12),tissue_plot)
# add missing tissues
new_levels <- c("ovaries", "testis")
# Define the new colors:
new_colors <- c("#008B8B", "#4682B4")
# Add these new levels and colors to your plots and palette
tissue_plot <- c(tissue_plot, new_levels)
extended_palette <- c(extended_palette, setNames(new_colors, new_levels))

# Read in data
file_name <- file.path(base_path, "imputed", "pmeth_imp_methyLimp_14.rds")
pmeth <- readRDS(file_name)

# Convert matrices to data frames for merging and create a new column with row names
list_of_dfs <- lapply(pmeth, function(x) {
  df <- as.data.frame(x)
  df$rownames <- row.names(df)
  return(df)
})

# Define function for joining two data frames
join_two_dfs <- function(df1, df2) {
  full_join(df1, df2, by = "rownames")  # Now there's a column named 'row.names'
}

# Reduce list of data frames to a single data frame by joining
merged_df <- Reduce(join_two_dfs, list_of_dfs)

# Keep only rows with no missing values
pmeth <- merged_df[complete.cases(merged_df),]
rownames(pmeth) <- pmeth$rownames
pmeth <- pmeth %>% select(-rownames)
pmeth <- t(as.matrix(pmeth)) #work with rows = samples

# Load metadata
metadata_path <- file.path(base_path, "metadata", "multitissue_metadata.txt")
metadata <- read.delim(metadata_path, stringsAsFactors = FALSE) %>%
  filter(lid_pid != "LID_109490_PID_10416")

pmeth <- pmeth[rownames(pmeth) != "LID_109490_PID_10416",]
metadata <- metadata[match(rownames(pmeth), metadata$lid_pid),]

set.seed(123)

# Compute Euclidean distance between samples (columns)
dist_matrix <- dist(pmeth, method = "euclidean")

# Perform hierarchical clustering using complete linkage
hc <- hclust(dist_matrix, method = "ward.D2") #complete

# convert to dendrogram
dend <- as.dendrogram(hc)

# Prepare annotation colors

# tissue type
tissue_fac <- factor(metadata$grantparent_tissueType)
tissue_cols <- extended_palette[metadata$grantparent_tissueType]

# reorder colors to match dendrogram leaf order
leaf_order <- labels(dend)
annot_cols <- tissue_cols[match(leaf_order, metadata$lid_pid)]

# Add colored squares ABOVE the leaves
dend2 <- dend %>%
  set("labels", rep("", length(leaf_order))) %>%      # remove sample names
  set("leaves_pch", 15) %>%                           # square = 15, circle = 19
  set("leaves_cex", 1.5) %>%                          # square size
  set("leaves_col", annot_cols)                       # square colors

# Save plot
svglite::svglite(file.path(output_path, "tissue_hierarchical_tree_ward.svg"),
                 width = 10, height = 6)

# Plot
plot(dend2, main = "")
dev.off()
