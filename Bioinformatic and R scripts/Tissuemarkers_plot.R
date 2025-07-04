rm(list = ls())

library_list <- c("tidyverse","dplyr","tidyr","purrr,"parallel","stringr","MatrixGenerics","svglite","gplots","RColorBrewer","pheatmap")
lapply(library_list, require, character.only = TRUE)

# Prepare renaming
recode_map <- c(
  "omental_at" = "omental adipose",
  "skeletal_muscle" = "skeletal muscle",
  "whole_blood" = "whole blood"
)

tissue_oi <- c("whole_blood","spleen","omental_at","heart","testis","ovaries","kidney","lung","adrenal","thymus","thyroid","pituitary","liver","skeletal_muscle")

# Define ggplot theme upfront
tissue_plot <- sort(c("whole blood","spleen","omental adipose","heart","kidney","lung","adrenal","thymus","thyroid", "pituitary", "liver", "skeletal muscle"))
extended_palette <- setNames(colorRampPalette(brewer.pal(8, "Dark2"))(12),tissue_plot)

# Load tissue markers
tissue_markers = read.table("/path/to/tissuespecific_methylation.txt", sep="\t", header = TRUE)

tissue_markers$marker <- ifelse(tissue_markers$mean_beta>0,"hyper","hypo")

# Select top20 markers per tissue
top_markers <- tissue_markers %>% group_by(sig_tissue) %>% slice_max(abs(mean_beta), n=20, with_ties = TRUE) %>% ungroup()

top_markers <- top_markers %>%
  mutate(tissue = recode(sig_tissue, !!!recode_map))

# Load imputed percent methylation data
pmeth_imp <- readRDS("/path/to/Regions_imputed_methyLimp_pmeth14.rds")

pmeth_imp = pmeth_imp[!names(pmeth_imp) %in% c("ovaries","testis")]

# Create a list of data frames with samples percent meth
top_mats <- lapply(pmeth_imp, function(x){
  x[rownames(x) %in% top_markers$sites,]
})

names(top_mats) <- recode(
  names(top_mats),
  !!!recode_map
)

top_mats <- top_mats[order(names(top_mats))]

# Iterate over the matrices in the list, append list element name to column names
# then convert each matrix to a tibble with rownames as a new column
top_mats_df <- map2(top_mats, names(top_mats),
                    function(mat, name) {
                      colnames(mat) <- paste(name, colnames(mat), sep = "_")
                      mat %>%
                        as_tibble(rownames = "loci")    # Convert each matrix to a df with rownames as a column
                    })

# Combine all tibbles into one, filling missing column combinations with NA
top_mat_combined <- purrr::reduce(top_mats_df, full_join, by = "loci")

# Convert back to a matrix everything except the 'loci' column to numeric
top_mat_combined <- top_mat_combined %>%
  mutate(across(-loci, as.numeric)) %>% 
  column_to_rownames("loci") %>%     
  as.matrix()

### Calculate the average percent methylation for each locus
rowMeans_meth <- rowMeans(top_mat_combined, na.rm = TRUE)

### Combine with the `tissue_focal` and `marker` information
top_markersnoX <- top_markers %>% filter(!grepl("Region_X_*",sites)) 
loci_to_tissue_focal <- top_markersnoX[, c("sites", "tissue", "marker")] %>%
  mutate(avg_meth = rowMeans_meth[match(sites, rownames(top_mat_combined))])

### Order the rows by `tissue_focal` and decreasing average methylation within those levels
loci_to_tissue_focal <- loci_to_tissue_focal %>%
  arrange(tissue, -avg_meth)

### Get the ordered row names
ordered_rows <- loci_to_tissue_focal$sites

### Reorder the rows of `top_mat_combined`
top_mat_combined_ordered <- top_mat_combined[ordered_rows, ]

### Reorder the columns based on tissue levels (already done but included for completeness)
tissue_levels <- unique(sub("_.*", "", colnames(top_mat_combined_ordered)))
sorted_col_indices <- unlist(sapply(tissue_levels, function(tl) which(startsWith(colnames(top_mat_combined_ordered), tl))))
top_mat_combined_ordered <- top_mat_combined_ordered[, sorted_col_indices]

# Create a data frame for the column annotation based on tissue identity
tissue_annotations <- data.frame(
  Tissue = sub("_.*", "", colnames(top_mat_combined_ordered)),
  row.names = colnames(top_mat_combined_ordered)
)

# Make sure the row names match the column names of the matrix
rownames(tissue_annotations) <- colnames(top_mat_combined_ordered)

# Create a data frame for row annotations
row_annotations <- top_markers %>%
  dplyr::select(sites, sig_tissue, marker) %>%
  dplyr::rename(Loci = sites) %>%
  filter(Loci %in% rownames(top_mat_combined_ordered)) %>%
  arrange(match(Loci, rownames(top_mat_combined_ordered))) %>%
  column_to_rownames("Loci") %>% 
  mutate(sig_tissue=recode(sig_tissue,!!!recode_map)) %>% 
  rename(Marker=marker)

# Fig. 1E
Cairo::CairoPNG(filename = "/path/to/Figures/Fig1.png",
                bg = "transparent", width = 6.5, height = 5, units="in")
pheatmap::pheatmap(top_mat_combined_ordered,
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   # treeheight_col = 0,
                   clustering_distance_cols = "euclidean",
                   clustering_method = "complete",
                   show_rownames = FALSE,
                   show_colnames = FALSE,
                   scale = "none",
                   # main = "Heatmap of Tissue x Tissue_Focal",
                   annotation_col = tissue_annotations,
                   annotation_row = row_annotations,
                   annotation_names_row = FALSE,
                   annotation_names_col = FALSE,
                   # legend = FALSE,
                   annotation_colors = list(
                     Tissue=extended_palette,
                     sig_tissue=extended_palette,
                     Marker = c(
                       "hyper" = "darkred",
                       "hypo" = "lightblue"
                     )
                   ),
                   annotation_legend = F #mute to show legend
)
dev.off()
