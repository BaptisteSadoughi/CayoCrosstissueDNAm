# ==============================================================================
# Plot top tissue-specific markers
# ------------------------------------------------------------------------------
# This script loads top tissue-specific methylation markers, 
# extracts methylation values from imputed data, 
# and generates a heatmap for visualization.
# ==============================================================================

rm(list = ls())

library_list <- c("tidyverse","dplyr","tidyr","purrr","parallel","stringr","MatrixGenerics","svglite","gplots","RColorBrewer","pheatmap")
lapply(library_list, require, character.only = TRUE)

# === Recode tissue names for plotting ===
recode_map <- c(
  "omental_at" = "omental adipose",
  "skeletal_muscle" = "skeletal muscle",
  "whole_blood" = "whole blood"
)

# === Paths ===
base_path <- "/path/to/your/directory" # <-- set the path once

figure_path <- file.path(base_path, "Figures")
tissuemarker_path <- file.path(base_path, "tissue_comparisons", "tissuespecific_methylation.txt")
imputed_pmeth <- file.path(base_path, "imputed", "Regions_imputed_methyLimp_pmeth14.rds")

# === Tissues ===
tissue_oi <- c("whole_blood","spleen","omental_at","heart","testis","ovaries","kidney","lung","adrenal","thymus","thyroid","pituitary","liver","skeletal_muscle")

# === Plot palette ===
tissue_plot <- sort(c("whole blood","spleen","omental adipose","heart","kidney","lung","adrenal","thymus","thyroid", "pituitary", "liver", "skeletal muscle"))
extended_palette <- setNames(colorRampPalette(brewer.pal(8, "Dark2"))(12),tissue_plot)

# -----------------------------------
# === LOAD DATA ===
# -----------------------------------

# === Load tissue marker ===
tissue_markers = read.table(tissuemarker_path, sep="\t", header = TRUE)
tissue_markers$marker <- ifelse(tissue_markers$mean_beta>0, "hyper", "hypo")

# === Select top20 markers per tissue ===
top_markers <- tissue_markers %>% group_by(sig_tissue) %>% slice_max(abs(mean_beta), n=20, with_ties = TRUE) %>% ungroup() %>%
  mutate(tissue = recode(sig_tissue, !!!recode_map))

# === Load and format imputed percent methylation data ===
files <- list.files(path = file.path(base_path, "imputed"), pattern = "methyLimp", full.names = TRUE)

pmeth_imp <- list()
for(file_n in files){
   name <- basename(file_n)
   tissue_name <- gsub("pmeth_imp_methyLimp_","",name)
   tissue_name <- gsub(".rds","",tissue_name)
   pmeth_imp[[tissue_name]] <- readRDS(file = file_n)
 }

# Exclude gonadal tissues
pmeth_imp = pmeth_imp[!names(pmeth_imp) %in% c("ovaries","testis")]

# === Extract top marker values from each tissue matrix ===
top_mats <- lapply(pmeth_imp, function(x){
  x[rownames(x) %in% top_markers$sites,]
})

names(top_mats) <- recode(
  names(top_mats),
  !!!recode_map
)

top_mats <- top_mats[order(names(top_mats))]

# === Combine into one matrix ===
top_mats_df <- map2(top_mats, names(top_mats), function(mat, name) {
  colnames(mat) <- paste(name, colnames(mat), sep = "_")
  mat %>%
  as_tibble(rownames = "loci") 
})

top_mat_combined <- purrr::reduce(top_mats_df, full_join, by = "loci") %>%
  mutate(across(-loci, as.numeric)) %>% 
  column_to_rownames("loci") %>%     
  as.matrix()

# -----------------------------------
# === PREPARE ANNOTATIONS FOR PLOTTING ===
# -----------------------------------

# === Row ordering ===
rowMeans_meth <- rowMeans(top_mat_combined, na.rm = TRUE)

top_markersnoX <- top_markers %>% filter(!grepl("Region_X_*",sites)) 
loci_to_tissue_focal <- top_markersnoX[, c("sites", "tissue", "marker")] %>%
  mutate(avg_meth = rowMeans_meth[match(sites, rownames(top_mat_combined))]) %>%
  arrange(tissue, -avg_meth)

ordered_rows <- loci_to_tissue_focal$sites
top_mat_combined_ordered <- top_mat_combined[ordered_rows, ]

# === Column ordering ===
tissue_levels <- unique(sub("_.*", "", colnames(top_mat_combined_ordered)))
sorted_col_indices <- unlist(sapply(tissue_levels, function(tl) which(startsWith(colnames(top_mat_combined_ordered), tl))))
top_mat_combined_ordered <- top_mat_combined_ordered[, sorted_col_indices]

# === Column (sample) annotations ===
tissue_annotations <- data.frame(
  Tissue = sub("_.*", "", colnames(top_mat_combined_ordered)),
  row.names = colnames(top_mat_combined_ordered)
)

# === Row (locus) annotations ===
row_annotations <- top_markers %>%
  dplyr::select(sites, sig_tissue, marker) %>%
  dplyr::rename(Loci = sites) %>%
  filter(Loci %in% rownames(top_mat_combined_ordered)) %>%
  arrange(match(Loci, rownames(top_mat_combined_ordered))) %>%
  column_to_rownames("Loci") %>% 
  mutate(sig_tissue=recode(sig_tissue,!!!recode_map)) %>% 
  rename(Marker=marker)

# === Fig. 1E ===
Cairo::CairoPNG(filename = file.path(figure_path, "Fig1E.png"), bg = "transparent", width = 6.5, height = 5, units="in")
                
pheatmap::pheatmap(top_mat_combined_ordered,
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   clustering_distance_cols = "euclidean",
                   clustering_method = "complete",
                   show_rownames = FALSE,
                   show_colnames = FALSE,
                   scale = "none",
                   annotation_col = tissue_annotations,
                   annotation_row = row_annotations,
                   annotation_names_row = FALSE,
                   annotation_names_col = FALSE,
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
