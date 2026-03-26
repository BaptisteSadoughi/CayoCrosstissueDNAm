##############################
# MULTI-TISSUE AGING ANALYSIS
##############################

rm(list = ls())

##############################
# LIBRARIES
##############################

library_list <- c("corrplot","dplyr","purrr","parallel","tidyverse", "DHARMa", "performance","mice","useful","NbClust","factoextra","RColorBrewer","ComplexHeatmap","reshape2","Cairo","circlize")

lapply(library_list, require, character.only=TRUE)

##############################
# GLOBAL SETTINGS
##############################

tissue_plot <- sort(c("whole blood","spleen","omental adipose","heart","kidney","lung","adrenal","thymus","thyroid", "pituitary", "liver", "skeletal muscle"))
extended_palette <- setNames(colorRampPalette(brewer.pal(8, "Dark2"))(12),tissue_plot)
new_levels <- c("ovaries", "testis")
new_colors <- c("#008B8B", "#4682B4")
tissue_plot <- c(tissue_plot, new_levels)
extended_palette <- c(extended_palette, setNames(new_colors, new_levels))

# === Paths ===

base_path   <- "/path/to/project" # <-- Define this path only once
metadata_path <- file.path(base_path, "metadata", "multitissue_metadata.txt")
age_deviations_path <- file.path(base_path,"DNAm_deviation_data.txt")
output_path <- file.path(base_path, "Figures")

##############################
# LOAD DATA
##############################

metadata = read.table(metadata_path, sep = "\t", header = TRUE) %>% filter(lid_pid != "LID_109490_PID_10416") %>% mutate(percent_unique = unique/reads)

combined_data <- read.table(age_deviations_path, sep="\t", header=TRUE) %>%
  filter(monkey_id != "22H" & !is.na(Residual_adult)) %>% dplyr::select(-Residual) %>%
  mutate(tissue = recode(tissue, "omental fat" = "omental adipose", "blood" = "whole blood")) %>%
  filter(!tissue %in%  c("testis","ovaries"))

##############################
# MATRIX FORMAT
##############################
combined_data_mat = dcast(combined_data, monkey_id ~ tissue, value.var = "Residual_adult")
rownames(combined_data_mat) = combined_data_mat$monkey_id
combined_data_mat = combined_data_mat[,-1]

# Replace spaces for mice
names(combined_data_mat) <- gsub(" ", "_", names(combined_data_mat))

##############################
# IMPUTATION (MICE)
##############################
set.seed(3500)
imp <- mice::mice(combined_data_mat, m = 5, method = "pmm", maxit = 35, seed = 500)

summary(imp)
mice::densityplot(imp)
mice::stripplot(imp)
plot(imp)
completed_data <- complete(imp, 1)

# Restore names
names(completed_data) <- gsub("_", " ", names(completed_data))
completed_data$monkey_id <- rownames(completed_data)

##############################
# FILTER BASED ON MISSINGNESS
##############################
# Drop monkeys for which more than 1/3 (i.e.4) of the tissues were imputed
missing_per_monkey = apply(combined_data_mat,1,function(x) sum(is.na(x)))
completed_data_f = completed_data %>% filter(monkey_id %in% names(missing_per_monkey)[missing_per_monkey<=4])


##############################
# CORRELATION ANALYSIS
##############################
#### Pair-wise correlations of tissue-specific age deviations (the focus is on the tissues)

cor_matrix <- cor(select(completed_data_f, -monkey_id))
          
# --- Correlation heatmap
Cairo::CairoPNG(file.path(output_path, "Fig3E.png"))
corrplot::corrplot(cor_matrix,type = "lower",method="color", diag = F,tl.col="black",
                   order = 'hclust', hclust.method = "complete",addCoef.col = "black",
                   tl.cex = 1.2,tl.srt=45,number.cex=0.9,cl.cex=1.5, cl.length = 5)
dev.off()

# --- Histogram of correlations
cor_matrix2 <- cor_matrix[upper.tri(cor_matrix)]
ggplot(as.data.frame(cor_matrix2), aes(x=cor_matrix2))+
  geom_histogram(bins = 10, fill="lightblue",col="black", boundary=0)+
  geom_vline(xintercept = mean(cor_matrix2), linetype="dashed", linewidth = 0.9)+
  labs(x="Pearson r")+
  theme_classic()+
  theme(axis.text = element_text(size=16, color="black"),
        axis.title = element_text(size=16, color="black"))
ggsave(file.path(output_path,"Fig3F.png"), width = 3.01, height = 2.28)
                       
##############################
# HEATMAP CLUSTERING
##############################

cluster_data <- select(completed_data_f, -monkey_id)

hc_rows <- hclust(dist(cluster_data, method = 'euclidean'), method = 'complete')
hc_cols <- hclust(dist(t(cluster_data), method = 'euclidean'), method = 'complete')

# Annotation (sex)
annotation_df <- data.frame(
  sex = metadata$individual_sex[
    match(rownames(cluster_data), metadata$monkey_id)
  ]
)
rownames(annotation_df) <- rownames(cluster_data)

sex_colors <- c("M" = "purple", "F" = "limegreen")

heat_colors <- circlize::colorRamp2(
  c(-5, 5),  # or other values if you want
  c("yellow", "red")
)

Cairo::CairoPNG(file.path(output_path,"Fig3G.png"), width = 8, height = 7, units="in",dpi=300)
Heatmap(
  as.matrix(cluster_data),
  name = "age deviation",
  col = heat_colors,
  cluster_rows = as.dendrogram(hc_rows),
  cluster_columns = as.dendrogram(hc_cols),
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_rot = 45,
  column_names_gp = grid::gpar(fontsize = 18, fontface = "plain"),
  row_names_side = "right",
  heatmap_legend_param = list(
    title_gp = grid::gpar(fontsize = 14, fontface = "plain"),
    labels_gp = grid::gpar(fontsize=14)
  ),
  right_annotation = rowAnnotation(
    df = annotation_df,
    col = list(sex = sex_colors),
    show_annotation_name=FALSE,
    annotation_legend_param = list(
      title_gp = grid::gpar(fontsize = 14, fontface = "plain"),
      labels_gp = grid::gpar(fontsize=14)
    )
  )
)
dev.off()

##############################
# CLUSTERING (AGEOTYPES)
##############################

# Elbow graph on within-cluster sum of square
factoextra::fviz_nbclust(cluster_data, kmeans, method = "wss", k.max = 12) + theme_minimal() + ggtitle("Elbow Method")
factoextra::fviz_nbclust(cluster_data, hcut, method = "wss", k.max = 12) + theme_minimal() + ggtitle("Elbow Method")

# Silhouette plot
factoextra::fviz_nbclust(cluster_data, kmeans, method = "silhouette", k.max = 12) + theme_minimal() + ggtitle("Elbow Method")
factoextra::fviz_nbclust(cluster_data, hcut, method = "silhouette", k.max = 12) + theme_minimal() + ggtitle("Elbow Method")
                       
res.eclust_kmean <- factoextra::eclust(cluster_data,
                             FUNcluster = "kmeans",
                             k = 2,
                             stand = FALSE,
                             graph = TRUE,
                             gap_maxSE = list(method = "firstSEmax", SE.factor = 1),
                             nboot = 100,
                             verbose = interactive())

factoextra::fviz_cluster(res.eclust_kmean, data = cluster_data) + theme_bw() + ggtitle("k = 2")
factoextra::fviz_silhouette(res.eclust_kmean)+theme_bw()

# --> Poor accuracy of group assignment.
