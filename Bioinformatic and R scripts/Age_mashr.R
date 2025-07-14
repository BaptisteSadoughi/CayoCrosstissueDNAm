# ==============================================================================
# Mashr Correction of Age-associated Effects
# ==============================================================================

rm(list=ls())

# === Libraries ===
required_libraries <- c("ashr", "mashr", "flashier", "svglite", "tidyverse", "corrplot", "MatrixGenerics","ggpubr","ComplexUpset","scales")
lapply(required_libraries, require, character.only = TRUE)

# === Paths ===
base_path <- "/path/to/project"  # <-- Define this path only once
pqlseq_path <- file.path(base_path, "PQLSEQ")
metadata_path <- file.path(base_path,"metadata", "multitissue_metadata.txt")
mash_path <- file.path(base_path, "MASH")
bed_path <- file.path(mash_path, "bedfiles")
figure_path <- file.path(base_path, "Figures")

if (!dir.exists(paths = mash_path)) dir.create(mash_path, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(paths = bed_path)) dir.create(bed_path, recursive = TRUE, showWarnings = FALSE)

# === Tissues of interest ===
tissue_oi <- c("whole_blood", "spleen", "omental_at", "heart", "testis", "ovaries", "kidney", "lung", "adrenal", "thymus", "thyroid", "pituitary", "liver", "skeletal_muscle")

# --- Significance threshold ---
lfsr_threshold <- 0.05

vars_to_keep <- c("base_path", "pqlseq_path", "metadata_path", "mash_path", "bed_path", "lfsr_threshold", "tissue_oi")

# -----------------------------------
# === CONSTRUCT MATRICES FOR MASHr ===
# -----------------------------------

# === Load PQLSEQ data into environment ===

read_table_and_rownames <- function(tissue) {
  pqlseq <- read.table(file.path(pqlseq_path, paste0("bed_", tissue, ".txt")), sep = "\t", header = TRUE)
  rownames(pqlseq) <- pqlseq$sites
  return(pqlseq)
}

data_list <- lapply(tissue_oi, read_table_and_rownames)
names(data_list) <- paste("pqlseq", tissue_oi, sep = "_")
list2env(data_list, envir = .GlobalEnv)

# === Build matrices ===

fill_matrix <- function(column) {
  # get all rownames from all dataframes
  all_rownames <- unique(unlist(lapply(tissue_oi, function(x) row.names(get(paste("pqlseq", x, sep = "_"))))))

  # initialize matrix with NA values
  matrix_data <- matrix(NA, nrow = length(all_rownames), ncol = length(tissue_oi))
  rownames(matrix_data) <- all_rownames
  colnames(matrix_data) <- tissue_oi

  # loop through each dataframe and fill matrix_data
  for (tissue in tissue_oi) {
    df <- get(paste("pqlseq", tissue, sep = "_"))
    matrix_data[rownames(df), tissue] <- df[[column]]
  }
  return(matrix_data)
}

matrix_beta <- fill_matrix("beta")
matrix_SE <- fill_matrix("se_beta")

saveRDS(list(matrix_beta=matrix_beta, matrix_SE = matrix_SE), file.path(pqlseq_path, "pre_MASH_effects.RDS"))

# -----------------------------------
# === MASHr ===
# -----------------------------------

rm(list = setdiff(ls(), vars_to_keep))

# Follows the flow of tutorial vignettes available at https://stephenslab.github.io/mashr/articles/index.html

pqlseq_effects <- readRDS(file.path(pqlseq_path, "pre_MASH_effects.RDS"))
pqlseq_effects$matrix_beta <- subset(pqlseq_effects$matrix_beta, select = -c(testis, ovaries))
pqlseq_effects$matrix_SE <- subset(pqlseq_effects$matrix_SE, select = -c(testis, ovaries))
pqlseq_effects <- lapply(pqlseq_effects, function(x) x[complete.cases(x), ])
  
data = mash_set_data(pqlseq_effects$matrix_beta, pqlseq_effects$matrix_SE)

# === Strong set ===
                         
tissue_files <- list.files(pqlseq_path, pattern = ".txt", full.names = TRUE)
strong_subset_list <- list()
for(file in tissue_files) {
  bed_data <- read.table(file, header = TRUE)
  filtered_data <- bed_data[bed_data$qval < quantile(bed_data$qval, 0.01),] # Filter rows for qval < 0.001
  tissue_name <- gsub("\\.txt$", "", basename(file))
  strong_subset_list[[tissue_name]] <- filtered_data$site
}

strong_subset_list <- strong_subset_list[!names(strong_subset_list) %in% c("bed_testis", "bed_ovaries")]
strong_subset <- Reduce(union, strong_subset_list)
strong_subset <- intersect(rownames(pqlseq_effects$matrix_beta),strong_subset)

# === Estimate background correlation ===

Vhat = estimate_null_correlation_simple(data)

# === Update mash objects with correlation structure ===
                         
data.cor = mash_update_data(data, V=Vhat)
data.strong.cor = mash_set_data(pqlseq_effects$matrix_beta[strong_subset,], pqlseq_effects$matrix_SE[strong_subset,], V=Vhat)

# === Investigate data-driven covariances ===

U.pca = cov_pca(data.strong.cor, 5)
U.f = cov_flash(data.strong.cor)
U.ed = cov_ed(data.strong.cor, c(U.pca, U.f))
U.c = cov_canonical(data.cor)

# === Estimate mixture proportions ===

m.all = mash(data.cor, Ulist = c(U.ed, U.c))
saveRDS(m.all, file.path(mash_path, "mash_object.RDS"))

mash_results <- with(m.all$result, list(beta = PosteriorMean, SD = PosteriorSD, LFSR = lfsr))
saveRDS(mash_results, file.path(mash_path, "mash_estimates.RDS"))

# -----------------------------------
# === DIRECTION OF CHANGE ANALYSIS ===
# -----------------------------------

# Load metadata
metadata <- read.table(metadata_path, sep = "\t", header = TRUE) %>%
  filter(lid_pid != "LID_109490_PID_10416") %>%
  mutate(percent_unique = unique / reads) %>%
  filter(age_at_sampling>2.9) # remove infants

# --- Load tissues as list ---
load_tissue_data <- function(tissue) {
  message("Loading tissue data for: ", tissue)
  file_path <- gsub("XXX", tissue, file.path(base_path, "tissues_meth", "XXX_meth", "Regions_pmeth_full_XXX_1000_14T.rds"))
  r <- readRDS(file_path)
  sample_ids <- intersect(colnames(r$coverage), metadata$lid_pid)
  r$metadata <- metadata %>% filter(lid_pid %in% sample_ids)
  
  # Optional: remove flawed sample for skeletal muscle
  if (tissue == "skeletal_muscle") {
    idx <- which(colnames(r$coverage) == "LID_109490_PID_10416")
    if (length(idx) > 0) {
      r$coverage <- r$coverage[, -idx]
      r$methylation <- r$methylation[, -idx]
      r$pmeth <- r$pmeth[, -idx]
    }
  }
  return(r)
}

# Load tissues into a named list
tissue_data_list <- setNames(lapply(tissue_oi, load_tissue_data), tissue_oi)

# --- Extract LFSR and PM ---

lfsr=as.data.frame(get_lfsr(m.all))
pm=as.data.frame(get_pm(m.all))
lfsr$region<-rownames(lfsr)

# --- Intersect with methylation levels ---
                         
coeff_intercept_list <- list()

for (tissue in setdiff(tissue_oi, c("testis", "ovaries"))) {
  sig_regions <- lfsr %>% filter(!!sym(tissue) < lfsr_threshold) %>% pull(region)
  tissue_betas <- pm[sig_regions, tissue, drop = FALSE]
  
  if (length(sig_regions) == 0) next  # skip empty
  
  tissue_beta <- tibble(
    region = sig_regions,
    beta = tissue_betas[[1]],
    sign_beta = ifelse(tissue_betas[[1]] > 0, "pos", "neg"),
    tissue = tissue
  )
  
  tissue_pmeth <- tissue_data_list[[tissue]]$pmeth
  sample_ids <- intersect(colnames(tissue_pmeth), metadata$lid_pid)
  tissue_pmeth <- tissue_pmeth[, sample_ids]
  
  pop_means <- rowMeans2(tissue_pmeth[tissue_beta$region, ], na.rm = TRUE)
  young_ids <- metadata %>%
    filter(tissue == tissue, age <= 6) %>%
    pull(lid_pid)
  
  young_means <- rowMeans2(tissue_pmeth[tissue_beta$region, young_ids], na.rm = TRUE)
  
  # Sanity check
  stopifnot(length(young_means) == nrow(tissue_beta), length(pop_means) == nrow(tissue_beta))
  
  tissue_beta <- tissue_beta %>%
    mutate(
      young_intercept = young_means,
      pop_intercept = pop_means,
      intercept_high_low_50 = ifelse(young_intercept < 0.5, "lower", "higher"),
      nonvariable = ifelse(pop_intercept < 0.1 | pop_intercept > 0.9, "nonvar", "var")
    )
  
  coeff_intercept_list[[tissue]] <- tissue_beta
}

# Save results
saveRDS(coeff_intercept_list, file.path(bed_path,"tissue_age_associated_sites_list.rds"))

# --- Pearson correlation (variable sites only) ---
apply_cor_test <- function(df) {
  df <- df %>% filter(nonvariable == "var")
  test <- cor.test(df$beta, df$young_intercept, method = "pearson")
  tibble(pearson_corr = round(test$estimate, 3), pvalue = ifelse(test$p.value < 0.001, "<0.001", test$p.value))
}

# --- Table S5 ---
results_corr_df_var <- map_dfr(coeff_intercept_list, apply_cor_test, .id = "tissue") %>%
  mutate(
    tissue = recode(tissue,
                    "omental_at" = "omental adipose",
                    "skeletal_muscle" = "skeletal muscle",
                    "whole_blood" = "whole blood"
    )
  ) %>%
  arrange(tissue)

# --- Fig S3 ---
all_coeff <- bind_rows(coeff_intercept_list)

cat_young_plot <- all_coeff %>%
  mutate(
    interaction = recode(paste(sign_beta, intercept_high_low_50),
                         "neg higher" = "high decreasing",
                         "neg lower"  = "low decreasing",
                         "pos higher" = "high increasing",
                         "pos lower"  = "low increasing"),
    tissue = recode(tissue,
                    "omental_at" = "omental adipose",
                    "skeletal_muscle" = "skeletal muscle",
                    "whole_blood" = "whole blood")
  ) %>%
  ggplot(aes(x = interaction, fill = interaction)) +
  geom_bar() +
  labs(x = "", title = "A") +
  facet_wrap(~tissue, scales = "free_y") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", plot.title = element_text(size = 28),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.text = element_text(face = "bold"))

corr_plot <- all_coeff %>%
  mutate(tissue = recode(tissue,
                         "omental_at" = "omental adipose",
                         "skeletal_muscle" = "skeletal muscle",
                         "whole_blood" = "whole blood")) %>%
  ggplot(aes(x = young_intercept, y = beta)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  labs(x = "Avg. methylation in young subjects", y = "Age-associated effect size", title = "B") +
  facet_wrap(~tissue, scales = "free_y") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", plot.title = element_text(size = 28))

# Combine and save plot
ggpubr::ggarrange(cat_young_plot, corr_plot, nrow = 2)
ggsave(file.path(figure_path,"FigS3.png"), width=9.5, height=12.5, dpi=300)

# -----------------------------------
# === TISSUE SHARING OF AGE EFFECTS ===
# -----------------------------------

rm(list = setdiff(ls(), vars_to_keep))

mash_results <- readRDS(file.path(mash_path, "mash_estimates.RDS"))
m.all <- readRDS(file.path(mash_path, "mash_object.RDS"))

# Extract lfsr and posterior means
lfsr=as.data.frame(get_lfsr(m.all))
pm=as.data.frame(get_pm(m.all))
lfsr$id<-1:dim(lfsr)[1]
pm$id<-1:dim(pm)[1]

 # Count number of tissues (columns) below threshold for each row
lfsr$sig_tissues <- apply(lfsr[, !(names(lfsr) == "id")] < lfsr_threshold, 1, sum)

# paste the names of the tissues (columns) below threshold for each row
lfsr <- lfsr %>%
  rowwise() %>%
  mutate(sig_which_tissues = paste(names(.)[1:12][c_across(1:12) < lfsr_threshold], collapse = ",")) %>%
  ungroup() %>%
  as.data.frame()

rownames(lfsr) <- rownames(pm)

# === Fig 2F ===          
x = get_pairwise_sharing(m.all, factor=0.5, lfsr_thresh = 0.05)

# Rename tissues for clarity
tissue_names <- c("whole_blood" = "whole blood", 
                  "omental_at" = "omental adipose", 
                  "skeletal_muscle" = "skeletal muscle")

rownames(x) <- colnames(x) <- recode(rownames(x), !!!tissue_names)
                 
Cairo::CairoPNG(file.path(figure_path,"Fig2F.png"))
corrplot(x, method='color', col.lim=c(0,1), type='lower', addCoef.col = "black", tl.col="black", tl.srt=45,number.cex=0.65,cl.cex=1.3, cl.length = 5,tl.cex = 1.3)
dev.off()

# === Fig 2A === 
beta_res <- mash_results$beta
lfsr_res <- mash_results$LFSR

beta_res <- beta_res %>%
  as_tibble(rownames = "Regions") %>%
  pivot_longer(-Regions, names_to = "tissue", values_to = "beta") %>%
  mutate(tissue = recode(tissue, !!!tissue_names))

lfsr_res <- lfsr_res %>%
  as_tibble(rownames = "Regions") %>%
  pivot_longer(-Regions, names_to = "tissue", values_to = "lfsr") %>%
  mutate(tissue = recode(tissue, !!!tissue_names))

lfsr_sig <- lfsr_res %>% filter(lfsr < lfsr_threshold)
beta_sig <- beta_res %>% semi_join(lfsr_sig, by=c("Regions", "tissue"))

summary_count <- beta_sig %>%
                 group_by(tissue) %>%
                 summarize(pos = sum(beta>0), neg = -sum(beta<0)) %>%
                 tibble() %>%
                 pivot_longer(cols = -tissue, names_to = "sign", values_to = "count")

ggplot(summary_count, aes(x=tissue, y=count, fill = tissue, alpha=sign))+
  geom_col()+
  labs(y="Number of differentially methylated regions")+
  scale_y_continuous(labels = label_number(scale = 1e-3, suffix = "k"))+
  scale_fill_manual(values = extended_palette)+
  scale_alpha_discrete(range = c(0.7,1))+
  theme_bw()+theme(
    axis.text.x=element_text(angle = 45, hjust=1, vjust = 1, color="black",size=17),
    axis.title.x = element_blank(),
    axis.title.y = element_text(color="black",size=17, hjust = 1),
    axis.text.y = element_text(color="black",size=17),
    legend.position = "n")
ggsave(file.path(figure_path,"Fig2A.png"))

# === Investigating shared vs. tissue-specific DMRs ===

lfsr_sig1 <- subset(lfsr,sig_tissues>0)
lfsr_sig1$within2_v2 <- NA
lfsr_sig1$focal_pm <- NA
lfsr_sig1$other_pm <- NA

pm_sig1<-subset(pm, id %in% lfsr_sig1$id)
rownames_sig1=rownames(lfsr_sig1)
                 
for (i in 1:dim(lfsr_sig1)[1]){
  
  z <- min(lfsr_sig1[i, tissue_oi])[1]
  focal <- pm[lfsr_sig1$id[i], which(lfsr_sig1[i, tissue_oi] == z)][1]
  not_focal <- pm[lfsr_sig1$id[i], which(lfsr_sig1[i, tissue_oi] != z)]
  not_focal2 <- log2(as.numeric(not_focal) / as.numeric(focal)) #negative numbers returning NaN
  
  lfsr_sig1$within2_v2[i] <- length(which(not_focal2 > -1 & not_focal2<1))
  not_focal_df <- t(as.data.frame(not_focal2))
  colnames(not_focal_df) <- colnames(not_focal)
    
  lfsr_sig1$tissues_within2[i] <- paste(colnames(not_focal_df)[1:12][which(not_focal2 > -1 & not_focal2<1)], collapse = ",")
  lfsr_sig1$focal_pm[i] <- mean(focal)
  lfsr_sig1$other_pm[i] <- mean(t(not_focal))
}

lfsr_sig1 <- lfsr_sig1 %>%
  rowwise() %>%
  mutate(focal_tissue = names(.)[1:12][which.min(c_across(1:12))]) %>%
  ungroup() %>%
  as.data.frame()

rownames(lfsr_sig1) <- rownames_sig1
lfsr_sig1$site <- rownames(lfsr_sig1)

write.table(as.data.frame(lfsr_sig1),file.path(bed_path,"age_sharing_sig005.txt"),row.names=F,sep='\t')

# === Export bed files ===
bed = data.frame(do.call("rbind",lapply(strsplit(gsub("Region_", "", rownames(lfsr)), split = "_", fixed=TRUE), function (x) x[1:3])))
colnames(bed) <- c("chr","start","end")
bed[c("start","end")]<- lapply(bed[c("start","end")], as.integer)
bed <- bed %>% mutate(chr = paste0("chr",chr))

write.table(bed,file.path(bed_path,"all_sites_age.bed"), col.names = FALSE, sep="\t", row.names = FALSE, quote=FALSE)
write.table(data.frame(regions=rownames(lfsr)),file.path(bed_path,"all_sites_age.txt"), sep="\t",row.names = FALSE,quote=FALSE)

# CpG-level bed export
regions_to_cpg <- read.table(file.path(base_path,"Regions","regions_to_cpgs_mapping.bed"))

bed_cpg <- regions_to_cpg %>%
           filter(paste(V1,V2,V3,sep="_") %in% paste(bed$chr,bed$start,bed$end, sep = "_")) %>%
           dplyr::select(V4,V5,V6,V1,V2,V3)

write.table(bed_cpg, file.path(bed_path, "all_sites_age_CpG.bed"), col.names = FALSE, sep="\t", row.names = FALSE, quote=FALSE)

# === Upset Plots ===

# Split the tissues by comma and create a list
sigtissues_listall <- strsplit(as.character(lfsr_sig1$sig_which_tissues), ",")

# Create a binary presence/absence data.frame
tissue_matrixall <- data.frame(matrix(0, nrow = nrow(lfsr_sig1), ncol = 0))
row.names(tissue_matrixall) <- lfsr_sig1$site

for (tissue in unique(unlist(sigtissues_listall))) {
  tissue_matrixall[[tissue]] <- as.integer(sapply(sigtissues_listall, function(x) tissue %in% x))
}

ComplexUpset::upset(tissue_matrixall, colnames(tissue_matrixall), n_intersections = 20,height_ratio = 0.9, name="", set_sizes=F,
                    themes=upset_modify_themes(list(
                      'intersections_matrix'=theme(text=element_text(size=14)),
                      "Intersection size"=theme(text=element_text(size=14), axis.title.y = element_text(margin = margin(t = 0, r = -10, b = 0, l = 0))))),
                    labeller=ggplot2::as_labeller(c(
                      "thymus"="thymus", "pituitary"="pituitary", "spleen"="spleen", 'whole_blood' = 'whole blood',
                      'omental_at'= "omental adipose", "liver"="liver", "thyroid"="thyroid", "adrenal"="adrenal",
                      'skeletal_muscle'= "skeletal muscle", "kidney"="kidney", "lung"="lung", "heart"="heart")),
                    base_annotations = list(
                      'Intersection size'= (intersection_size(counts=F) +
                                            scale_y_continuous(labels = label_number(scale = 1e-3, suffix = "k")))
                    )
)

ggsave(file.path(figure_path,"Fig.2E.pdf"))
