#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript

### Differential Methylation Analysis with PQLSEQ

# Submit with:
# sbatch --cpus-per-task=40 --mem=100G -p general -q public -t 0-4 --array=1-14 /path/to/this_script.R

rm(list = ls())

# ====================
# Library dependencies
# ====================
library_list <- c("bsseq", "BiocGenerics", "GenomicRanges", "GenomicFeatures",
                  "tidyverse", "PQLseq", "purrr", "svglite", "qvalue")
lapply(library_list, require, character.only = TRUE)

# ================
# CONFIGURABLE PATHS
# ================
base_path <- "/your/base/path"
metadata_path <- file.path(base_path, "metadata.txt")
meth_dir <- file.path(base_path, "tissues_meth")
kinship_path <- file.path(base_path, "wgs_kinmat.rds")
output_dir <- file.path(base_path, "PQLSEQ")
diagnostics_dir <- file.path(output_dir, "diagnostics")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(diagnostics_dir)) dir.create(diagnostics_dir, recursive = TRUE)

# =============
# Parallel cores
# =============
n.cores <- 40
print(paste("Using", n.cores, "cores"))

# ====================
# Tissue selection
# ====================
tissue_oi <- c("liver", "whole_blood", "spleen", "omental_at", "heart",
               "ovaries", "testis", "kidney", "lung", "adrenal",
               "thymus", "thyroid", "pituitary", "skeletal_muscle")

SAMP <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

if (is.na(SAMP) || SAMP < 1 || SAMP > length(tissue_oi)) {
  stop("Invalid or missing SLURM_ARRAY_TASK_ID")
}

tissue <- tissue_oi[SAMP]
cat("Analyzing tissue:", tissue, "\n")

# =====================
# Load data and metadata
# =====================
metadata <- read.table(metadata_path, sep = "\t", header = TRUE) %>%
  filter(lid_pid != "LID_109490_PID_10416") %>%
  mutate(percent_unique = unique / reads)

rds_path <- file.path(meth_dir, paste0(tissue, "_meth"), paste0("Regions_pmeth_full_", tissue, "_1000_14T.rds"))
if (!file.exists(rds_path)) stop("Missing RDS file: ", rds_path)

r <- readRDS(rds_path)

# Match metadata to available samples
subset_metadata <- metadata[metadata$lid_pid %in% colnames(r[["coverage"]]), ]
r$metadata <- subset_metadata
assign(tissue, r)
rm(r, subset_metadata)

# Remove flawed skeletal muscle sample
if (tissue == "skeletal_muscle") {
  sm <- get(tissue)
  bad_id <- "LID_109490_PID_10416"
  for (layer in c("coverage", "methylation", "pmeth")) {
    sm[[layer]] <- sm[[layer]][, !(colnames(sm[[layer]]) == bad_id)]
  }
  assign(tissue, sm)
}

# ===============
# Load kinship matrix
# ===============
kinmat <- readRDS(kinship_path)

# =========================
# Prepare data for PQLseq
# =========================
tissue_data <- get(tissue)
methcount <- tissue_data$methylation
coverage  <- tissue_data$coverage
metadata_tissue <- tissue_data$metadata

# Filter metadata for age and matching IDs
metadata_tissue <- metadata_tissue %>%
  filter(age_at_sampling >= 2.9) %>%
  column_to_rownames("lid_pid")

coverage <- coverage[, colnames(coverage) %in% rownames(metadata_tissue)]
methcount <- methcount[, colnames(methcount) %in% rownames(metadata_tissue)]

# Match ordering
metadata_tissue <- metadata_tissue[colnames(coverage), ]

# Subset kinship
monkey_ids <- metadata_tissue$monkey_id
dnam_kinship <- kinmat[monkey_ids, monkey_ids]

# Check for consistency
if (!identical(colnames(coverage), rownames(metadata_tissue)) ||
    !identical(colnames(methcount), rownames(metadata_tissue)) ||
    !identical(colnames(dnam_kinship), monkey_ids) ||
    !identical(rownames(dnam_kinship), monkey_ids)) {
  stop("Mismatch in data matrices")
}

# ================
# Define covariates
# ================
age <- metadata_tissue$age
sex_bino <- as.numeric(metadata_tissue$sex == "M")
group_bino <- as.numeric(metadata_tissue$group == "KK")
percent_unique <- metadata_tissue$percent_unique

covariates <- if (tissue %in% c("ovaries", "testis")) {
  as.matrix(cbind(percent_unique, group_bino))
} else {
  as.matrix(cbind(percent_unique, sex_bino, group_bino))
}

# ==============================================
# Remove sites with 0 coverage in either group
# ==============================================
remove_zero_coverage_sites <- function(group_label) {
  sample_ids <- rownames(metadata_tissue)[metadata_tissue$group == group_label]
  zero_rows <- rowSums(coverage[, sample_ids]) == 0
  rownames(coverage)[zero_rows]
}

sites_to_remove <- unique(c(remove_zero_coverage_sites("HH"), remove_zero_coverage_sites("KK")))

coverage <- coverage[!rownames(coverage) %in% sites_to_remove, ]
methcount <- methcount[!rownames(methcount) %in% sites_to_remove, ]

# =============
# Run PQLseq
# =============
fit <- PQLseq::pqlseq(
  RawCountDataSet = methcount,
  Phenotypes = age,
  RelatednessMatrix = dnam_kinship,
  LibSize = coverage,
  Covariates = covariates,
  fit.model = "BMM",
  numCore = n.cores
)

# =========================
# Post-processing and output
# =========================
fit$chr <- str_split_i(rownames(fit), 2, "_")
fit$start <- str_split_i(rownames(fit), 3, "_")
fit$end <- str_split_i(rownames(fit), 4, "_")
fit$sites <- rownames(fit)

fit_conv <- fit %>% filter(converged == TRUE)

# Save p-value histogram
ggplot(fit_conv, aes(x = pvalue)) +
  geom_histogram(bins = 20) +
  theme_minimal() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14))

ggsave(filename = file.path(diagnostics_dir, paste0(tissue, "_pvalhist.png")),
       width = 6, height = 4, dpi = 300)

# Compute q-values and BH-adjusted p-values
qresults <- qvalue::qvalue(fit_conv$pvalue)
fit_conv$qval <- qresults$qvalues
fit_conv$padj <- p.adjust(fit_conv$pvalue, method = "BH")

# Save BED-style output
write.table(fit_conv,
            file = file.path(output_dir, paste0("bed_", tissue, ".txt")),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("Finished PQLseq for", tissue, "\n")
