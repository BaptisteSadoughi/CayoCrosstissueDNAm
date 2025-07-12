#!/usr/bin/env /packages/apps/spack/21/opt/spack/linux-rocky8-zen3/gcc-12.1.0/r-4.4.0-4yi4nm4foi7jsbczjxvv77uq7adnzb67/bin/Rscript

# -------------------------------------------------------------
# SLURM submission command (run as array over tissues):
# sbatch --cpus-per-task=22 --mem-per-cpu=1G -p general -q public -t 0-02:00:00 --array=1-(number of tissue) base_path/Bioinformatic and R scripts/thiscript
# -------------------------------------------------------------

# - Load tissue-specific percent methylation matrices
# - Filter for regions with full coverage
# - Impute missing values

rm(list = ls())

library_list <- c("parallel", "dplyr", "tidyr","methyLImp2","BiocParallel")

lapply(library_list, require, character.only = TRUE)

# === USER DEFINED BASE PATH ===

base_path <- "/path/to/project"  # <-- Define this path only once

# === Configuration ===

tissue_oi <- c("lung", "kidney", "heart", "liver", "spleen", "skeletal_muscle","adrenal",
               "thyroid", "thymus", "whole_blood", "omental_at","pituitary", "testis", "ovaries")

pattern_ <- "_1000_14T"

# === Paths ===

main_folder <- file.path(base_path,"tissues_meth")
output_full_path <- file.path(base_path, "full_matrices")
output_imputed_path <- file.path(base_path, "imputed")

dir.create(output_full_path, recursive = TRUE, showWarnings = FALSE)
dir.create(output_imputed_path, recursive = TRUE, showWarnings = FALSE)

sub_folders <- paste0(tissue_oi,"_meth")

# Construct full paths for additional folders
full_folders <- file.path(main_folder, sub_folders)

# === Extract pmeth from RDS files ===

extract_pmeth <- function(file_path) {
  readRDS(file_path)$pmeth
}

extract_pmeth_from_files <- function(folder_path, pattern) {
  files_ <- list.files(folder_path, full.names = TRUE)
  file_ <- files_[grep(pattern, files_)]
  extract_pmeth(file_)
}

# === Assemble list of pmeth matrices across tissues ===

pmeth_list <- setNames(lapply(tissue_oi, function(tissue) {
  extract_pmeth_from_files(file.path(main_folder, paste0(tissue, "_meth")), pattern_)
}), tissue_oi)

# === Save matrices with complete coverage (no missing values) ===

pmeth_list_full <- lapply(pmeth_list, function(mat) mat[complete.cases(mat), , drop = FALSE])
saveRDS(pmeth_list_full, file = file.path(output_full_path, paste0("Regions_full_pmeth", length(tissue_oi), ".rds")))

# ------------------------------------------
# === DATA IMPUTATION ===
# ------------------------------------------

# For details on methyLimp2 see
# https://scholar.google.com/citations?view_op=view_citation&hl=en&user=YXzVOS8AAAAJ&citation_for_view=YXzVOS8AAAAJ:qjMakFHDy7sC
# https://academic.oup.com/bioinformatics/article/35/19/3786/5364016?login=false

# === One SLURM task per tissue ===
SAMP <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

tissue <- tissue_oi[SAMP]
message("Processing tissue: ", tissue)

pmeth <- pmeth_list[[tissue]]
rm(pmeth_list)

# Omit X which displays strong sex-differences
pmeth <- pmeth[!grepl("Region_X",rownames(pmeth)), ]

# === Build annotation ===
annotation <- data.frame("cpg"=rownames(pmeth)) %>%
                          mutate(site=cpg, strand="+") %>%
                          tidyr::separate(site, into = c("site","chr","start","end"), sep = "_", fill = "right") %>%
                          select(-site) %>%
                          mutate(across(c("chr", "start", "end"), as.numeric)) %>%
                          arrange(chr, start)

# Reorder matrix to match annotation
annotation <- annotation[order(annotation$chr, annotation$start), ]
pmeth <- pmeth[order(annotation$chr, annotation$start), ]

# === Impute missing data ===

imputed_data <- methyLImp2(input=t(pmeth),
                           type = "user",
                           annotation = annotation,
                           BPPARAM = SnowParam(workers = 20, exportglobals = FALSE),
                           overwrite_res = TRUE)

# === Save imputed matrix ===

saveRDS(t(imputed_data), file = file.path(base_path, "imputed", paste0("pmeth_imp_methyLimp_", tissue, ".rds")))

print(paste("COMPLETED",tissue))
