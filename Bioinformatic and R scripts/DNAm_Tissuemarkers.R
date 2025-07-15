# ==============================================================================
# Define tissue-specific differentially methylated regions
# -----------------------------------
# The aim is to identify loci exhibiting consistent differences in DNAm levels with other tissues.
# ==============================================================================

rm(list=ls())

library_list <- c("dplyr","qvalue","PQLseq")
lapply(library_list, require, character.only = TRUE)

# === Paths ===

base_path <- "/path/to/your/directory" # <-- set the path once

metadata_path <- file.path(base_path,"metadata", "multitissue_metadata.txt")
kinship_path <- file.path(base_path, "metadata", "wgs_kinmat.rds")
alltissues <- c("liver", "omental_at", "spleen", "kidney", "lung", "heart", "skeletal_muscle", "adrenal", "pituitary", "thymus", "thyroid", "whole_blood")

# === Load data ===

metadata <- read.table(metadata_path, sep = "\t", header = TRUE) %>%
  mutate(AnimalID.text= paste("'", monkey_id, sep=""),
        percent_unique = unique / reads) %>%
  filter(lid_pid != "LID_109490_PID_10416")

kinmat <- readRDS(kinship_path)

# === Load methylation and coverage data ===

cov_list <- list()
for (i in alltissues) {
  filename <- gsub("XXX", i, file.path(base_path, "tissues_meth", "XXX_meth", 
                                       "Regions_pmeth_full_XXX_1000_14T.rds"))
  r <- readRDS(filename)
  cov_list[[i]] <- r$coverage
}

# Filter shared rows across tissues
shared_rows <- Reduce(intersect, lapply(cov_list, rownames))

rm(cov_list)
cov_list <- list()
meth_list <- list()

# Reload and filter to shared rows
for (i in alltissues) {
  filename <- gsub("XXX", i, file.path(base_path, "tissues_meth", "XXX_meth", 
                                       "Regions_pmeth_full_XXX_1000_14T.rds"))
  r <- readRDS(filename)
  cov_list[[i]] <- r$coverage[shared_rows, ]
  meth_list[[i]] <- r$methylation[shared_rows, ]
}

message("Binding methylation and coverage matrices...")

meth_all <- do.call(cbind, meth_list)
cov_all <- do.call(cbind, cov_list)

message("Finished binding.")

# === Focal tissue from command line ===

in_args <- commandArgs(trailingOnly = TRUE)
focaltissue = in_args
message("Focal tissue: ", focaltissue)

# Match methylation and metadata
countsfiltered <- meth_all[, colnames(meth_all) %in% metadata$lid_pid]
covfiltered <- cov_all[, colnames(cov_all) %in% metadata$lid_pid]
covordered <- covfiltered[, colnames(countsfiltered)]
stopifnot(identical(colnames(covfiltered), colnames(countsfiltered)))

# === Prepare Kinship Matrix ===

metadata_lid <- metadata[, c("monkey_id", "lid_pid")]
kinmatdf <- as.data.frame(kinmat)
kinmatdf$monkey_id <- rownames(kinmat)

kinmat_expand1 <- left_join(kinmatdf, metadata_lid)
rownames(kinmat_expand1) <- kinmat_expand1$lid_pid
kinmat_expand1_tidy <- kinmat_expand1[, !(names(kinmat_expand1) %in% c("monkey_id", "lid_pid"))]
kinmat_expand1_tidy_t <- t(kinmat_expand1_tidy)

kinmat_t1 <- as.data.frame(kinmat_expand1_tidy_t)
kinmat_t1$monkey_id <- rownames(kinmat_t1)
kinmat_expand2 <- left_join(kinmat_t1, metadata_lid)
rownames(kinmat_expand2) <- kinmat_expand2$lid_pid
kinmat_expand2_tidy <- kinmat_expand2[, !(names(kinmat_expand2) %in% c("monkey_id", "lid_pid"))]

dnam_kinship_new <- kinmat_expand2_tidy[colnames(countsfiltered), colnames(countsfiltered)]

# === Build phenotype and covariates ===

rownames(metadata) <- metadata$lid_pid
metadata_ordered <- metadata[colnames(countsfiltered), ]

metadata_tissuebin <- metadata_ordered %>%
  mutate(tissue_bin = ifelse(tissue == focaltissue, "1", "0"))

predictor <- metadata_tissuebin$tissue_bin
  
age <- metadata_tissuebin$age
sex_bino <- as.numeric(ifelse(metadata$sex == "F", 0, 1))
percent_unique <- metadata_tissuebin$percent_unique
covariates <- as.matrix(cbind(age, percent_unique, sex_bino))

print(length(which((colnames(dnam_kinship_new)==colnames(covordered))==F)))
print(length(which((rownames(dnam_kinship_new)==colnames(countsfiltered))==F)))

# === Subset to 10,000 sites ===
countsfiltered10 <- countsfiltered[1:10000, ]
covordered10 <- covordered[1:10000, ]

# === Fit PQLseq Model ===

fit <- PQLseq::pqlseq(RawCountDataSet=countsfiltered10,
                      Phenotypes=predictor,
                      RelatednessMatrix=dnam_kinship_new,
                      LibSize=covordered10,
                      Covariates = covariates,
                      fit.model="BMM",
                      numCore = 40)

# === Post-processing ===

converged <- subset(fit, converged == "TRUE")
converged$qvalue <- qvalue(converged[,"pvalue"])$qvalues
out_path <- gsub("XXX", focaltissue, file.path(base_path, "tissue_comparisons", "XXX_v_all_pqlseq.rds"))
saveRDS(converged, file = out_path)
