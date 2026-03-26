#=====================================================================
# Tissue Classification using Multinomial GLMNET
# --------------------------------------------------------------
# submit using: sbatch --cpus-per-task=10 --mem=100G -p htc -q public -t 0-4 --array=1-14  base_path/Bioinformatic and R scripts/Tissue_classification_glmnet.R
#=====================================================================

# === Clear workspace ===
rm(list = ls())

# === Load libraries ===
library_list <- c("glmnet","tidyverse","dplyr","parallel","stringr")
lapply(library_list, require, character.only = TRUE)

# === Paths ===
base_path <- "/path/to/project"  # <-- Define this path only once

metadata_path <- file.path(base_path, "metadata", "multitissue_metadata.txt")
output_path <- file.path(base_path, "TissueClassificationGLMNET")
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

# === Optional: Restrict to specific tissues (set to NULL for all tissues) ===
tissue_oi <- NULL
# tissue_oi <- c("whole_blood", "spleen", "omental_at", "heart", "ovaries", "kidney", 
#                "lung", "adrenal", "thymus", "thyroid", "liver", "skeletal_muscle")

# === Cross-validation parameters ===
n_batches     <- 10
max_per_batch <- 250

# Set seed for reproducibility
set.seed(1003)

# Get SLURM array task ID (with local fallback)
SAMP <- Sys.getenv("SLURM_ARRAY_TASK_ID")
if (SAMP == "") {
  SAMP <- 1  # Default for local testing
  cat("Running in local test mode (batch 1)\n")
}
SAMP <- as.integer(SAMP)

# === Load data ===

file_name <- file.path(base_path, "imputed", "pmeth_imp_methyLimp_14.rds")
pmeth <- readRDS(file_name)

# Restrict methylation matrices to tissues of interest
if (!is.null(tissue_oi) && length(tissue_oi) > 0) {
  pmeth <- pmeth[names(pmeth) %in% tissue_oi]
}

# Keep autosomes only
pmeth <- lapply(pmeth, function(x) x[!grepl("Region_X|CpG_X", rownames(x)), , drop = FALSE])

# === Load metadadata ===
meta <- read.delim(metadata_path, stringsAsFactors = FALSE) %>%
  filter(lid_pid != "LID_109490_PID_10416") %>%
  select(lid_pid, monkey_id, tissue)

# Restrict to tissues of interest
if (!is.null(tissue_oi) && length(tissue_oi) > 0) {
  meta <- meta[meta$tissue %in% tissue_oi, ]
}

# Get sample names from matrices in pmeth
all_colnames <- unique(unlist(lapply(pmeth, colnames)))
meta <- meta[meta$lid_pid %in% all_colnames, ]

# === Format data ===

# Remove incomplete rows
pmeth <- lapply(pmeth, function(x){
  x[complete.cases(x), , drop = FALSE]
})

# Assemble the tissue-specific matrices into one
combine_matrices <- function(mat_list) {
  # Get the union of all row names across all matrices
  all_rows <- unique(unlist(lapply(mat_list, rownames)))
  
  # Initialize a list to hold reindexed matrices
  reindexed_matrices <- list()
  
  # Reindex each matrix to the union of row names
  for (i in seq_along(mat_list)) {
    mat <- mat_list[[i]]
    
    # Create a new matrix filled with NA for missing rows
    full_mat <- matrix(NA, nrow = length(all_rows), ncol = ncol(mat), 
                       dimnames = list(all_rows, colnames(mat)))
    
    # Fill in the values from the original matrix
    full_mat[rownames(mat), ] <- mat
    
    # Add the reindexed matrix to the list
    reindexed_matrices[[i]] <- full_mat
  }
  
  # Combine the matrices by columns (cbind)
  combined_matrix <- do.call(cbind, reindexed_matrices)
  
  return(combined_matrix)
}

# Combine matrices
combined <- combine_matrices(pmeth)
rm(pmeth)  # Free memory

# Keep rows with full coverage
combined <- combined[complete.cases(combined), ]
combined <- combined[, colnames(combined) %in% meta$lid_pid]

# Order metadata to match combined
meta <- meta[match(colnames(combined), meta$lid_pid), ]

stopifnot(identical(colnames(combined), meta$lid_pid))

# === Initialize empty batches ===
                
batches <- vector("list", n_batches)
for (i in seq_len(n_batches)) batches[[i]] <- character(0)

# Split samples by tissue and shuffle (reproducible)
tissue_samples <- lapply(
  split(meta$lid_pid, meta$tissue),
  function(x) sample(x)   # shuffles samples within tissue
)

# Fill batches in round-robin manner
b <- 1
repeat {
  for (tissue in names(tissue_samples)) {
    if (length(tissue_samples[[tissue]]) > 0 &&
        length(batches[[b]]) < max_per_batch) {
      
      # Assign one sample to the current batch
      batches[[b]] <- c(batches[[b]], tissue_samples[[tissue]][1])
      tissue_samples[[tissue]] <- tissue_samples[[tissue]][-1]
    }
  }
  
  # Move to next batch
  b <- b + 1
  if (b > n_batches) b <- 1
  
  # Stop when all batches full or all samples used
  if (all(sapply(batches, length) >= max_per_batch) ||
      all(sapply(tissue_samples, length) == 0)) break
}

# Define batch processing function
run_batch <- function(batch_idx, data_matrix) {
  test_samples <- batches[[batch_idx]]
  train_samples <- unlist(batches[-batch_idx])
  
  x_train <- t(data_matrix[, train_samples, drop = FALSE])
  y_train <- meta$tissue[match(train_samples, meta$lid_pid)]
  
  x_test <- t(data_matrix[, test_samples, drop = FALSE])
  y_test <- meta$tissue[match(test_samples, meta$lid_pid)]
  y_names <- meta$lid_pid[match(test_samples, meta$lid_pid)]
  
  fit <- cv.glmnet(x_train, y_train, family = "multinomial", type.multinomial = "grouped")
  
  preds <- predict(fit, x_test, s = "lambda.min", type = "class")
  
  res <- data.frame(sample = y_names, tissue = y_test, predicted = preds, batch = batch_idx)
  return(res)
}

# --------------------------------
# === Process current batch ===
# --------------------------------
                
if (SAMP >= 1 && SAMP <= n_batches) {
  # Process the current batch
  res <- run_batch(batch_idx = SAMP, data_matrix = combined)
  
  # Write results for this batch
  output_file <- paste0(output_path, "/glmnet_tissue_assignment_batch", SAMP, ".csv")
  write.csv(res, file = output_file, quote = FALSE, row.names = FALSE)
  
  cat("Batch", SAMP, "completed. Results saved to:", output_file, "\n")
  
  # After processing all batches, run aggregation
  # This will be triggered by the last batch (n_batches)
  if (SAMP == n_batches) {
    cat("All batches processed. Starting aggregation...\n")
    
    # List all batch files
    batch_files <- list.files(path = output_path, pattern = "glmnet_tissue_assignment_batch.*\\.csv$", full.names = TRUE)
    
    if (length(batch_files) > 0) {
      # Read and combine all batch results
      list_batches <- lapply(batch_files, read.csv)
      res_ <- do.call(rbind, list_batches)
      
      # Calculate accuracy
      accuracy <- sum(res_$tissue == res_$predicted) / nrow(res_) * 100
      cat("Overall accuracy:", round(accuracy, 2), "%\n")
      
      # Save aggregated results
      output_file <- paste0(output_path, "/glmnet_tissue_assignments.csv")
      write.csv(res_, file = output_file, quote = FALSE, row.names = FALSE)
      cat("Aggregated results saved to:", output_file, "\n")
    } else {
      cat("No batch files found for aggregation.\n")
    }
  }
} else {
  cat("Invalid batch index. SAMP =", SAMP, "but valid range is 1-", n_batches, "\n")
}

cat("Script completed.\n")
