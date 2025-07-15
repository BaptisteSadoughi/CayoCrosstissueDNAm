#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript

# ==============================================================================
# Leave One Out Age Prediction
# SLURM submission command (run as array over samples):
# submit using: sbatch --cpus-per-task=10 --mem=100G -p general -q public -t 0-4 --array=1-(max sample size within one tissue) base_path/Bioinformatic and R scripts/thiscript.R
# Produces: single_sample-single_alpha files saved in tissue specific subfolders.
# ==============================================================================

rm(list = ls())

library_list <- c("glmnet","tidyverse","parallel","stringr")
lapply(library_list, require, character.only = TRUE)

num_cores <- 15 #number of cores to use for parallel tasks (tissue+1 is good)

maturity <- 5 #set age at maturity for non-linear relationship

base_path <- "/path/to/your/directory" # <- set this path once
metadata_path <- file.path(base_path,"metadata","multitissue_metadata.txt")

# -------------------
# === LOAD METHYLATION DATA AND METADATA ===
# -------------------

# Load methylation data
file_name <- file.path(base_path, "full_matrices","Regions_full_pmeth14.rds")
pmeth <- readRDS(file_name)

# Remove one flawed sample
pmeth$skeletal_muscle <- pmeth$skeletal_muscle[,-which(colnames(pmeth$skeletal_muscle) == "LID_109490_PID_10416")]

# keep autosomes only (run call) or all chr (mute call)
pmeth <- lapply(pmeth,function(x) x[!grepl("Region_X|CpG_X",rownames(x)), ,drop = FALSE])

# Load metadata
meta <- read.table(metadata_path, sep = "\t", header = TRUE) %>% 
                filter(lid_pid !="LID_109490_PID_10416") %>% 
                dplyr::select(lid_pid, age, monkey_id, tissue)

# Match metadata and methylation data
meta <- meta[meta$grantparent_tissueType %in% names(pmeth), ]
all_colnames <- unique(unlist(lapply(pmeth, colnames)))
meta <- meta[meta$lid_pid %in% all_colnames, ]

# Check that metadata and methylation data match
for(tissue_value in names(pmeth)){
  if (length(meta$lid_pid[meta$grantparent_tissueType == tissue_value]) != length(colnames(pmeth[[as.character(tissue_value)]])) ||
      !all.equal(meta$lid_pid[meta$grantparent_tissueType == tissue_value], colnames(pmeth[[tissue_value]]))) {
    stop(paste("STOP metadata and samples info are not equal for", tissue_value))
  }
}

# Confirm absence of NAs for the elastic net regression.
pmeth <- lapply(pmeth, function(x) x[complete.cases(mat), , drop = FALSE])

# -------------------
# === AGE TRANSFORMATION ===
# -------------------

# Horvath (2013)-style transformation https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-10-r115#MOESM2
Fage = Vectorize(function(x){
  if (is.na(x) | is.na(maturity)) {return(NA)}
  k <- 0
  y <- 0
  if (x < maturity) {y = log(x+k)-log(maturity+k)}
  else {y = (x-maturity)/(maturity+k)}
  return(y)
})
                
## Inverse transformation
F.inverse= Vectorize(function(y) {
  if (is.na(y) | is.na(maturity)) {return(NA)}
  k <- 0
  x <- 0
  if (y < 0) {x = (maturity+k)*exp(y)-k}
  else {x = (maturity+k)*y+maturity}
  return(x)
})

# -------------------
# === ELASTIC NET MODEL ===
# -------------------

tissues <- names(pmeth)

# run tissues in parallel independently
parallel::mclapply(tissues, function(tissue) {
  
  epi <- pmeth[[tissue]]
  meta_s <- meta[meta$tissue == tissue, ]
  
  SAMP <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
    
  # Define training and test sets
  train <- epi[, -SAMP]
  test <- epi[, SAMP]
  trainage <- meta_s$age[-SAMP]
  testage <- meta_s$age[SAMP]
  test_id <- meta_s$monkey_id[SAMP]
  test_lid <- meta_s$lid_pid[SAMP]
  
  # Alphas values (0=ridge coeff are shrunk towards each other, 1=lasso coeff are shrunk to 0)
  alphas <- seq(0, 1, by = 0.1)

  for (alph in alphas) {
    
    # IMPORTANT matrices must be samples (rows) x loci (columns)
    
    if(exists("Fage", envir = .GlobalEnv) & exists("F.inverse", envir = .GlobalEnv)){
    model <- cv.glmnet(t(train), Fage(trainage), type.measure = "mse", nfolds = 10, alpha = alph, standardize = FALSE) # with age transformation
    } else {
    model <- cv.glmnet(t(train), trainage, type.measure = "mse", nfolds = 10, alpha = alph, standardize = FALSE)
    }
    
    # Predict age of the test sample with parameters that minimized MSE during internal CV
    if(exists("Fage", envir = .GlobalEnv) & exists("F.inverse", envir = .GlobalEnv)){
    predicted <- F.inverse(predict(model, newx = t(test), type = "response", s = "lambda.min")) # with age transformation
    } else {
    predicted <- predict(model, newx = t(test), type = "response", s = "lambda.min")
    }
    
    predicted_out <- t(c(test_id, test_lid, testage, predicted))
        
    path_p <- file.path(base_path, "predicted_age", tissue)

    if (!dir.exists(paths = path_p)) dir.create(path_p, recursive = TRUE)

    write.table(predicted_out, paste0(path_p,"/AgePred_",nrow(meta_s),"_alpha", alph,"_", tissue, SAMP, ".txt"), quote = F, row.names = F, col.names = F)
    }
  }, mc.cores = num_cores)
