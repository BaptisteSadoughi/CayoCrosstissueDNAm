#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript

#submit using: sbatch --cpus-per-task=10 --mem=100G -p general -q public -t 0-4 --array=1-(max sample size within one tissue) /path/to/thiscript

# The aim is to split a list of percent methylation matrices per tissue of interest
# and run a LOOCV for each tissue.
# The output are single_sample-single_alpha files saved in subfolders for each tissue.

rm(list = ls())

library_list <- c("glmnet","tidyverse","parallel","stringr")
lapply(library_list, require, character.only = TRUE)

num_cores <- 15 #number of cores to use for parallel tasks (tissue+1 is good)

maturity <- 5 #set age at maturity for non-linear relationship

#---------------------------------------------------------------------
# LOAD METHYLATION DATA AND METADATA
#---------------------------------------------------------------------
file_name<-("/scratch/sbaptis7/Cayotissue_CpG_coverage/combined_CpG_Regions/full_matrices/Regions_full_pmeth14.rds")

# Read in data
pmeth <- readRDS(file_name)

# Remove one flawed sample
pmeth$skeletal_muscle <- pmeth$skeletal_muscle[,-which(colnames(pmeth$skeletal_muscle) == "LID_109490_PID_10416")]

# keep autosomes only (run call) or all chr (mute call)
pmeth <- lapply(pmeth,function(x) x[!grepl("Region_X|CpG_X",rownames(x)), ,drop = FALSE])

# Read in metadata info with known chronological ages/sex and technical variables.
metadata_lid = read.table("/scratch/sbaptis7/Cayo_meth_metadata/metadata_final_lidpids_Nov24.txt", sep = "\t", header = TRUE) %>% filter(lid_pid !="LID_109490_PID_10416")

# Simplify metadata for glmnet
meta <- meta %>% dplyr::select(lid_pid,age_at_sampling,monkey_id,grantparent_tissueType)

# Keep in metadata only tissues and lid_pid matching the methylation data
meta <- meta[meta$grantparent_tissueType %in% names(pmeth), ]

# Get unique colnames from all matrices in pmeth
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
pmeth <- lapply(pmeth, function(x){
  mat <- x
  mat <- mat[complete.cases(mat),,drop = FALSE]
  return(mat)
})

#----------------------------------------------
# Transformation of age
#----------------------------------------------
# Transform age according to age at sexual maturity
# based on Horvath 2013 https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-10-r115#MOESM2
Fage = Vectorize(function(x){
  if (is.na(x) | is.na(maturity)) {return(NA)}
  k <- 0
  y <- 0
  if (x < maturity) {y = log(x+k)-log(maturity+k)}
  else {y = (x-maturity)/(maturity+k)}
  return(y)
})
## Inverse log linear transformation
F.inverse= Vectorize(function(y) {
  if (is.na(y) | is.na(maturity)) {return(NA)}

  k <- 0
  x <- 0
  if (y < 0) {x = (maturity+k)*exp(y)-k}
  else {x = (maturity+k)*y+maturity}
  return(x)
})

#--------------------------------------------------
# Elastic net model
#--------------------------------------------------

# Store the different tissues
tissues <- names(pmeth)

# run tissues in parallel independently
parallel::mclapply(tissues,function(tissue){
  
  # extract methylation data and lid_pid for one tissue
  epi <- pmeth[[tissue]]
  meta_s <- meta[meta$grantparent_tissueType == tissue,]
  
  SAMP <- Sys.getenv("SLURM_ARRAY_TASK_ID")
  SAMP <- as.integer(SAMP)
 
    # for (SAMP in 1:ncol(epi)) { #loop can be used to run the code interactively
  # Remove test subject(s)
  # SAMP indexes from 1 to N samples
  train <- epi[, -SAMP]
  test <- epi[, SAMP]
  
  # Create a vector of training and test ages for elastic net model construction
  trainage <- meta_s$age_at_sampling[-SAMP]
  testage <- meta_s$age_at_sampling[SAMP]
  test_id <- meta_s$monkey_id[SAMP]
  test_lid <- meta_s$lid_pid[SAMP]
  
  # Create a vector of alphas to test (0=ridge coeff are shrunk towards each other, 1=lasso coeff are shrunk to 0)
  alphas <- seq(0, 1, by = 0.1)

  for (alph in alphas) {
    
    # Using N-fold internal CV, train the elastic net model using the training data

    ## Transpose the train matrix to be samples x features
    
    if(exists("Fage", envir = .GlobalEnv) & exists("F.inverse", envir = .GlobalEnv)){
      # if the age transformation is loaded in the environment use it
    model <- cv.glmnet(t(train), Fage(trainage), type.measure = "mse", nfolds = 10, alpha = alph, standardize = FALSE) # with age transformation
    } else {
      # otherwise set the model without
    model <- cv.glmnet(t(train), trainage, type.measure = "mse", nfolds = 10, alpha = alph, standardize = FALSE)
    }
    
    # Predict age using the test sample from parameters that minimized MSE during internal CV
    if(exists("Fage", envir = .GlobalEnv) & exists("F.inverse", envir = .GlobalEnv)){
    predicted <- F.inverse(predict(model, newx = t(test), type = "response", s = "lambda.min")) # with age transformation
    } else {
    predicted <- predict(model, newx = t(test), type = "response", s = "lambda.min")
    }
    
    ## Add monkey_id and LID info
    predicted_out <- c(test_id, test_lid, testage, predicted)
    
    ## Transpose so that it writes out as one row of information
    predicted_out <- t(predicted_out)
        
    path_p <- file.path(".", "predicted_age", tissue)

    if (!dir.exists(paths = path_p)) dir.create(path_p, recursive = TRUE)

    write.table(predicted_out, paste0(path_p,"/AgePred_",nrow(meta_s),"_alpha", alph,"_", tissue, SAMP, ".txt"), quote = F, row.names = F, col.names = F)
    }
  }, mc.cores = num_cores)
##################################################################
