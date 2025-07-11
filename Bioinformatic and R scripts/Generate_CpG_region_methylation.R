#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript

# ==============================================================================
# Generate CpG Region Methylation Data
# -----------------------------------
# The aim is to extract methylation and coverage at regions of interest.
# The procedure involves a number of variables which should all be set at the start of this script.
#
# The resulting tissue-level matrices are exported as separate objects.
# ==============================================================================

rm(list = ls())

# sbatch --cpus-per-task=22 --mem=300G -p general -q public -t 0-08:00:00 --mail-user=sbaptis7@asu.edu --mail-type=ALL /path/to/Generate_CpG_region_methylation.R

library_list <- c("bsseq","BiocGenerics","GenomicRanges","GenomicFeatures","tidyverse","comethyl","parallel","purrr","DelayedMatrixStats")

lapply(library_list, require, character.only = TRUE)

# -------------------
# === PARAMETERS ===
# -------------------

base_path <- "base_path" #adjust accordingly
chrs <- paste0("",c(1:20,"X")) #adapt to genome of the species
num_cores <- 22 #number of cores to use for parallel tasks (chrs+1 is advised)
gap_region <- 1000 #maximum gap between two consecutive CpG sites in a region
hypomethylated <- 0.1 #upper threshold to consider a region as hypomethylated
hypermethylated <- 0.9 #lower threshold to consider a region as hypermethylated
in_at_least_Xperc_samples <- 0.25 # minimum fraction of samples with coverage above threshold in a tissue
med_cov_filter <- 0 #median coverage filter within tissue (0 means no filter)
path_bsseq_file <- file.path(base_path, "Rhesus_bismarkBSseq.rds") #define the path to the bsseq formatted methylation data
path_to_metadata <- file.path(base_path, "metadata", "multitissue_metadata.txt")
min_cov <- 1       # minimum coverage per CpG in a sample to be counted
p_sample <- 0.25   # minimum fraction of samples with min_cov coverage per CpG
med_cov <- 5       # median coverage threshold across samples per CpG
path_to_output <- file.path(base_path, "Regions")

# Source external helper functions
source(file.path(base_path, "Bioinformatic and R scripts", "SupportFunctions_generate_region_methylation.R"))

#------------------------------------------
# ########## READ IN BSSEQ RAW DATA ####
#------------------------------------------
DNAm_data <- readRDS(path_bsseq_file)

# Simplify sample names
colnames(DNAm_data) <- gsub(".CpG_report.merged_CpG_evidence.cov.gz", "",
                       str_split_i(colnames(DNAm_data), "/", 6))

# Read prepared metadata
metadata_lid <- read.table(path_to_metadata, sep = "\t", header = TRUE)

# Subset bsseq object by colnames present in the metadata lid_pid
DNAm_data_red <- DNAm_data[,metadata_lid$lid_pid]

# Check consistency between metadata and bsseq sample names
if (length(metadata_lid$lid_pid) != length(colnames(DNAm_data_red)) || !all.equal(metadata_lid$lid_pid, colnames(DNAm_data_red))) {
  stop("STOP: Metadata and samples info are not equal")
}

# Remove original full bsseq object to free memory
rm(DNAm_data)

# --------------------------------------
# === SPLIT BSSEQ OBJECT BY CHROMOSOME ===
# --------------------------------------

DNAm_data_red_list <- parallel::mclapply(chrs,function(x){
  chr <- chrSelectBSseq(DNAm_data_red, seqnames = x, order = TRUE)
  return(chr)
  }, mc.cores = num_cores)

# ---------------------------------
# === FILTER CpGs BASED ON COVERAGE ===
# ---------------------------------

# Get unique tissue types from metadata
tissues = c(unique(metadata_lid$grantparent_tissueType))
tissue_oi = tissues

DNAm_data_red_filtered_list=mclapply(chrs, filter_CpGs,
                                bsseq_list=DNAm_data_red_list,
                                metadata=metadata_lid,
                                tissues=tissues,
                                min_cov=min_cov,
                                p_sample=p_sample,
                                med_cov=med_cov,
                                mc.cores= num_cores)

# ----------------------------------------------
# === ADD PERCENT METHYLATION TO FILTERED DATA ===
# ----------------------------------------------

DNAm_data_red_filtered_list=parallel::mclapply(
  DNAm_data_red_filtered_list, function(x){
    pmeth = getCoverage(x, type = "M") / getCoverage(x, type = "Cov")
    obj = x
    obj@assays@data@listData[["pmeth"]] = pmeth
    return(obj)
  },mc.cores = num_cores
)

# Create output directory if needed
if (!dir.exists(paths = path_to_output)) dir.create(path_to_output)

# Save the list format
save(DNAm_data_red_filtered_list, file = paste0(path_to_output,"/DNAm_data_red_CpG_filtered.RData"))

rm(DNAm_data_red)

# -----------------------------------
# === DEFINE REGIONS OF INTEREST ===
# -----------------------------------

regions_list <- parallel::mclapply(DNAm_data_red_filtered_list,
                                   getting_regions,
                                   gap_region = gap_region,
                                   mc.cores = num_cores)

# Create output directory for batches of regions processed sequentially
regions_batch_dir <- file.path(path_to_output, "regions_batches", as.character(gap_region))
if (!dir.exists(regions_batch_dir)) {
  dir.create(regions_batch_dir, recursive = TRUE)
}

# --------------------------------------------------------
# === WRITE COVERAGE AND METHYLATION FOR REGION BATCHES ===
# --------------------------------------------------------

parallel::mclapply(1:length(DNAm_data_red_filtered_list),
                   function(x) {
                     getting_regions_coverage(
                       DNAm_data_meth = DNAm_data_red_filtered_list[[x]],
                       regions = regions_list[[x]]
                     )
                   },
                   mc.cores = 22)

# -----------------------------------------
# === GROUP REGION BATCHES BY CHROMOSOME AND MERGE ===
# -----------------------------------------

# List all files created above
file_list <- list.files(path = regions_batch_dir, pattern = "_chr(X|\\d+)_", full.names = TRUE)

# Group files by the common _chr_ identifier
grouped_files <- split(file_list, gsub(".*_chr(X|\\d+)_batch.*", "\\1", file_list))

# Load back and assemble per chromosome
result_list <- map(grouped_files, function(files) {
  data_list <- map(files, readRDS)

  # Ensure number of columns match for all components
  num_cols_regions <- map_dbl(data_list, ~ ncol(.x$regions))
  num_cols_coverage <- map_dbl(data_list, ~ ncol(.x$coverage))
  num_cols_methylation <- map_dbl(data_list, ~ ncol(.x$methylation))
  num_cols_pmeth <- map_dbl(data_list, ~ ncol(.x$pmeth))

  if (length(unique(num_cols_regions)) > 1 ||
      length(unique(num_cols_coverage)) > 1 ||
      length(unique(num_cols_methylation)) > 1 ||
      length(unique(num_cols_pmeth)) > 1) {
    stop("Number of columns is not consistent across files.")
  }

  # Combine each data type across batches
  regions_combined <- do.call(rbind, lapply(data_list, function(x) x$regions))
  coverage_combined <- do.call(rbind, lapply(data_list, function(x) x$coverage))
  methylation_combined = do.call(rbind, lapply(data_list, function(x) x$methylation))
  pmeth_combined = do.call(rbind, lapply(data_list, function(x) x$pmeth))

  # Combine all components into a list
  combined_data <- list(
    regions = regions_combined,
    coverage = coverage_combined,
    methylation = methylation_combined,
    pmeth = pmeth_combined
  )

  return(combined_data)
})

# ----------------------------
# === CLEANUP TEMPORARY FILES ===
# ----------------------------

unlink(regions_batch_dir, recursive = TRUE)

# -----------------------------------
# === REORDER RESULT LIST BY CHROMOSOME ===
# -----------------------------------
                           
# Identify numeric chr and sort
numeric_names <- grep("^\\d+$", names(result_list), value = TRUE)
sorted_numeric_names <- sort(as.numeric(numeric_names))

# Ensure "X" is in the list
if ("X" %in% names(result_list)) {
  sorted_names <- c(as.character(sorted_numeric_names), "X")
} else {
  sorted_names <- as.character(sorted_numeric_names)
}

# Reorder the list
result_list <- result_list[sorted_names]

# combine the results
final_result <- list(
  regions = do.call(rbind, lapply(result_list, function(x) x$regions)),
  coverage = do.call(rbind, lapply(result_list, function(x) x$coverage)),
  methylation = do.call(rbind, lapply(result_list, function(x) x$methylation)),
  pmeth = do.call(rbind, lapply(result_list, function(x) x$pmeth))
)

# ------------------------------------------
# === CALCULATE DESCRIPTIVE STATS ON REGIONS ===
# ------------------------------------------
                  
final_result <- amend_regions(final_result)
final_result <- get_summary_regions(final_result)

# Export the data generated
save(final_result, result_list, file = paste0(path_to_ouput,"/Regions_",gap_region,"_",length(tissue_oi),"T.RData"))

rm(final_result)

# ------------------------------------------
# === FILTER HYPO/HYPER-METHYLATED REGIONS ===
# ------------------------------------------
                  
# Filter regions consistently hypo/hyper methylated across tissues
DNAm_data_red_list_methfiltered <- mclapply(chrs, compute_and_filter_regions,
                                       bsseq_list = result_list,
                                       metadata = metadata_lid,
                                       tissues = tissue_oi,
                                       hypomethylated = hypomethylated,
                                       hypermethylated = hypermethylated,
                                       mc.cores = num_cores)

# ------------------------------------------
# === GENERATE TISSUE-SPECIFIC REGION MATRICES ===
# ------------------------------------------

for (tissue_level in tissue_oi){

 # Apply tissue_specific_coverage to each element of DNAm_data_red_list_methfiltered
 result_tissue <- mclapply(1:length(DNAm_data_red_list_methfiltered),
                           function(i) tissue_specific_coverage(
                             matrix_list=DNAm_data_red_list_methfiltered[[i]],
                             tissue = tissue_level,
                             metadata_sample = metadata_lid,
                             in_at_least_Xperc_samples = in_at_least_Xperc_samples,
                             med_cov_filter = med_cov_filter
                           ),
                           mc.cores = num_cores)

 # Combine all chr
 all_chr_result <- list(
   regions = do.call(rbind, lapply(result_tissue, function(x) x$regions)),
   coverage = do.call(rbind, lapply(result_tissue, function(x) x$coverage)),
   methylation = do.call(rbind, lapply(result_tissue, function(x) x$methylation)),
   pmeth = do.call(rbind, lapply(result_tissue, function(x) x$pmeth))
 )

 # Calculate descriptive stats on the regions for that tissue
 all_chr_result <- amend_regions(all_chr_result)
 all_chr_result <- get_summary_regions(all_chr_result)

 # Create tissue specific output directory
 tissue_dir <- file.path(path_to_output, "tissues_meth", paste0(tissue_level, "_meth"))
 if (!dir.exists(tissue_dir)) dir.create(tissue_dir, recursive = TRUE)

 # Save tissue-level results
  saveRDS(all_chr_result, file = file.path(
    tissue_dir,
    paste0("Regions_pmeth_full_", tissue_level, "_", gap_region, "_", length(tissue_oi), "T.rds")
  ))
}

# ------------------------------------------
# === CLEAN WORKSPACE ===
# ------------------------------------------

rm(all_chr_result,grouped_files,regions_list,result_list,result_tissue) #clean workspace
