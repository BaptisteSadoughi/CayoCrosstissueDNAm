#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript

# The aim is to extract methylation and coverage at regions of interest.
# The procedure involves a number of variables which should all be set at the start of this script.
#
# The resulting tissue-level matrices are exported as separate objects.

rm(list = ls())

# sbatch --cpus-per-task=22 --mem=300G -p general -q public -t 0-08:00:00 --mail-user=sbaptis7@asu.edu --mail-type=ALL /path/to/Generate_CpG_region_methylation.R

library_list <- c("bsseq","BiocGenerics","GenomicRanges","GenomicFeatures","tidyverse","comethyl","parallel","purrr","DelayedMatrixStats")

lapply(library_list, require, character.only = TRUE)

# Set the following variable !!!
chrs=paste0("",c(1:20,"X")) #adapt to genome of the species
num_cores <- 22 #number of cores to use for parallel tasks (chrs+1 is advised)
gap_region <- 1000 #maximum gap between two consecutive CpG sites in a region
hypomethylated <- 0.1 #upper threshold to consider a region as hypomethylated
hypermethylated <- 0.9 #lower threshold to consider a region as hypermethylated
in_at_least_Xperc_samples=0.25 #coverage above threshold in at least X% of the samples from a single tissue
med_cov_filter <- 0 #median coverage ????
path_bsseq_file = "/path/to/bsseq_file.rds" #define the path to the bsseq formatted methylation data
path_to_metadata = "/path/to/multitissue_metadata.txt.txt"
min_cov=1 #minimum coverage for a single CpG site in p_sample
p_sample=0.25 #percentage of samples in which a single CpG site meets the min_cov threshold
med_cov=5 #threshold of median coverage of a single CpG site across samples
path_to_output = "/path/to/output"

#------------------------------------------
# ########## READ IN BSSEQ RAW DATA ####
#------------------------------------------
cayo=readRDS(path_bsseq_file)

# Simplify sample names
colnames(cayo) <- gsub(".CpG_report.merged_CpG_evidence.cov.gz", "",
                       str_split_i(colnames(cayo), "/", 6))

# Read prepared metadata
metadata_lid <- read.table(path_to_metadata, sep = "\t", header = TRUE)

# Subset bsseq object by colnames present in the metadata lid_pid
cayo_red <- cayo[,metadata_lid$lid_pid]

# Check that metadata and methylation data match. Otherwise if condition is FALSE stop the script
if (length(metadata_lid$lid_pid) != length(colnames(cayo_red)) || !all.equal(metadata_lid$lid_pid, colnames(cayo_red))) {
  stop("STOP metadata and samples info are not equal")
}

rm(cayo)

# Split bsseq per chr
cayo_red_list <- parallel::mclapply(chrs,function(x){
  #select one chr at a time
  chr <- chrSelectBSseq(cayo_red, seqnames = x, order = TRUE)
  return(chr)
  }, mc.cores =22)




tissues = c(unique(metadata_lid$grantparent_tissueType))
tissue_oi = tissues

cayo_red_filtered_list=mclapply(chrs, filter_CpGs,
                                bsseq_list=cayo_red_list,metadata=metadata_lid,
                                tissues=tissues, min_cov=min_cov,
                                p_sample=p_sample, med_cov=med_cov,
                                mc.cores=22)

# Add a percent methylation matrix on CpG sites
cayo_red_filtered_list=parallel::mclapply(
  cayo_red_filtered_list,function(x){
    pmeth=getCoverage(x,type="M")/getCoverage(x,type="Cov")
    obj=x
    obj@assays@data@listData[["pmeth"]]=pmeth
    return(obj)
  },mc.cores = 22)

# Create output directory if not already done
if (!dir.exists(paths = path_to_output)) dir.create(path_to_output)

# Save the list format
save(cayo_red_filtered_list, file = paste0(path_to_output,"/cayo_red_CpG_filtered.RData"))

rm(cayo_red)

# Below we define regions from the BSSeq objects based on a maximum gap between consecutive CpGs

# Create output directory for batches of regions processed sequentially
if (!dir.exists(paths = paste0(path_to_output,"/region_batches"))) dir.create(paste0(path_to_output,"/region_batches"))

#-----------------------------------------------------------------------------
# Generating the data
#-----------------------------------------------------------------------------

# Define the regions on each chr across all samples
regions_list <- parallel::mclapply(cayo_red_filtered_list,
                                   getting_regions,
                                   gap_region = gap_region,
                                   mc.cores = num_cores)

# number of regions generated
my_list = lapply(regions_list,nrow)
sum_tot = 0
for(i in my_list){
  sum_tot <- sum_tot + i[[1]]
}

if (!dir.exists(paths = paste0(path_to_output,"/regions_batches/",gap_region))) dir.create(paste0(path_to_output,"/regions_batches/",gap_region))

# Write out files with coverage, methylation, and %meth for batches of regions
parallel::mclapply(1:length(cayo_red_filtered_list),
                   function(x) {
                     getting_regions_coverage(
                       cayo_meth = cayo_red_filtered_list[[x]],
                       regions = regions_list[[x]]
                     )
                   },
                   mc.cores = 22)

# List all files created above in the folder at the current gap_region
file_list <- list.files(path = paste0(path_to_output,"/regions_batches/",gap_region), pattern = "_chr(X|\\d+)_", full.names = TRUE)

# Group files by the common _chr_ identifier
grouped_files <- split(file_list, gsub(".*_chr(X|\\d+)_batch.*", "\\1", file_list))

# Load and combine data for each group of files to get back one entry per chromosome in a list
result_list <- map(grouped_files, function(files) {
  # Load data from each file
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

  # Combine specific columns from each file
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

# Erase the folders containing the temporary batch files to always start fresh.
unlink(paste0(path_to_output,"/regions_batches/", gap_region), recursive = TRUE)


# Reorder result_list by chromosome
# Identify numeric names
numeric_names <- grep("^\\d+$", names(result_list), value = TRUE)
# Sort numeric names
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

#Calculate descriptive stats on the regions
final_result <- amend_regions(final_result)
final_result <- get_summary_regions(final_result)

# Export the data generated
save(final_result,result_list, file = paste0(path_to_ouput,"/Regions_",gap_region,"_",length(tissue_oi),"T.RData"))

rm(final_result)

# Filter regions which are consistently hypo/hyper methylated across tissues
cayo_red_list_methfiltered <- mclapply(chrs, compute_and_filter_regions,
                                       bsseq_list = result_list,
                                       metadata=metadata_lid, tissues = tissue_oi,
                                       hypomethylated = hypomethylated, hypermethylated = hypermethylated,
                                       mc.cores = num_cores)

# Sample filtering for region coverage QC

# Calculate the number of regions covered (x1) for each sample
CpGs_covered_in_sample <- mclapply(cayo_red_list_methfiltered,
                                   function(x) {
                                     cov_mat <- x$coverage
                                     col_tot <- colSums(cov_mat>0)
                                   }, mc.cores = 22)
                  
# Combine across chr
total_covered_perLID <- Reduce("+", CpGs_covered_in_sample)
total_covered_perLID <- left_join(data.frame(lid_pid=names(total_covered_perLID),
                                             nb_covered=total_covered_perLID),
                                  metadata_lid[,c("lid_pid","grantparent_tissueType")])

# Based on inspection, samples with missing coverage at >30% regions were discared

# Extract total number of regions
all_regions <- do.call(rbind, lapply(cayo_red_list_methfiltered, function(x) x$regions))
print(nrow(all_regions))

flawed_lids <- total_covered_perLID %>% filter(nb_covered < nrow(all_regions)*0.70) %>% pull(lid_pid)

metadata_lid = metadata_lid %>% filter(!lid_pid %in% flawed_lids)

# Create tissue specific datasets
for (tissue_level in tissue_oi){

 # Apply tissue_specific_coverage to each element of cayo_red_list_methfiltered
 result_tissue <- mclapply(1:length(cayo_red_list_methfiltered),
                           function(i) tissue_specific_coverage(
                             matrix_list=cayo_red_list_methfiltered[[i]],
                             tissue = tissue_level,
                             metadata_sample = metadata_lid,
                             in_at_least_Xperc_samples = in_at_least_Xperc_samples,
                             med_cov_filter = med_cov_filter
                           ),
                           mc.cores = num_cores)

 # for this tissue combine all chr
 all_chr_result <- list(
   regions = do.call(rbind, lapply(result_tissue, function(x) x$regions)),
   coverage = do.call(rbind, lapply(result_tissue, function(x) x$coverage)),
   methylation = do.call(rbind, lapply(result_tissue, function(x) x$methylation)),
   pmeth = do.call(rbind, lapply(result_tissue, function(x) x$pmeth))
 )

 # Calculate descriptive stats on the regions for that tissue
 all_chr_result <- amend_regions(all_chr_result)
 all_chr_result <- get_summary_regions(all_chr_result)

 if (!dir.exists(paths = paste0(path_to_output,"/tissues_meth/",tissue_level,"_meth"))) dir.create(paste0(path_to_output,"/tissues_meth/",tissue_level,"_meth"), recursive = TRUE)

 # Save separately the different tissues for memory
 saveRDS(all_chr_result, file = paste0(path_to_output,"/tissues_meth/",tissue_level,"_meth/Regions_pmeth_full_",tissue_level,"_",gap_region,"_",length(tissue_oi),"T.rds"))
}

rm(all_chr_result,grouped_files,regions_list,result_list,result_tissue) #clean workspace

write.table(metadata_lid, paste0(path_to_metadata,"/metadata_final.txt", sep = "\t", row.names = F, quote = F)
#### END OF THE SCRIPT
