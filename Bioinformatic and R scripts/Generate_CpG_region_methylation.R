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
in_at_least_Xperc_samples=0.25 #coverage above threshold in at least X% of the samples
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


##### FILTER CpGs

# FUNCTION to filter CpGs based on mininum coverage and median coverage.
filter_CpGs <- function(x, bsseq_list, metadata, tissues, min_cov, p_sample, med_cov) {

  # If x is "X", change it to "21"
  if (x == "X") {
    x <- "21"
  }

  # Get the corresponding BSseq object for the chromosome
  bsseq_object <- bsseq_list[[as.numeric(x)]]

  # Initialize a logical vector to track rows to keep (same length as number of CpGs in bsseq_object)
  overall_rows_to_keep <- rep(FALSE, nrow(bsseq_object))

  # Iterate over each tissue type to filter the BSseq object
  for (subgroup in tissues) {
    # now<-Sys.time()
    # Subset metadata to get samples for the current tissue
    group_samples <- metadata$lid_pid[metadata$grantparent_tissueType == subgroup]

    # Subset the BSseq object to the current tissue samples
    subgroupbsseq <- bsseq_object[, group_samples]

    # Extract coverage matrix for the current tissue subgroup
    coverage_matrix <- getCoverage(subgroupbsseq)

    # Calculate the number of samples where each CpG is covered with min_cov
    temp_cov <- (coverage_matrix >= min_cov) %>% DelayedMatrixStats::rowSums2(useNames=TRUE)

    # Identify CpGs that meet the sample coverage requirement
    rows_meeting_sample_coverage <- temp_cov >= (p_sample * ncol(subgroupbsseq))

    # Calculate the median coverage across samples for each CpG row (using rowMedians from DelayedMatrixStats significantly faster than apply)
    row_median <- DelayedMatrixStats::rowMedians(coverage_matrix)

    # Identify rows that meet the median coverage requirement
    rows_meeting_median_coverage <- row_median >= med_cov

    # Combine the sample and median coverage filters for the current subgroup
    final_rows_to_keep <- rows_meeting_sample_coverage & rows_meeting_median_coverage

    # Update the overall logical vector to track rows that should be kept (across tissues)
    overall_rows_to_keep <- overall_rows_to_keep | final_rows_to_keep
    # print(Sys.time()-now)
  }

  # After iterating through tissues, subset the original bsseq_object using the combined logical filter
  final_bsseq <- bsseq_object[overall_rows_to_keep, ]

  return(final_bsseq)
}

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

#-------------------------------------------------------------------------------
# FUNCTIONS
#-------------------------------------------------------------------------------
# FUNCTION to extract regions
getting_regions <- function(cayo_meth,
                            gap_region){
  
  # define regions with regionFinder3
  regions_oi <- bsseq:::regionFinder3(x = as.integer(rep(1,length(cayo_meth))),
                                      chr = as.character(GenomeInfoDb::seqnames(cayo_meth)),
                                      positions = BiocGenerics::start(cayo_meth), maxGap = gap_region, verbose = FALSE)[["up"]]
  
  # name regions
  regions_oi$region_name <- paste(regions_oi$chr, regions_oi$start,
                                  regions_oi$end, sep = "_")
  
  # when run in parallel, the automatically generated "Region_ID" are no longer unique because each chr restarts at Region_1, so replace the Region_ID with the unique region_name
  regions_oi$RegionID <- paste0("Region_",regions_oi$region_name)
  return(regions_oi)
}

# FUNCTION splits the regions in smaller batches and writes out coverage, methylation count, and percent methylation for each batch.
getting_regions_coverage <- function(cayo_meth,
                                     regions){
  
  # cut the region list into 5 batches
  batch_list <- split(regions, cut(seq(nrow(regions)), breaks = 5, labels = FALSE))
  
  # get coverage and save each batch
  for(i in seq_along(batch_list)){
    
    regions_oi <- batch_list[[i]]
    
    coverage_regions_oi <- bsseq::getCoverage(cayo_meth, type = "Cov", regions = regions_oi,
                                              what="perRegionTotal", withDimnames = TRUE)
    
    rownames(coverage_regions_oi) <- regions_oi$RegionID
    
    methylation_regions_oi <- bsseq::getCoverage(cayo_meth, type = "M", regions = regions_oi,
                                                 what="perRegionTotal", withDimnames = TRUE)
    
    rownames(methylation_regions_oi) <- regions_oi$RegionID
    
    perc_meth <- methylation_regions_oi / coverage_regions_oi
    
    return_ <- list(regions_oi, coverage_regions_oi, methylation_regions_oi, perc_meth)
    names(return_)=c("regions", "coverage", "methylation", "pmeth")
    
    write_rds(return_,file = paste0(path_to_output,"/regions_batches/",gap_region,"/regions_batch_chr",unique(regions_oi$chr),"_batch",i,".rds"))
    
    # clean workspace for memory
    rm(regions_oi,coverage_regions_oi,methylation_regions_oi,perc_meth,return_)
  }
}

# FUNCTION to filter across tissues on tissue-specific average percent methylation
compute_and_filter_regions <- function(x, bsseq_list, metadata, tissues, hypomethylated, hypermethylated) {
  
  if(x == "X"){
    x <- "21"
  }
  
  bsseq_object <- bsseq_list[[as.numeric(x)]]
  
  # Compute row means for each grouping factor tissue
  each_tissue_means <- sapply(tissues, function(subgroup) {
    
    group_samples <- metadata$lid_pid[metadata$grantparent_tissueType == subgroup]
    temp_pmeth <- bsseq_object$pmeth[,colnames(bsseq_object$pmeth) %in% group_samples]
    rowMeans(temp_pmeth, na.rm = TRUE)
    
  })
  
  # Remove sites which are not covered in any tissues and sites constantly hypo/hyper methylated across tissues
  row_remove <- apply(each_tissue_means, 1, function(x) {
    # If all values are NA, return TRUE
    if (all(is.na(x))) {
      return(TRUE)
    }
    # Otherwise ignore NA values when comparing with the thresholds
    decent_values <- x[!is.na(x)]
    all(decent_values < hypomethylated) | all(decent_values > hypermethylated)
  })
  
  # Remove the identified rows from BSseq object
  pmeth <- bsseq_object$pmeth[!row_remove,]
  coverage <- bsseq_object$coverage[!row_remove,]
  methylation <- bsseq_object$methylation[!row_remove,]
  regions <- bsseq_object$region %>% filter(!RegionID %in% names(row_remove[row_remove]))
  bsseq_filtered <- list(regions=regions, coverage=coverage, methylation=methylation, pmeth=pmeth)
  
  return(bsseq_filtered)
}

# FUNCTION to calculate descriptive statistics on the regions
amend_regions <- function(x){
  # Calculate width
  x$regions$width <- x$regions$end - x$regions$start +1
  
  # Calculate covMean (Mean coverage for each row)
  covMean <- rowMeans(x$coverage, na.rm = TRUE)
  
  # Calculate covMin (Minimum coverage for each row)
  covMin <- apply(x$coverage, 1, min, na.rm = TRUE)
  
  # Calculate covSD (Standard deviation of coverage for each row)
  covSD <- apply(x$coverage, 1, sd, na.rm = TRUE)
  
  # Calculate methMean (Mean %methylation for each row)
  methMean <- rowMeans(x$pmeth, na.rm = TRUE)
  
  # Calculate methSD (Standard deviation of %methylation for each row)
  methSD <- apply(x$pmeth, 1, sd, na.rm = TRUE)
  
  x$regions$covMean <- covMean
  x$regions$covMin <- covMin
  x$regions$covSD <- covSD
  x$regions$methMean <- methMean
  x$regions$methSD <- methSD
  
  return(x)
}

# FUNCTION to calculate descriptive stats across the regions stats (i.e., overall summary)
get_summary_regions <- function(x){
  
  # Apply summary function to each vector
  summary_regions <- sapply(list(
    "width" = x$regions$width,
    "n" = x$regions$n,
    "covMean" =  x$regions$covMean,
    "covMin" =  x$regions$covMin,
    "covSD" =  x$regions$covSD,
    "methMean" =  x$regions$methMean,
    "methSD" =  x$regions$methSD
  ), summary)
  summary_regions <- round(summary_regions, 2)
  summary_regions <- as.data.frame(t(summary_regions))
  summary_regions$parameter <- rownames(summary_regions)
  
  x$summary_regions <- summary_regions 
  return(x)
}

# FUNCTION to obtain tissue-specific coverage data on regions
tissue_specific_coverage <- function(matrix_list, tissue, metadata_sample, in_at_least_Xperc_samples, med_cov_filter) {

  # Extract 'lid_pid' values for the current tissue level
  tissue_lids <- metadata_sample$lid_pid[metadata_sample$grantparent_tissueType %in% tissue]

  # Filter cov matrix based on tissue_lids
  subset_cov <- matrix_list$coverage[, colnames(matrix_list$coverage) %in% tissue_lids]

  # Count number of samples with coverage per row
  nonzero_counts <- rowSums(subset_cov !=0)

  # Calculate median coverage per row
  row_median <- apply(subset_cov, 1, median)

  # Filter to keep rows with coverage in >=25% of samples AND median coverage >= med_cov_filter
  rows_to_keep <- nonzero_counts >= in_at_least_Xperc_samples*length(tissue_lids) & row_median >= med_cov_filter

  subset_cov <- subset_cov[rows_to_keep, ]

  # Get row names and column names after filtering pmeth
  rows_cov <- rownames(subset_cov)
  cols_cov <- colnames(subset_cov)

  # Filter pmeth and methylation matrices based on coverage
  subset_result <- lapply(matrix_list[c("methylation","pmeth")], function(df) {
    df[rows_cov, cols_cov, drop = FALSE]
  })

  # Now subset the regions based on coverage row names
  regions <- matrix_list$regions[matrix_list$regions$RegionID %in% rows_cov, ]

  # Return the filtered data
  return(list(coverage = subset_cov,
              methylation = subset_result$methylation,
              pmeth = subset_result$pmeth,
              regions = regions))
}
##############################################################################

#-----------------------------------------------------------------------------
# Generating the data
#-----------------------------------------------------------------------------

# Define the regions on each chr across all samples
regions_list <- parallel::mclapply(cayo_red_filtered_list,
                                   getting_regions,
                                   gap_region = gap_region,
                                   mc.cores = num_cores)

# To have a quick look at the number of regions generated
my_list = lapply(regions_list,nrow)
sum_tot = 0
for(i in my_list){
  sum_tot <- sum_tot + i[[1]]
}

if (!dir.exists(paths = paste0("/scratch/sbaptis7/Cayotissue_CpG_coverage/regions_batches/",gap_region))) dir.create(paste0("/scratch/sbaptis7/Cayotissue_CpG_coverage/regions_batches/",gap_region))

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
file_list <- list.files(path = paste0("/scratch/sbaptis7/Cayotissue_CpG_coverage/regions_batches/",gap_region), pattern = "_chr(X|\\d+)_", full.names = TRUE)

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
unlink(paste0("/scratch/sbaptis7/Cayotissue_CpG_coverage/regions_batches/", gap_region), recursive = TRUE)


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


# Use do.call and rbind to combine the results
final_result <- list(
  regions = do.call(rbind, lapply(result_list, function(x) x$regions)),
  coverage = do.call(rbind, lapply(result_list, function(x) x$coverage)),
  methylation = do.call(rbind, lapply(result_list, function(x) x$methylation)),
  pmeth = do.call(rbind, lapply(result_list, function(x) x$pmeth))
)


#Calculate descriptive stats on the regions
final_result <- amend_regions(final_result)
final_result <- get_summary_regions(final_result)

save(final_result,result_list, file = paste0("/scratch/sbaptis7/Cayotissue_CpG_coverage/Regions_",gap_region,"_",length(tissue_oi),"T.RData"))
# load("/scratch/sbaptis7/Cayotissue_CpG_coverage/Regions_1000_14T.RData")

# Comethyl plot on regions statistics
plotRegionStats(
  final_result$regions %>% select(-c("idxStart", "idxEnd","cluster")),
  # apply the function
  maxQuantile = 0.99, # this is useful to remove outlier from the visual
  bins = 50,
  histCol = "#132B43",
  lineCol = "red",
  nBreaks = 4,
  save = TRUE, #change to save/update the plot
  file = paste0("/scratch/sbaptis7/Cayotissue_CpG_coverage/RegionPlots/Region_Plots_",gap_region,"_",length(tissue_oi),"T.pdf"),
  width = 11,
  height = 8.5,
  verbose = TRUE
)

rm(final_result)


# Filter regions which are consistently hypo/hyper methylated across tissues
cayo_red_list_methfiltered <- mclapply(chrs, compute_and_filter_regions,
                                       bsseq_list = result_list,
                                       metadata=metadata_lid, tissues = tissue_oi,
                                       hypomethylated = hypomethylated, hypermethylated = hypermethylated,
                                       mc.cores = num_cores)

saveRDS(cayo_red_list_methfiltered,"/scratch/sbaptis7/Cayotissue_CpG_coverage/temp_RegionsData.RDS")

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

# After inspection (cf script	3z2_Plot_CpG_regions_tissue_specific_cov.R)
# we decided to remove samples with missing coverage in >30% regions and this regardless of tissue.

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

 if (!dir.exists(paths = paste0("/scratch/sbaptis7/Cayotissue_CpG_coverage/tissues_meth/",tissue_level,"_meth"))) dir.create(paste0("/scratch/sbaptis7/Cayotissue_CpG_coverage/tissues_meth/",tissue_level,"_meth"), recursive = TRUE)

 # Save separately the different tissues for memory
 saveRDS(all_chr_result, file = paste0("/scratch/sbaptis7/Cayotissue_CpG_coverage/tissues_meth/",tissue_level,"_meth/Regions_pmeth_full_",tissue_level,"_",gap_region,"_",length(tissue_oi),"T.rds"))
}

rm(all_chr_result,grouped_files,regions_list,result_list,result_tissue) #clean workspace

write.table(metadata_lid, "/scratch/sbaptis7/Cayo_meth_metadata/metadata_final_lidpids_Nov24.txt", sep = "\t", row.names = F, quote = F)
#### END OF THE SCRIPT
