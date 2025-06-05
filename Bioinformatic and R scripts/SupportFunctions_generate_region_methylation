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
