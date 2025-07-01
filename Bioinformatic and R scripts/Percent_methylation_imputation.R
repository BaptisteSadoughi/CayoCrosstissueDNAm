#!/usr/bin/env /packages/apps/spack/21/opt/spack/linux-rocky8-zen3/gcc-12.1.0/r-4.4.0-4yi4nm4foi7jsbczjxvv77uq7adnzb67/bin/Rscript

# Imputation section performed as a job submitted with :
# sbatch --cpus-per-task=22 --mem-per-cpu=1G -p general -q public -t 0-02:00:00 --array=1-(number of tissue) /path/to/thiscrtipt

# Load the pmeth matrices from all the tissues_meth folders RDS objects
# Assembles the matrices in a list, and applied a filtering to only
# keep regions with full coverage for the clocks.
# Imputate missing values for analyses and graphs that cannot deal with missing data.

rm(list = ls())

library_list <- c("parallel", "dplyr", "tidyr","methyLImp2","BiocParallel")

lapply(library_list, require, character.only = TRUE)

# Set the directory where your folders are located
main_folder <- "/path/to/tissues_meth"

# Define subfolders to consider
tissue_oi <- c("lung","kidney","heart","liver","spleen","skeletal_muscle","adrenal", "thyroid", "thymus","whole_blood","omental_at","pituitary","testis","ovaries")

sub_folders <- paste0(tissue_oi,"_meth")

# Construct full paths for additional folders
full_folders <- file.path(main_folder, sub_folders)

# Define a function to extract pmeth from each file
extract_pmeth <- function(file_path) {
  # Read the file
  data <- readRDS(file_path)
  
  # Extract pmeth
  pmeth <- data$pmeth
  
  return(pmeth)
}

# Define a function to filter files based on the pattern and extract pmeth
extract_pmeth_from_files <- function(folder_path, pattern) {
  
  # Get the list of files in the folder
  files_ <- list.files(folder_path, full.names = TRUE)
  
  # Filter files based on the pattern
  file_ <- files_[grep(pattern, files_)]
  
  # Extract pmeth from the target file
  pmeth <- extract_pmeth(file_)
  
  return(pmeth)
}

# Initialize a list to store pmeth from each folder
pmeth_list <- list()

pattern_ <- "_1000_14T"

# Loop through each folder

for (folder in full_folders) {

    # Get the folder name
  folder_name <- gsub("_meth$", "", basename(folder))
  
  # Extract pmeth from files named with the pattern
  pmeth_list[[folder_name]] <- extract_pmeth_from_files(folder, pattern = pattern_)
}

# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------

# #### SAVE FULL COVERAGE MATRIX WITH NO IMPUTATION
# 
# Keep regions with full coverage for the elastic net regression.
pmeth_list_full <- lapply(pmeth_list, function(mat){
  mat <- mat[complete.cases(mat),,drop = FALSE]
  return(mat)
})
# 
saveRDS(pmeth_list_full, file = paste0("/path/to/full_matrices/Regions_full_pmeth",length(tissue_oi),".rds"))

#-------------------------
##### DATA IMPUTATION
#-------------------------

# library(methyLImp2)
# library(BiocParallel)
# For details on methyLimp2 see
# https://scholar.google.com/citations?view_op=view_citation&hl=en&user=YXzVOS8AAAAJ&citation_for_view=YXzVOS8AAAAJ:qjMakFHDy7sC
# https://academic.oup.com/bioinformatics/article/35/19/3786/5364016?login=false

SAMP <- Sys.getenv("SLURM_ARRAY_TASK_ID")
SAMP <- as.integer(SAMP)

# Define tissue for the array
tissue <- tissue_oi[SAMP]
print(tissue)

# Isolate tissue matrix
pmeth <- pmeth_list[[tissue]]
rm(pmeth_list)

# Omit X which displays strong sex-differences
pmeth=pmeth[!grepl("Region_X",rownames(pmeth)),]

# Save proportion of missing sites per sample --> sites with >30% missing values are probably best discarded.
p_sample_missing_sites <- round(apply(pmeth, 2, function(x) sum(is.na(x))/length(x))*100,2)
p_uncomplete_sites <- round(sum(apply(pmeth, 1, function(x) any(is.na(x))))/nrow(pmeth),3)
if (!dir.exists("/path/to/combined_CpG_Regions/imputed/prop_missing_sites")) dir.create("/path/to/combined_CpG_Regions/imputed/prop_missing_sites", recursive = TRUE)
write.table(data.frame(lid_pid=names(p_sample_missing_sites),prop_missing_sites=p_sample_missing_sites,p_uncomplete_sites=p_uncomplete_sites), file=paste0("/path/to/combined_CpG_Regions/imputed/prop_missing_sites/prop_missing_sites_",tissue), sep = "\t", row.names = FALSE, quote = FALSE)

# Create the annotation file to split by chromosome
annotation <- data.frame("cpg"=rownames(pmeth))
annotation <- annotation %>% mutate(site=cpg, strand="+") %>% separate(., site, into = c("site","chr","start","end"), sep = "_", fill = "right") %>% select(-site)

# Convert to numeric
annotation[c("chr","start","end")] <- lapply(annotation[c("chr","start","end")], as.numeric)

# Order the rows based on chromosome number and start position
order <- order(annotation$chr, annotation$start)
annotation <- annotation[order,]
pmeth=pmeth[order,]

imputed_data <- methyLImp2(input=t(pmeth),
                           type = "user",
                           annotation = annotation,
                           # groups = tissue_groups
                           BPPARAM = SnowParam(workers = 20,exportglobals = FALSE),
                           overwrite_res = TRUE)

saveRDS(t(imputed_data), paste0("/path/to/combined_CpG_Regions/imputed/pmeth_imp_methyLimp_",tissue,".rds"))

print(paste("COMPLETED",tissue))
