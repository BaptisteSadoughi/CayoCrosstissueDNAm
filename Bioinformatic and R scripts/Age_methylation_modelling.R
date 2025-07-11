#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript

### Differential Methylation Analysis with PQLSEQ

#submit as a job using: sbatch --cpus-per-task=40 --mem=100G -p general -q public -t 0-4 --array=1-14 /path/to/thiscript

rm(list = ls())

library_list <- c("bsseq","BiocGenerics","GenomicRanges",
                  "GenomicFeatures","tidyverse","PQLseq","purrr","svglite","qvalue")

lapply(library_list, require, character.only = TRUE)

n.cores <- 40
print(n.cores)

# Define tissues of interest => determines files to load
tissue_oi <- c("liver","whole_blood","spleen","omental_at","heart","ovaries","testis","kidney","lung","adrenal","thymus","thyroid","pituitary","skeletal_muscle")

SAMP <- Sys.getenv("SLURM_ARRAY_TASK_ID")
SAMP <- as.integer(SAMP)

tissue = tissue_oi[SAMP]

#------------------------------
#### LOADING DATA AND METADATA
#------------------------------

# Load metadata
metadata_lid = read.table("/path/to/metadata.txt", sep = "\t", header = TRUE) %>% filter(lid_pid != "LID_109490_PID_10416")
metadata_lid$percent_unique<-metadata_lid$unique/metadata_lid$reads

# define pattern
filename <- gsub("XXX",tissue,"/path/to/tissues_meth/XXX_meth/Regions_pmeth_full_XXX_1000_14T.rds")
  
# load data
r <- readRDS(filename)
  
# subset metadata to matching samples
subset_metadata <- metadata_lid[metadata_lid$lid_pid %in% colnames(r[["coverage"]]),]
r$metadata <- subset_metadata
  
# name the object in the environment
assign(tissue, r)

rm(r,subset_metadata)

# Remove flawed skeletal muscle
if(exists("skeletal_muscle")) {
skeletal_muscle$coverage <- skeletal_muscle$coverage[,-which(colnames(skeletal_muscle$coverage) == "LID_109490_PID_10416")]
skeletal_muscle$methylation <- skeletal_muscle$methylation[,-which(colnames(skeletal_muscle$methylation) == "LID_109490_PID_10416")]
skeletal_muscle$pmeth <- skeletal_muscle$pmeth[,-which(colnames(skeletal_muscle$pmeth) == "LID_109490_PID_10416")]
}

# load kinship
kinmat<-readRDS("/path/to/kinmats/wgs_kinmat.rds")

#------------------
#### PQLSEQ
#------------------

print(tissue)
  
# load tissue
tissue_data <- get(tissue)
  
methcount = tissue_data$methylation
coverage = tissue_data$coverage

# subset metadata to the coverage matrix
metadata_tissue <- tissue_data$metadata 
rownames(metadata_tissue) <- metadata_tissue$lid_pid
subset_metadata_tissue <- metadata_tissue[colnames(coverage),]
  
subset_metadata_tissue <- subset_metadata_tissue %>% filter(age_at_sampling>=2.9)
  
# subset coverage and methylation matrix
coverage <- coverage[,colnames(coverage) %in% subset_metadata_tissue$lid_pid]
methcount <- methcount[,colnames(methcount) %in% subset_metadata_tissue$lid_pid]
  
# match columns order accross all data
subset_metadata_tissue <- subset_metadata_tissue[match(subset_metadata_tissue$lid_pid, colnames(coverage)),]
  
# subset kinship matrix
dnam_kinship_new <- kinmat[subset_metadata_tissue$monkey_id,subset_metadata_tissue$monkey_id]
  
# Check that metadata and methylation data match. Otherwise if condition is FALSE stop the script
  if (!identical(rownames(subset_metadata_tissue),colnames(coverage))| !identical(rownames(subset_metadata_tissue),colnames(methcount)) |
      !identical(colnames(dnam_kinship_new),subset_metadata_tissue$monkey_id) | !identical(rownames(dnam_kinship_new),subset_metadata_tissue$monkey_id)){
    stop("STOP there are some mismatches across matrices")
  }
  
  # extract predictors and covariates
  age <- subset_metadata_tissue$age_at_sampling
  sex_bino <- subset_metadata_tissue %>% mutate(sex_bino = as.numeric(ifelse(individual_sex == "F",0,1))) %>% pull(sex_bino)
  percent_unique <- subset_metadata_tissue$percent_unique
  group_bino <- subset_metadata_tissue %>% mutate(group_bino = as.numeric(ifelse(subset_metadata_tissue$group=="HH", 0, 1))) %>% pull(group_bino)

  # build covariates
  if (tissue == "testis" | tissue == "ovaries") {
    covariates<- as.matrix(cbind(percent_unique, group_bino))
  } else {
    covariates<- as.matrix(cbind(percent_unique, sex_bino, group_bino))
  }
  
  # remove sites with 0 coverage in a group (happens for gonads)

  group0 <- subset_metadata_tissue[subset_metadata_tissue[,"group"] == "HH", ]
  cov_group <- as.data.frame(coverage[,colnames(coverage) %in% group0$lid_pid])
  cov_group0 <- coverage[rowSums(cov_group) == 0,]
   
  group1 <- subset_metadata_tissue[subset_metadata_tissue[,"group"] == "KK", ]
  cov_group<-as.data.frame(coverage[,colnames(coverage) %in% group1$lid_pid])
  cov_group1<-coverage[rowSums(cov_group) == 0,]
  
  sites_to_remove <-c(rownames(cov_group0), rownames(cov_group1)); length(sites_to_remove)
  
  methcount <- methcount[!rownames(methcount) %in% sites_to_remove, ]
  coverage <- coverage[!rownames(coverage) %in% sites_to_remove, ]
  
  # run PQLseq
  fit = PQLseq::pqlseq(RawCountDataSet=methcount,
                     Phenotypes=age,
                     RelatednessMatrix=dnam_kinship_new,
                     LibSize=coverage,
                     Covariates = covariates,
                     fit.model="BMM",
                     numCore = 40)
  
  # Prep bed file
  fit$chr <- str_split_i(rownames(fit), i=2, pattern = "_")
  fit$start <- str_split_i(rownames(fit), i=3, pattern = "_")
  fit$end <- str_split_i(rownames(fit), i=4, pattern = "_")
  fit <- fit %>% mutate(sites = rownames(.))

  # keep regions which converged
  fit_conv <- fit %>% filter(converged==TRUE)

  # # Check pvalue uniform distribution before using qvalue
  # # https://www.bioconductor.org/packages/3.15/bioc/vignettes/qvalue/inst/doc/qvalue.pdf
  if(!dir.exists("/path/to/PQLSEQ/diagnostics")) dir.create("/path/to/PQLSEQ/diagnostics")
  
  ggplot(data = fit_conv, aes(x=pvalue))+geom_histogram(bins = 20)+theme_minimal()+theme(axis.title = element_text(size=14),axis.text = element_text(size=14))
  ggsave(filename = paste0("/path/to/PQLSEQ/diagnostics/",tissue,"_pvalhist.png"), width = 6, height = 4, dpi = 300)
  
  # Calculate qvalues and export diagnostic plots
  qresults <- qvalue::qvalue(fit_conv$pvalue)
  fit_conv$qval <- qresults$qvalues
  
  # Calculate adjusted BH pvalues
  fit_conv$padj <- stats::p.adjust(fit_conv$pvalue, method = "BH")
  
  # name the output
  assign(paste0("pqlseq_",tissue),fit_conv)
  
  # export the bed files for converged = TRUE regions
  write.table(fit_conv,paste0("/path/to/PQLSEQ/bed_",tissue,".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
