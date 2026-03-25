# ==============================================================================
# Define tissue-specific differentially methylated regions
# -----------------------------------
# The aim is to identify loci exhibiting consistent differences in DNAm levels with other tissues.
# ==============================================================================

rm(list=ls())

# === Load libraries

library_list <- c("dplyr","qvalue","PQLseq")
lapply(library_list, require, character.only = TRUE)

# === Paths ===

base_path <- "/path/to/your/directory" # <-- set the path once

metadata_path <- file.path(base_path,"metadata", "multitissue_metadata.txt")
kinship_path <- file.path(base_path, "metadata", "wgs_kinmat.rds")
alltissues <- c("liver", "omental_at", "spleen", "kidney", "lung", "heart", "skeletal_muscle", "adrenal", "pituitary", "thymus", "thyroid", "whole_blood")

# === Load data ===

metadata = read.table(metadata_path, sep = "\t", header = TRUE) %>%
  mutate(AnimalID.text= paste("'", monkey_id, sep="")) %>%
  filter(age_at_sampling > 2.9)
metadata$percent_unique<-metadata$unique/metadata$reads
kinmat<-readRDS(kinship_path)
kinmat_adults<-kinmat[rownames(kinmat) %in% metadata$monkey_id, colnames(kinmat) %in% metadata$monkey_id]

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
focaltissue=in_args
message("Focal tissue: ", focaltissue)

compare_tissues <- alltissues[alltissues != focaltissue]

for(compare_tissue in compare_tissues){

message("Compare tissue: ", compare_tissue)

#import cov and meth files and filter for shared sites
  
 focal <- readRDS(gsub("XXX",focaltissue, file.path(base_path, "tissues_meth", "XXX_meth", 
                                                     "Regions_pmeth_full_XXX_1000_14T.rds")))
  meth<-focal$methylation
  meth_focal<-meth[shared_rows,]

  cov<-focal$coverage
  cov_focal<-cov[shared_rows,]

  
 compare <- readRDS(gsub("XXX",compare_tissue, file.path(base_path, "tissues_meth", "XXX_meth", 
                                                          "Regions_pmeth_full_XXX_1000_14T.rds")))
  meth<-compare$methylation
  meth_compare<-meth[shared_rows,]

  cov<-compare$coverage
  cov_compare<-cov[shared_rows,]
  
  meth_both <-cbind(meth_focal, meth_compare)
  cov_both <- cbind(cov_focal, cov_compare)
  
  metadata_tissuespecific <- subset(metadata, tissue == focaltissue | tissue == compare_tissue)
  
  countsfiltered<-meth_both[,colnames(meth_both) %in% metadata_tissuespecific$lid_pid]
  covfiltered<-cov_both[,colnames(cov_both) %in% metadata_tissuespecific$lid_pid]

 #order the same
covordered<-covfiltered[,colnames(countsfiltered)]
stopifnot(identical(colnames(covfiltered), colnames(countsfiltered)))

# === Prepare Kinship Matrix ===

metadata_lid<-metadata[,c("monkey_id", "lid_pid")]
kinmatdf<- as.data.frame(kinmat_adults)
kinmatdf$monkey_id<-rownames(kinmat_adults)

kinmat_expand1<-left_join(kinmatdf, metadata_lid)
rownames(kinmat_expand1)<-kinmat_expand1$lid_pid
kinmat_expand1_tidy <- kinmat_expand1[, !(names(kinmat_expand1) %in% c("monkey_id", "lid_pid"))]
kinmat_expand1_tidy_t <- t(kinmat_expand1_tidy)


kinmat_t1<- as.data.frame(kinmat_expand1_tidy_t)
kinmat_t1$monkey_id<-rownames(kinmat_t1)
kinmat_expand2<-left_join(kinmat_t1, metadata_lid)
rownames(kinmat_expand2)<-kinmat_expand2$lid_pid
kinmat_expand2_tidy <- kinmat_expand2[, !(names(kinmat_expand2) %in% c("monkey_id", "lid_pid"))]

dnam_kinship_new<-kinmat_expand2_tidy[colnames(countsfiltered), colnames(countsfiltered)]

# === Build phenotype and covariates ===

rownames(metadata_tissuespecific)<-metadata_tissuespecific$lid_pid
metadata_tissuespecific_ordered<- metadata_tissuespecific[colnames(countsfiltered),]

metadata_lid_tissuebin<- metadata_tissuespecific_ordered %>%
  mutate(tissue_bin= ifelse(tissue == focaltissue, "1", "0"))

predictor <- metadata_lid_tissuebin$tissue_bin
  
  age <- metadata_lid_tissuebin$age_at_sampling
  sex <- metadata_lid_tissuebin$individual_sex
  sex_bino <-as.numeric(ifelse(metadata_lid_tissuebin$individual_sex=="F", 0, 1))
  percent_unique <- metadata_lid_tissuebin$percent_unique
  covariates<- as.matrix(cbind(age, percent_unique, sex_bino))

  print(length(which((colnames(dnam_kinship_new)==colnames(covordered))==F)))
  print(length(which((rownames(dnam_kinship_new)==colnames(countsfiltered))==F)))

# === Fit PQLseq Model ===
  
   fit = PQLseq::pqlseq(RawCountDataSet=countsfiltered,
                     Phenotypes=predictor,
                     RelatednessMatrix=dnam_kinship_new,
                     LibSize=covordered,
                     Covariates = covariates,
                     fit.model="BMM",
                     numCore = 40)

# === Post-processing === 
  
   converged<- subset(fit, converged == "TRUE")
   converged$qvalue <-qvalue(converged[,"pvalue"])$qvalues
  if (!dir.exists(tissue_dir <- file.path(base_path, "tissue_comparisons"))) {
  dir.create(tissue_dir)
}
   string1 <- gsub("XXX", focaltissue, file.path(base_path, "tissue_comparisons", "XXX_v_000_pairwise_pqlseq.rds"))
   string2 <- gsub("000", compare_tissue, string1)
   saveRDS(converged, file= string2)
}


# === Identify tissue-specific regions === 

results<-list()

for (focaltissue in alltissues){
  
  compare_tissues <- alltissues[alltissues != focaltissue]
  
  for (compare_tissue in compare_tissues){
    
    string1 <- gsub("XXX", focaltissue, "XXX_v_000_pairwise_pqlseq.rds")
    string2 <- gsub("000", compare_tissue, string1)
    result<-readRDS(string2)
    sig_results<-subset(result, qvalue < 0.1)
    betas<-as.data.frame(sig_results[,"beta"])
    rownames(betas)<-rownames(sig_results)
    colnames(betas)<-compare_tissue
    results[[compare_tissue]]<- betas 
  }
  
    #create dataframe of sig regions across all comparisons
    merged_df <- results %>%
      lapply(rownames_to_column, var = "rowname") %>%
      reduce(inner_join, by = "rowname") %>%
      column_to_rownames("rowname")
  
    filtered_df <- merged_df[apply(merged_df, 1, function(x) all(x < -0.2231436) | all(x > 0.1823216)),]
    
    filtered_df$mean<-apply(filtered_df, 1, mean)
    filtered_df$Region<-rownames(filtered_df)
    format<-filtered_df %>%
      separate(Region, c("reg", "chr", "start", "end")) %>%
      mutate(chr= paste0("chr", chr)) %>%
      mutate(region= paste(chr, start, end, sep="_")) %>%
      mutate(tissue= "testis") %>%
      mutate(marker= ifelse(mean < 0, "hypo", "hyper")) %>%
      dplyr::select(chr, start, end, region, tissue, mean, marker)
    
    write.csv(format, gsub("XXX", focaltissue, file.path(base_path, "tissue_comparisons","XXX_tissuespecific_regions.csv")), quote=F)
}
