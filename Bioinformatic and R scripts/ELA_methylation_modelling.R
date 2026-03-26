# ==============================================================================
# Differential Methylation Analysis with PQLSEQ
# -------------------------------------------------------------
# SLURM submission command (run as array over tissues):
# ==============================================================================

# === Clear workspace ===
rm(list = ls())

# === Load libraries ===
library_list <- c("tidyverse", "PQLseq", "qvalue")
lapply(library_list, require, character.only = TRUE)

# === Paths ===

base_path <- "/path/to/project"  # <-- Define this path only once

metadata_path <- file.path(base_path,"metadata","multitissue_metadata.txt")
meth_dir <- file.path(base_path, "tissues_meth")
kinship_path <- file.path(base_path, "metadata", "wgs_kinmat.rds")
output_dir <- file.path(base_path, "PQLSEQ")

rds_path <- file.path(meth_dir, paste0(tissue, "_meth"), paste0("Regions_pmeth_full_", tissue, "_1000_14T.rds"))

# === Load data ===

metadata_lid = read.csv(file.path(metadata_path,"ELA_metadata.csv"), header = TRUE) %>%
  mutate(AnimalID.text= paste("'", monkey_id, sep="")) %>%
  filter(lid_pid != "LID_109490_PID_10416") %>%
  filter(age_at_sampling > 2.9) %>%
  filter(monkey_id != "22H") %>%
  mutate(percent_unique = unique/reads)

ELA_data = read.csv(file.path(metadata_path, "ELA_data.csv"))

meta_ELA<-merge(metadata_lid, ELA_data, by="AnimalID.text")

kinmat<-readRDS(kinship_path)

# === Settings ===

in_args <- commandArgs(trailingOnly = TRUE)
tissue=in_args
print(tissue)

plots <- list()

ELAvariables <- c("RankIndex01", "q_kinBirth_sub", "MomAllLossType", "q_grp_sub", "PrimpIndex", "CompetingSibIndex", "CumlQuartIndex6_4_sub")

# load data
r <- readRDS(rds_path)

# Meth count
counts <- r$methylation
  
# Total count (cov)
cov <- r$coverage

# filter counts for those in metadata
countsfiltered <- counts[,colnames(counts) %in% metadata_lid$lid_pid]
covfiltered <- cov[,colnames(cov) %in% metadata_lid$lid_pid]
  
# order the same
covordered<-covfiltered[,colnames(countsfiltered)]

# filter metadata for IDs in that tissue & included in kin matrix (all should be included in kin matrix)
rownames(meta_ELA)<- meta_ELA$lid_pid
subset_metadata_tissue<- meta_ELA[colnames(covordered),]
subset_metadata_kin<-subset_metadata_tissue[subset_metadata_tissue$monkey_id %in% rownames(kinmat),]

# ---------------------------------
# === RUN PQLSEQ MODELLING ===
# ---------------------------------

for(ELA in ELAvariables){
 tryCatch({
    #Filter metadata for IDs with that ELA variable
    ELAsubsetmeta <- subset_metadata_kin[complete.cases(subset_metadata_kin[,ELA]), ]
    predictor<-ELAsubsetmeta[,ELA]
    
    #Filter counts for IDs with that ELA variable
    ELAsubsetcounts<- countsfiltered[,rownames(ELAsubsetmeta)]
    ELAsubsetcov <- covordered[,rownames(ELAsubsetmeta)]
    
  # build covariates
  age <- ELAsubsetmeta$age_at_sampling
  sex <- ELAsubsetmeta$individual_sex
  sex_bino <-as.numeric(ifelse(sex=="F", 0, 1))
  percent_unique <- ELAsubsetmeta$percent_unique
  rank <- ELAsubsetmeta$compiled_rank
  rank_num <-as.numeric(ifelse(rank=="H", 0, ifelse(rank == "M", 0.5, 1)))
  group <- ELAsubsetmeta$group.x
  group_bino <-as.numeric(ifelse(group=="HH", 0, 1))

  if (tissue == "testis" | tissue == "ovaries") {
    covariates<- as.matrix(cbind(age, percent_unique, rank_num, group_bino))
} else {
  covariates<- as.matrix(cbind(age, percent_unique, sex_bino, rank_num, group_bino))
}
  
  # subset kinship matrix, order same as metadata
dnam_kinship_new <- kinmat[ELAsubsetmeta$monkey_id,ELAsubsetmeta$monkey_id]

  #remove sites without coverage for either level of ELA or group

  sites_to_remove <- vector("character", length = 0)
  
  ELA0 <- ELAsubsetmeta[ELAsubsetmeta[[ELA]] == "0", ]
cov_ELA<-as.data.frame(ELAsubsetcov[,colnames(ELAsubsetcov) %in% ELA0$lid_pid])
cov_ELA0<-cov_ELA[rowSums(cov_ELA) == 0,]

  ELA1 <- ELAsubsetmeta[ELAsubsetmeta[[ELA]] == "1", ]
cov_ELA<-as.data.frame(ELAsubsetcov[,colnames(ELAsubsetcov) %in% ELA1$lid_pid])
cov_ELA1<-cov_ELA[rowSums(cov_ELA) == 0,]

group0 <- ELAsubsetmeta[ELAsubsetmeta[,"group.x"] == "HH", ]
cov_group<-as.data.frame(ELAsubsetcov[,colnames(ELAsubsetcov) %in% group0$lid_pid])
cov_group0<-cov_ELA[rowSums(cov_group) == 0,]

group1 <- ELAsubsetmeta[ELAsubsetmeta[,"group.x"] == "KK", ]
cov_group<-as.data.frame(ELAsubsetcov[,colnames(ELAsubsetcov) %in% group1$lid_pid])
cov_group1<-cov_ELA[rowSums(cov_group) == 0,]


sites_to_remove <-c(rownames(cov_ELA0), rownames(cov_ELA1), rownames(cov_group0), rownames(cov_group1)); length(sites_to_remove)
  
filtered_counts <- ELAsubsetcounts[!rownames(ELAsubsetcounts) %in% sites_to_remove, ]
filtered_cov <- ELAsubsetcov[!rownames(ELAsubsetcov) %in% sites_to_remove, ]

# model the data

print("starting modeling")

   fit = PQLseq::pqlseq(RawCountDataSet=filtered_counts,
                     Phenotypes=predictor,
                     RelatednessMatrix=dnam_kinship_new,
                     LibSize=filtered_cov,
                     Covariates = covariates,
                     fit.model="BMM",
                     numCore = 40)

print("finished modeling")

# filter for converged sites and calculate qvalue
converged<- subset(fit, converged == "TRUE")
converged$qvalue <-qvalue(converged[,"pvalue"])$qvalues
   
# save resultss
string1 <- gsub("XXX", tissue, file.path(output_dir,"XXX_000_pqlseq.rds"))
string2 <- gsub("000", ELA, string1)
saveRDS(converged, file= string2)
   
   print(paste("finished",tissue, ELA, sep=" "))
 
  }, error = function(e) {
  print(paste("Error occurred:", e$message))
  }) 
}
