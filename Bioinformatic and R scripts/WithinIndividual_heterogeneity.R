#######################################################
#
# WITHIN INDIVIDUAL HETEROGENEITY IN PREDICTED AGES
#
#######################################################

rm(list = ls())

library_list <- c("corrplot","dplyr","purrr","parallel","tidyverse", "DHARMa", "performance")

lapply(library_list, require, character.only=TRUE)

# Define ggplot theme upfront
palette <- c(RColorBrewer::brewer.pal(12, "Set3")[-2], RColorBrewer::brewer.pal(8, "Set2")[8], RColorBrewer::brewer.pal(12,"Paired")[1],RColorBrewer::brewer.pal(12, "Set3")[2])
my_theme <- theme_classic() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.background = element_blank())

combined_data = read.table("/scratch/sbaptis7/Cayo_clocks/tissue_predicted_agetransfo/DNAm_deviation_data.txt", sep="\t", header=TRUE)

# Load metadata
metadata_lid = read.table("/scratch/sbaptis7/Cayo_meth_metadata/metadata_final_lidpids_Nov24.txt", sep = "\t", header = TRUE)
metadata_lid = metadata_lid %>% filter(lid_pid != "LID_109490_PID_10416")
metadata_lid$percent_unique<-metadata_lid$unique/metadata_lid$reads

# Add metadata to the age deviations datasets
combined_data = merge(combined_data, metadata_lid[,c("lid_pid", "monkey_id", "individual_sex")])
