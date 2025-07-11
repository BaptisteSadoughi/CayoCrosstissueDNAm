## This script performs MASHr correction on estimated effect sizes for age and determines tissue shared and tissue-specific effects.

### Load the bed files from PQLSEQ

read_table_and_rownames <- function(tissue) {
  pqlseq <- read.table(paste("/path/to/PQLSEQ/bed_", tissue, ".txt", sep = ""),
                       sep = "\t",
                       header = TRUE)

  rownames(pqlseq) <- pqlseq$sites
  return(pqlseq)
}

data_list <- lapply(tissue_oi, read_table_and_rownames)

names(data_list) <- paste("pqlseq", tissue_oi, sep = "_")

list2env(data_list, envir = .GlobalEnv)

#### BUILD BETA AND SE MATRICES (for MASHr)
####
# Function to fill in matrix
fill_matrix <- function(column) {
  # get all rownames from all dataframes
  all_rownames <- unique(unlist(lapply(tissue_oi, function(x) row.names(get(paste("pqlseq", x, sep = "_"))))))

  # initialize matrix with NA values
  matrix_data <- matrix(NA, nrow = length(all_rownames), ncol = length(tissue_oi))
  rownames(matrix_data) <- all_rownames
  colnames(matrix_data) <- tissue_oi

  # loop through each dataframe and fill matrix_data
  for (tissue in tissue_oi) {
    df <- get(paste("pqlseq", tissue, sep = "_"))
    matrix_data[rownames(df), tissue] <- df[[column]]
  }

  return(matrix_data)
}

# Export matrices of beta estimates from PQLSEQ for MASH
matrix_beta <- fill_matrix("beta")
matrix_SE <- fill_matrix("se_beta")

saveRDS(list(matrix_beta=matrix_beta, matrix_SE = matrix_SE),"/path/to/PQLSEQ/pre_MASH_effects.RDS")

############################
#
# MASHr
#
############################

# Permorm MASH correction on estimated effect sizes calculated from PQLseq.
# The script follows the flow of tutorial vignettes available at https://stephenslab.github.io/mashr/articles/index.html

rm(list = ls())

library_list <- c("ashr", "mashr", "flashier", "svglite", "tidyverse")
lapply(library_list, require, character.only=TRUE)

# Load effect sizes and SE calculated with PQLSEQ
pqlseq_effects = readRDS("/path/to/PQLSEQ/pre_MASH_effects.RDS")

# Split gonads from other tissues because values are only available for one of the two sexes
pqlseq_effects$matrix_beta <- subset(pqlseq_effects$matrix_beta,select=-c(testis,ovaries))
pqlseq_effects$matrix_SE <- subset(pqlseq_effects$matrix_SE,select=-c(testis,ovaries))

# Work on fully covered sites
pqlseq_effects = lapply(pqlseq_effects, function(x) { x[complete.cases(x),]})

# Turn to mash format
data = mash_set_data(pqlseq_effects$matrix_beta, pqlseq_effects$matrix_SE)

##### Identify a subset of strong tests
tissue_files <- list.files("/path/to/PQLSEQ", pattern=".txt", full.names=TRUE)
strong_subset_list <- list()
for(file in tissue_files) {
  # Load the data from the .txt file into a data frame
  bed_data <- read.table(file, header = TRUE)
  
  # Filter rows for qval < 0.001
  filtered_data <- bed_data[bed_data$qval < quantile(bed_data$qval, 0.01),]
  
  # Extract tissue type from the file name
  tissue_name <- gsub("\\.txt$", "", basename(file))
  
  # Add the 'site' column of the filtered data for this tissue to the list
  strong_subset_list[[tissue_name]] <- filtered_data$site
}

strong_subset_list <- strong_subset_list[!names(strong_subset_list) %in% c("bed_testis", "bed_ovaries")]
strong_subset <- Reduce(union, strong_subset_list)

# we can only keep the sites which were fully covered
strong_subset <- intersect(rownames(pqlseq_effects$matrix_beta),strong_subset)

# Estimate background correlation in the data
Vhat = estimate_null_correlation_simple(data)

# Update mash objects with correlation structure
data.cor = mash_update_data(data, V=Vhat)
data.strong.cor = mash_set_data(pqlseq_effects$matrix_beta[strong_subset,], pqlseq_effects$matrix_SE[strong_subset,], V=Vhat)

# Investigate data-driven covariances on the strong set
U.pca = cov_pca(data.strong.cor, 5)
U.f = cov_flash(data.strong.cor)

# Run ED (Extreme Deconvolution) with data-driven covariances
U.ed = cov_ed(data.strong.cor, c(U.pca, U.f))

# canonical covariance
U.c = cov_canonical(data.cor)

# Estimate mixture proportions with MASH (we use both canonical and data-driven covariances)
m.all = mash(data.cor, Ulist = c(U.ed, U.c))
saveRDS(m.all, "/path/to/MASH/mash_object.RDS")

mash_results <- with(m.all$result, list(beta = PosteriorMean, SD = PosteriorSD, LFSR = lfsr))
saveRDS(mash_results, "/path/to/MASH/mash_estimates.RDS")

############################################################################
#
# Direction of change according to average percent methylation levels
#
#############################################################################
rm(list=ls())
library_list <- c("ashr", "mashr", "flashier", "corrplot","svglite","tidyverse","MatrixGenerics")
lapply(library_list, require, character.only=TRUE)

tissue_oi <- c("whole_blood","spleen","omental_at","heart","testis","ovaries","kidney","lung","adrenal","thymus","thyroid","pituitary","liver","skeletal_muscle")

#----------------------------------------------------------
# INSPECTING AND PLOTTING RESULTS
#----------------------------------------------------------
mash_results <- readRDS("/path/to/MASH/mash_estimates_June25.RDS")
m.all <- readRDS("/path/to/MASH/mash_object_June25.RDS")

# Load metadata
metadata_lid = read.table("/path/to/metadata.txt", sep = "\t", header = TRUE)
metadata_lid = metadata_lid %>% filter(lid_pid != "LID_109490_PID_10416")
metadata_lid$percent_unique<-metadata_lid$unique/metadata_lid$reads

# Remove infants
metadata_lid_ad <- metadata_lid %>% filter(age_at_sampling>2.9)

# Loop through the folders to load the tissue files and associated metadata
now <- Sys.time()
for(tissue in tissue_oi){
  
  # define pattern
  filename <- gsub("XXX",tissue,"/path/to/tissues_meth/XXX_meth/Regions_pmeth_full_XXX_1000_14T.rds")
  
  # load data
  r <- readRDS(filename)
  
  # subset metadata to matching samples #not that I use the full metadata because infants were included in the data generation process
  subset_metadata <- metadata_lid[metadata_lid$lid_pid %in% colnames(r[["coverage"]]),]
  r$metadata <- subset_metadata
  
  # name the object in the environment
  assign(tissue, r)
  
  #rm
  rm(r,subset_metadata)
}
print(Sys.time()-now)

# Remove flawed skeletal muscle
if(exists("skeletal_muscle")) {
  skeletal_muscle$coverage <- skeletal_muscle$coverage[,-which(colnames(skeletal_muscle$coverage) == "LID_109490_PID_10416")]
  skeletal_muscle$methylation <- skeletal_muscle$methylation[,-which(colnames(skeletal_muscle$methylation) == "LID_109490_PID_10416")]
  skeletal_muscle$pmeth <- skeletal_muscle$pmeth[,-which(colnames(skeletal_muscle$pmeth) == "LID_109490_PID_10416")]
}

# Extract results for further inspection (note that infants are not included in these results)
lfsr=as.data.frame(get_lfsr(m.all))
pm=as.data.frame(get_pm(m.all))
lfsr$region<-rownames(lfsr)

#set sig threshold 
threshold <- 0.05

coeff_intercept_list <- list()

for(i in tissue_oi[tissue_oi!=c("testis","ovaries")]){
  
  # get significant regions in the tissue
  tissue_lfsr <- lfsr[lfsr[,i] <threshold,]$region
  
  # extract their coefficients in the tissue
  tissue_beta <- pm[tissue_lfsr,i]
  
  # built dataframe
  tissue_beta <- data.frame(region=tissue_lfsr, beta = tissue_beta, tissue=i)
  
  tissue_beta$sign_beta <- ifelse(tissue_beta$beta >0, "pos", "neg")
  
  # Extract percent methylation for samples in that tissue
  tissue_meth <- get(i)
  
  # Extract percent methylation matrix and Remove infants
  tissue_pmeth <- tissue_meth$pmeth
  tissue_pmeth <- tissue_pmeth[,colnames(tissue_pmeth) %in% metadata_lid_ad$lid_pid]
  
  # calculate study sample mean percent methylation to flag hypo hyper sites
  pop_pmeth_intercept <- rowMeans2(tissue_pmeth[tissue_beta$region,], na.rm = TRUE)

  # calculate mean percent methylation in the subset of young individuals
  meta_young_lids <- metadata_lid_ad %>% filter(grantparent_tissueType == i & age_at_sampling <= 6) %>% pull(lid_pid)
  
  tissue_pmeth_young <- tissue_pmeth[tissue_beta$region,meta_young_lids]
  
  young_pmeth_intercept <- rowMeans2(tissue_pmeth_young, na.rm = TRUE)

  # check that vectors' length is still equal before merging
  if(length(young_pmeth_intercept)!=nrow(tissue_beta) | length(pop_pmeth_intercept)!=nrow(tissue_beta)) stop(paste0("length of regions means and coefficients differ for ",name))
  
  tissue_beta$young_intercept <- young_pmeth_intercept
  tissue_beta$pop_intercept <- pop_pmeth_intercept
  tissue_beta$intercept_high_low_50 <- ifelse(tissue_beta$young_intercept<0.5,"lower","higher")
  tissue_beta$nonvariable <- ifelse(tissue_beta$pop_intercep <0.1 | tissue_beta$pop_intercep >0.9 ,"nonvar","var")

  # save in list
  coeff_intercept_list[[i]] <- tissue_beta
}

# Export this list which can be useful for many other analysis
# saveRDS(coeff_intercept_list, "/path/to/MASH/bedfiles/tissue_age_associated_sites_list.rds") #exported with lfsr<0.05
                                       
############################ Test whether sites tend to trend towards 50 with age

coeff_intercept_list=readRDS("/path/to/MASH/bedfiles/tissue_age_associated_sites_list.rds")

### NOTE that the categories for intercept "higher" and "lower" (than 50% pmeth) were defined
### based on individuals <6 years of age to test whether sites starting in early life above or below
### threshold tend to move towards 50% with advancing age. Results using the average calculated on all individuals are qualitatively the same.

# Fisher's exact test per tissue
apply_fisher_test <- function(df){
  contingency_table <- with(df, table(sign_beta, intercept_high_low_50))
  test_result <- fisher.test(contingency_table)
  return(test_result)
}

# apply function to each dataframe in the list
test_results_list <- lapply(coeff_intercept_list, apply_fisher_test)

### Note the table is neg - pos x higher - lower

# function to extract components
extract_components <- function(test_result){
  # Fisher's test returns a list, we will extract the components we need
  estimate <- test_result$estimate
  lower_CI <- test_result$conf.int[1]
  upper_CI <- test_result$conf.int[2]
  pvalue <- test_result$p.value
  
  # return a data frame
  return(data.frame(estimate = estimate, lower_CI = lower_CI, upper_CI = upper_CI, pvalue = pvalue))
}

# apply the function to each Fisher's test result
components <- lapply(test_results_list, extract_components)

results_tests_df <- do.call(rbind, components)
results_tests_df$tissue <- rownames(results_tests_df)

# Repeat the process on variable sites for a given tissue to test the robustness of the relationship

# Fisher's exact test per tissue
apply_fisher_test_variable <- function(df){
  contingency_table <- with(df %>% filter(nonvariable == "var"), table(sign_beta, intercept_high_low_50))
  test_result <- fisher.test(contingency_table)
  return(test_result)
}

# apply function to each dataframe in the list
test_results_list_var <- lapply(coeff_intercept_list, apply_fisher_test_variable)

components_var <- lapply(test_results_list_var, extract_components)
results_tests_df_var <- do.call(rbind, components_var)
results_tests_df_var$tissue <- rownames(results_tests_df_var)
results_tests_df_var <- results_tests_df_var %>% dplyr::select(tissue, everything())

# Export results
results_tests_df_var <- results_tests_df_var %>%
  mutate_at(vars(-pvalue), ~if(is.numeric(.)) round(., digits=3) else .) %>%
  mutate(pvalue = ifelse(pvalue < 0.001, "<0.001", pvalue))

results_tests_df_var <- results_tests_df_var %>% 
  mutate(tissue = recode(tissue, "omental_at" = "omental adipose")) %>%
  mutate(tissue = recode(tissue, "skeletal_muscle" = "skeletal muscle")) %>% 
  mutate(tissue = recode(tissue, "whole_blood" = "whole blood"))

results_tests_df_var = results_tests_df_var %>% arrange(tissue)
rownames(results_tests_df_var)<-NULL

all_coeff_af <- do.call(rbind,coeff_intercept_list)

cat_young_plot <- ggplot(all_coeff_af %>% mutate(interaction = paste(sign_beta,intercept_high_low_50,"_")) %>% 
         mutate(interaction = recode(interaction,"neg higher _"="high decreasing",
                                            "neg lower _"="low decreasing",
                                            "pos higher _"="high increasing",
                                            "pos lower _"="low increasing")) %>% 
         mutate(tissue = recode(tissue, "omental_at" = "omental adipose",
                                "skeletal_muscle" = "skeletal muscle",
                                "whole_blood" = "whole blood")),
       aes(x=interaction, fill=interaction))+
  geom_bar()+
  labs(x="",title="A")+
  facet_wrap(~tissue, scales = "free_y")+
  theme_bw()+theme(legend.position = "none",
                   plot.title = element_text(size=28),
                   axis.title = element_text(size=12),
                   axis.text.x = element_text(size=12,angle=45, color="black",hjust=1,vjust=1),
                   axis.text.y = element_text(size=12, color="black"),
                   strip.text = element_text(size = 12, face = "bold"),
                   strip.background = element_rect(fill = "white"))
ggsave("/path/to/MASH/Figures/FigS7.png")

# Pearson's correlation per tissue
apply_cor_test_variable <- function(df){
  data_cor <- df %>% filter(nonvariable == "var") %>% dplyr::select(beta,young_intercept)
  test_result <- cor.test(data_cor$beta,data_cor$young_intercept, method="pearson")
  return(test_result)
}

# function to extract components
extract_corr_components <- function(test_result){
  # Fisher's test returns a list, we will extract the components we need
  pearson_corr <- test_result$estimate
  pvalue <- test_result$p.value
  # return a data frame
  return(data.frame(pearson_corr = pearson_corr, pvalue = pvalue))
}
# apply function to each dataframe in the list
cor_results_list_var <- lapply(coeff_intercept_list, apply_cor_test_variable)
components_corr <- lapply(cor_results_list_var, extract_corr_components)
results_corr_df_var <- do.call(rbind, components_corr)
results_corr_df_var$tissue <- rownames(results_corr_df_var)
results_corr_df_var <- results_corr_df_var %>% dplyr::select(tissue, everything())

results_corr_df_var <- results_corr_df_var %>%
  mutate(pearson_corr=round(pearson_corr, digits=3),pvalue = ifelse(pvalue < 0.001, "<0.001", pvalue))

results_corr_df_var <- results_corr_df_var %>% 
  mutate(tissue = recode(tissue, "omental_at" = "omental adipose")) %>%
  mutate(tissue = recode(tissue, "skeletal_muscle" = "skeletal muscle")) %>% 
  mutate(tissue = recode(tissue, "whole_blood" = "whole blood"))

results_corr_df_var = results_corr_df_var %>% arrange(tissue)

# write.table(results_corr_df_var %>% rename(`Pearson correlation`=pearson_corr,`p-value`=pvalue),
#             "/path/to/TableS6.csv",row.names = F,quote = F)

# Plot the three info: change with age, direction of tissue-specificity, and average pmeth in the population 
coeff_intercept_list= lapply(names(coeff_intercept_list),function(x) {
  tissue <- x
  data <- coeff_intercept_list[[x]]
  data$tissue <- tissue
  return(data)
  })

testy=do.call(rbind,coeff_intercept_list)

corr_plot<-ggplot(testy %>% 
                    mutate(tissue = recode(tissue, "omental_at" = "omental adipose",
                                           "skeletal_muscle" = "skeletal muscle",
                                           "whole_blood" = "whole blood")),
       aes(x=young_intercept,y=beta))+
  geom_point(alpha=0.3)+
  labs(x="Average percent methylation in young subjects", y = "Age-associated effect sizes", title="B")+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0.5)+
  facet_wrap(~tissue,scales = "free_y")+
  geom_smooth(method = "lm")+
  theme_bw()+theme(legend.position = "none",
                   plot.title = element_text(size=28),
                   axis.title = element_text(size=12),
                   axis.text.x = element_text(size=12,color="black",hjust=1,vjust=1),
                   axis.text.y = element_text(size=12, color="black"),
                   strip.text = element_text(size = 12, face = "bold"),
                   strip.background = element_rect(fill = "white"))
# ggsave("/path/to/Figures/FigS9.png", width=9.5, height=7.5, dpi=300)

ggpubr::ggarrange(cat_young_plot,corr_plot,nrow=2)
ggsave("/scratch/sbaptis7/MASH/Figures/Direction_of_agechange_meanYoung_June25.png", width=9.5, height=12.5, dpi=300)
ggsave("/scratch/sbaptis7/MASH/Figures/Direction_of_agechange_meanYoung_June25.pdf", width=9.5, height=12.5, dpi=300)
