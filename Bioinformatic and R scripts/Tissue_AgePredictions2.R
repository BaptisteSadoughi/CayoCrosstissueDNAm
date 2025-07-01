#### After the Elastic Net Regression is done this script calculates the Mean Squared error for each values of Alpha.
#### The calculation is performed independently for each tissue and saved as a dataframe.

rm(list = ls())
library_list <- c("glmnet","tidyverse","parallel","svglite","dyplr","ggpubr")
lapply(library_list, require, character.only = TRUE)

tissue_oi <- c("whole_blood","pituitary","thymus","thyroid","testis","ovaries","lung","kidney","heart","liver","spleen","skeletal_muscle","adrenal","omental_at")

# Function to calculate mean squared error (MSE) for a given tissue and alpha
calculate_mse <- function(folder, alpha) {
  files <- list.files(path = folder, pattern = paste0(".*Combined_\\d+_alpha", alpha, "_.*\\.txt"), full.names = TRUE)
  
  if (length(files) == 0) {
    return(NA)
  }
  
  mse_values <- sapply(files, function(file) {
    data <- read.table(file, header = FALSE)
    
    # Calculate (predicted - testage)^2 for each file (i.e., sample)
    return((data[,4] - data[,3])^2)
  })
  # The mean is calculated across samples at a given alpha for a specific tissue (i.e., folder)
  return(mean(mse_values, na.rm = TRUE))
}

folders <- paste0("./predicted_agetransfo/",tissue_oi,"/")

## Loop through folders and alphas

for (i in seq_along(folders)) {
  folder <- folders[i]

  files <- list.files(path = folder, pattern = "^AgePred_\\d+_alpha\\d+\\.?\\d*_.+\\.txt$", full.names = TRUE)
  
  # Extract alpha values from filenames
  alphas <- unique(as.numeric(sub(".*_alpha(\\d+\\.?\\d*)_.*", "\\1", files)))
  
  # Create an empty dataframe to store results
  MSE_result <- data.frame(Alpha = numeric(), MSE = numeric())
  
  for (alpha in alphas) {
    mse <- calculate_mse(folder, alpha)
    MSE_result <- bind_rows(MSE_result, data.frame(Alpha = as.numeric(alpha), MSE = mse))
  }
  
  # Save MSE_result as a text file in the corresponding folder
  write.table(MSE_result, file.path(folder, paste0("/MSE_agepred_result.txt")), row.names = FALSE)
}

################################################################################
#### PREDICTED AGE
################################################################################

setwd("/path/to/predicted_agetransfo")

# List all subfolders
subfolders <- list.dirs(full.names = FALSE)
subfolders <- subfolders[-1]

# Initialize a list to store results
result_list <- list() #store min alpha value to inspect
combined_data_list <- list() #list of one dataframe per tissue

# Loop through each subfolder to load the data corresponding to the alpha minimizing MSE
# and concatenate into a list of dataframes

for(subfolder in subfolders){
  # Create the file path to MSE_result.txt in the current subfolder
  file_path <- file.path(subfolder, paste0("MSE_result.txt"))
  
  # Read the data into a data frame
  data <- read.table(file_path, header = TRUE)
  
  # Find the row where MSE is minimum
  min_row <- which.min(data$MSE)
  
  # Extract the Alpha value associated with the minimum MSE
  min_alpha <- data$Alpha[min_row]
  
  # Store the result in the list (simply to print at the end)
  result_list[[subfolder]] <- min_alpha
  
  # Create the pattern for files names with "_alpha" and the min_alpha
  # pattern <- paste(".*R", ".*_alpha",min_alpha,".*",sep="")
  if(floor(min_alpha) == min_alpha){
    pattern <- paste(".*Combined",".*_alpha", min_alpha, "_.*", sep="") #for min alpha 0 or 1
  } else {
    pattern <- paste(".*Combined",".*_alpha", min_alpha,".*",sep="")
  }
  
  # list all the files in the subfolder matching the pattern
  matching_files <- list.files(path = subfolder, pattern = pattern, full.names = TRUE)
  
  # Initalize empty dataframe to concatenate single files
  combined_data <- data.frame(monkey_id=as.character(),lid_pid=as.character(),
                              predicted_age=as.numeric(),chronological_age=as.numeric())
  
  for (file in matching_files) {
    data <- read.table(file, sep = " ", header = FALSE, #careful with the sep
                       col.names = c("monkey_id","lid_pid","chronological_age","predicted_age"),
                       colClasses = c("character", "character", "numeric", "numeric"))
    data$alpha <- min_alpha
    
    # Concatenate data
    combined_data <- rbind(combined_data,data)
  }
  combined_data_list[[subfolder]] <- combined_data
}

# Print the results
for (i in seq_along(subfolders)) {
  cat("Subfolder:", subfolders[i], "Min_MSE:", result_list[[i]], "\n")
}

### Save the min_alpha for later
subfolder_list <- list()
mse_list <- list()

for (i in seq_along(subfolders)) {
  subfolder_list[[i]] <- subfolders[i]
  mse_list[[i]] <- result_list[[i]]
}

# Create a dataframe
min_alphas <- data.frame(Subfolder = unlist(subfolder_list),
                         Min_MSE = unlist(mse_list))
# Save min alphas
write.table(min_alphas,"/path/to/tissue_predicted_agetransfo/min_alphas_age_transfo.txt",sep = "\t", quote = FALSE, row.names = FALSE)

# Concatenate data for ggplot
combined_data <- do.call(rbind, Map(cbind, combined_data_list, subfolder = names(combined_data_list)))

combined_data = combined_data %>% rename(tissue=subfolder)

#------------------------------
# CLOCK PERFORMANCE
#------------------------------

# Calculate correlation and MAE (median average error)
corr_and_mae <- function(df){
  coef <- round(cor.test(df$chronological_age,df$predicted_age, method = "pearson")$estimate,2)
  MAE <- round(median(abs(df$predicted_age - df$chronological_age)), digits = 2)
  return(c(coef=coef, MAE= MAE))}

# Apply to the list of prediction to calculate summary statistics     
summ.stats <- lapply(combined_data_list,corr_and_mae)

# Convert to data frame
summ.stats <- data.frame(
  tissue = names(summ.stats),
  coef.cor = sapply(summ.stats, function(x) x["coef.cor"]),
  MAE = sapply(summ.stats, function(x) x["MAE"])
)

## Performance calculated on adults only
adult_combined_data_list <- lapply(combined_data_list, function(x) x[x$chronological_age > 2.9,])

# Apply to the list of prediction to calculate summary statistics     
summ.stats.adult <- lapply(adult_combined_data_list,corr_and_mae)

# Convert to data frame
summ.stats.adult <- data.frame(
  tissue = names(summ.stats.adult),
  coef.cor.adult = sapply(summ.stats.adult, function(x) x["coef.cor"]),
  MAE.adult = sapply(summ.stats.adult, function(x) x["MAE"])
)

# Combine the model performance for all datapoints and >2.9 only
summ.stats = merge(summ.stats,summ.stats.adult)

# Save the output
saveRDS(list(combined_data,summ.stats), file = paste0("/path/to/tissue_predicted_agetransfo/AgePredictions.rds"))

combined_data <- combined_data %>% 
  mutate(tissue = recode(tissue, "omental_at" = "omental adipose")) %>%
  mutate(tissue = recode(tissue, "skeletal_muscle" = "skeletal muscle")) %>% 
  mutate(tissue = recode(tissue, "whole_blood" = "whole blood"))
summ.stats <- summ.stats %>% mutate(tissue = recode(tissue, "omental_at" = "omental adipose")) %>%
  mutate(tissue = recode(tissue, "skeletal_muscle" = "skeletal muscle")) %>% 
  mutate(tissue = recode(tissue, "whole_blood" = "whole blood"))

# Necessary for labeller below
summ.stats <- summ.stats[order(summ.stats$tissue),]

#-------------------------------------
##### PLOT CLOCK PREDICTIONS - Fig. 3A
#-------------------------------------

# Necessary for labeller below
summ.stats <- summ.stats[order(summ.stats$tissue),]

clocks_plot<-ggplot(data = combined_data, aes(x=chronological_age, y=predicted_age, col=tissue)) +
  geom_point(size=0.7) +
  scale_color_manual(values = extended_palette)+
  # adjust limits to the data
  xlim(min(cbind(min(combined_data$chronological_age),min(combined_data$predicted_age))),max(cbind(max(combined_data$chronological_age),max(combined_data$predicted_age))))+
  ylim(min(cbind(min(combined_data$chronological_age),min(combined_data$predicted_age))),max(cbind(max(combined_data$chronological_age),max(combined_data$predicted_age))))+
  labs(x="Chronological age", y="Predicted age")+
  theme_bw()+
  # facet per tissue
  facet_wrap(~tissue, ncol=7, labeller = as_labeller(function(label) {
    index <- which(summ.stats$tissue == label)
    paste(label)

  })) +
  theme(strip.text = element_text(size = 6, face = "bold"),
        strip.background = element_rect(fill="white"),
        axis.text = element_text(color="black"),
        axis.title = element_text(color="black"))+
  # add reference line chronological=predicted
  geom_abline()+
  
  # add linear regression line
  stat_smooth(method = "lm", span = 1.1, se = FALSE)+
  
  # aesthetics
  theme(plot.title = element_text(size=16, hjust=0.5,color = "black"),
        strip.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 16, color = "black"),
        axis.text = element_text(size = 16, color = "black"),legend.position = "none")
print(clocks_plot)
ggsave("/path/to/Figures/Panels_Figure3/MultitissueClocks.pdf",width = 12.5, height = 4.0)

#-------------------------------------
##### PLOT CLOCK PERFORMANCE - Fig. 3B
#-------------------------------------
                          
all_performance = summ.stats

# Order tissue levels by performance of Pearson's correlation
all_performance <- all_performance %>% 
  arrange(-desc(MAE)) %>%
  mutate(tissue=factor(tissue,levels=unique(tissue)))
  
# Long format for plot
all_performance_long <- all_performance %>% pivot_longer(cols = c(coef.cor,MAE), names_to = "metric", values_to = "performance")

# Define a named vector for mapping original labels to custom labels
label_mapping <- c("coef.cor" = "Pearson's correlation",
                   "MAE" = "MAE")

# Define custom labeling function using the named vector
custom_labels <- function(labels) {
  # Replace the original labels with custom labels using the named vector
  unname(label_mapping[as.character(labels)])
}

all_performance_long$metric <- factor(all_performance_long$metric, levels = c("MAE","coef.cor"))

ggplot(all_performance_long, aes(x = tissue, y = performance, fill = tissue)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = extended_palette) +
  labs(y = "performance", x = "") +
  facet_wrap(~ metric, ncol = 1, scales = "free_y",
             labeller = as_labeller(custom_labels)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, size = 16, colour = "black", vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.title = element_text(size = 20, colour = "black"),
        strip.text = element_text(size = 20),
        strip.background = element_rect(fill = "white"),
        legend.position = "none")
ggsave("/path/to/Figures/Panels_Figure3/Performance_clocks.pdf",width = 5,height = 5)
                          
#------------------------------------
##### AGE DEVIATIONS
#------------------------------------

# Function to extract lm coefficients for each facet
get_lm_coefficients <- function(df) {
  lm_fit <- lm(predicted_age ~ chronological_age, data = df)
  coefficients <- coef(summary(lm_fit))
  return(coefficients)
}

# Extract lm coefficients for each facet
facet_coefficients <- lapply(split(combined_data, combined_data$tissue), get_lm_coefficients)
print(facet_coefficients)

# Residual resulting from regressing DNAmAge on chronological age
lm_residuals <- list()
lm_mse <- list()
qq_plots <- list()

# load tissue types
tissue_list <- unique(combined_data$tissue)

# Loop through all tissues
for (tissues in tissue_list) {
  # Subset the data for the current tissue
  subset_data <- combined_data %>% filter(tissue == tissues)
  
  # Fit a linear model
  lm_model <- lm(predicted_age ~ chronological_age, data = subset_data)

  # Store the model in the list with species as the key
  lm_residuals[[tissues]] <- data.frame(lid_pid = subset_data$lid_pid, Residual = residuals(lm_model), tissue = tissues)
  
  # Create the ggplot QQ plot
  p <- ggqqplot(residuals(lm_model)) + ggtitle(tissues)
  
  qq_plots[[tissues]] <- p
  
}

# Plot in one panel
nrow <- ceiling(length(qq_plots) / 4)
gridExtra::grid.arrange(grobs = qq_plots, ncol = 4, nrow = nrow)

# unlist and combined data
lm_residuals <- imap_dfr(lm_residuals, ~mutate(.x, tissue = .y))
combined_data <- left_join(combined_data, lm_residuals)

### Same for adults only
lm_residuals.ad <- list()
lm_mse.ad <- list()
qq_plots.ad <- list()

# load tissue types
tissue_list.ad <- unique(combined_data_adult$tissue)

# Loop through all tissues
for (tissues in tissue_list.ad) {
  # Subset the data for the current tissue
  subset_data <- combined_data_adult %>% filter(tissue == tissues)
  
  # Fit a linear model for weight ~ size
  lm_model <- lm(predicted_age ~ chronological_age, data = subset_data)
  
  # Store the model in the list with species as the key
  lm_residuals.ad[[tissues]] <- data.frame(lid_pid = subset_data$lid_pid, Residual = residuals(lm_model), tissue = tissues)
  
  # Store the model precision (MSE)
  lm_mse.ad[[tissues]] <- mean(lm_model$residuals^2)
  
  # Create the ggplot QQ plot
  p <- ggqqplot(residuals(lm_model)) + ggtitle(tissues)
  
  qq_plots.ad[[tissues]] <- p
  
}

# Plot in one panel
nrow <- ceiling(length(qq_plots.ad) / 4)
gridExtra::grid.arrange(grobs = qq_plots.ad, ncol = 4, nrow = nrow)

# unlist and combined data
lm_residuals.ad <- imap_dfr(lm_residuals.ad, ~mutate(.x, tissue = .y))
lm_residuals.ad <- lm_residuals.ad %>% rename(Residual_adult = Residual)
combined_data <- merge(combined_data, lm_residuals.ad, all.x = TRUE)

# Export the dataset with age deviation
write.table(combined_data ,"/path/to/DNAm_deviation_data.txt", sep = "\t", row.names = F, quote = F)
                          
### Table S10
summary_residual <- summary(combined_data$Residual)
summary_residual.ad <- round(summary(combined_data$Residual_adult)[-7], digits = 4)
summary_corr <- summary(summ.stats$coef.cor)
summary_corr.ad <- summary(summ.stats$coef.cor.adult)
summary_mae <- summary(summ.stats$MAE)
summary_mae.ad <- summary(summ.stats$MAE.ad)
summary_info <- as.data.frame(rbind(summary_corr,rbind(summary_mae,rbind(summary_residual,rbind(summary_corr.ad,rbind(summary_mae.ad,summary_residual.ad))))))
summary_info$parameter <- rownames(summary_info)
summary_info <- mutate_all(summary_info,~if(is.numeric(.)) round(., digits = 3) else .)
summary_info <- summary_info %>%
  rename("Parameter" = parameter) %>%
  mutate(Parameter = recode(Parameter, "summary_corr" = "Pearson correlation",
                            "summary_mae" = "MAE", "summary_residual" = "Linear Regression Residual",
                            "summary_corr.ad" = "Pearson correlation (adult)",
                            "summary_mae.ad" = 'MAE (adult)', "summary_residual.ad" = "Linear Regression Residual (adult)")) %>%
  select(Parameter, everything())

write.csv(summary_info, "/path/to/ClocksPerformance.csv", row.names = F, quote = F)

#---------------------------------------------------------------
###### QC: TECHNICAL CORRELATES OF AGE DEVIATION
#---------------------------------------------------------------
                           
# Non-biological sources of variation in epigenetic age prediction
# Take inspiration from https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-10-r115#Fig3
# and from https://www.nature.com/articles/s43587-021-00134-3

metadata_lid = read.table("/path/to/metadata_final.txt", sep = "\t", header = TRUE) %>%  filter(lid_pid != "LID_109490_PID_10416")

qc_technical_df <- combined_data %>% left_join(., metadata_lid, by="lid_pid") %>% select(lid_pid, tissue, reads, unique, meth_cpg, unmeth_cpg, Residual, predicted_age, chronological_age)

## Fig. S10
p1<-ggplot(qc_technical_df, aes(x=reads, y=Residual))+geom_point()+labs(title="A")+theme_bw()+theme_plot
p2<-ggplot(qc_technical_df, aes(x=reads, y=predicted_age))+geom_point()+theme_bw()+theme_plot+labs(y="DNAm age")
p3<-ggplot(qc_technical_df, aes(x=unique, y=Residual))+geom_point()+labs(title="B")+theme_bw()+theme_plot+labs(x="uniquely mapped reads")
p4<-ggplot(qc_technical_df, aes(x=unique, y=predicted_age))+geom_point()+theme_bw()+theme_plot+labs(x="uniquely mapped reads",y="DNAm age")
p5<-ggplot(qc_technical_df %>% mutate(p_unique=unique/reads), aes(x=p_unique, y=Residual))+geom_point()+labs(title="C")+theme_bw()+theme_plot+labs(x="% uniquely mapped reads")
p6<-ggplot(qc_technical_df %>% mutate(p_unique=unique/reads), aes(x=p_unique, y=predicted_age))+geom_point()+theme_bw()+theme_plot+labs(x="% uniquely mapped reads", y="DNAm age")

Technical_QC_AgeDev1 <- ggpubr::ggarrange(p1,p3,p5,p2,p4,p6, nrow=2,ncol=3)
Technical_QC_AgeDev2 <- ggplot(qc_technical_df %>% mutate(p_unique=unique/reads), aes(x=p_unique, y=Residual))+geom_point()+labs(title="D")+facet_wrap(~tissue, nrow=2)+theme_bw()+theme_plot+labs(x="% uniquely mapped reads")+theme(strip.text = element_text(size=14))
ggpubr::ggarrange(Technical_QC_AgeDev1,Technical_QC_AgeDev2,nrow=2)
ggsave("/path/to/Figures/Technical_QC_AgeDev.pdf",width=12,height = 13)

### Technical correlates of mean squared error

# Loop through each subfolder to load the minimimal MSE and concatenate into a list of dataframes
result_mse <- data.frame(tissue=as.character(),minMSE = as.numeric(), stringsAsFactors = FALSE)
for(subfolder in subfolders){
  # Create the file path to MSE_result.txt in the current subfolder
  file_path <- file.path(paste0("/path/to/predicted_agetransfo/",subfolder), paste0("MSE_combined_result.txt"))
  
  # Read the data into a data frame
  data <- read.table(file_path, header = TRUE)
  
  # Find the min MSE
  min_mse <- min(data$MSE)
  
  # Store the result in the list (simply to print at the end)
  result_mse <- rbind(result_mse,data.frame(tissue = subfolder, minMSE = min_mse))
}

result_mse <- result_mse %>% mutate(tissue = forcats::fct_recode(tissue, "skeletal muscle" = "skeletal_muscle")) %>%
  mutate(tissue = forcats::fct_recode(tissue, "omental adipose" = "omental_at")) %>%
  mutate(tissue = forcats::fct_recode(tissue, "whole blood" = "whole_blood"))

# Rbind tissue specific summary statistics
result_mse=result_mse %>% left_join(.,combined_data %>% group_by(tissue) %>%
                                      summarize(sample_size=n(),age_range = (max(chronological_age)-min(chronological_age)), age_variation = sd(chronological_age)),
                                    join_by(tissue))

# MSE vs sample size
sample_size<-ggplot(result_mse,aes(x=sample_size,y=minMSE))+geom_point(size=2.5)+labs(
  x="sample size", y="minimal MSE",title="A"
  )+theme_bw()+theme_plot

# MSE vs age range
age_range <- ggplot(result_mse,aes(x=age_range,y=minMSE))+geom_point(size=2.5)+
  labs(
    x="age range (max-min) in years", y="minimal MSE",title="B"
  )+theme_bw()+theme_plot

# MSE vs age variation
age_var <- ggplot(result_mse,aes(x=age_variation,y=minMSE))+geom_point(size=2.5)+labs(
  x="chronological age SD", y="minimal MSE",title="C"
)+theme_bw()+theme_plot

ggpubr::ggarrange(sample_size,age_range, age_var,nrow = 1)
ggsave("/path/to/Figures/tissue_predicted_agetransfo/Technical_QC_mse.pdf",width = 9, height = 3)
