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

combined_data = read.table("/path/to/DNAm_deviation_data.txt", sep="\t", header=TRUE)

# Load metadata
metadata_lid = read.table("/path/to/metadata.txt", sep = "\t", header = TRUE)
metadata_lid = metadata_lid %>% filter(lid_pid != "LID_109490_PID_10416")
metadata_lid$percent_unique<-metadata_lid$unique/metadata_lid$reads

# Add metadata to the age deviations datasets
combined_data = merge(combined_data, metadata_lid[,c("lid_pid", "monkey_id", "individual_sex")])

## Add z-transformed tissue-specific age deviation to make the residuals comparable across tissues
combined_data_ad <- combined_data %>% 
  filter(!is.na(Residual_adult)) %>%
  group_by(tissue) %>% 
  mutate(Zresidual_ad = as.vector(scale(Residual_adult))) %>% 
  ungroup()

# Calculate observed mean withinID variance
observed_mean_ad <- combined_data_ad %>%
  group_by(monkey_id) %>%
  summarise(withinID_var = var(Zresidual_ad, na.rm=TRUE)) %>%
  ungroup() %>%
  summarise(mean_Obs_withinID_var = mean(withinID_var, na.rm=TRUE)) %>%
  pull()

shuffle_residuals_ad <- function(df) {
  df %>%
    mutate(Zresidual_ad = sample(Zresidual_ad))
}

# Function to calculate variance within monkey_id
calculate_variance_ad <- function(df) {
  df %>%
    group_by(monkey_id) %>%
    summarise(variance = var(Zresidual_ad, na.rm=TRUE)) %>%
    ungroup() %>% 
    summarise(mean_variance = mean(variance, na.rm=TRUE)) %>%
    pull()
}

# Run permutations to extract randomized withinID variance
set.seed(123)  # for reproducibility

simulations <- 10000

MeanRandomized_withinID_var_ad <- replicate(simulations, {
  permuted_df <- combined_data_ad %>%
    group_by(tissue) %>%
    nest() %>%
    mutate(data = map(data, ~ shuffle_residuals_ad(.x))) %>%
    unnest(cols = data)
  
  variance_result <- calculate_variance_ad(permuted_df)
  variance_result
})

# Calculate p-value
num_extreme_values_ad <- sum(abs(MeanRandomized_withinID_var_ad - mean(MeanRandomized_withinID_var_ad)) >= abs(observed_mean_ad - mean(MeanRandomized_withinID_var_ad)))
p_value_ad <- num_extreme_values_ad / length(MeanRandomized_withinID_var_ad)
p_value_ad

## Plot observed vs randomized withinID variance

# Extract withinID variance on observed data
withinID_var_ad <- combined_data_ad %>%
  group_by(monkey_id) %>%
  summarise(withinID_var_ad = var(Zresidual_ad, na.rm=TRUE)) %>% pull(withinID_var_ad)

# Combine randomized and observed values
df_plot_withinID_var_ad = data.frame(source = c("observed_mean",rep("randomized", length(MeanRandomized_withinID_var_ad)),rep("observed",length(withinID_var_ad))),
                                  mean_var = c(observed_mean_ad,MeanRandomized_withinID_var_ad,withinID_var_ad))

# Fig. 3D #one outlier was kept in the analysis (which is considered conservative with respect to the null hypothesis) but is omitted from the plot.
withinID_consistency = ggplot(df_plot_withinID_var_ad %>% filter(source=="observed",mean_var<5),
                              aes(x=mean_var,fill=source))+
  geom_density(alpha=0.8,fill=NA)+
  geom_vline(data = subset(df_plot_withinID_var_ad, source == "observed_mean"),
             aes(xintercept = mean_var), color= "darkblue", linetype="dashed", linewidth=1)+
  geom_vline(data = subset(df_plot_withinID_var_ad, source == "randomized"),
             aes(xintercept = min(mean_var)), color= "black", linetype="dashed", linewidth=1)+
  labs(x="Within-individual variance\nin age deviation", y="Density")+
  theme_bw()+theme(axis.text = element_text(size=24, color="black"),
                   axis.title.y = element_text(size=24,color="black"),
                   axis.title.x = element_text(size=24, color="black",hjust=0.5), legend.position = "none")
ggsave("/path/to/Figures/3E.png")

######### Correlation of within-individual variance with age

# This time the outlier is removed to avoid influencing the relationship between age and variance (which is again conservative with respect to the null hypothesis)
combined_data <- combined_data %>% left_join(., metadata_lid[,c("lid_pid","group")])

df_heterogeneity_withinID_noOutlier_ad = combined_data %>%
  filter(monkey_id != "22H") %>%
  group_by(tissue) %>% 
  mutate(Zresidual_ad = as.vector(scale(Residual_adult))) %>% 
  ungroup() %>% 
  group_by(monkey_id, age, sex, group) %>% 
  summarize(withinID_var = var(Zresidual_ad),
            nb_tissue = n()) %>%
  filter(complete.cases(withinID_var))

## Expressing the degree of within individual repeatability with 
df_heterogeneity_withinID_noOutlier_ad$z.age <- as.vector(scale(df_heterogeneity_withinID_noOutlier_ad$chronological_age))
df_heterogeneity_withinID_noOutlier_ad$z.nb_tissue <- as.vector(scale(df_heterogeneity_withinID_noOutlier_ad$nb_tissue))

model_heterogeneity_aging_noOutlier_ad <- lm(withinID_var ~ z.age + individual_sex + group + z.nb_tissue,
                                          data = df_heterogeneity_withinID_noOutlier_ad)
summary(model_heterogeneity_aging_noOutlier_ad)

# Extract and tidy results "z.age" term only
forestplot_model_heterogeneity_aging_noOutlier_ad = broom::tidy(model_heterogeneity_aging_noOutlier_ad, conf.int = TRUE)
forestplot_model_heterogeneity_aging_noOutlier_ad = forestplot_model_heterogeneity_aging_noOutlier_ad %>% 
  mutate(term = recode(term, "(Intercept)"="intercept",
                       "individual_sexM"="sex(m)",
                       "z.age"="age",
                       "groupKK"="group(KK)",
                       "z.nb_tissue"="n° tissue")) %>% 
  mutate(term = factor(term, levels = rev(c("intercept","age","sex(m)","group(KK)","n° tissue"))))

# Fig. S12
ggplot(forestplot_model_heterogeneity_aging_noOutlier_ad, aes(x = estimate, y = term, xmin = conf.low, xmax = conf.high)) +
  geom_pointrange(size=1.1, color="lightblue", linewidth = 1.5) +  # Draw the points and whiskers
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # Zero line
  labs(x = "Effect Size (Estimate ± 95% CI)", y = NULL) +
  theme_bw() +
  theme(axis.text = element_text(size=14, colour = "black"),
        axis.title = element_text(size=14, colour = "black"))
ggsave("/path/to/Figures/FigS12.png")

####### Linear mixed models for the investigation of within-individual consistency

## Estimating the degree of within individual consistency with lmer - allowing to account for age, sex, group, and number of tissues
df_heterogeneity_withinID_noOutlier_ad_glm = combined_data %>%
  filter(monkey_id != "22H") %>% #note that this one is still an influential case among adult residuals
  group_by(tissue) %>% 
  mutate(Zresidual_ad = as.vector(scale(Residual_adult))) %>% 
  ungroup()

df_heterogeneity_withinID_noOutlier_ad_glm$z.age <- as.vector(scale(df_heterogeneity_withinID_noOutlier_ad_glm$chronological_age))

model_heterogeneity_aging_noOutlier_ad_lmer <- lme4::lmer(Zresidual_ad ~ z.age + individual_sex + group + tissue + (1|monkey_id),
                                                          data = df_heterogeneity_withinID_noOutlier_ad_glm, REML = F)
summary(model_heterogeneity_aging_noOutlier_ad_lmer)
qqnorm(residuals(model_heterogeneity_aging_noOutlier_ad_lmer))
qqline(residuals(model_heterogeneity_aging_noOutlier_ad_lmer))
hist(residuals(model_heterogeneity_aging_noOutlier_ad_lmer))

model_heterogeneity_aging_noOutlier_ad_lm <- lm(Zresidual_ad ~ z.age + individual_sex + group + tissue,
                                                data = df_heterogeneity_withinID_noOutlier_ad_glm)
anova(model_heterogeneity_aging_noOutlier_ad_lmer, model_heterogeneity_aging_noOutlier_ad_lm)

RLRsim::exactLRT(model_heterogeneity_aging_noOutlier_ad_lmer, model_heterogeneity_aging_noOutlier_ad_lm)

variance_components <- lme4::VarCorr(model_heterogeneity_aging_noOutlier_ad_lmer)

individual_variance <- as.numeric(variance_components$monkey_id[1])  # variance for monkey_id
residual_variance <- attr(variance_components, "sc")^2  # this looks weird but is the way to extract the residual variance

# Calculate repeatability
R <- individual_variance / (individual_variance + residual_variance)

# Print R
print(paste("Repeatability (R):", R))
