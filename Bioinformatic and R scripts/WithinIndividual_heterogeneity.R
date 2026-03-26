# ==============================================================================
# WITHIN INDIVIDUAL HETEROGENEITY IN PREDICTED AGES
# --------------------------------------------------------------------
# The aim is to investigate the level of within-individual consistency in DNAm
# age deviations across tissues.
# ==============================================================================

# === Clear workspace ===
rm(list = ls())

# === Load libraries ===
library_list <- c("corrplot","purrr","parallel","tidyverse", "DHARMa", "performance")
lapply(library_list, require, character.only=TRUE)

# === Paths ===

base_path <- "/path/to/project"  # <-- Define this path only once

dnam_devation_path <- file.path(base_path,"output","DNAm_deviation_data.txt")
metadata_path <- file.path(base_path,"metadata","multitissue_metadata.txt")
figure_path <- file.path(base_path, "Figures")

# === Color theme ===
palette <- c(RColorBrewer::brewer.pal(12, "Set3")[-2], RColorBrewer::brewer.pal(8, "Set2")[8], RColorBrewer::brewer.pal(12,"Paired")[1],RColorBrewer::brewer.pal(12, "Set3")[2])
my_theme <- theme_classic() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.background = element_blank())

# === Load data ===
combined_data = read.table(dnam_devation_path, sep="\t", header=TRUE)

# Load metadata
metadata <- read.table(metadata_path, sep="\t", header=TRUE) %>% 
  filter(lid_pid != "LID_109490_PID_10416") %>%
  mutate(percent_unique = unique / reads)

# Merge data
combined_data = merge(combined_data, metadata[,c("lid_pid", "monkey_id", "individual_sex", "group")])

# ----------------------------------------------------------------------------
# === ESTIMATING WITHIN-INDIVIDUAL CONSISTENCY VIA PERMUTATIONS ===
# ----------------------------------------------------------------------------

# Scaling tissue deviations
combined_data <- combined_data %>% 
  filter(!is.na(Residual_adult)) %>%
  group_by(tissue) %>% 
  mutate(Zresidual_ad = as.vector(scale(Residual_adult))) %>% 
  ungroup()

# Calculate observed mean withinID variance
observed_mean <- combined_data %>%
  group_by(monkey_id) %>%
  summarise(withinID_var = var(Zresidual_ad, na.rm=TRUE)) %>%
  ungroup() %>%
  summarise(mean_Obs_withinID_var = mean(withinID_var, na.rm=TRUE)) %>%
  pull()

# Permutation function
shuffle_residuals <- function(df) mutate(df, Zresidual_ad = sample(Zresidual_ad))

calculate_mean_variance <- function(df) {
  df %>%
    group_by(monkey_id) %>%
    summarise(var_id = var(Zresidual_ad, na.rm=TRUE)) %>%
    summarise(mean(var_id, na.rm = TRUE)) %>%
    pull()
}

# Run permutations to extract randomized withinID variance
set.seed(123)  # for reproducibility

simulations <- 10000

mean_randomized_var <- replicate(simulations, {
  permuted <- combined_data %>%
    group_by(tissue) %>% nest() %>%
    mutate(data = map(data, shuffle_residuals)) %>%
    unnest(cols = data)
  calculate_mean_variance(permuted)
})

# Calculate p-value
p_value_ad <- mean(abs(mean_randomized_var - mean(mean_randomized_var)) >= abs(observed_mean - mean(mean_randomized_var)))
p_value_ad

# === Plot observed vs randomized withinID variance ===

withinID_var <- combined_data %>%
  group_by(monkey_id) %>%
  summarise(withinID_var = var(Zresidual_ad, na.rm=TRUE)) %>% pull()

df_plot_withinID_var = data.frame(source = c("observed_mean",rep("randomized", length(mean_randomized_var)),rep("observed",length(withinID_var))),
                                  mean_var = c(observed_mean,mean_randomized_var,withinID_var))

# === Fig. 3E === #one outlier was kept in the analysis (which is considered conservative with respect to the null hypothesis) but is omitted from the plot.

withinID_consistency = ggplot(df_plot_withinID_var %>% filter(source=="observed",mean_var<5),
                              aes(x=mean_var,fill=source))+
  geom_density(alpha=0.8,fill=NA)+
  geom_vline(data = subset(df_plot_withinID_var, source == "observed_mean"),
             aes(xintercept = mean_var), color= "darkblue", linetype="dashed", linewidth=1)+
  geom_vline(data = subset(df_plot_withinID_var, source == "randomized"),
             aes(xintercept = min(mean_var)), color= "black", linetype="dashed", linewidth=1)+
  labs(x="Within-individual variance\nin age deviation", y="Density")+
  theme_bw()+theme(axis.text = element_text(size=24, color="black"),
                   axis.title.y = element_text(size=24,color="black"),
                   axis.title.x = element_text(size=24, color="black",hjust=0.5), legend.position = "none")

ggsave(file.path(figure_path,"3E.png"))

# -------------------------------------------------------
# === WITHIN-INDIVIDUAL CONSISTENCY VIA LINEAR MODELS ===
# -------------------------------------------------------

df_no_outlier <- combined_data %>%
  filter(!is.na(Residual_adult), monkey_id != "22H") %>%
  group_by(tissue) %>%
  mutate(Zresidual_ad = as.vector(scale(Residual_adult))) %>%
  ungroup() %>%
  mutate(z_age = as.vector(scale(chronological_age)))

# Mixed model
model_heterogeneity_aging_noOutlier_lmer <- lme4::lmer(Zresidual_ad ~ z_age + individual_sex + group + tissue + (1|monkey_id),
                                                          data = df_no_outlier, REML = F)
summary(model_heterogeneity_aging_noOutlier_lmer)

## Diagnostics
qqnorm(residuals(model_heterogeneity_aging_noOutlier_lmer))
qqline(residuals(model_heterogeneity_aging_noOutlier_lmer))
hist(residuals(model_heterogeneity_aging_noOutlier_lmer))

# Repeatability (R)
variance_components <- lme4::VarCorr(model_heterogeneity_aging_noOutlier_lmer)
individual_variance <- as.numeric(variance_components$monkey_id[1])  # variance for monkey_id
residual_variance <- attr(variance_components, "sc")^2  # this looks weird but is the way to extract the residual variance
R <- individual_variance / (individual_variance + residual_variance)
print(paste("Repeatability (R):", R))

performance::icc(
  model = model_heterogeneity_aging_noOutlier_lmer,
  by_group = FALSE, tolerance = 1e-05, ci = 0.95,
  iterations = 100, ci_method = NULL,  null_model = NULL,
  approximation = "lognormal", model_component = NULL, verbose = TRUE
)

# linear model
model_heterogeneity_aging_noOutlier_lm <- lm(Zresidual_ad ~ z_age + individual_sex + group + tissue,
                                                data = df_no_outlier)

# Testing the explanatory power of individual ID
anova(model_heterogeneity_aging_noOutlier_lmer, model_heterogeneity_aging_noOutlier_lm)

RLRsim::exactLRT(model_heterogeneity_aging_noOutlier_lmer, model_heterogeneity_aging_noOutlier_lm)

# ------------------------------------------------------
# === LINK BETWEEN ICC AND AGE SPREAD ===
# ------------------------------------------------------

# Test whether excluding certain ages impact ICC

col_df_icc <- c("age_thres","age_spread","age_means","sample_size","icc","ci_low","ci_high")
df_icc <- data.frame(matrix(nrow=0, ncol=length(col_df_icc)))
colnames(df_icc) <- col_df_icc

set.seed(589)
for(i in 3:16){
  df_subset <- combined_data %>%
    filter(monkey_id != "22H" & chronological_age > i) %>%
    group_by(tissue) %>% mutate(Zresidual_ad = as.vector(scale(Residual_adult))) %>% ungroup() %>%
    mutate(z_age = as.vector(scale(chronological_age)))
  
  model_i <- lme4::lmer(Zresidual_ad ~ z_age + individual_sex + group + tissue + (1|monkey_id), data = df_subset, REML = FALSE)
  icc_i <- performance::icc(model_i, by_group = FALSE, tolerance = 1e-05, ci = 0.95, iterations = 100, approximation = "lognormal", verbose = TRUE)
  
  df_icc <- rbind(df_icc, data.frame(
    age_thres = i,
    age_spread = sd(df_subset$chronological_age),
    age_means = mean(df_subset$chronological_age),
    sample_size = length(unique(df_subset$monkey_id)),
    icc = icc_i$ICC_adjusted[1],
    ci_low = icc_i$ICC_adjusted[2],
    ci_high = icc_i$ICC_adjusted[3]
  ))
}

# Plot ICC vs age threshold / sample size / mean age / age spread
plot_icc <- function(df, x, y, xlab){
  ggplot(df, aes_string(x = x, y = y)) +
    geom_point() +
    labs(x = xlab, y = "ICC") +
    theme_bw() +
    theme(axis.text = element_text(size=14, color="black"),
          axis.title = element_text(size=14, color="black"))
}

p1_obs <- plot_icc(df_icc, "age_thres", "icc", "Minimum age included")
p2_obs <- plot_icc(df_icc, "sample_size", "icc", "Number of individuals")
p3_obs <- plot_icc(df_icc, "age_means", "icc", "Average age")
p4_obs <- plot_icc(df_icc, "age_spread", "icc", "SD age")

# Is the drop in ICC due to low sample size or narrower age spread? 

col_df_icc_rand <- c("age_spread","age_means","sample_size","icc","ci_low","ci_high")
df_icc_rand <- data.frame(matrix(nrow=0, ncol=length(col_df_icc_rand)))
colnames(df_icc_rand) <- col_df_icc_rand

set.seed(789)
sample_sizes <- df_icc$sample_size[9:21]
for(i in rep(sample_sizes, times=5)){
  
  subject <- sample(unique(
    combined_data$monkey_id[combined_data$monkey_id != "22H" &
                              !is.na(Residual_adult)]),
    size = i)
  
  df_subset <- combined_data %>%
    filter(monkey_id %in% subject) %>%
    group_by(tissue) %>% mutate(Zresidual_ad = as.vector(scale(Residual_adult))) %>% ungroup() %>%
    mutate(z_age = as.vector(scale(chronological_age)))

  model_i <- lme4::lmer(Zresidual_ad ~ z_age + individual_sex + group + tissue + (1|monkey_id), data = df_subset, REML = FALSE)
  icc_i <- performance::icc(model_i, by_group = FALSE, tolerance = 1e-05, ci = 0.95, iterations = 100, approximation = "lognormal", verbose = TRUE)

    df_icc_rand <- rbind(df_icc_rand, data.frame(
    age_thres = i,
    age_spread = sd(df_subset$chronological_age),
    age_means = mean(df_subset$chronological_age),
    sample_size = length(unique(df_subset$monkey_id)),
    icc = icc_i$ICC_adjusted[1],
    ci_low = icc_i$ICC_adjusted[2],
    ci_high = icc_i$ICC_adjusted[3]
  ))
}

p1_subset <- plot_icc(df_icc_rand, "age_thres", "icc", "Minimum age included")
p2_subset <- plot_icc(df_icc_rand, "sample_size", "icc", "Number of individuals")
p3_subset <- plot_icc(df_icc_rand, "age_means", "icc", "Average age")
p4_subset <- plot_icc(df_icc_rand, "age_spread", "icc", "SD age")

ggpubr::ggarrange(p1_obs,p1_subset,p2_obs,p2_subset,p3_obs,p3_subset,p4_obs,p4_subset,ncol=2,nrow=4)
ggsave(file.path(figure_path,"FigS16.png"))

# --------------------------------------------------------------------
# === WITHIN-INDIVIDUAL VARIANCE AND CORRELATES ===
# --------------------------------------------------------------------

df_var <- df_no_outlier %>%
  group_by(monkey_id, chronological_age, individual_sex, group) %>%
  summarise(withinID_var = var(Zresidual_ad, na.rm = TRUE),
            nb_tissue = n()) %>%
  ungroup() %>%
  filter(complete.cases(withinID_var)) %>%
  mutate(z_age = scale(chronological_age),
         z_nb_tissue = scale(nb_tissue))

model_var <- lm(withinID_var ~ z_age + individual_sex + group + z_nb_tissue, data = df_var)
summary(model_var)

# Model diagnostic
{
  # Residuals
  hist(residuals(model_var), probability=T)
  qqnorm(residuals(model_var))
  qqline(residuals(model_var))
  plot(x=fitted(model_var), y=residuals(model_var),
       pch=19)
  abline(h=0, lty=3)
  # DFFITs
  max(abs(dffits(model_var))) # should be <2
  # DFBETAs
  xx=cbind(coef(model_var),
           coef(model_var)+t(apply(X=dfbeta(model_var), MARGIN=2, FUN=range)))
  colnames(xx)=c("orig", "min", "max")
  round(xx, 5)
  # Cook's distance
  max(cooks.distance(model_var)) #should be <1, >4/n
}

# Extract coefficients for reporting
coef_table <- coef(summary(model_var)) %>%
  as.data.frame() %>%
  rownames_to_column("Variable") %>%
  rename(Estimate = Estimate, `Std. Error` = `Std. Error`, `t-value` = `t value`, `p-value` = `Pr(>|t|)`) %>%
  mutate(across(c(Estimate, `Std. Error`, `t-value`), round, 3),
         `p-value` = ifelse(`p-value` < 0.001, "<0.001", round(`p-value`,3)))

# Recode variables for clarity
coef_table <- coef_table %>%
  mutate(Variable = recode(Variable,
                           "(Intercept)" = "intercept",
                           "individual_sexM" = "sex(m)",
                           "z_age" = "age",
                           "groupKK" = "group(KK)",
                           "z_nb_tissue" = "n° tissue"))

# Forest plot
forest_df <- broom::tidy(model_var, conf.int = TRUE) %>%
  mutate(term = recode(term,
                       "(Intercept)" = "intercept",
                       "individual_sexM" = "sex(m)",
                       "z_age" = "age",
                       "groupKK" = "group(KK)",
                       "z_nb_tissue" = "n° tissue"),
         term = factor(term, levels = rev(c("intercept","age","sex(m)","group(KK)","n° tissue"))))

ggplot(forest_df, aes(x=estimate, y=term, xmin=conf.low, xmax=conf.high)) +
  geom_pointrange(size=1.1, color="lightblue", linewidth=1.5) +
  geom_vline(xintercept = 0, linetype="dashed", color="black") +
  labs(x="Effect Size (Estimate ± 95% CI)", y=NULL) +
  theme_bw() +
  theme(axis.text = element_text(size=14, color="black"),
        axis.title = element_text(size=14, color="black"))

ggsave(file.path(figure_path,"FigS12.png"))
