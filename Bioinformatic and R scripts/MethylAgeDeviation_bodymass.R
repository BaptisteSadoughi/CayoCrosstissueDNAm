# ==============================================================================
# Does weight adjusted for age predict DNAm clock deviations?
# ==============================================================================

# === Clear workspace ===
rm(list = ls())

# === Load libraries ===
library_list <- c("lme4","car","DHARMa","RColorBrewer","ggpubr","dplyr","splines","jtools","interactions")
lapply(library_list, require, character.only=TRUE)

# === Paths ===
base_path <- "/path/to/project"  # <-- Define this path only once

body_mass_data <- file.path(base_path,"metadata","MethylAgeDeviation_bodymass.txt")
output_path <- file.path(base_path, "output")

# === Color theme ===
tissue_plot <- sort(c("whole blood","spleen","omental adipose","heart","kidney","lung","adrenal","thymus","thyroid", "pituitary", "liver", "skeletal muscle"))
extended_palette <- setNames(colorRampPalette(brewer.pal(8, "Dark2"))(12),tissue_plot)
new_levels <- c("ovaries", "testis")
new_colors <- c("#008B8B", "#4682B4")
tissue_plot <- c(tissue_plot, new_levels)
extended_palette <- c(extended_palette, setNames(new_colors, new_levels))

# --------------------------------------------------------------------------
# === LINEAR REGRESSION OF BODY MASS ON AGE DEVIATIONS CORRECTED FOR AGE ===
# --------------------------------------------------------------------------

# === Load data ===
complete_data <- read.table(body_mass_data, sep = "txt", header = TRUE)

## Test association with age deviation (sex-specific mixed-effects models)

# null model
null <- lmer(Residual_adult ~ individual_sex + tissue + chronological_age + group + (1|monkey_id), data = complete_data, REML = F)

# full model
mod_sex_int <- lmer(Residual_adult ~ bw_adj * (individual_sex + tissue) +
                      chronological_age + group +
                      (1 | monkey_id),
                    data = complete_data,
                    REML = FALSE
                    )

anova(mod_sex_int, null, test = "Chisq") #sig full-null comparison

drop1(mod_sex_int) #ns interaction bw_adj * sex & sig bw_adj * tissue

# simplified model
mod_int_tissue <- lmer(Residual_adult ~ (bw_adj * tissue) + individual_sex + 
                         chronological_age + group +
                         (1 | monkey_id),
                       data = complete_data,
                       REML = FALSE
                       )

summary(mod_int_tissue)

# model without interaction terms
mod_bw <- lmer(Residual_adult ~ bw_adj + tissue + individual_sex +
                 chronological_age + group +
                 (1 | monkey_id),
               data = complete_data,
               REML = FALSE
               )

summary(mod_bw)

# === Model comparison ===
anova(mod_int_tissue, mod_bw, test = "Chisq") #sig interaction
anova(mod_bw, null, test = "Chisq") #sig for the main effect

# === Model diagnostics ===
vif(lm(Residual_adult ~ bw_adj + tissue + chronological_age + group + individual_sex, data = complete_data))

simulationoutput <- DHARMa::simulateResiduals(fittedModel = mod_int_tissue, plot = F, re.form = NULL)
residuals(simulationoutput, quantileFunction = qnorm, outlierValues = c(-7,7))
plot(simulationoutput)

# === Plotting predicted values ===
effect_plot(mod_bw, pred=bw_adj, interval=TRUE, plot.points=TRUE)
interact_plot(mod_int_tissue, pred=bw_adj, modx= tissue, interval=T, plot.points=F, colors = extended_palette)

# === Export model results (Tables S21) ===
library(flextable)
# some renaming needed for the function
boot.res.mod_int_tissue <- boot.res.int
boot.res.mod_bw <- boot.res.bw

format_model_result <- function(m){

  model_name <- deparse(substitute(m))
  conf_data <- get(paste0("boot.res.",model_name))

  table_model_int <- as.data.frame(coef(summary(m)))
  table_model_int <- cbind(table_model_int, conf_data$ci.estimates[,2:3])
  table_model_int <- cbind(Variable = rownames(table_model_int), table_model_int)
  rownames(table_model_int) <- NULL
  names(table_model_int) <- c("Variable", "Estimate", "Std. Error", "df", "t-value", "p-value", "lowCI","highCI")
  table_model_int = table_model_int %>%
    mutate(Variable= dplyr::recode(Variable, "(Intercept)"="intercept",
                                   "individual_sexM"="sex(m)",
                                   "groupKK"="group(KK)"),
           Variable = str_replace(Variable, "^bw_adj", "age-adjusted body mass")
           )

  # Round numeric columns for better formatting
  table_model_int <- table_model_int %>%
    mutate(across(c(Estimate, `Std. Error`, df, `t-value`, lowCI, highCI), ~ round(., 2)),
           `p-value` = ifelse(`p-value` < 0.001, "<0.001", round(`p-value`, 2)))

  # Create the table
  ft_int <- flextable(table_model_int) %>%
    theme_vanilla() %>%
    autofit() %>%
    bold(part = "header")

  flextable::save_as_image(ft_int, paste0(output_path,"/bodymass_agedev_",model_name,"_.png"), zoom = 2)
}

format_model_result(mod_int_tissue)
format_model_result(mod_bw)
