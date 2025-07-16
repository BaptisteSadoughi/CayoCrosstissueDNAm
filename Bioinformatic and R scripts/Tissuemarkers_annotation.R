# ==============================================================================
# Testing for enrichment for chromHMM annotations
# ==============================================================================

rm(list = ls())

library_list <- c("corrplot","svglite","tidyverse","RColorBrewer","patchwork")
lapply(library_list, require, character.only=TRUE)

# === Paths ===

base_path <- "/path/to/project"  # <-- Define this path only once
bed_path <- file.path(base_path, "bedfiles")
figure_path <- file.path(base_path, "Figures")
tissuemarker_path <- file.path(base_path, "tissue_comparisons", "tissuespecific_methylation.txt")

# === Tissues of interest ===

tissue_oi <- c("whole_blood","spleen","omental_at","heart","kidney","lung","adrenal","thymus","thyroid","pituitary","liver","skeletal_muscle")

# === Plot palette ===

tissue_plot <- sort(c("whole blood","spleen","omental adipose","heart","kidney","lung","adrenal","thymus","thyroid", "pituitary", "liver", "skeletal muscle"))
extended_palette <- setNames(colorRampPalette(brewer.pal(8, "Dark2"))(12),tissue_plot)

# -----------------------------------
# === LOAD AND FORMAT DATA ===
# -----------------------------------

# === Load chromHMM annotated bed files ===

temp <- list.files(path = bed_path, pattern = "all_sites_age_.*chromHMM\\.bed")
                      
chromFiles = lapply(temp, function(x){
  tmp <- read.delim(file = file.path(bed_path, x), sep = '\t', header = FALSE)
  colnames(tmp)[1:11] <- c("chr", "start", "end","chr_region","start_region","end_region","chr_state","start_state","end_state","chrom","bp_overlap")
  tmp <- tmp %>%
    mutate(site = paste(chr,start, end, sep="_"))
  return(tmp)
})
names(chromFiles) <- gsub(c(temp), pattern = ".bed", replacement = "")
rm(temp)

# === Define 7 states from chrommHMM ===

chromFiles <- lapply(chromFiles, function(x){
  x <- x %>% mutate(chrom = case_when(
    chrom %in% c("1_TssA", "2_TssAFlnk", "3_TxFlnk") ~ "promoter",
    chrom %in% c("4_Tx", "5_TxWk") ~ "active_tx",
    chrom %in% c("6_EnhG", "7_Enh") ~ "enhancer",
    chrom %in% c("8_ZNF/Rpts") ~ "znf",
    chrom %in% c("9_Het") ~ "heterochromatin",
    chrom %in% c("10_TssBiv","11_BivFlnk","12_EnhBiv") ~ "bivalent",
    chrom %in% c("13_ReprPC", "14_ReprPCWk") ~ "repressed_pc",
    chrom %in% c("15_Quies") ~ "quiescent",
    TRUE~chrom)
  )
}
)

chrom_anno <- c("promoter", "active_tx","enhancer","znf","heterochromatin","bivalent","repressed_pc","quiescent","unassigned")

# === Load tissue marker ===

tDMR = read.table(tissuemarker_path, sep="\t", header = TRUE)
tDMR <- tDMR %>% mutate(marker = ifelse(mean_beta>0, "hyper", "hypo"),
                        sites = gsub(pattern = "Region_","",sites)) %>% mutate(sites = paste0("chr",sites))

# Format site names to match pqlseq output
regions_to_cpg <- read.table(file.path(base_path,"metadata","regions_to_cpgs_mapping.bed"), header = FALSE)
regions_to_cpg = regions_to_cpg %>% mutate(region=paste(V1,V2,V3,sep = "_"))

regions_to_cpg_short <- regions_to_cpg[, c("V4", "V5", "V6")] %>% mutate(site = paste(V4,V5, V6, sep="_"))

# Initialize annotation columns
for (i in chrom_anno) { 
  regions_to_cpg_short[, paste0(i)] <- 0 
} 

# === Load regions tested as tissue marker ===

dfs <- lapply(tissue_oi, function(tissue) {
  file_path <-  gsub("XXX", tissue, file.path(base_path, "tissues_meth", "XXX_meth", "Regions_pmeth_full_XXX_1000_14T.rds"))
  readRDS(file_path)$coverage
})

names(dfs) <- tissue_oi
shared_rows <- Reduce(intersect, lapply(dfs, rownames))

# Standardize region names to match CpG format
shared_rows <- gsub("Region_", "chr", shared_rows)

# ==============================================================================
# MAP CHROMHMM STATES TO CpG-LEVEL SITES PER TISSUE
# ==============================================================================

all_sites_alltissues <- list()

for(level in names(chromFiles)){
   if(grepl("adipose",level)){
    tissue_level <- "omental_at"
  }else{
    tissue_level <- tissue_oi[which(stringr::str_detect(level, fixed(tissue_oi)))]
  }
  
  regions_to_cpg_level <- regions_to_cpg %>% filter(region %in% shared_rows)  
  regions_to_cpg_short_level <- regions_to_cpg_short %>% filter(site %in% with(regions_to_cpg_level,paste(V4,V5,V6,sep="_")))
  all_sites_alltissues[[level]] <- regions_to_cpg_short_level
  
  chromFile <- chromFiles[[level]]
  
  for (column_name in chrom_anno) { 
    all_sites_alltissues[[level]][[column_name]] <- ifelse(
      all_sites_alltissues[[level]]$site %in% chromFile$site[chromFile$chrom == column_name],
      1,
      all_sites_alltissues[[level]][[column_name]]
    ) 
  }
}

# Mark 'unassigned' sites where all other annotations are 0
all_sites_alltissues <- lapply(all_sites_alltissues, function(df) {
  df$unassigned <- ifelse(rowSums(df[, chrom_anno], na.rm = TRUE) == 0, 1, 0)
  df
})

# ==============================================================================
# ADD SIGNIFICANTLY TISSUE-SPECIFIC SITES TO THE DATASET
# ==============================================================================

#intersect the annotation datasets with significance for tissue-specificity
for(level in names(all_sites_alltissues)){
  
 if(grepl("adipose",level)){
    tissue_level <- "omental_at"
  }else{
    tissue_level <- tissue_oi[which(stringr::str_detect(level, fixed(tissue_oi)))]
  }

  hypo_sites <- tDMR %>% filter(sig_tissue == tissue_level, sign == "hypo") %>% pull(sites)
  hyper_sites <- tDMR %>% filter(sig_tissue == tissue_level, sign == "hyper") %>% pull(sites)

  cpg_hypo <- regions_to_cpg %>% filter(region %in% hypo_sites) %>% 
    mutate(cpg = paste(V4, V5, V6, sep = "_")) %>% pull(cpg)
  cpg_hyper <- regions_to_cpg %>% filter(region %in% hyper_sites) %>% 
    mutate(cpg = paste(V4, V5, V6, sep = "_")) %>% pull(cpg)

  all_sites_alltissues[[level]] <- all_sites_alltissues[[level]] %>%
    mutate(sig_hypo = ifelse(site %in% cpg_hypo, 1, 0),
           sig_hyper = ifelse(site %in% cpg_hyper, 1, 0))
}

# -----------------------------------
# === ENRICHMENT TESTS ===
# -----------------------------------

# Function to perform Fisher test for each annotation per tissue
run_fisher_tests <- function(sites_list, sig_column, chrom_anno){
  fisher_results <- list()
  
  for (i in seq_along(sites_list)) {
    df <- sites_list[[i]]
    tissue <- names(sites_list)[i]
    fisher_results[[tissue]] <- lapply(chrom_anno, function(annotation) {
      res <- fisher.test(table(df[[annotation]], df[[sig_column]]))
      df_res <- as.data.frame(t(unlist(res)), stringsAsFactors = FALSE)
      colnames(df_res) <- c("pvalue", "conf.low", "conf.high", "OR", "null", "sided", "method", "table")
      df_res$annotation <- annotation
      df_res$chromHMM <- tissue
      return(df_res)
    })
    names(fisher_results[[tissue]]) <- chrom_anno
  }

  # Combine results
  combined <- do.call(rbind, lapply(fisher_results, function(tissue) do.call(rbind, tissue)))
  combined <- combined %>%
    mutate(across(c(pvalue, OR, conf.low, conf.high), as.numeric),
           sig05 = ifelse(pvalue < 0.05, "sig", "nonsig"),
           FDR = p.adjust(pvalue, method = "BH"),
           annotation = factor(annotation, levels = chrom_anno))
  
  return(combined)
}

# Run for hypo and hyper
fisher_results_hypo <- run_fisher_tests(all_sites_alltissues, "sig_hypo", chrom_anno)
fisher_results_hyper <- run_fisher_tests(all_sites_alltissues, "sig_hyper", chrom_anno)                                    

# ==============================================================================
# ANNOTATION FORMATTING
# ==============================================================================

map_chromHMM_to_tissuename <- function(level) {
  if (grepl("adipose", level)) {
    return("omental_at")
  } else {
    tissue_level <- tissue_oi[sapply(tissue_oi, function(tissue) grepl(paste0(".*", tissue, ".*"), level))]
  if (length(tissue_level) > 0) return(tissue_level[1]) else return(NA)
  }
}

apply_common_formatting <- function(df) {
  df %>%
    mutate(tissue = sapply(chromHMM, map_chromHMM_to_tissuename)) %>%
    mutate(tissue = recode(tissue,
                           "skeletal_muscle" = "skeletal muscle",
                           "whole_blood" = "whole blood",
                           "omental_at" = "omental adipose")) %>%
    mutate(annotation = recode(annotation,
                               "active_tx" = "Active transcription",
                               "znf" = "ZNF genes & repeats",
                               "bivalent" = "Bivalent marks",
                               "repressed_pc" = "Repressed polycomb",
                               "enhancer" = "Enhancer",
                               "heterochromatin" = "Heterochromatin",
                               "quiescent" = "Quiescent",
                               "promoter" = "Promoters"))
}

fisher_results_hypo <- apply_common_formatting(fisher_results_hypo)
fisher_results_hyper <- apply_common_formatting(fisher_results_hyper)

# ==============================================================================
# HEATMAP PLOT FUNCTION
# ==============================================================================

plot_enrichment_heatmap <- function(df, filename, width = 5.5, height = 6.5, include_strip = TRUE) {
  
  df <- df %>%
    group_by(tissue) %>%
    mutate(scale_OR = scale(OR)) %>%
    ungroup()
  
  df$tissue <- factor(df$tissue, levels = rev(sort(unique(df$tissue))))
  df$annotation <- factor(df$annotation, levels = sort(unique(df$annotation)))

  present_tissues <- sort(unique(df$tissue))
  present_tissues <- intersect(tissue_plot, present_tissues)
  filtered_palette <- extended_palette[present_tissues]

  heatmap <- ggplot(df, aes(x = annotation, y = tissue, fill = scale_OR)) +
    geom_tile() +
    scale_fill_gradient(low = "red", high = "yellow") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, vjust=1, hjust = 1, size = 16, color = "black"),
      axis.text.y = element_blank(),
      axis.title = element_blank(),
      plot.title = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "top",
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16),
      plot.margin = margin(t = 10, r = 0, b = 10, l = 0)
    ) +
    labs(title = NULL, fill = "Standardized odds ratio")
  
  if (include_strip) {
    color_df <- data.frame(tissue = present_tissues)
    color_df$tissue <- factor(color_df$tissue, levels = levels(df$tissue))
    
    strip <- ggplot(color_df, aes(x = 1, y = tissue, fill = tissue)) +
      geom_tile() +
      scale_fill_manual(values = filtered_palette) +
      theme_void() +
      theme(
        legend.position = "none",
        axis.text = element_blank(),
        plot.margin = margin(t = 10, r = 0, b = 10, l = 0)
      )
    
    combined_plot <- strip + heatmap +
      patchwork::plot_layout(widths = c(0.3, 10)) &
      theme(plot.margin = margin(r = 0, l = 5))
    
    ggsave(file.path(figure_path, filename), combined_plot, width = width, height = height)
  } else {
    ggsave(file.path(figure_path, filename), heatmap, width = width, height = height)
  }
}

# ==============================================================================
# GENERATE AND SAVE FIGURES
# ==============================================================================

plot_enrichment_heatmap(fisher_results_hypo, "Fig1F.pdf", width = 5.5, height = 6.5, include_strip = TRUE)
plot_enrichment_heatmap(fisher_results_hyper, "FigS1.pdf", width = 5, height = 5, include_strip = FALSE)
