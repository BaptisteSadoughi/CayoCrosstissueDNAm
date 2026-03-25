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

tissue_oi <- c("whole_blood","spleen","omental_at","heart","kidney","lung","adrenal","thymus","thyroid","pituitary","liver","skeletal_muscle","testis","ovaries")

# === Plot palette ===

tissue_plot <- sort(c("whole blood","spleen","omental adipose","heart","kidney","lung","adrenal","thymus","thyroid", "pituitary", "liver", "skeletal muscle"))
extended_palette <- setNames(colorRampPalette(brewer.pal(8, "Dark2"))(12),tissue_plot)

new_levels <- c("ovaries", "testis")
new_colors <- c("#008B8B", "#4682B4")

tissue_plot <- c(tissue_plot, new_levels)
extended_palette <- c(extended_palette, setNames(new_colors, new_levels))

# -----------------------------------
# === LOAD AND FORMAT DATA ===
# -----------------------------------

# === Load chromHMM annotated bed files ===

temp <- list.files(path = bed_path, pattern = "all_sites_age_.*chromHMM\\.bed")
                      
chromFiles = lapply(temp, function(x){
  tmp <- read.delim(file = file.path(bed_path, x), sep = '\t', header = FALSE)
 colnames(tmp)[1:8] <- c("chr", "start", "end","chr_state","start_state","end_state","chrom","bp_overlap")
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

# Load the concatenated list of tissue markers from Table S4
tDMR <- readxl::read_excel(paste0(base_path, "/SupplementaryTables.xlsx"), sheet = "TableS4")

tDMR <- tDMR %>% rename(sites = region)

# Format site names to match pqlseq output
regions_to_cpg <- read.table(file.path(base_path,"metadata","regions_to_cpgs_mapping.bed"), header = FALSE)
regions_to_cpg <- regions_to_cpg %>%
  mutate(site = paste(V4,V5, V6, sep="_"),
         region=paste(V1,V2,V3,sep = "_"))

# Initialize annotation columns
for (i in chrom_anno) { 
  regions_to_cpg[, paste0(i)] <- 0 
} 

# ==============================================================================
# MAP CHROMHMM STATES TO CpG-LEVEL SITES PER TISSUE
# ==============================================================================

## IMPORTANT NOT ALL REGIONS WERE TESTED
alltissues<-c("liver", "omental_at", "spleen", "kidney", "lung", "heart", "skeletal_muscle", "adrenal", "pituitary", "thymus", "thyroid", "whole_blood")

for(i in alltissues){
  filename <- gsub("XXX",i, file.path(base_path, "tissues_meth", "XXX_meth", "Regions_pmeth_full_XXX_1000_14T.rds"))
  r <- readRDS(filename)
  cov<-r$coverage
  assign(paste(i, "cov", sep="_"), cov)
}

dfs <- list(liver_cov, kidney_cov, lung_cov, heart_cov, omental_at_cov, spleen_cov, adrenal_cov, thymus_cov, thyroid_cov, pituitary_cov, whole_blood_cov, skeletal_muscle_cov)

# Find shared row names across all dataframes
shared_rows <- Reduce(intersect, lapply(dfs, rownames))

# format
shared_rows <-sapply(shared_rows, function(x) gsub(pattern = "Region_","",paste0("chr",x)),USE.NAMES = FALSE)

## FOR GONADS - NOT ALL REGIONS WERE TESTED
gonads <-c("ovaries", "testis")

for(i in gonads){
  filename <-  gsub("XXX", i, file.path(base_path, "tissues_meth", "XXX_meth", "Regions_pmeth_full_XXX_1000_14T.rds"))
  r <- readRDS(filename)
  cov<-r$coverage
  assign(paste(i, "cov", sep="_"), cov)
}

dfs_gonads <- list(liver_cov, kidney_cov, lung_cov, heart_cov, omental_at_cov, spleen_cov, adrenal_cov, thymus_cov, thyroid_cov, pituitary_cov, whole_blood_cov, skeletal_muscle_cov, ovaries_cov, testis_cov)
shared_rows_gonads <- Reduce(intersect, lapply(dfs_gonads, rownames))
                     
# format
shared_rows_gonads <-sapply(shared_rows_gonads, function(x) gsub(pattern = "Region_","",paste0("chr",x)),USE.NAMES = FALSE)

## Annotate and filter sites for each tissue
all_sites_alltissues <- list()

# Fill the levels of the list with dataframe for each tissue and intersecting the corresponding marks.
for(level in names(chromFiles)){
  
  #isolate tissue level for that 
  if(grepl("adipose",level)){
    tissue_level <- "omental_at"
  } else if (grepl("ovary",level)){
      tissue_level <- "ovaries"
    } else {
      tissue_level <- tissue_oi[which(sapply(tissue_oi,function(tissue)grepl(paste0(".*",tissue,".*"), level)))]
  }
  
  # restrict to tested tissue markers
  region_set <- if(tissue_level == "ovaries"){
    region_set <- shared_rows_gonads
    }else{
      region_set <- shared_rows
    }
  
  regions_to_cpg_level <- regions_to_cpg %>% filter(region %in% region_set)

  all_sites_alltissues[[level]] <- regions_to_cpg_level
  
  chromFile <- chromFiles[[level]]
  
  for (column_name in chrom_anno) { 
    all_sites_alltissues[[level]][[column_name]] <- ifelse(
      all_sites_alltissues[[level]]$site %in% chromFile$site[chromFile$chrom == column_name],
      1,
      all_sites_alltissues[[level]][[column_name]]
    ) 
  }
}

# Fill unassigned (check rows where all elements are 0 using rowSums)
all_sites_alltissues <- lapply(all_sites_alltissues, function(x) {
  df_subset <- x[,chrom_anno]
  x[rowSums(df_subset, na.rm=TRUE) == 0, "unassigned"] <- 1
  return(x)
})

#intersect the annotation datasets with significance for tissue-specificity
for(level in names(all_sites_alltissues)){
  
  #isolate tissue level for that 
  if(grepl("adipose",level)){
    tissue_level <- "omental_at"
  } else if (grepl("ovary",level)){
    tissue_level <- "ovaries"
  } else {
    tissue_level <- tissue_oi[which(sapply(tissue_oi,function(tissue)grepl(paste0(".*",tissue,".*"), level)))]
  }
  
  #extract significantly tissue-associated sites in the tissue_level
  marker_in_tissue_hypo <- tDMR[tDMR$sig_tissue == tissue_level & tDMR$marker=="hypo","sites"]
  marker_in_tissue_hyper <- tDMR[tDMR$sig_tissue == tissue_level & tDMR$marker=="hyper","sites"]
  
  # convert to the corresponding CpG level data
  cpg_sig_in_tissue_hypo <- regions_to_cpg %>% filter(region %in% marker_in_tissue_hypo)
  cpg_sig_in_tissue_hyper <- regions_to_cpg %>% filter(region %in% marker_in_tissue_hyper)
  
  all_sites_alltissues[[level]]$sig_hypo <- ifelse(all_sites_alltissues[[level]]$site %in%
                                                     cpg_sig_in_tissue_hypo$site, 1, 0)
  all_sites_alltissues[[level]]$sig_hyper <- ifelse(all_sites_alltissues[[level]]$site %in%
                                                      cpg_sig_in_tissue_hyper$site, 1, 0)
}

# -----------------------------------
# === ENRICHMENT TESTS ===
# -----------------------------------

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

  combined <- do.call(rbind, lapply(fisher_results, function(tissue) do.call(rbind, tissue)))
  
  combined <- combined %>%
    mutate(across(c(pvalue, OR, conf.low, conf.high), as.numeric),
           sig05 = ifelse(pvalue < 0.05, "sig", "nonsig"),
           FDR = p.adjust(pvalue, method = "BH"),
           annotation = factor(annotation, levels = chrom_anno))  # keep original order
  
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
  } else if (grepl("ovary", level)) {   # NEW
    return("ovaries")
  } else {
    tissue_level <- tissue_oi[sapply(tissue_oi, function(tissue) 
      grepl(paste0(".*", tissue, ".*"), level))]
    
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

fisher_results_hypo <- fisher_results_hypo %>%
  select(tissue, annotation, chromHMM, method, OR, conf.low, conf.high, pvalue, sig05, FDR)
rownames(fisher_results_hypo) <- NULL
fisher_results_hyper <- fisher_results_hyper %>%
  select(tissue, annotation, chromHMM, method, OR, conf.low, conf.high, pvalue, sig05, FDR)
rownames(fisher_results_hyper) <- NULL

# ==============================================================================
# HEATMAP PLOT FUNCTION
# ==============================================================================

plot_enrichment_heatmap <- function(df, filename, width = 5.5, height = 6.5, include_strip = TRUE) {
  
  cols <- colorRampPalette(c(
    "#3B0F70", "#8C2981", "#DE4968", "#FE9F6D", "#FEEB9C"
  ))(256)
  
  df <- df %>%
    group_by(tissue) %>%
    mutate(scale_OR = scale(OR)) %>%
    ungroup()
  
  df$tissue <- factor(df$tissue, levels = rev(sort(unique(df$tissue))))
  df$annotation <- factor(df$annotation, levels = sort(unique(df$annotation)))

  present_tissues <- sort(unique(df$tissue))
  present_tissues <- intersect(tissue_plot, present_tissues)
  filtered_palette <- extended_palette[present_tissues]

  missing_colors <- setdiff(present_tissues, names(filtered_palette))
  if(length(missing_colors) > 0){
    warning("Missing colors for: ", paste(missing_colors, collapse = ", "))
  }

  df <- df %>%
  mutate(
    annotation = stringr::str_to_sentence(annotation),
    annotation = recode(annotation,
                        "Znf genes & repeats" = "ZNF genes & repeats")
  )

  df$annotation <- factor(df$annotation, levels = unique(df$annotation))
  
  heatmap <- ggplot(df, aes(x = annotation, y = tissue, fill = scale_OR)) +
    geom_tile() +
    scale_fill_gradientn(
      colours = cols,
      values = scales::rescale(seq(-1.6, 2.6, length.out = 256))
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, vjust=1, hjust = 1, size = 16, color = "black"),
      axis.text.y = element_blank(),
      axis.title = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "top",
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16),
      plot.margin = margin(t = 10, r = 0, b = 10, l = 0)
    ) +
    labs(fill = "Standardized odds ratio")

  if (include_strip) {
    color_df <- data.frame(tissue = present_tissues)
    color_df$tissue <- factor(color_df$tissue, levels = levels(df$tissue))
    
    strip <- ggplot(color_df, aes(x = 1, y = tissue, fill = tissue)) +
      geom_tile() +
      scale_fill_manual(values = filtered_palette) +
      theme_void() +
      theme(legend.position = "none",
            plot.margin = margin(t = 10, r = 0, b = 10, l = 0))
    
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
