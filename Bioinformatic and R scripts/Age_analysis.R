# -------------------------------------------------------------
# Enrichment analysis of age-associated differentially methylated sites 
# -------------------------------------------------------------

rm(list = ls())

library_list <- c("corrplot","svglite","tidyverse","RColorBrewer","ggpubr")
lapply(library_list, require, character.only=TRUE)

# === Paths ===
base_path <- "/path/to/project"  # <-- Define this path only once

bed_path <- file.path(base_path, "bedfiles")
figure_path <- file.path(base_path, "Figures")

# === Color theme ===
tissue_plot <- sort(c("whole blood","spleen","omental adipose","heart","kidney","lung","adrenal","thymus","thyroid", "pituitary", "liver", "skeletal muscle"))
extended_palette <- setNames(colorRampPalette(brewer.pal(8, "Dark2"))(12),tissue_plot)

# === Tissues of interest ===
tissue_oi <- c("whole_blood","spleen","omental_at","heart","kidney","lung","adrenal","thymus","thyroid","pituitary","liver","skeletal_muscle")

# -----------------------------------
# === ENRICHEMENT OF CGISLANDS, GENE BODIES, AND PROMOTERS ===
# -----------------------------------

# === Helper function to read region files, select and rename cols, add site ===

read_region_file <- function(filepath) {
  read.table(filepath, header = FALSE, sep = "\t") %>%
    dplyr::select(V1, V2, V3) %>%
    rename(chr = V1, start = V2, end = V3) %>%
    mutate(site = paste(chr, start, end, sep = "_"))
}

# === Load annotated region files excluding chromHMM ===

region_files <- list.files(path = bed_path, pattern = "all_sites_age") %>%
  grep(pattern = "chromHMM", value = TRUE, invert = TRUE)

regionFiles <- lapply(region_files, function(f) read_region_file(file.path(bed_path, f)))
names(regionFiles) <- gsub(".bed", "", region_files)

# === Load CpG data and prepare annotation columns ===

regionToCpG <- read.table(file.path(bed_path, "all_sites_age_CpG.bed"), header = FALSE) %>%
  mutate(site = paste(V1, V2, V3, sep = "_"))

# Add blank columns for each region type and then add 1's where the site is in that region (both DMsites and non-DMsites)
reg_regions <- c("gene_body", "promoters", "cgIslands", "unassigned")
regionToCpG <- regionToCpG %>% mutate(!!!setNames(as.list(rep(0, length(reg_regions))), reg_regions))

# Assign 1 if site falls within region
regionToCpG <- regionToCpG %>%
  mutate(
    cgIslands = if_else(site %in% regionFiles[[1]]$site, 1, cgIslands),
    gene_body = if_else(site %in% regionFiles[[2]]$site, 1, gene_body),
    promoters = if_else(site %in% regionFiles[[3]]$site, 1, promoters),
    unassigned = if_else(gene_body == 0 & promoters == 0 & cgIslands == 0, 1, 0)
  )

# === Load age-associated sites ===

sigsites <- read.table(file.path(bed_path, "age_sharing.txt"), header = TRUE) %>%
  separate(site, into = c("type", "chr", "start", "stop")) %>%
  mutate(chrnew = paste0("chr", chr),
         site = paste(chrnew, start, stop, sep = "_"))

# Split sites by methylation direction of change with age
sigsites_pos <- sigsites %>% filter(focal_pm >0)
sigsites_neg <- sigsites %>% filter(focal_pm <0)

# Format                      
regionToCpG <- regionToCpG %>%
  mutate(region = paste(V4,V5,V6,sep="_")) %>%
  mutate(sig_pos = ifelse(region %in% sigsites_pos$site,1,0),
         sig_neg = ifelse(region %in% sigsites_neg$site,1,0))

# === Enrichment tests ===

# Function to perform enrichment test and tidy output
perform_enrichment <- function(region_col, sig_col) {
  res <- fisher.test(table(region_col, sig_col))
  data.frame(
    pvalue = res$p.value,
    conf.low = res$conf.int[1],
    conf.high = res$conf.int[2],
    OR = res$estimate,
    method = res$method
  )
}

# Run enrichment tests for each region and direction
enrichments <- function(sig_col) {
  map_dfr(reg_regions, ~ {
    perform_enrichment(regionToCpG[[.x]], regionToCpG[[sig_col]]) %>%
      mutate(type = .x)
  })
}

allenrich_pos <- enrichments("sig_pos") %>% mutate(direction_of_change = "hypermethylating")
allenrich_neg <- enrichments("sig_neg") %>% mutate(direction_of_change = "hypomethylating")

# === Plot Fig S4 ===

# Plot enrichment results helper
plot_enrichment <- function(df, title) {
  ggplot(df %>% mutate(type = recode(type, gene_body = "gene body")),
         aes(x = type, y = log2(OR), color = type)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = log2(conf.low), ymax = log2(conf.high))) +
    scale_color_brewer(palette = "Dark2") +
    ylim(-2.1, 2.8) +
    labs(y = ifelse(title == "A", "log2-transformed odds ratio", ""),
         x = "", title = title) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, color = "black", size = 12),
      axis.text.y = element_text(color = "black", size = 12),
      axis.title.y = element_text(color = "black", size = 12),
      plot.title = element_text(color = "black", size = 16),
      legend.position = "none"
    ) +
    geom_hline(yintercept = 0, linetype = "dashed")
}

enrich_hyper <- plot_enrichment(allenrich_pos, "A")
enrich_hypo <- plot_enrichment(allenrich_neg, "B")

ggpubr::ggarrange(enrich_hyper, enrich_hypo, align = "v")
ggsave(file.path(figure_path,"FigS4.pdf"), enrich_combined_plot, width=7.5, height=7.5)

# Export combined results Table S7
all_res <- bind_rows(allenrich_neg, allenrich_pos) %>%
  select(direction_of_change, type, OR, conf.low, conf.high, pvalue)

# -----------------------------------
# === ENRICHEMENT OF CHROMMHMM STATES ===
# -----------------------------------

rm(list = setdiff(ls(), c("base_path", "bed_path", "figure_path", "reg_regions", "perform_enrichment", "enrichments", "plot_enrichment","sigsites",
                          "tissue_plot","extended_palette", "tissue_oi")))

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

# === Load CpG data and prepare annotation columns ===

all_site_bed <- read.table(file.path(bed_path, "all_sites_age_CpG.bed"), header = FALSE) %>%
  mutate(site = paste(V1, V2, V3, sep = "_"),
         region = paste(V4, V5, V6, sep = "_"))

# Add blank columns for each region type and then add 1's where the site is in that region (both DMsites and non-DMsites)
chrom_anno <- c("1_TssA", "2_TssAFlnk","3_TxFlnk","4_Tx","5_TxWk","6_EnhG","7_Enh","8_ZNF/Rpts","9_Het",
                 "10_TssBiv","11_BivFlnk","12_EnhBiv","13_ReprPC","14_ReprPCWk","15_Quies","unassigned")

all_site_bed <- all_site_bed %>% mutate(!!!setNames(as.list(rep(0, length(chrom_anno))), chrom_anno))

# Initialize list to hold data frames for each tissue/chromHMM level
all_sites_alltissues <- list()

# Assign annotation 1 if site is in chromFile corresponding to the annotation
for(level in names(chromFiles)){
  all_sites_alltissues[[level]] <- all_site_bed
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

# Intersect the annotation datasets with significance for age-association
for(level in names(all_sites_alltissues)){
    if(grepl("adipose",level)){
    tissue_level <- "omental_at"
  }else{
    tissue_level <- tissue_oi[which(stringr::str_detect(level, fixed(tissue_oi)))]
  }
  
  sig_in_tissue <- sigsites[sigsites[[tissue_level]]<0.05,"site"]
  cpg_sig_in_tissue <- all_site_bed %>% filter(region %in% sig_in_tissue)
  
  all_sites_alltissues[[level]]$sig <- ifelse(all_sites_alltissues[[level]]$site %in% cpg_sig_in_tissue$site, 1, 0)
}

# === Enrichment tests ===

tissue_chromHmm_enrich <- list()

for(i in seq_along(all_sites_alltissues)) {
  x <- all_sites_alltissues[[i]]
  tissue <- names(all_sites_alltissues)[i]
  tissue_chromHmm_enrich[[tissue]] <- list()
  
  for(annotation in c(chrom_anno)){
    annotation_enrich <- fisher.test(table(x[[annotation]], x$sig))
    annotation_df <- as.data.frame(t(unlist(annotation_enrich)), stringsAsFactors = FALSE)
    
    colnames(annotation_df)<-c("pvalue", "conf.low", "conf.high", "OR", "null", "sided", "method", "table")
    annotation_df$annotation <- annotation
    annotation_df$chromHMM <- tissue
    
    tissue_chromHmm_enrich[[names(all_sites_alltissues)[i]]][[annotation]] <- annotation_df
  }
}

chrommHmm_fischer_results <- do.call(rbind, lapply(tissue_chromHmm_enrich,
                                                   function(tissue_results){
                                                     do.call(rbind,tissue_results)
                                                   }))

# Convert to numeric
chrommHmm_fischer_results <- chrommHmm_fischer_results %>%
  mutate(
    OR = as.numeric(OR),
    pvalue = as.numeric(pvalue),
    conf.low = as.numeric(conf.low),
    conf.high = as.numeric(conf.high),
    sig05 = ifelse(pvalue < 0.05, "Y", "N"),
    FDR = p.adjust(pvalue, method = "BH"),
    FDR05 = ifelse(FDR < 0.05, "Y", "N")
  )

chrommHmm_fischer_results$annotation <- factor(chrommHmm_fischer_results$annotation, levels=chrom_anno)

# Map original names to shortened tissue names
map_chromHMM_to_tissuename <- function(level) {
  if (grepl("adipose", level)) {
    return("omental_at")
  } else {
    # Find matching tissue_oi
    tissue_level <- tissue_oi[sapply(tissue_oi, function(tissue) grepl(paste0(".*", tissue, ".*"), level))]
    if (length(tissue_level) > 0) {
      return(tissue_level[1])  # Return first match if there are any
    } else {
      return(NA)  # Return NA if no match is found
    }
  }
}

chrommHmm_fischer_results <- chrommHmm_fischer_results %>%
  mutate(tissue = sapply(chromHMM, map_chromHMM_to_tissuename)) %>%
  mutate(tissue = recode(tissue,
                         "skeletal_muscle"="skeletal muscle",
                         "whole_blood"="whole blood",
                         "omental_at"="omental adipose"),
         annotation = recode(annotation,
                             "1_TssA"="Active TSS",
                             "2_TssAFlnk"="Flanking Active TSS",
                             "3_TxFlnk"="Transcr. at gene 5' and 3'",
                             "4_Tx"="Strong transcription",
                             "5_TxWk"="Weak transcription",
                             "6_EnhG"="Genic enhancers",
                             "7_Enh"="Enhancers",
                             "8_ZNF/Rpts"="ZNF genes & repeats",
                             "9_Het"="Heterochromatin",
                             "10_TssBiv"="Bivalent/Poised TSS",
                             "11_BivFlnk"="Flanking Bivalent TSS/Enh",
                             "12_EnhBiv"="Bivalent Enhancer",
                             "13_ReprPC"="Repressed PolyComb",
                             "14_ReprPCWk"="Weak Repressed PolyComb",
                             "15_Quies"="Quiescent",))

# === Table S6 ===
                                     
chrommHmm_fischer_results %>% select(tissue, annotation, OR, conf.low, conf.high, pvalue, sided, method, chromHMM)

# === Plot Fig 2B ===

ggplot(data=chrommHmm_fischer_results, aes(x=annotation, y=log2(OR), group=tissue, color=tissue, fill=tissue)) +
  geom_line(linewidth=0.8) +
  geom_point(aes(alpha=FDR05),size=3, pch=21, color="black") +
  scale_alpha_manual(values = c("Y" = 1, "N" = 0.2)) +
  scale_fill_manual(values= extended_palette, guide="none") +
  scale_color_manual(values= extended_palette, guide="none") +
  theme_bw() +
  geom_hline(yintercept=0, linetype="dashed") +
  labs(x="",y="log2-transformed odds ratio",alpha="FDR < 0.05") +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust = 1, color="black",size=14),
        axis.text.y=element_text(color="black",size=20),
        axis.title = element_text(color="black", size=20),
        legend.position = 'top',
        legend.key = element_rect(size=20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.box.margin = margin(-10, 0, 0, 0),
        legend.margin = margin(t=5, 0, 0, 0))
ggsave(file.path(figure_path,"Fig2B.pdf"), width=7.5,height=7.5)

# -----------------------------------
# === INTERSECTION OF TISSUE AGE ASSOCIATED SITES AND TISSUE-SPECIFIC MARKERS ===
# -----------------------------------

# === Load tissue markers and age DMRs ===
                                     
# Load tissue markers                            
tDMR <- read.table(file.path(base_path,"tissue_markers","tissuespecific_methylation.txt"),header=TRUE)

# Load the list of age-associated sites for all tissues (i.e., any sites that passed the lfsr in that tissue)                             
coeff_intercept_list <- readRDS(file.path(bed_path,"tissue_age_associated_sites_list.rds"))

# Load all sites tested for age                             
allsites_age <- read.table(file.path(bed_path,"all_sites_age.txt"), header=TRUE)

# Load all sites tested for tissue markers                               
dfs <- lapply(tissue_oi, function(tissue) {
  file_path <-  gsub("XXX", i, file.path(base_path, "tissues_meth", "XXX_meth", "Regions_pmeth_full_XXX_1000_14T.rds"))
  readRDS(file_path)$coverage
})

names(dfs) <- tissue_oi

# Find shared row names across all dataframes                              
allsites_markers <- Reduce(intersect, lapply(dfs, rownames)) #179,969 sites measured across all 12 tissues
markers_tested_for_age <- intersect(allsites_markers, allsites_age$regions)
tDMR_filtered <- tDMR %>% filter(sites %in% markers_tested_for_age)
                                     
# === Build data with merged tissue markers and aDMRs ===

# Combine all aDMRs into one dataframe
aDMRs <- do.call(rbind, coeff_intercept_list)

# Merge tissue markers with age markers
tDMR_aDMR <- merge(tDMR_filtered, aDMRs, by.x=c("sites", "sig_tissue"),
                   by.y=c("region", "tissue"), all.x=TRUE)

# Annotate aDMR
tDMR_aDMR <- tDMR_aDMR %>%
  mutate(aDMR = ifelse(is.na(beta), 0, 1),
         beta_age = beta,
         beta_marker = mean_beta,
         tissue_marker = sig_tissue,
         sign_beta_age = sign_beta) %>%
  select(-beta, -sign_beta, -mean_beta, -sig_tissue, -young_intercept, -pop_intercept, -intercept_high_low_50, -nonvariable)

# === Fisher test to determine enrichment ===

res_list <- lapply(coeff_intercept_list, function(x){
  result <- allsites_age %>% left_join(x %>%
                                  dplyr::select(region, beta, sign_beta),
                                join_by(regions==region)) %>% 
  dplyr::rename(beta_age=beta,sign_beta_age=sign_beta) %>%
  mutate(tissue=unique(x$tissue)) %>% 
  merge(.,tDMR_filtered, by.x=c("regions", "tissue"), by.y=c("sites","sig_tissue"),
        all.x=TRUE) %>% 
  mutate(aDMR = ifelse(is.na(beta_age),0,1),
         tDMR = ifelse(is.na(mean_beta),0,1))
  return(result)
})

fisher_list <- lapply(res_list, function(x){
  cbinding <- table(x$aDMR, x$tDMR)  
  spe_enrich <- fisher.test(cbinding)  
  spe_enrich <- as.data.frame(do.call(rbind, list(unlist(spe_enrich))))  
  tissue_value <- unique(x$tissue)
  spe_enrich$tissue <- tissue_value
  return(spe_enrich)
})

tissues_spe_enrich<- as.data.frame(do.call(rbind, fisher_list))
colnames(tissues_spe_enrich)<-c("pvalue", "conf.low", "conf.high", "OR", "null", "sided", "method", "table", "tissue")

# === Table S8 ===
                                     
tissues_spe_enrich <- tissues_spe_enrich %>%
  mutate(across(c(OR, pvalue, conf.low, conf.high), as.numeric))

# === Recode helper function ===
recode_tissue_names <- function(df, column) {
  df %>% mutate({{ column }} := recode({{ column }},
    "skeletal_muscle" = "skeletal muscle",
    "whole_blood" = "whole blood",
    "omental_at" = "omental adipose"
  ))
}

# Apply recoding
tissues_spe_enrich <- recode_tissue_names(tissues_spe_enrich, tissue)
tDMR_aDMR <- recode_tissue_names(tDMR_aDMR, tissue_marker)
                                     
# === Plot Fig 2G ===
                                     
ggplot(data=tissues_spe_enrich, aes(x=tissue, y=log2(OR), color=tissue)) +
  geom_point() +
  geom_errorbar(aes(ymin=log2(conf.low), ymax=log2(conf.high)))+
  scale_color_manual(values = extended_palette)+
  geom_text(aes(label = ifelse(pvalue < 0.05, "*", ""), # Add asterisk conditionally
                y = ifelse(log2(conf.low) < 0, log2(conf.low) - 0.2, log2(conf.high) + 0.05), # Position above the bar
                group = tissue), # Ensure proper grouping
    color="black",
    size=4,
    vjust = 0) +
  geom_hline(yintercept=0, linetype="dashed")+
  labs(y="log2-transformed odds ratio")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1,size=16, color="black"),
        axis.text.y=element_text(size=16, color="black"),
        axis.title.y=element_text(size=16, color="black"),
        legend.position = "none",
        axis.title.x = element_blank())
ggsave(file.path(figure_path, "Fig2G.pdf"), width=6,height=5)

# === FigS7 ===
                                     
ggplot(tDMR_aDMR, aes(x = beta_marker, y = beta_age)) +
  geom_point(alpha = 0.2) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  labs(x = "effect sizes for tissue marker", y = "effect sizes for age") +
  facet_wrap(~tissue_marker, scales = "free") +
  theme_bw() +
  theme(axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12, face = "bold"))
ggsave(file.path(figure_path,"FigS7.png"), width = 10, height = 6)
                                     
# === Fig2H ===

tDMR_aDMR = tDMR_aDMR %>%
    mutate(direction_of_change = case_when(
      beta_age<0 & beta_marker>0 ~ "loss_specificity",
      beta_age>0 & beta_marker<0 ~ "loss_specificity",
      beta_age<0 & beta_marker<0 ~ "gain_specificity",
      beta_age>0 & beta_marker>0 ~ "gain_specificity"
    ),
    direction_of_change_bino = ifelse(direction_of_change=="loss_specificity",1,0))

# Tissue summary
tissue_summ_direction_of_change <- tDMR_aDMR %>% 
  filter(aDMR==1) %>% 
  group_by(tissue_marker) %>% 
  summarize(prop_loss_specificity = sum(direction_of_change=="loss_specificity")/n(),
            n_sites = n()) %>% 
  arrange(-prop_loss_specificity) %>% 
  mutate(tissue_marker = factor(tissue_marker, levels=tissue_marker, ordered = TRUE))

tissue_summ_direction_of_change <- recode_tissue_names(tissue_summ_direction_of_change, tissue_marker)

ggplot(tissue_summ_direction_of_change, aes(x=tissue_marker, y=prop_loss_specificity, col=tissue_marker,size=n_sites))+
  geom_count()+
  scale_size_area(max_size = 20)+
  geom_hline(yintercept = 0.5, linetype="dashed")+
  scale_color_manual(values=extended_palette)+
  labs(x="", y="Prop. of marker sites losing specificity", size="Sample size")+
  theme_bw()+
  guides(size = guide_legend(override.aes = list(shape = c(21))))+
  theme(axis.text.x=element_text(size=22, angle=45, vjust=1, hjust = 1, color="black"),
        axis.text.y=element_text(size=22, color="black"),
        axis.title=element_text(size=22, color="black",hjust=1),
        legend.text = element_text(size=18),
        legend.title = element_text(size=18),
        # legend.position = "bottom",
        # legend.justification = c("right"),
        legend.box.margin = margin(-10, 0, 5, 5),
        legend.margin = margin(-10, 0, 0, 0))+
  guides(color="none")
ggsave(file.path(figure_path,"Fig2H.pdf"), width=8.6,height=6.5,dpi=300)
