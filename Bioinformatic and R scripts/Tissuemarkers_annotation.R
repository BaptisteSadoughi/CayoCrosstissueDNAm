#### Testing for enrichment for chromHMM annotations

rm(list = ls())

library_list <- c("corrplot","svglite","tidyverse","RColorBrewer")
lapply(library_list, require, character.only=TRUE)

# Define ggplot theme upfront
tissue_plot <- sort(c("whole blood","spleen","omental adipose","heart","kidney","lung","adrenal","thymus","thyroid", "pituitary", "liver", "skeletal muscle"))
extended_palette <- setNames(colorRampPalette(brewer.pal(8, "Dark2"))(12),tissue_plot)

## Load output from bedtools intersect
temp <- list.files(path = "/path/to/MASH/",
                   pattern="all_sites_age_.*chromHMM\\.bed")

# read in files into a list
chromFiles = lapply(temp, function(x){
  tmp <- read.delim(file = paste0("/path/to/MASH/", x),
                    sep = '\t', header = F)
  colnames(tmp)[1:11] <- c("chr", "start", "end","chr_region","start_region","end_region","chr_state","start_state","end_state","chrom","bp_overlap")
  tmp <- tmp %>%
    mutate(site = paste(chr,start, end, sep="_"))
  return(tmp)
})
names(chromFiles) <- gsub(c(temp), pattern = ".bed", replacement = "")
rm(temp)

# 7 categories annotation of chromHMM
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

# Load the list of sites tested for tissue-specificity
tissue_oi <- c("whole_blood","spleen","omental_at","heart","kidney","lung","adrenal","thymus","thyroid","pituitary","liver","skeletal_muscle")

# Load differentially methylated regions in one tissue compared to all other tissues
tDMR <- read.table("/path/to/tissuespecific_methylation.txt",header=TRUE)
tDMR$sign <- ifelse(tDMR$mean_beta>0,"hyper","hypo")

# Format site names to match pqlseq output
regions_to_cpg<-read.table("/path/to/regions_to_cpgs_mapping.bed", header =F)
regions_to_cpg_short <- regions_to_cpg[, c("V4", "V5", "V6")] %>%
  mutate(site = paste(V4,V5, V6, sep="_"))
#format regions in regions_to_cpg
regions_to_cpg = regions_to_cpg %>% mutate(region=paste(V1,V2,V3,sep = "_"))

#7 levels
chrom_anno <- c("promoter", "active_tx","enhancer","znf","heterochromatin","bivalent","repressed_pc","quiescent","unassigned")

for (i in chrom_anno) { 
  regions_to_cpg_short[, paste0(i)] <- 0 
} 

## IMPORTANT NOT ALL REGIONS WERE TESTED
alltissues<-c("liver", "omental_at", "spleen", "kidney", "lung", "heart", "skeletal_muscle", "adrenal", "pituitary", "thymus", "thyroid", "whole_blood")

for(i in alltissues){
  filename <- gsub("XXX",i,"/path/to/tissues_meth/XXX_meth/Regions_pmeth_full_XXX_1000_14T.rds")
  # load data
  r <- readRDS(filename)
  cov<-r$coverage
  assign(paste(i, "cov", sep="_"), cov)
}

dfs <- list(liver_cov, kidney_cov, lung_cov, heart_cov, omental_at_cov, spleen_cov, adrenal_cov, thymus_cov, thyroid_cov, pituitary_cov, whole_blood_cov, skeletal_muscle_cov)  # Add all 12 dataframes to this list

# Find shared row names across all dataframes
shared_rows <- Reduce(intersect, lapply(dfs, rownames)) #179,969 sites measured across all 12 tissues

# format
shared_rows<-sapply(shared_rows, function(x) gsub(pattern = "Region_","",paste0("chr",x)),USE.NAMES = FALSE)

## Annotate and filter sites for each tissue
all_sites_alltissues <- list()

# Fill the levels of the list with dataframe for each tissue and intersecting the corresponding marks.
for(level in names(chromFiles)){
  
  #isolate tissue level for that 
  if(grepl("adipose",level)){
    tissue_level <- "omental_at"
  }else{
    tissue_level <- tissue_oi[which(sapply(tissue_oi,function(tissue)grepl(paste0(".*",tissue,"*"), level)))]
  }
  
  # regions_to_cpg_level <- regions_to_cpg %>% filter(site %in% all_sites_list[[tissue_level]]$sites)
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

# Fill unassigned (check rows where all elements are 0 using rowSums)
all_sites_alltissues <- lapply(all_sites_alltissues, function(x) {
  df_subset <- x[,chrom_anno]
  x[rowSums(df_subset, na.rm=TRUE) == 0, "unassigned"] <- 1
  return(x)
})

# Harmonize format in tDMRs
tDMR <- tDMR %>% mutate(sites = gsub(pattern = "Region_","",sites)) %>% mutate(sites = paste0("chr",sites))

#intersect the annotation datasets with significance for tissue-specificity
for(level in names(all_sites_alltissues)){
  
  #isolate tissue level for that 
  if(grepl("adipose",level)){
    tissue_level <- "omental_at"
  }else{
    tissue_level <- tissue_oi[which(sapply(tissue_oi,function(tissue)grepl(paste0(".*",tissue,"*"), level)))]
  }
  
  #extract significantly tissue-associated sites in the tissue_level
  marker_in_tissue_hypo <- tDMR[tDMR$sig_tissue == tissue_level & tDMR$sign=="hypo","sites"]
  marker_in_tissue_hyper <- tDMR[tDMR$sig_tissue == tissue_level & tDMR$sign=="hyper","sites"]
  
  # convert to the corresponding CpG level data
  cpg_sig_in_tissue_hypo <- regions_to_cpg %>% filter(region %in% marker_in_tissue_hypo)
  cpg_sig_in_tissue_hyper <- regions_to_cpg %>% filter(region %in% marker_in_tissue_hyper)
  
  cpg_sig_in_tissue_hypo_formatted <- cpg_sig_in_tissue_hypo %>% mutate(cpg=paste(V4,V5,V6,sep = "_")) %>% pull(cpg)
  cpg_sig_in_tissue_hyper_formatted <- cpg_sig_in_tissue_hyper %>% mutate(cpg=paste(V4,V5,V6,sep = "_")) %>% pull(cpg)
  
  all_sites_alltissues[[level]]$sig_hypo <- ifelse(all_sites_alltissues[[level]]$site %in%
                                                cpg_sig_in_tissue_hypo_formatted, 1, 0)
  all_sites_alltissues[[level]]$sig_hyper <- ifelse(all_sites_alltissues[[level]]$site %in%
                                                     cpg_sig_in_tissue_hyper_formatted, 1, 0)
}

# Export the annotation of tissue specific markers with chromHMM
saveRDS(all_sites_alltissues,"/path/to/tissue_markers_chromHMM_annotated.rds")

###### ENRICHMENT FOR CHROMHMM ANNOTATIONS AMONG TISSUE MARKERS

# Perform enrichment test with Fisher's test 

### Hypomarker enrichment
tissue_chromHmm_enrich_hypo <- list()

for(i in seq_along(all_sites_alltissues)) {
  
  x <- all_sites_alltissues[[i]]
  tissue <- names(all_sites_alltissues)[i]
  tissue_chromHmm_enrich_hypo[[tissue]] <- list()
  
  for(annotation in c(chrom_anno, "unassigned")){
    
    annotation_enrich <- fisher.test(table(x[[annotation]], x$sig_hypo))
    
    # Convert the result to an unlisted vector
    annotation_df <- as.data.frame(t(unlist(annotation_enrich)), stringsAsFactors = FALSE)
    
    colnames(annotation_df)<-c("pvalue", "conf.low", "conf.high", "OR", "null", "sided", "method", "table")
    
    annotation_df$annotation <- annotation
    
    annotation_df$chromHMM <- tissue
    
    # Add the result to the list
    tissue_chromHmm_enrich_hypo[[names(all_sites_alltissues)[i]]][[annotation]] <- annotation_df
  }
}

chrommHmm_fischer_results_hypo <- do.call(rbind, lapply(tissue_chromHmm_enrich_hypo,
                                                        function(tissue_results){
                                                          do.call(rbind,tissue_results)
                                                        }))
chrommHmm_fischer_results_hypo$OR<- as.numeric(chrommHmm_fischer_results_hypo$OR)
chrommHmm_fischer_results_hypo$pvalue<- as.numeric(chrommHmm_fischer_results_hypo$pvalue)
chrommHmm_fischer_results_hypo$conf.low<- as.numeric(chrommHmm_fischer_results_hypo$conf.low)
chrommHmm_fischer_results_hypo$conf.high<- as.numeric(chrommHmm_fischer_results_hypo$conf.high)
chrommHmm_fischer_results_hypo$sig05 <- ifelse(chrommHmm_fischer_results_hypo$pvalue<0.05,"sig","nonsig")
chrommHmm_fischer_results_hypo$annotation <- factor(chrommHmm_fischer_results_hypo$annotation, levels=chrom_anno)
# Add BH corrected p-values
chrommHmm_fischer_results_hypo$FDR <- p.adjust(chrommHmm_fischer_results_hypo$pvalue, method = "BH")

# Add back the tissue names for the plot
# Function to map original levels to short tissue names
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

# Add a new column 'tissue_short' to the data frame
chrommHmm_fischer_results_hypo <- chrommHmm_fischer_results_hypo %>%
  mutate(tissue = sapply(chromHMM, map_chromHMM_to_tissuename))

chrommHmm_fischer_results_hypo = chrommHmm_fischer_results_hypo %>% mutate(tissue = recode(tissue,
                                                                                 "skeletal_muscle"="skeletal muscle",
                                                                                 "whole_blood"="whole blood",
                                                                                 "omental_at"="omental adipose"))
## Export Table S2
write.csv(chrommHmm_fischer_results_hypo,"/path/to/TableS2.csv", row.names = F, quote = F)

## Fig. 1F
                                     
library(RColorBrewer)
library(patchwork)
                                     
data_heatmap=chrommHmm_fischer_results_hypo %>% 
  group_by(tissue) %>% 
  mutate(scale_OR = scale(OR)) %>% 
  ungroup()

data_heatmap$tissue <- factor(data_heatmap$tissue, levels = rev(sort(unique(data_heatmap$tissue))))
data_heatmap$annotation <- factor(data_heatmap$annotation, levels = sort(unique(data_heatmap$annotation)))

data_heatmap <- data_heatmap %>% mutate(annotation=recode(annotation, "active_tx"="Active transcription",
                                                          "znf" = "ZNF genes & repeats",
                                                          "bivalent" = "Bivalent marks",
                                                          "repressed_pc" = "Repressed polycomb",
                                                          "enhancer"="Enhancer",
                                                          "heterochromatin"="heterochromatin",
                                                          "quiescent"="Quiescent",
                                                          "promoter"="Promoters"
                                                          ))

data_heatmap$tissue <- factor(data_heatmap$tissue, levels = sort(unique(data_heatmap$tissue)))

present_tissues <- sort(unique(data_heatmap$tissue))
present_tissues <- intersect(tissue_plot, present_tissues)

filtered_palette <- extended_palette[present_tissues]

missing_colors <- setdiff(present_tissues, names(filtered_palette))
if(length(missing_colors) > 0){
  warning("The following tissues do not have assigned colors: ", paste(missing_colors, collapse = ", "))
  # Optionally assign default colors or handle accordingly
}

heatmap_plot <- ggplot(data_heatmap, aes(x = annotation, y = tissue, fill = scale_OR)) +
  geom_tile() +
  scale_fill_gradient(low = "red", high = "yellow") +  # Red to Yellow gradient
  labs(title = "Organ-specific markers", fill="Standardized odds ratio")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust = 1, color="black",size=16),  # Rotate x-axis labels
        # axis.text.y = element_text(color="black",size=14),  # Adjust y-axis text size
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        # plot.title = element_text(face = "bold",colour = "black",size=10),
        panel.grid.major.y = element_blank(),
        plot.title = element_blank(),
        legend.position = "top",
        legend.title = element_text(size=16),
        legend.text = element_text(size=16),
        plot.margin = margin(t = 10, r = 0, b = 10, l = 0)
  )

color_key_df <- data.frame(tissue = present_tissues)

color_key_df$tissue <- factor(color_key_df$tissue, levels = levels(data_heatmap$tissue))

color_key_plot <- ggplot(color_key_df, aes(x = 1, y = tissue, fill = tissue)) +
  geom_tile() +
  scale_fill_manual(values = filtered_palette) +
  theme_void() +  # Removes all background, gridlines, and axes
  theme(
    legend.position = "none",  # Hide the legend
    axis.text = element_blank(),  # Add tissue names next to tiles
    plot.margin = margin(t = 10, r = 0, b = 10, l = 0)  # Adjust margins to align with heatmap
  )

combined_plot <- color_key_plot + 
  heatmap_plot  + 
  patchwork::plot_layout(widths = c(0.3, 10)) & theme(plot.margin = margin(r=0,l=5)) # Adjust widths to allocate space appropriately

ggsave("/path/to/Figures/Fig1F.pdf", width=5.5, height = 6.5)

######### Hyper marker enrichment

tissue_chromHmm_enrich_hyper <- list()

for(i in seq_along(all_sites_alltissues)) {
  
  x <- all_sites_alltissues[[i]]
  tissue <- names(all_sites_alltissues)[i]
  tissue_chromHmm_enrich_hyper[[tissue]] <- list()
  
  for(annotation in c(chrom_anno, "unassigned")){
    
    annotation_enrich <- fisher.test(table(x[[annotation]], x$sig_hyper))
    
    # Convert the result to an unlisted vector
    annotation_df <- as.data.frame(t(unlist(annotation_enrich)), stringsAsFactors = FALSE)
    
    colnames(annotation_df)<-c("pvalue", "conf.low", "conf.high", "OR", "null", "sided", "method", "table")
    
    annotation_df$annotation <- annotation
    
    annotation_df$chromHMM <- tissue
    
    # Add the result to the list
    tissue_chromHmm_enrich_hyper[[names(all_sites_alltissues)[i]]][[annotation]] <- annotation_df
  }
}

chrommHmm_fischer_results_hyper <- do.call(rbind, lapply(tissue_chromHmm_enrich_hyper,
                                                        function(tissue_results){
                                                          do.call(rbind,tissue_results)
                                                        }))

chrommHmm_fischer_results_hyper$OR<- as.numeric(chrommHmm_fischer_results_hyper$OR)
chrommHmm_fischer_results_hyper$pvalue<- as.numeric(chrommHmm_fischer_results_hyper$pvalue)
chrommHmm_fischer_results_hyper$conf.low<- as.numeric(chrommHmm_fischer_results_hyper$conf.low)
chrommHmm_fischer_results_hyper$conf.high<- as.numeric(chrommHmm_fischer_results_hyper$conf.high)
chrommHmm_fischer_results_hyper$sig05 <- ifelse(chrommHmm_fischer_results_hyper$pvalue<0.05,"sig","nonsig")
chrommHmm_fischer_results_hyper$annotation <- factor(chrommHmm_fischer_results_hyper$annotation, levels=chrom_anno)
# Add BH corrected p-values
chrommHmm_fischer_results_hyper$FDR <- p.adjust(chrommHmm_fischer_results_hyper$pvalue, method = "BH")

# Add a new column 'tissue_short' to the data frame
chrommHmm_fischer_results_hyper <- chrommHmm_fischer_results_hyper %>%
  mutate(tissue = sapply(chromHMM, map_chromHMM_to_tissuename))

chrommHmm_fischer_results_hyper = chrommHmm_fischer_results_hyper %>% mutate(tissue = recode(tissue,
                                                                                 "skeletal_muscle"="skeletal muscle",
                                                                                 "whole_blood"="whole blood",
                                                                                 "omental_at"="omental adipose"))
## Export Table S3
write.csv(chrommHmm_fischer_results_hyper,"/path/to/TableS3.csv", row.names = F, quote = F)

## Heatmap version
data_heatmap_hyper=chrommHmm_fischer_results_hyper %>% 
  group_by(tissue) %>% 
  mutate(scale_OR = scale(OR)) %>% 
  ungroup()

data_heatmap_hyper$tissue <- factor(data_heatmap_hyper$tissue, levels = rev(sort(unique(data_heatmap_hyper$tissue))))
data_heatmap_hyper$annotation <- factor(data_heatmap_hyper$annotation, levels = sort(unique(data_heatmap_hyper$annotation)))

data_heatmap_hyper <- data_heatmap_hyper %>% mutate(annotation=recode(annotation, "active_tx"="Active transcription",
                                                          "znf" = "ZNF genes & repeats",
                                                          "bivalent" = "Bivalent marks",
                                                          "repressed_pc" = "Repressed polycomb",
                                                          "enhancer"="Enhancer",
                                                          "heterochromatin"="heterochromatin",
                                                          "quiescent"="Quiescent",
                                                          "promoter"="Promoters"
))

ggplot(data_heatmap_hyper, aes(x = annotation, y = tissue, fill = scale_OR)) +
  geom_tile() +
  scale_fill_gradient(low = "red", high = "yellow") +  # Red to Yellow gradient
  # scale_x_discrete(position="top")+
  labs(title = "Organ-specific markers", fill="Standardized odds ratio")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust = 1, color="black",size=12),  # Rotate x-axis labels
        axis.text.y = element_text(color="black",size=12),  # Adjust y-axis text size
        axis.title = element_blank(),
        text = element_text(size=12),
        # plot.title = element_text(face = "bold",colour = "black",size=10),
        plot.title = element_blank(),
        legend.position = "top",
      legend.justification = "right"
  )
ggsave("/path/to/Figures/FigS1.pdf", width=5, height = 5)
