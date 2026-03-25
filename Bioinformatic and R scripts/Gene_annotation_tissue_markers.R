####### GENE & ENHANCER ANNOTATION OF TISSUE-SPECIFIC MARKERS
rm(list = ls())

library_list <- c("corrplot","svglite","tidyverse","RColorBrewer","patchwork","GenomicRanges","GenomicFeatures","rtracklayer","msigdbr","fgsea")
lapply(library_list, require, character.only=TRUE)

# === Paths ===
base_path <- "/path/to/project"
tissuemarker_path <- file.path(base_path, "tissue_comparisons", "tissuespecific_methylation.txt")
gtf_path <- file.path(base_path, "mmulatta_gtf")

# === Plot palette ===
tissue_plot <- sort(c("whole blood","spleen","omental adipose","heart","kidney","lung","adrenal","thymus","thyroid", "pituitary", "liver", "skeletal muscle"))
extended_palette <- setNames(colorRampPalette(brewer.pal(8, "Dark2"))(12), tissue_plot)
new_levels <- c("ovaries", "testis")
new_colors <- c("#008B8B", "#4682B4")
tissue_plot <- c(tissue_plot, new_levels)
extended_palette <- c(extended_palette, setNames(new_colors, new_levels))

# Optional recoding for tissues
recode_map <- c(
  "omental_at" = "omental adipose",
  "skeletal_muscle" = "skeletal muscle",
  "whole_blood" = "whole blood"
)

# === Load tissue marker table ===
tDMR <- readxl::read_excel(file.path(base_path, "SupplementaryTables.xlsx"), sheet = "TableS4")
tDMR <- tDMR %>% rename(site = region)
tDMR$start <- as.numeric(tDMR$start)
tDMR$end <- as.numeric(tDMR$end)

# Tissues of interest
tissue_oi <- c("whole_blood","spleen","omental_at","heart","testis","ovaries","kidney","lung","adrenal","thymus","thyroid","pituitary","liver","skeletal_muscle")

# =============================================================================
# GENE ANNOTATION
# =============================================================================

# Load GTF
gtf_file <- file.path(gtf_path, "Macaca_mulatta.Mmul_10.110.gtf")
gtf_macaque <- import(gtf_file)

# Extract gene info
gene_data <- data.frame(
  gene_id = elementMetadata(gtf_macaque)$gene_id,
  Gene    = elementMetadata(gtf_macaque)$gene_name,
  stringsAsFactors = FALSE
) %>% filter(!is.na(Gene)) %>% distinct(gene_id, Gene)

# TxDb object
txdb <- makeTxDbFromGFF(gtf_file)
genes_ <- genes(txdb)
promoters_ <- promoters(genes_, upstream = 2000, downstream = 200)
genomic_annotations <- list(genes = genes_, promoters = promoters_)

# GRanges for tissue markers
markers_gr <- with(tDMR, GRanges(seqnames = chr, IRanges(start, end), beta = mean_beta, tissue = sig_tissue))

# Overlaps
overlap_results_list <- list()
for (annot in names(genomic_annotations)) {
  genome_feat <- genomic_annotations[[annot]]
  hits <- findOverlaps(markers_gr, genome_feat)
  df <- cbind(tDMR[queryHits(hits), ],
              genome_feat[subjectHits(hits), ])
  names(df)[7:9] <- c("seqnames_annotation", "start_annotation", "end_annotation")
  df <- merge(df, gene_data, all.x = TRUE)
  overlap_results_list[[annot]] <- df
}

gene_results <- overlap_results_list$genes %>%
  select(chr, start, end, sig_tissue, mean_beta, gene_id, Gene,
         seqnames_annotation, start_annotation, end_annotation, width, strand) %>%
  rename(marker_tissue = sig_tissue) %>%
  mutate(region = paste(chr, start, end, sep = "_"))

promoter_results <- overlap_results_list$promoters %>%
  select(chr, start, end, sig_tissue, mean_beta, gene_id, Gene,
         seqnames_annotation, start_annotation, end_annotation, width, strand) %>%
  rename(marker_tissue = sig_tissue) %>%
  mutate(region = paste(chr, start, end, sep = "_"))

datasets <- list(gene_results = gene_results, promoter_results = promoter_results)

# =============================================================================
# FGSEA FUNCTION
# =============================================================================

compute_gene_stats <- function(df) {
  df %>%
    filter(!is.na(Gene)) %>%
    group_by(Gene, marker_tissue) %>%
    summarise(mean_beta = mean(mean_beta), .groups = "drop") %>%
    split(.$marker_tissue) %>%
    lapply(function(df_t) {
      stats <- df_t$mean_beta
      names(stats) <- df_t$Gene
      stats <- stats[!duplicated(names(stats))]
      sort(stats, decreasing = TRUE)
    })
}

# Load MSigDB once
gene_sets <- list(
  BP   = list(data = msigdbr(species = "Macaca mulatta", subcollection = "GO:BP"),
              label = "biological_processes", clean = "GOBP_"),
  CC   = list(data = msigdbr(species = "Macaca mulatta", subcollection = "GO:CC"),
              label = "cellular_component", clean = "GOCC_"),
  MF   = list(data = msigdbr(species = "Macaca mulatta", subcollection = "GO:MF"),
              label = "molecular_function", clean = "GOMF_"),
  KEGG = list(data = msigdbr(species = "Macaca mulatta", subcollection = "CP:KEGG_LEGACY"),
              label = "kegg", clean = "GOKEGG_"),
  HALL = list(data = msigdbr(species = "Macaca mulatta", collection = "H"),
              label = "hallmarks", clean = "HALLMARK_")
)

gene_sets <- lapply(gene_sets, function(x) { x$pathways <- split(x$data$gene_symbol, x$data$gs_name); x })

# FGSEA runner
run_fgsea_all <- function(datasets, gene_sets) {
  results <- list()
  for (ds_name in names(datasets)) {
    stats_list <- compute_gene_stats(datasets[[ds_name]])
    for (gs_name in names(gene_sets)) {
      gs <- gene_sets[[gs_name]]
      res <- lapply(names(stats_list), function(tissue) {
        fg <- fgseaMultilevel(pathways = gs$pathways,
                              stats    = stats_list[[tissue]],
                              minSize  = 15, maxSize = 500, eps = 0)
        fg$tissue <- tissue
        fg
      })
      res <- bind_rows(res) %>%
        filter(padj < 0.05) %>%
        mutate(annotation = ds_name,
               gene_set   = gs$label,
               pathway    = gsub(gs$clean, "", pathway),
               pathway    = gsub("_", " ", pathway),
               pathway    = tolower(pathway))
      results[[paste(ds_name, gs_name, sep = "_")]] <- res
    }
  }
  bind_rows(results)
}

final_gene_set_enrichment <- run_fgsea_all(datasets, gene_sets) %>%
  mutate(annotation = recode(annotation, gene_results = "gene_body", promoter_results = "promoter"))

# =============================================================================
# ENHANCER DMR FGSEA
# =============================================================================

# Load tissue-marker chromHMM annotation
tissue_marker_chrom <- readRDS(file.path(base_path, "revised_tissue_markers_chromHMM_annotated.rds"))
names_chrom <- names(tissue_marker_chrom)
for(i in seq_along(tissue_marker_chrom)) tissue_marker_chrom[[i]]$chromHMM <- names_chrom[i]

# Filter enhancer hypo
enhancer_markers <- lapply(tissue_marker_chrom, function(x) x %>% filter(enhancer == 1 & sig_hypo == 1))

map_chromHMM_to_tissuename_fast <- function(levels) {
  sapply(levels, function(lvl) {
    if (grepl("adipose", lvl)) return("omental_at")
    if (grepl("ovary", lvl)) return("ovaries")
    matches <- tissue_oi[str_detect(lvl, paste0(".*", tissue_oi, ".*"))]
    if(length(matches) > 0) matches[1] else NA
  })
}

enhancer_markers <- lapply(seq_along(enhancer_markers), function(i) {
  df <- enhancer_markers[[i]]
  df$chromHMM <- names(enhancer_markers)[i]
  df <- df %>% mutate(tissue = map_chromHMM_to_tissuename_fast(chromHMM))
})

names(enhancer_markers) <- names(tissue_marker_chrom)
rm(tissue_marker_chrom)

enhancer_markers_betas <- lapply(enhancer_markers, function(df) {
  merge(df %>% select(V1:V6, site, region, enhancer, sig_hypo, sig_hyper, chromHMM, tissue),
        tDMR %>% mutate(chr = paste0("chr", chr)),
        by.x = c("V1","V2","V3","tissue"), by.y = c("chr","start","end","sig_tissue"))
})

enhancer_markers_gr_list <- lapply(enhancer_markers_betas, function(df) {
  GRanges(seqnames = gsub("chr","",df$V4),
          ranges   = IRanges(start = df$V5, end = df$V6),
          strand   = if("strand" %in% colnames(df)) df$strand else "*",
          df[, intersect(colnames(df), c("sig_hypo","sig_hyper","mean_beta","marker","tissue")), drop = FALSE])
})

# FGSEA for enhancer markers
go_categories <- c(BP = "biological_processes", CC = "cellular_component", MF = "molecular_function")
go_sets <- lapply(names(go_categories), function(cat) {
  gs <- msigdbr(species = "Macaca mulatta", subcollection = paste0("GO:", cat))
  split(gs$gene_symbol, gs$gs_name)
})
names(go_sets) <- names(go_categories)

all_enhancer_results <- list()

for(subset_name in names(enhancer_markers_gr_list)) {
  gr <- enhancer_markers_gr_list[[subset_name]]
  nearest_hits <- GenomicRanges::nearest(gr, genes_)
  gr$gene_id <- mcols(genes_)$gene_id[nearest_hits]
  gene_map <- setNames(gene_data$Gene, gene_data$gene_id)
  gr$Gene <- gene_map[gr$gene_id]
  d <- distanceToNearest(gr, genes_)
  gr$distance <- mcols(d)$distance
  gr$weight <- 1 / (gr$distance + 1)
  
  df <- as.data.frame(gr)
  gene_scores <- df %>% filter(!is.na(Gene)) %>%
    group_by(Gene) %>%
    summarise(weighted_mean_beta = sum(mean_beta * weight)/sum(weight), .groups="drop")
  
  stats <- sort(setNames(gene_scores$weighted_mean_beta, gene_scores$Gene), decreasing=TRUE)
  
  for(cat in names(go_categories)) {
    pathways <- go_sets[[cat]]
    set.seed(2340)
    fg <- fgseaMultilevel(pathways = pathways, stats = stats, minSize = 15, maxSize = 500, eps = 0)
    fg$tissue <- unique(gr$tissue)
    fg$gene_set <- go_categories[cat]
    fg <- fg %>% filter(pval < 0.05) %>%
      mutate(pathway = gsub(paste0("GO", cat), "", pathway),
             pathway = gsub("_", " ", pathway),
             pathway = tolower(pathway))
    all_enhancer_results[[paste(subset_name, cat, sep="_")]] <- fg
  }
}

final_enhancer_enrichment <- bind_rows(all_enhancer_results)

# =============================================================================
# SAVE RESULTS
# =============================================================================

dir.create(file.path(base_path, "Tissue_spe_methylation"), showWarnings = FALSE)

write.csv(final_gene_set_enrichment %>%
            mutate(leadingEdge = sapply(leadingEdge, function(x) paste(x, collapse=";"))),
          file = file.path(base_path, "Tissue_spe_methylation", "gene_set_enrichment_tissuemarkers.csv"),
          quote = FALSE, row.names = FALSE)

write.csv(final_enhancer_enrichment %>%
            mutate(leadingEdge = sapply(leadingEdge, function(x) paste(x, collapse=";"))),
          file = file.path(base_path, "Tissue_spe_methylation", "enhancer_set_enrichment_tissuemarkers.csv"),
          quote = FALSE, row.names = FALSE)
