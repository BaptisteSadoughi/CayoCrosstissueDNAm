#### Testing for enrichment for cgIslands, Gene bodies, and Promoters

rm(list = ls())
library_list <- c("corrplot","svglite","tidyverse")
lapply(library_list, require, character.only=TRUE)

### Assign regulatory regions -----------------
## Load output from bedtools intersect
temp <- list.files(path = "/path/to/MASH/bedfiles/",
                   pattern="all_sites_age.bed")
# drop chroHMM annotations and TEs
temp <- grep("chromHMM|repeatmasker", temp, value = TRUE, invert = TRUE)

# read in files into a list
regionFiles = lapply(temp, function(x){
  tmp <- read.delim(file = paste0("/path/to/MASH/bedfiles/", x),
                    sep = '\t', header = F) %>%
    dplyr::select(c("V1", "V2", "V3"))
  colnames(tmp)[1:3] <- c("chr", "start", "end")
  tmp <- tmp %>%
    mutate(site = paste(chr,start, end, sep="_"))
  return(tmp)
})
names(regionFiles) <- gsub(c(temp), pattern = ".bed", replacement = "")
rm(temp)


#Load CpG level data for all sites tested for age
regionToCpG=read.table("/path/to/MASH/bedfiles/all_sites_age_CpG.bed")

regionToCpG = regionToCpG %>% mutate(site=paste(V1,V2,V3,sep="_"))

# Add blank columns for each region type and then add 1's where the site is in that region (both DMsites and non-DMsites)
reg_regions <- c("gene_body", "promoters", "cgIslands", "unassigned")
for (i in c(reg_regions)) { # make new column for each regulatory region, fill with 0
  regionToCpG[, paste0(i)] <- 0
}

# Fill with 1 if site is in the region

cpgislands<-regionFiles[[1]]
genes<-regionFiles[[2]]
promoters<-regionFiles[[3]]

regionToCpG[regionToCpG$site %in% cpgislands$site,]$cgIslands <- 1
regionToCpG[regionToCpG$site %in% genes$site,]$gene_body <- 1
regionToCpG[regionToCpG$site %in% promoters$site,]$promoters <- 1

regionToCpG[regionToCpG$gene_body == 0 & regionToCpG$promoters == 0 &
              regionToCpG$cgIslands == 0,]$unassigned <- 1

### Sites significant in at least one tissue ###

#match the annotations of all sites with which sites are significant 
sigsites<- read.table('/path/to/MASH/bedfiles/age_sharing.txt', header=T)

sigsites<-sigsites %>% 
  separate(site, c("type", "chr", "start", "stop")) %>%
  mutate(chrnew = paste0("chr",chr)) %>%
  mutate(site= paste(chrnew, start, stop, sep="_"))


## Differentiate between sites increasing/decreasing methylation with age
sigsites_pos <- sigsites %>% filter(focal_pm >0)
sigsites_neg <- sigsites %>% filter(focal_pm <0)

regionToCpG <- regionToCpG %>% mutate(region = paste(V4,V5,V6,sep="_"))

regionToCpG <- regionToCpG %>% mutate(sig_pos = ifelse(region %in% sigsites_pos$site,1,0),
                                      sig_neg = ifelse(region %in% sigsites_neg$site,1,0))

gene_enrich_pos<-fisher.test(table(regionToCpG$gene_body, regionToCpG$sig_pos))
promoter_enrich_pos<-fisher.test(table(regionToCpG$promoters, regionToCpG$sig_pos))
cgIslands_enrich_pos<-fisher.test(table(regionToCpG$cgIslands, regionToCpG$sig_pos)) #all sites enriched for cpIslands
unassigned_enrich_pos<-fisher.test(table(regionToCpG$unassigned, regionToCpG$sig_pos))

allenrich_pos<-as.data.frame(do.call(rbind, list(unlist(gene_enrich_pos), unlist(promoter_enrich_pos), unlist(cgIslands_enrich_pos),unlist(unassigned_enrich_pos))))

allenrich_pos$type<- c("gene_body", "promoters", "cgIslands","unassigned")

colnames(allenrich_pos)<-c("pvalue", "conf.low", "conf.high", "OR", "null", "sided", "method", "table", "type")

allenrich_pos$OR<- as.numeric(allenrich_pos$OR)
allenrich_pos$pvalue<- as.numeric(allenrich_pos$pvalue)
allenrich_pos$conf.low<- as.numeric(allenrich_pos$conf.low)
allenrich_pos$conf.high<- as.numeric(allenrich_pos$conf.high)

enrich_hyper <- ggplot(data=allenrich_pos %>% mutate(type = recode(type, "gene_body"="gene body")),
                       aes(x=type, y=log2(OR), color=type)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=log2(conf.low), ymax=log2(conf.high))) +
  # geom_segment(aes(x=type, xend = type, y=log2(conf.low),yend=log2(conf.high)),linewidth=2)+
  scale_color_brewer(palette = "Dark2")+
  ylim(-2.1,2.8)+
  labs(y="log2-transformed odds ratio",x="",title="A")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, vjust=1,hjust=1, color="black",size=12),
        axis.text.y=element_text(color="black",size=12),
        axis.title.y=element_text(color="black",size=12),
        plot.title = element_text(color = "black",size=16),
        legend.position = "none") +
  geom_hline(yintercept=0,linetype="dashed")

gene_enrich_neg<-fisher.test(table(regionToCpG$gene_body, regionToCpG$sig_neg))
promoter_enrich_neg<-fisher.test(table(regionToCpG$promoters, regionToCpG$sig_neg))
cgIslands_enrich_neg<-fisher.test(table(regionToCpG$cgIslands, regionToCpG$sig_neg))
unassigned_enrich_neg<-fisher.test(table(regionToCpG$unassigned, regionToCpG$sig_neg))

allenrich_neg<-as.data.frame(do.call(rbind, list(unlist(gene_enrich_neg), unlist(promoter_enrich_neg), unlist(cgIslands_enrich_neg),unlist(unassigned_enrich_neg))))

allenrich_neg$type<- c("gene_body", "promoters", "cgIslands","unassigned")

colnames(allenrich_neg)<-c("pvalue", "conf.low", "conf.high", "OR", "null", "sided", "method", "table", "type")

allenrich_neg$OR<- as.numeric(allenrich_neg$OR)
allenrich_neg$pvalue<- as.numeric(allenrich_neg$pvalue)
allenrich_neg$conf.low<- as.numeric(allenrich_neg$conf.low)
allenrich_neg$conf.high<- as.numeric(allenrich_neg$conf.high)

enrich_hypo <- ggplot(data=allenrich_neg%>% mutate(type = recode(type, "gene_body"="gene body")),
                      aes(x=type, y=log2(OR), color=type)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=log2(conf.low), ymax=log2(conf.high))) +
  scale_color_brewer(palette = "Dark2")+
  ylim(-2.1,2.8)+
  labs(y="",x="",title="B")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, vjust=1,hjust=1, color="black",size=12),
        axis.text.y=element_text(color="black",size=12),
        axis.title.y=element_text(color="black",size=12),
        plot.title = element_text(color = "black",size=16),
        legend.position = "none") +
  geom_hline(yintercept=0,linetype="dashed")

ggpubr::ggarrange(enrich_hyper,enrich_hypo, align = "v")
# ggsave("/path/to/Figures/FigS4.png", width=6,height=4)

#### Export result
all_res = rbind(allenrich_neg %>% mutate(direction_of_change="hypomethylating"),
      allenrich_pos %>% mutate(direction_of_change="hypermethylating"))
all_res = all_res %>% select(direction_of_change,type,OR,conf.low,conf.high,pvalue)
# write.csv(all_res, "/path/to/SupMat/TableS5.csv", row.names = F, quote = F)

# First, select only the necessary columns then pivot
transposed_data <- all_res %>%
  select(type, direction_of_change, OR, conf.low, conf.high) %>%
  pivot_wider(names_from = direction_of_change, values_from = c(OR, conf.low, conf.high))

ggplot(transposed_data, aes(x=log2(OR_hypomethylating), y=log2(OR_hypermethylating), color=type)) +
  geom_point() +
  # x-axis thin segments
  geom_segment(aes(xend = log2(OR_hypomethylating), y = log2(OR_hypomethylating), yend = log2(OR_hypomethylating)), size = 0.5) +
  # y-axis thin segments
  geom_segment(aes(x = log2(OR_hypomethylating), xend = log2(OR_hypomethylating), yend = log2(OR_hypermethylating)), size = 0.5) +
  # labels with white background
  geom_label(aes(label = type), size = 5, color = "black", fill = "white", hjust = -0.3) +
  geom_vline(xintercept = 0, linetype="dashed")+geom_hline(yintercept = 0, linetype="dashed")+
  scale_color_brewer(palette = "Dark2") +
  theme_bw() +
  theme(axis.text=element_text(color="black",size=12),
        axis.title=element_text(color="black",size=12),
        legend.position = "none") +
  labs(x="log2-transformed odds ratio hypomethylating", y="log2-transformed odds ratio hypermethylating")

#### Testing for enrichment for chromHMM annotations

# The aim here is to intersect the correct annotations for the correct tissue, and each time
# identify which sites where differentially methylated for that tissue.

rm(list = ls())

library_list <- c("corrplot","svglite","tidyverse","RColorBrewer")
lapply(library_list, require, character.only=TRUE)

# Define ggplot theme upfront
tissue_plot <- sort(c("whole blood","spleen","omental adipose","heart","kidney","lung","adrenal","thymus","thyroid", "pituitary", "liver", "skeletal muscle"))
extended_palette <- setNames(colorRampPalette(brewer.pal(8, "Dark2"))(12),tissue_plot)

## Load output from bedtools intersect
temp <- list.files(path = "/path/to/MASH/bedfiles/",
                   pattern="all_sites_age_.*chromHMM\\.bed")

# read in files into a list
chromFiles = lapply(temp, function(x){
  tmp <- read.delim(file = paste0("/path/to/MASH/bedfiles/", x),
                    sep = '\t', header = F)
  colnames(tmp)[1:11] <- c("chr", "start", "end","chr_region","start_region","end_region","chr_state","start_state","end_state","chrom","bp_overlap")
  tmp <- tmp %>%
    mutate(site = paste(chr,start, end, sep="_"))
  return(tmp)
})
names(chromFiles) <- gsub(c(temp), pattern = ".bed", replacement = "")
rm(temp)

# Format so the site names match pqlseq output so can add region information directly to pqlseq output'

all_site_bed<-read.table("/path/to/MASH/bedfiles/all_sites_age_CpG.bed",header=F)

all_site_bed <- all_site_bed %>%
  mutate(site = paste(V1,V2, V3, sep="_"),
         region = paste(V4,V5, V6, sep="_"))

# Add blank columns for each region type and then add 1's where the site is in that region (both DMsites and non-DMsites)
chrom_anno <- c("1_TssA", "2_TssAFlnk","3_TxFlnk","4_Tx","5_TxWk","6_EnhG","7_Enh","8_ZNF/Rpts","9_Het",
                 "10_TssBiv","11_BivFlnk","12_EnhBiv","13_ReprPC","14_ReprPCWk","15_Quies","unassigned")

for (i in chrom_anno) { 
  all_site_bed[, paste0(i)] <- 0 
} 

all_sites_alltissues <- list()

# Fill the levels of the list with dataframe for each tissue and intersecting the corresponding marks.
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

# Fill unassigned (check rows where all elements are 0 using rowSums)
all_sites_alltissues <- lapply(all_sites_alltissues, function(x) {
  df_subset <- x[,chrom_anno]
  x[rowSums(df_subset, na.rm=TRUE) == 0, "unassigned"] <- 1
  return(x)
})

# match the annotations of all sites with which sites are significantly associated with age

sigsites <- read.table('/path/to/MASH/bedfiles/age_sharing.txt', header=T)

hist(sigsites$sig_tissues)

#change to similar site format
sigsites_format<- sigsites %>%
  separate(site, c("type", "chr", "start", "stop")) %>%
  mutate(chrnew = paste0("chr",chr)) %>%
  mutate(site = paste(chrnew, start, stop, sep="_"))

tissue_oi <- c("whole_blood","spleen","omental_at","heart","kidney","lung","adrenal","thymus","thyroid","pituitary","liver","skeletal_muscle")

#intersect the annotation datasets with significance for age-association
for(level in names(all_sites_alltissues)){
  
  #isolate tissue level for that 
  if(grepl("adipose",level)){
    tissue_level <- "omental_at"
  }else{
    tissue_level <- tissue_oi[which(sapply(tissue_oi,function(tissue)grepl(paste0(".*",tissue,"*"), level)))]
  }
  
  #extract significantly age-associated sites in the tissue_level
  sig_in_tissue <- sigsites_format[sigsites_format[[tissue_level]]<0.05,"site"]
  # sign_in_tissues <- ifelse(sigsites_format[sigsites_format[[focal_pm]]<0,"site"],"neg","pos")
  
  # convert to the corresponding CpG level data
  cpg_sig_in_tissue <- all_site_bed %>% filter(region %in% sig_in_tissue)
  
  all_sites_alltissues[[level]]$sig <- ifelse(all_sites_alltissues[[level]]$site %in%
                                                cpg_sig_in_tissue$site, 1, 0)
}

# Perform enrichment test with Fisher's test 
tissue_chromHmm_enrich <- list()

for(i in seq_along(all_sites_alltissues)) {
  
  x <- all_sites_alltissues[[i]]
  tissue <- names(all_sites_alltissues)[i]
  tissue_chromHmm_enrich[[tissue]] <- list()
  
  for(annotation in c(chrom_anno, "unassigned")){
    
    annotation_enrich <- fisher.test(table(x[[annotation]], x$sig))
    
    # Convert the result to an unlisted vector
    annotation_df <- as.data.frame(t(unlist(annotation_enrich)), stringsAsFactors = FALSE)
    
    colnames(annotation_df)<-c("pvalue", "conf.low", "conf.high", "OR", "null", "sided", "method", "table")
    
    annotation_df$annotation <- annotation
    
    annotation_df$chromHMM <- tissue
    
    # Add the result to the list
    tissue_chromHmm_enrich[[names(all_sites_alltissues)[i]]][[annotation]] <- annotation_df
  }
}

chrommHmm_fischer_results <- do.call(rbind, lapply(tissue_chromHmm_enrich,
                                                   function(tissue_results){
                                                     do.call(rbind,tissue_results)
                                                   }))


chrommHmm_fischer_results$OR<- as.numeric(chrommHmm_fischer_results$OR)
chrommHmm_fischer_results$pvalue<- as.numeric(chrommHmm_fischer_results$pvalue)
chrommHmm_fischer_results$conf.low<- as.numeric(chrommHmm_fischer_results$conf.low)
chrommHmm_fischer_results$conf.high<- as.numeric(chrommHmm_fischer_results$conf.high)
chrommHmm_fischer_results$sig05 <- ifelse(chrommHmm_fischer_results$pvalue<0.05,"Y","N")
chrommHmm_fischer_results$FDR <- as.numeric(p.adjust(chrommHmm_fischer_results$pvalue, method="BH"))
chrommHmm_fischer_results$FDR05 <- ifelse(chrommHmm_fischer_results$FDR<0.05,"Y","N")

chrommHmm_fischer_results$annotation <- factor(chrommHmm_fischer_results$annotation, levels=chrom_anno)

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
chrommHmm_fischer_results <- chrommHmm_fischer_results %>%
  mutate(tissue = sapply(chromHMM, map_chromHMM_to_tissuename))

chrommHmm_fischer_results = chrommHmm_fischer_results %>%
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

export_res=chrommHmm_fischer_results %>% select(tissue, annotation,OR,conf.low,conf.high,pvalue,sided,method,chromHMM)
rownames(export_res) <- NULL
# write.csv(export_res %>% select(-c(method,sided,chromHMM)), "/path/to/SupMat/TableS6.csv",row.names=F, quote = F)

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
ggsave("/path/to/Figures/Fig2B.pdf", width=7.5,height=7.5)


#### Are tissue specific regions over represented within age-associated sites?

rm(list=ls())

library_list <- c("corrplot","svglite","tidyverse","RColorBrewer","parallel")
lapply(library_list, require, character.only=TRUE)

# Define ggplot theme upfront
tissue_plot <- sort(c("whole blood","spleen","omental adipose","heart","kidney","lung","adrenal","thymus","thyroid", "pituitary", "liver", "skeletal muscle"))
extended_palette <- setNames(colorRampPalette(brewer.pal(8, "Dark2"))(12),tissue_plot)

##############################################################################################################################
#
# SECTION 1: INTERSECTION OF TISSUE AGE ASSOCIATED SITES AND TISSUE-SPECIFIC MARKERS
#
##############################################################################################################################

# Load differentially methylated regions in one tissue compared to all other tissues
tDMR <- read.table("/path/tissue_markers/pairwise_methcomparisons/tissuespecific_methylation.txt",header=TRUE)

tissue_type <- c("whole_blood","spleen","omental_at","heart","kidney","lung","adrenal","thymus","thyroid","pituitary","liver","skeletal_muscle")

# Load the list of age-associated sites for all tissues (i.e., any sites that passed the lfsr in that tissue)
coeff_intercept_list <- readRDS("/path/to/MASH/bedfiles/tissue_age_associated_sites_list.rds")

# Load all sites tested for age
allsites_age <- read.table('/path/to/MASH/bedfiles/all_sites_age.txt',header=TRUE)

# Load all sites tested for tissue markers
for(i in tissue_type){
  filename <- gsub("XXX",i,"/path/to/tissues_meth/XXX_meth/Regions_pmeth_full_XXX_1000_14T.rds")
  # load data
  r <- readRDS(filename)
  cov<-r$coverage
  assign(paste(i, "cov", sep="_"), cov)
}
dfs <- list(liver_cov, kidney_cov, lung_cov, heart_cov, omental_at_cov, spleen_cov, adrenal_cov, thymus_cov, thyroid_cov, pituitary_cov, whole_blood_cov, skeletal_muscle_cov)  # Add all 12 dataframes to this list

# Find shared row names across all dataframes
allsites_markers <- Reduce(intersect, lapply(dfs, rownames)) #179,969 sites measured across all 12 tissues

# We can only test for the effect of age on tissue markers that were included in allsites_markers in the first place
markers_tested_for_age <- intersect(allsites_markers,allsites_age$regions)

tDMR_aDMR <- tDMR %>% filter(sites %in% markers_tested_for_age)

setdiff(allsites_age$regions,allsites_markers) # All sites tested for age-association were included
## in the test for markers so we can merge the data
## taking allsites_age as the reference

res_list <- lapply(coeff_intercept_list, function(x){
  
  result <- allsites_age %>% left_join(x %>%
                                  dplyr::select(region, beta, sign_beta),
                                join_by(regions==region)) %>% 
  dplyr::rename(beta_age=beta,sign_beta_age=sign_beta) %>%
  mutate(tissue=unique(x$tissue)) %>% 
  merge(.,tDMR_aDMR, by.x=c("regions", "tissue"), by.y=c("sites","sig_tissue"),
        all.x=TRUE) %>% 
  mutate(aDMR = ifelse(is.na(beta_age),0,1),
         tDMR = ifelse(is.na(mean_beta),0,1))
  return(result)
}
)

fisher_list <- list()
fisher_list <- lapply(res_list, function(x){
  
  # Create the contingency table
  cbinding <- table(x$aDMR, x$tDMR)
  
  # # Check if Haldane-Anscombe correction should be applied
  # if (any(cbinding == 0)) {
  #   cbinding <- cbinding + 1
  # }
  
  # Perform the Fisher exact test 
  spe_enrich <- fisher.test(cbinding)
  
  spe_enrich<-as.data.frame(do.call(rbind, list(unlist(spe_enrich))))  
  
  tissue_value <- unique(x$tissue)
  spe_enrich$tissue <- tissue_value
  
  return(spe_enrich)
}
)

tissues_spe_enrich<- as.data.frame(do.call(rbind, fisher_list))

colnames(tissues_spe_enrich)<-c("pvalue", "conf.low", "conf.high", "OR", "null", "sided", "method", "table", "tissue")

tissues_spe_enrich$OR<- as.numeric(tissues_spe_enrich$OR)
tissues_spe_enrich$pvalue<- as.numeric(tissues_spe_enrich$pvalue)
tissues_spe_enrich$conf.low<- as.numeric(tissues_spe_enrich$conf.low)
tissues_spe_enrich$conf.high<- as.numeric(tissues_spe_enrich$conf.high)

tissues_spe_enrich = tissues_spe_enrich %>% mutate(tissue = recode(tissue,"omental_at"="omental adipose",
                                                                   "skeletal_muscle"="skeletal muscle",
                                                                   "whole_blood"="whole blood"))

ggplot(data=tissues_spe_enrich, aes(x=tissue, y=log2(OR), color=tissue)) +
  geom_point() +
  geom_errorbar(aes(ymin=log2(conf.low), ymax=log2(conf.high)))+
  scale_color_manual(values = extended_palette)+
  geom_text(
    aes(
      label = ifelse(pvalue < 0.05, "*", ""), # Add asterisk conditionally
      y = ifelse(log2(conf.low) < 0, log2(conf.low) - 0.2, log2(conf.high) + 0.05), # Position above the bar
      group = tissue # Ensure proper grouping
    ),
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
# ggsave("/path/to/MASH/Figures/Fig2G.pdf", width=6,height=5)


#### Are tissue specific regions over represented within tissue specific aDMRs?

rm(list=ls())

library_list <- c("corrplot","svglite","tidyverse","RColorBrewer")
lapply(library_list, require, character.only=TRUE)

# Define ggplot theme upfront
tissue_plot <- sort(c("whole blood","spleen","omental adipose","heart","kidney","lung","adrenal","thymus","thyroid", "pituitary", "liver", "skeletal muscle"))
extended_palette <- setNames(colorRampPalette(brewer.pal(8, "Dark2"))(12),tissue_plot)


# Load differentially methylated regions in one tissue compared to all other tissues
tDMR <- read.table("/path/to/tissue_markers/pairwise_methcomparisons/tissuespecific_methylation.txt",header=TRUE)

tissue_type <- c("whole_blood","spleen","omental_at","heart","kidney","lung","adrenal","thymus","thyroid","pituitary","liver","skeletal_muscle")

# Load the list of age-associated sites for all tissues (i.e., any sites that passed the lfsr in that tissue)
coeff_intercept_list <- readRDS("/path/to/MASH/bedfiles/tissue_age_associated_sites_list.rds")

# Load all sites tested for age
allsites_age <- read.table('/path/to/MASH/bedfiles/all_sites_age.txt',header=TRUE)

# Load all sites tested for tissue markers
for(i in tissue_type){
  filename <- gsub("XXX",i,"/path/to/tissues_meth/XXX_meth/Regions_pmeth_full_XXX_1000_14T.rds")
  # load data
  r <- readRDS(filename)
  cov<-r$coverage
  assign(paste(i, "cov", sep="_"), cov)
}
dfs <- list(liver_cov, kidney_cov, lung_cov, heart_cov, omental_at_cov, spleen_cov, adrenal_cov, thymus_cov, thyroid_cov, pituitary_cov, whole_blood_cov, skeletal_muscle_cov)  # Add all 12 dataframes to this list

# Find shared row names across all dataframes
allsites_markers <- Reduce(intersect, lapply(dfs, rownames)) #179,969 sites measured across all 12 tissues

# We can only test for the effect of age on tissue markers that were included in allsites_markers in the first place
markers_tested_for_age <- intersect(allsites_markers,allsites_age$regions)

tDMR_aDMR <- tDMR %>% filter(sites %in% markers_tested_for_age)

aDMRs<-do.call(rbind, coeff_intercept_list)

tDMR_aDMR <- merge(tDMR_aDMR, aDMRs, by.x=c("sites", "sig_tissue"),
                   by.y=c("region", "tissue"), all.x=TRUE)

tDMR_aDMR$aDMR <- ifelse(is.na(tDMR_aDMR$beta),0,1)

tDMR_aDMR = tDMR_aDMR %>% dplyr::rename(beta_age=beta,
                                 beta_marker=mean_beta,
                                 tissue_marker=sig_tissue,
                                 sign_beta_age=sign_beta)

# Proportion of markers affected by age changes
tDMR_aDMR %>% group_by(tissue_marker) %>% 
  summarize(prop_affected = sum(aDMR)/length(aDMR)) %>% 
  arrange(prop_affected)

# Summary across tissues
tDMR_aDMR %>% group_by(tissue_marker) %>% 
  summarize(prop_affected = sum(aDMR)/length(aDMR)) %>% 
  summarize(av=mean(prop_affected),sd=sd(prop_affected))

tDMR_aDMR <- tDMR_aDMR %>% mutate(tissue_marker = recode(tissue_marker,
                                                         "skeletal_muscle"="skeletal muscle",
                                                         "whole_blood"="whole blood",
                                                         "omental_at"="omental adipose"))

#### FigS7
ggplot(tDMR_aDMR, aes(x=beta_marker,y=beta_age))+
  geom_point(alpha=0.2)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  labs(x="effect sizes for tissue marker", y="effect sizes for age")+
  facet_wrap(~tissue_marker, scales="free")+
  theme_bw()+
  theme(axis.text = element_text(size=14, color="black"),
        axis.title = element_text(size=14, color="black"),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(size=12, face="bold"))
ggsave("/path/to/Figures/FigS7.png", width=10,height=6)

### Divergent vs convergent sites
tDMR_aDMR = tDMR_aDMR %>%
  mutate(direction_of_change = case_when(
    beta_age<0 & beta_marker>0 ~ "loss_specificity",
    beta_age>0 & beta_marker<0 ~ "loss_specificity",
    beta_age<0 & beta_marker<0 ~ "gain_specificity",
    beta_age>0 & beta_marker>0 ~ "gain_specificity")
  ) %>% 
  mutate(direction_of_change_bino = ifelse(direction_of_change=="loss_specificity",1,0))

tissue_summ_direction_of_change <- tDMR_aDMR %>% 
  filter(aDMR==1) %>% 
  group_by(tissue_marker) %>% 
  summarize(prop_loss_specificity = sum(direction_of_change=="loss_specificity")/n(),
            n_sites = n()) %>% 
  arrange(-prop_loss_specificity) %>% 
  mutate(tissue_marker = factor(tissue_marker, levels=tissue_marker, ordered = TRUE))

tissue_summ_direction_of_change = tissue_summ_direction_of_change %>%
  mutate(tissue_marker = recode(tissue_marker,
                             "skeletal_muscle"="skeletal muscle",
                             "whole_blood"="whole blood",
                             "omental_at"="omental adipose"))

#### Fig2H
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
ggsave("/path/to/Figures/Fig2H.pdf", width=8.6,height=6.5,dpi=300)
