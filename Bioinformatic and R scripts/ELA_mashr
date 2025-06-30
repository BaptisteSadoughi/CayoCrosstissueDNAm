library(dplyr)
library(mashr)

tissue_oi<-c("liver", "omental_at", "spleen", "kidney", "lung", "heart", "skeletal_muscle", "adrenal", "thyroid", "thymus", "pituitary", "whole_blood")

ELAlist <- c("MomAllLossType", "q_kinBirth_sub", "q_grp_sub", "PrimpIndex", "CompetingSibIndex", "RankIndex01", "CumlQuartIndex6_4_sub")

for(ELA in ELAlist){
  
#reads in all the model results
  for(tissue in tissue_oi){
    filename0 <- gsub("XXX",tissue,"/path/to/results/XXX_000_pqlseq.rds")
    filename <- gsub("000", ELA, filename0)
  
    # load data
  r <- readRDS(filename)
  
  assign(tissue, r)
  }

  print(ELA)

#find sites measured in all tissues
common_row_names <- Reduce(intersect, list(rownames(liver), rownames(lung), rownames(kidney), rownames(spleen), rownames(heart), rownames(adrenal), rownames(skeletal_muscle), rownames(omental_at), rownames(thyroid), rownames(thymus), rownames(pituitary), rownames(whole_blood))); length(common_row_names) #130,278 sites present in 12 tissues 

write.table(common_row_names, gsub("000", ELA,'/path/to/save/000_allsites.txt')) #130,278

# Filter data frames for sites tested in all tissues
filtered_liver <- liver[common_row_names, ]
filtered_lung <- lung[common_row_names, ]
filtered_kidney <- kidney[common_row_names, ]
filtered_spleen <- spleen[common_row_names, ]
filtered_heart <- heart[common_row_names, ]
filtered_adrenal <- adrenal[common_row_names, ]
filtered_skeletal_muscle <- skeletal_muscle[common_row_names, ]
filtered_omental_at <- omental_at[common_row_names, ]
filtered_thyroid <- thyroid[common_row_names, ]
filtered_thymus <- thymus[common_row_names, ]
filtered_pituitary <- pituitary[common_row_names, ]
filtered_whole_blood <- whole_blood[common_row_names, ]

# Find strong sites
strong_liver<-rownames(subset(filtered_liver, pvalue<0.01)) 
strong_lung<-rownames(subset(filtered_lung, pvalue<0.01)) 
strong_kidney<-rownames(subset(filtered_kidney, pvalue<0.01)) 
strong_spleen<-rownames(subset(filtered_spleen, pvalue<0.01)) 
strong_heart<-rownames(subset(filtered_heart, pvalue<0.01)) 
strong_adrenal<-rownames(subset(filtered_adrenal, pvalue<0.01)) 
strong_skeletal_muscle<-rownames(subset(filtered_skeletal_muscle, pvalue<0.01)) 
strong_omental_at<-rownames(subset(filtered_omental_at, pvalue<0.01))
strong_thyroid<-rownames(subset(filtered_thyroid, pvalue<0.01)) 
strong_thymus<-rownames(subset(filtered_thymus, pvalue<0.01)) 
strong_pituitary<-rownames(subset(filtered_pituitary, pvalue<0.01)) 
strong_whole_blood<-rownames(subset(filtered_whole_blood, pvalue<0.01)) 

strong.sites<- unique(c(strong_liver, strong_lung, strong_kidney, strong_spleen, strong_heart, strong_adrenal, strong_skeletal_muscle, strong_omental_at, strong_thyroid, strong_thymus, strong_pituitary, strong_whole_blood)); length(strong.sites)
write.table(strong.sites, gsub("000", ELA,'/path/to/save/000_strongsites.txt'))


# setup mashR
coef=as.data.frame(cbind(filtered_liver$beta,filtered_lung$beta, filtered_kidney$beta, filtered_spleen$beta, filtered_heart$beta, filtered_adrenal$beta, filtered_skeletal_muscle$beta, filtered_omental_at$beta, filtered_thyroid$beta, filtered_thymus$beta, filtered_pituitary$beta, filtered_whole_blood$beta))
rownames(coef)<- common_row_names

se=as.data.frame(cbind(filtered_liver$se_beta,filtered_lung$se_beta, filtered_kidney$se_beta, filtered_spleen$se_beta, filtered_heart$se_beta, filtered_adrenal$se_beta, filtered_skeletal_muscle$se_beta, filtered_omental_at$se_beta, filtered_thyroid$se_beta, filtered_thymus$se_beta, filtered_pituitary$se_beta, filtered_whole_blood$se_beta))
rownames(se)<- common_row_names

#make mash set for all sites
data.all = mash_set_data(as.matrix(coef), as.matrix(se))

#make subset of strong sites
data.strong = mash_set_data(as.matrix(coef[strong.sites,]), as.matrix(se[strong.sites,]))

# estimate correlations
V.simple = estimate_null_correlation_simple(data.all)
data.all.cor = mash_update_data(data.all, V=V.simple)
data.strong.cor = mash_update_data(data.strong, V=V.simple)

#run mashr
U.pca = cov_pca(data.strong.cor,5)
U.f = cov_flash(data.strong.cor)
U.ed = cov_ed(data.strong, c(U.f, U.pca))
U.c = cov_canonical(data.all.cor)
m   = mash(data.all.cor, c(U.c,U.ed))

# save output
write.table(get_lfsr(m), gsub("000", ELA,'/path/to/save/000_LFSR_mashr.txt'),row.names=F,sep='\t')
write.table(get_pm(m),gsub("000", ELA,'/path/to/save/000_pm_mashr.txt'),row.names=F,sep='\t')

lfsr=as.data.frame(get_lfsr(m))
pm=as.data.frame(get_pm(m))
lfsr$id<-1:dim(lfsr)[1]
pm$id<-1:dim(pm)[1]

names(lfsr)<-names(pm)<-c('liver','lung', 'kidney', 'spleen', 'heart', 'adrenal', 'skeletal_muscle', 'omental_at', 'thyroid', 'thymus', 'pituitary', 'whole_blood', 'id')

#set sig threshold 
threshold <- 0.10

# Count number of tissues (columns) below threshold for each row
lfsr$sig_tissues <- apply(lfsr < threshold, 1, sum)

table(lfsr$sig_tissues)

#calculate tissue sharing
lfsr_sig1<-subset(lfsr,sig_tissues>0)
lfsr_sig1$within2_v2<-NA
lfsr_sig1$focal_pm<-NA
lfsr_sig1$other_pm<-NA

pm_sig1<-subset(pm, id %in% lfsr_sig1$id)

for (i in 1:dim(lfsr_sig1)[1]){
  
	z<-min(lfsr_sig1[i,1:12])[1]
	focal<-pm[lfsr_sig1$id[i],which(lfsr_sig1[i,1:12] == z)][1]
	not_focal<-pm[lfsr_sig1$id[i],which(lfsr_sig1[i,1:12] != z)]
	not_focal2<-log2(as.numeric(not_focal)/as.numeric(focal)) #negative numbers returning NaN
	
	lfsr_sig1$within2_v2[i]<-length(which(not_focal2 > -1 & not_focal2<1))
	lfsr_sig1$focal_pm[i]<-mean(focal)
	lfsr_sig1$other_pm[i]<-mean(t(not_focal))
}

table(lfsr_sig1$within2_v2)

write.table(as.data.frame(lfsr_sig1), gsub("000", ELA,'/path/to/save/000_sharing.txt'),row.names=F,sep='\t')

# tissue-unique sites
lfsr_sig2<-subset(lfsr_sig1,sig_tissues==1 & within2_v2==0); dim(lfsr_sig2)
lfsr_sig2$tissue<-apply(lfsr_sig2[,1:12],1,function(x) names(lfsr_sig2)[which(x==min(x))] ) 
uniq_sites<-as.data.frame(table(lfsr_sig2$tissue))
colnames(uniq_sites)<-c("tissue", "sites")

write.table(as.data.frame(lfsr_sig2), gsub("000", ELA,'/path/to/save/000_unique.txt'),row.names=F,sep='\t')


}
