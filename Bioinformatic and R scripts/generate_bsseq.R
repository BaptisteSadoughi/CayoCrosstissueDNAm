#!/bin/bash
# sbatch --cpus-per-task=1 --mem=500G -p general -q public -t 1-00:00:00 /path/to/this/script/generate_bsseq.R

library(Biostrings)
library(bsseq)

# Define your base directory here
base_path="/path/to/base_path"

# load macaque genome
mmul10_fa = readDNAStringSet("/path/to/genomes/bismark/Macaca_mulatta.Mmul_10.chr_only.fa",format="fasta")

## filter for the first 21 chromosomes (chr1-chrX)
names(mmul10_fa)[1:21]->chr
mmul10_fa=mmul10_fa[chr]

## find CG sites in the genome
loci = findLoci("CG", subject=mmul10_fa)

files <- list.files("${base_path}/final_cov", full.names = TRUE, pattern = ".cov.gz")

bismarkBSseq <- read.bismark(files = files,
                             strandCollapse = FALSE,
                             verbose = TRUE,
                             BACKEND = "HDF5Array",
                             loci = loci,
                             replace=TRUE,
                             rmZeroCov = FALSE,
                             dir = file.path(base_path, "hdf5"))

saveRDS(bismarkBSseq, file = file.path(base_path, "Rhesus_bismarkBSseq.rds"))

q(save = "no")
