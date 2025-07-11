# CayoCrosstissueDNAm
**Repository for the project on multi-tissue DNA methylation in Cayo Santiago rhesus macaques (bulk-tissue rrbs).**

Analyses were performed using the Arizona State University SOL supercomputer (DOI: [10.1145/3569951.3597573](https://doi.org/10.1145/3569951.3597573)) on high-performance computing (HPC) clusters. We have aimed to generalize the code by removing system-specific references to installed software and modules. A list of required software and versions is provided below. On HPC systems, all required scripts and binaries must be accessible through the system's PATH.

All analyses were performed via the command line in a Bash environment, using the Slurm workload manager or through an RStudio interface for R version 4.4.0, with the following R package dependencies:

CRAN packages:
corrplot
RColorBrewer
pheatmap
svglite
ggpubr
grid
glmnet
jtools
lmerTest
mashr
ashr
flashier
MatrixGenerics
umap


<pre> ```r # Install BiocManager if not already installed if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager") # Install Bioconductor packages BiocManager::install(c( "bsseq", "BiocGenerics", "GenomicRanges", "GenomicFeatures", "PQLseq", "qvalue", "methyLImp2", "comethyl", "DelayedMatrixStats", "BiocParallel" )) # Install CRAN packages install.packages(c( "corrplot", "RColorBrewer", "pheatmap", "svglite", "ggpubr", "grid", "glmnet", "jtools", "lmerTest", "mashr", "ashr", "flashier", "MatrixGenerics", "umap" )) ``` </pre>
