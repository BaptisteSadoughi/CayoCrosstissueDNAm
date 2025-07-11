# CayoCrosstissueDNAm
**Repository for the project on multi-tissue DNA methylation in Cayo Santiago rhesus macaques (bulk-tissue rrbs).**

Analyses were performed using the Arizona State University SOL supercomputer (DOI: [10.1145/3569951.3597573](https://doi.org/10.1145/3569951.3597573)) on high-performance computing (HPC) clusters. We have aimed to generalize the code by removing system-specific references to installed software and modules. A list of required software and versions is provided below. On HPC systems, all required scripts and binaries must be accessible through the system's PATH.

All analyses were performed via the command line in a Bash environment, using the Slurm workload manager or through an RStudio interface for R version 4.4.0, with the following R package dependencies:

**Required R packages**

*CRAN packages:*
- corrplot  
- RColorBrewer  
- pheatmap  
- svglite  
- ggpubr  
- grid  
- glmnet  
- jtools  
- lmerTest  
- mashr  
- ashr  
- flashier  
- MatrixGenerics  
- umap  

*Bioconductor packages:*
- bsseq  
- BiocGenerics  
- GenomicRanges  
- GenomicFeatures  
- PQLseq  
- qvalue  
- methyLImp2  
- comethyl  
- DelayedMatrixStats  
- BiocParallel  

**External software:** (must be installed as a system-level tool or loaded via module)
- bedtools
- bowtie2  
- samtools  
- cutadapt  
- trim_galore  
- bismark 

---

### Installation instructions

To install the required R packages, run the following commands in your R session:

```r
# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c("bsseq","BiocGenerics","GenomicRanges","GenomicFeatures","PQLseq","qvalue","methyLImp2","comethyl","DelayedMatrixStats","BiocParallel"
))

# Install CRAN packages
install.packages(c("corrplot","RColorBrewer","pheatmap","svglite","ggpubr","grid","glmnet","jtools","lmerTest","mashr","ashr","flashier","MatrixGenerics","umap"
))
```

## Processing of fastq files

This pipeline script [fastq_processing.sh](https://github.com/BaptisteSadoughi/CayoCrosstissueDNAm/blob/main/Bioinformatic%20and%20R%20scripts/fastq_processing.sh) performs read trimming, alignment, and methylation extraction on RRBS data. It is designed to run on HPC clusters with SLURM.

---

### Prerequisites and Setup

Before running the script, please prepare the following:

1. Create a central directory to hold data and outputs:  
   `base_path/`

2. Download raw FASTQ files into a subdirectory:
   `base_path/fastq/`

3. Inside the `base_path/fastq/` directory, create a plain text file named `lids` that lists all sample IDs (one per line).  
   These sample IDs should correspond exactly to the prefixes of the FASTQ files.  
   For example, if a FASTQ file is named `LID_106879_PID_15915*R1*.fastq.gz`, the `lids` file should include:  
   LID_106879_PID_15915

4. Ensure that all required software and tools are installed and available in your system PATH or loaded as modules.  
   Required tools include: `bowtie2`, `samtools`, `cutadapt`, `trim_galore`, `bismark`, `bedtools`  

5. Modify the script to set these paths according to your directory structure:
   - genome_path: path to the Bismark genome folder for mmul10 (e.g., /path/to/bismark/mmul10/). To generate the Bismark folder, see [tutorial](https://felixkrueger.github.io/Bismark/bismark/genome_preparation/) 
   - base_path: newly created central directory  
   - PATH: access to the modules

### Running the Script

Submit the script as a Slurm array job, specifying the total number of samples listed in the `lids` file:

```r
sbatch --cpus-per-task=12 --array=1-N mapping_and_methylation.sh
```

- Replace N with the number of lines (samples) in your `lids` file.  

---

## Generating the BSseq data

Process methylation coverage files using [generate_bsseq.R](https://github.com/BaptisteSadoughi/CayoCrosstissueDNAm/blob/main/Bioinformatic%20and%20R%20scripts/generate_bsseq.R) to generate a BSseq object for downstream analysis. A HDF5 backend is used for efficient storage. The resulting BSseq object is saved in the  as `Rhesus_bismarkBSseq.rds` in the base_path directory.

### Prerequisites and Setup:
Before running the script:

1. Set the path to your central directory (`base_path`) inside the script.
2. Ensure the macaque genome FASTA file is available and the path to it is correctly set in the script (`mmul10_fa`).
3. Required R packages: `Biostrings`, `bsseq`.

### Running the script

```bash
sbatch --cpus-per-task=1 --mem=500G -p general -q public -t 1-00:00:00 /path/to/generate_bsseq.R

---

## Building genomic regions
