# CayoCrosstissueDNAm
**Repository for the project on multi-tissue DNA methylation in Cayo Santiago rhesus macaques (bulk-tissue rrbs).**

Analyses were performed using the Arizona State University SOL supercomputer (DOI: [10.1145/3569951.3597573](https://doi.org/10.1145/3569951.3597573)) on high-performance computing (HPC) clusters. We have aimed to generalize the code by removing system-specific references to installed software and modules. A list of required software and versions is provided below. On HPC systems, all required scripts and binaries must be accessible through the system's PATH. We suggest that users create a main directory `base_path` containing the `Bioinformatic and R scripts` and `metadata` folders. 

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
sbatch --cpus-per-task=1 --mem=500G -p general -q public -t 1-00:00:00 base_path/Bioinformatic and R scripts/generate_bsseq.R
```

## Building genomic regions

CpG density based regions are generated using [Generate_CpG_region_methylation.R](https://github.com/BaptisteSadoughi/CayoCrosstissueDNAm/blob/main/Bioinformatic%20and%20R%20scripts/Generate_CpG_region_methylation.R). Several path and parameters must be defined at the start of the script.  The script requires to source [SupportFunctions_generate_region_methylation.R](https://github.com/BaptisteSadoughi/CayoCrosstissueDNAm/blob/main/Bioinformatic%20and%20R%20scripts/SupportFunctions_generate_region_methylation.R) to run helper functions. The output are tissue-specific methylation data saved in separate subdirectories.

Matrices of fully covered regions, and imputed percent methylation values are generated using [Percent_methylation_imputation.R](https://github.com/BaptisteSadoughi/CayoCrosstissueDNAm/blob/main/Bioinformatic%20and%20R%20scripts/Percent_methylation_imputation.R).

## Visualization of sample clustering with UMAP.
[UMAP_dimensionality_reduction.R](https://github.com/BaptisteSadoughi/CayoCrosstissueDNAm/blob/main/Bioinformatic%20and%20R%20scripts/UMAP_dimensionality_reduction.R).

## Tissue specific DNAm markers
Top tissue-specific markers are visualized using [Tissuemarkers_plot.R](https://github.com/BaptisteSadoughi/CayoCrosstissueDNAm/blob/main/Bioinformatic%20and%20R%20scripts/Tissuemarkers_plot.R), and enrichment for chrommHMM states is plotted using [Tissuemarkers_annotation.R](https://github.com/BaptisteSadoughi/CayoCrosstissueDNAm/blob/main/Bioinformatic%20and%20R%20scripts/Tissuemarkers_annotation.R).

## Age-associated differential methylation
Before running the scripts:

1. Unzip regions_to_cpgs_mapping.zip in the folder /metadata
   
### Modelling steps
Binomial mixed models testing for the association between age and methylation levels are performed using [Age_methylation_modelling.R](https://github.com/BaptisteSadoughi/CayoCrosstissueDNAm/blob/main/Bioinformatic%20and%20R%20scripts/Age_methylation_modelling.R), and effect size estimates are refined using Multivariate Adaptive Shrinkage in [Age_mashr.R](https://github.com/BaptisteSadoughi/CayoCrosstissueDNAm/blob/main/Bioinformatic%20and%20R%20scripts/Age_mashr.R).

### Annotation and enrichment
Before running the script [script_annotating_sites.sh](https://github.com/BaptisteSadoughi/CayoCrosstissueDNAm/blob/main/Bioinformatic%20and%20R%20scripts/script_annotating_sites.sh):

1. Set the path to your central directory (`base_path`) inside the script.
2. Download annotation files into a subdirectory (`base_path/annotations`)
   
Age-associated sites are intersected with genomic annotations and chrommHMM states using [script_annotating_sites.sh](https://github.com/BaptisteSadoughi/CayoCrosstissueDNAm/blob/main/Bioinformatic%20and%20R%20scripts/script_annotating_sites.sh).
```bash
sbatch base_path/Bioinformatic and R scripts/run_intersect_annotations.sh
```

Analyses are then performed using [Age_analysis.R](https://github.com/BaptisteSadoughi/CayoCrosstissueDNAm/blob/main/Bioinformatic%20and%20R%20scripts/Age_analysis.R).

## Tissue specific DNAm age clocks

## ELA-associated differential methylation

### Modelling steps
Binomial mixed models testing for the association between each ELA and methylation levels are performed using [ELA_methylation_modelling.R](https://github.com/BaptisteSadoughi/CayoCrosstissueDNAm/blob/main/Bioinformatic%20and%20R%20scripts/ELA_methylation_modelling.R), and effect size estimates are refined using Multivariate Adaptive Shrinkage in [ELA_mashr.R](https://github.com/BaptisteSadoughi/CayoCrosstissueDNAm/blob/main/Bioinformatic%20and%20R%20scripts/ELA_mashr.R).

### ELA analyses
Annotation of ELA-associated sites is perfomed following the same procedure as for age by simply modifying the input bed file. Annotation, enrichment analysis, linear models testing for the effect of tissue types on effect sizes, and comparison of age and ELA effects are all performed with [ELA_analyses](https://github.com/BaptisteSadoughi/CayoCrosstissueDNAm/blob/main/Bioinformatic%20and%20R%20scripts/ELA_analyses). 
