#!/bin/bash

# Define your base directory here
base_path="/path/to/base_path"

# Define subdirectories relative to base_path
fastq_path="${base_path}/fastq"
out_path="${base_path}/mapped"
cov_out="${base_path}/cov"
trim_path="${base_path}/trimmed"
finalcov_out="${base_path}/final_cov"

# Define genome folder path (update this to your genome folder)
genome_path="/path/to/bismark/mmul10/"

# Load required modules
module load bowtie2-2.4.2-gcc-11.2.0
module load samtools-1.13-gcc-11.2.0
module load py-cutadapt-2.10-gcc-11.2.0

# Get the sample ID from the lids file based on SLURM_ARRAY_TASK_ID
sampleID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${fastq_path}/lids")

# Create directories if they do not exist
mkdir -p "${out_path}" "${cov_out}" "${trim_path}" "${finalcov_out}"

# Add custom program paths to PATH if needed
export PATH=$PATH:/path/to/programs/bin:/path/to/programs/TrimGalore-0.6.6:/path/to/programs/bismark2/Bismark-0.24.0

# Trim adapters with trim_galore, flexible matching for R1 fastq files
trim_galore --rrbs -j 8 -o "${trim_path}" \
    --basename "${sampleID}" --gzip \
    ${fastq_path}/${sampleID}*R1*fastq.gz

# Map reads with Bismark
bismark --genome_folder "${genome_path}" \
    -o "${out_path}" \
    --score_min L,0,-0.6 -R 10 \
    --parallel 12 "${trim_path}/${sampleID}*.gz"

# Extract methylation data
bismark_methylation_extractor -s -o "${cov_out}" \
    --gzip \
    --bedGraph --comprehensive --parallel 12 \
    --merge_non_CpG --genome_folder "${genome_path}" \
    "${out_path}/${sampleID}*.bam"

# Clean up unwanted files
rm "${cov_out}"/*"${sampleID}"*txt.gz
rm "${cov_out}"/*"${sampleID}"*bedGraph.gz

# Merge CpG coverage files into cytosine report
coverage2cytosine --merge_CpG --gzip --genome_folder "${genome_path}" \
    -o "${sampleID}" \
    --dir "${finalcov_out}/" \
    "${cov_out}/${sampleID}"*.cov.gz
