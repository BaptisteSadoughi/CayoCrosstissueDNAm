#!/bin/bash

module load bowtie2-2.4.2-gcc-11.2.0
module load samtools-1.13-gcc-11.2.0
module load py-cutadapt-2.10-gcc-11.2.0

fastq_path=${1}

sampleID=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${fastq_path}/lids)

genome_path=/path/to/bismark/mmul10/
out_path=/base_path/mapped
cov_out=/base_path/cov
trim_path=/base_path/trimmed
finalcov_out=/base_path/final_cov

# Create output directories if they don't exist
mkdir -p "${out_path}" "${cov_out}" "${trim_path}" "${finalcov_out}"

export PATH=$PATH:/path/to/programs/bin:/path/to/programs/TrimGalore-0.6.6:/path/to/programs/bismark2/Bismark-0.24.0

## trim adaptors
trim_galore --rrbs -j 8 -o ${trim_path} \
    --basename ${sampleID} --gzip \
    ${fastq_path}/${sampleID}*R1*fastq.gz

## map the reads
bismark --genome_folder ${genome_path} \
    -o ${out_path} \
    --score_min L,0,-0.6 -R 10 \
    --parallel 12 ${trim_path}/${sampleID}*.gz

## extract methylation data
bismark_methylation_extractor -s -o ${cov_out} \
        --gzip \
        --bedGraph --comprehensive --parallel 12 \
        --merge_non_CpG --genome_folder ${genome_path} \
        ${out_path}/${sampleID}*.bam

## clean up unwanted files
rm ${cov_out}/*${sampleID}*txt.gz
rm ${cov_out}/*${sampleID}*bedGraph.gz

coverage2cytosine --merge_CpG --gzip --genome_folder $genome_path \
        -o ${sampleID} \
        --dir ${finalcov_out}/ \
        ${cov_out}/${sampleID}*.cov.gz
