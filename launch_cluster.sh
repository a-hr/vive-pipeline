#!/usr/bin/env bash

#SBATCH --time=20:00:00
#SBATCH --mem=8G

module load Singularity Nextflow Go Mamba

# limit the resources of the JVM
export NXF_OPTS="-Xms500M -Xmx8G"

# launch the main process
nextflow run main.nf -profile cluster \
    --run_name "vive-jul23_cDNA-pipeline" \
    --input_fastq "inputs-0/*.fastq*" \
    --output_dir "output-jul23_cDNA" \
    --is_rna false

nextflow run main.nf -profile cluster \
    --run_name "vive-dec23_cDNA-pipeline" \
    --input_fastq "inputs-1/*.fastq*" \
    --output_dir "output-dec23_cDNA" \
    --is_rna false

nextflow run main.nf -profile cluster \
    --run_name "vive-jul23_dRNA-pipeline" \
    --input_fastq "inputs-2/*.fastq*" \
    --output_dir "output-jul23_dRNA" \
    --is_rna true

nextflow run main.nf -profile cluster \
    --run_name "vive-dec23_dRNA-pipeline" \
    --input_fastq "inputs-3/*.fastq*" \
    --output_dir "output-dec23_dRNA" \
    --is_rna true