#!/usr/bin/env bash

#SBATCH --time=1:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=nextflow-cDNA-fw_new
#SBATCH --output=logs/%x-%j.out

module load Singularity Go Nextflow 
mkdir -p logs

# limit the resources of the JVM
export NXF_OPTS="-Xms500M -Xmx4G"

# set the log directory
export NXF_LOG="$PWD/logs"

# launch the main process
nextflow run main.nf -profile cluster -params-file inputs/params/input_params_cDNA_fw.yaml -resume
