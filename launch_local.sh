#!/usr/bin/env bash

# make sure to have Nextflow and Docker/Singularity available
mkdir -p logs

# limit the resources of the JVM
export NXF_OPTS="-Xms500M -Xmx4G"

# set the log directory
export NXF_LOG="$PWD/logs"

# launch the main process (one of the options below, depending on your backend)

# docker
nextflow run main.nf -profile standard -params-file input_params.yaml -resume

# singularity
nextflow run main.nf -profile local_singularity -params-file input_params.yaml -resume
