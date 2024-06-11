# VIVE Biotech pipeline


## Introduction

A pipeline to process Nanopore sequencing data. It performs the following steps:

1. QC of the FASTQ files
2. Removal of adapters and primers
3. Alignment of the reads to the reference genome
4. Optional: alignment to a secondary reference genome (i.e. spike-ins)
5. Generation of BAM files
6. Splice-site scoring on the reference genome

## Installation

The pipeline is written in Nextflow, a workflow manager that allows to run the pipeline in a wide variety of systems. It is configured to be run either on a SLURM-managed HPC cluster or a local machine, though it can be run on a cloud instance or using other workload managers by editing the configuration file according to [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html#config-scopes).

## Requirements

If running in local, you will need a PC with:

- At least 32GB of RAM (ideally 64GB)
- A multi-core CPU (ideally i7 or better)
- At least 500GB of free disk space (recommended to have an SSD)
- Docker or Singularity installed
- Nextflow installed OR EPI2ME installed


### Installing in a local machine

In your local machine, there are two ways of running the pipeline:

1. Using EPI2ME: the easiest way for those without bioinformatics experience. It's just a graphical interface that allows you to run the pipeline. Below will be explained how to install it.
2. Using the command line: for those with bioinformatics experience. It follows the same procedure as running in a cluster, so it will be explained in the `Installing in a cluster` section.

#### Installing EPI2ME

[Install EPI2ME](https://labs.epi2me.io/installation/) on your system and follow the instructions on the app to install all the dependencies (Java, Docker and Nextflow). To add the pipeline to your saved workflows simply copy this repository's URL and paste it on the "Add workflow" section of the EPI2ME interface.

### Installing in a cluster

If using Nextflow/nf-core, clone the repository and install the basic dependencies (Nextflow). The easiest way to do so is using conda. The pipeline can be run on any system that supports Docker or Singularity.

```bash
git clone https://github.com/a-hr/vive-pipeline.git
```

The internal dependencies of the pipeline are managed by Nextflow, so you don't need to worry about them. If for some reason Nextflow fails to download them when using Singularity (they are provided as Docker containers), you can manually download them with the Makefile:

```bash
# make sure you have Singularity installed and available
make pull
```

The pipeline is especially tailored to be run on a HPC cluster, though it can seamlessly be run on a local machine and, with some configuration, on a cloud instance.

## Usage

### Running with EPI2ME

1. Open EPI2ME and go to the "Workflows" tab.
2. Select the workflow.
3. Fill in the parameters.
4. In the `profile` section, make sure its set to `standard` (the default), which runs the backend on top of Docker containers. If you are using Singularity, set it to `local_singularity`.
5. Run the pipeline.

### Running with the command line

1. Go to the directory where you cloned the repository.
2. Fill in the parameters in the `input_params.yaml` file.
3. Make sure your system has Docker/Singularity and Nextflow available.
4. Run the pipeline in the cluster/local machine with the following command:

```bash
sbatch launch_cluster.sh  # for SLURM-managed HPC clusters
bash launch_local.sh  # for local machines
```

> The cluster launch script is configured to run the pipeline in a SLURM-managed HPC cluster. If you are using another workload manager, you will need to edit the script accordingly.  

***If you are running the pipeline in a local machine, you can run it with the following command:***

```bash
# with Docker
nextflow run main.nf -profile standard -params-file input_params.yaml
# or with Singularity
nextflow run main.nf -profile local_singularity -params-file input_params.yaml
```

## Parameters

Below, a description of the parameters that can either be set in the `input_params.yaml` file or provided through the EPI2ME interface.

### Mandatory parameters

- `experiment_name`: all the output files will be prefixed with this name
- `input_fastqs`: folder containing the FASTQ files output by the basecaller. *Note that they will all be processed together, so make sure they are from the same run.*
- `plasmid_ref_fa`: reference genome of the sequenced target
- `output_dir`: folder where the output files will be saved (only available if running through the command line)
- `is_rna`: whether the input is RNA or DNA. If RNA, the FASTQ files will be accordingly processed.

### Quality control parameters

- `min_len`: expected minimum length of the reads (all reads below this length will be discarded)
- `max_len`: expected maximum length of the reads (all reads above this length will be discarded)

### Secondary target

- `assess_secondary`: whether to process a secondary target
- `secondary_ref_fa`: reference genome of the secondary target

# Optional arguments

- `get_bams`: whether to generate BAM files as output
