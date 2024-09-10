#!/bin/bash
#SBATCH --job-name=fastq_pipeline
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

source /cluster/apps/x86_64/packages/modules-tcl/init/bash

# Path to your configuration file
CONFIG_FILE="/home/lowbh-int/ipsc_differentiation/SRA_config.yaml"

# Run the Python script with srun
source /home/lowbh-int/anaconda3/bin/activate; srun python /home/lowbh-int/sra_data_pipeline/scripts/fastq_download_unstable.py -c ${CONFIG_FILE}