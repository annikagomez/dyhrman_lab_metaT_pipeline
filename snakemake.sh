#!/bin/bash
#SBATCH --job-name=snakemake
#SBATCH --partition=RM-shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=3-00:00:00

source /jet/home/agomez3/.bashrc
source activate /jet/home/agomez3/miniforge3/envs/snakemake

snakemake --use-conda --workflow-profile profile/
