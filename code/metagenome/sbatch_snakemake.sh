#!/bin/bash

#SBATCH --job-name=snakemake
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --account=class
#SBATCH --mail-user=csrkang@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --qos=nolimit

snakemake --cores 16 --use-conda
