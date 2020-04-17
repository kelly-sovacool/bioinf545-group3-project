#!/bin/bash

#SBATCH --job-name=smk
#SBATCH --cpus-per-task=4
#SBATCH --mem=39G
#SBATCH --time=24:00:00
#SBATCH --output=log/slurm-%j_%x.out
#SBATCH --account=class
#SBATCH --mail-user=YOUR_EMAIL
#SBATCH --mail-type=BEGIN,END,FAIL

time snakemake -n --use-conda --configfile config.yml
