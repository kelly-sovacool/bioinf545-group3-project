#!/bin/bash

#SBATCH --job-name=post_trim_fastQC
#SBATCH --cpus-per-task=4
#SBATCH --mem=39G
#SBATCH --time=24:00:00
#SBATCH --account=class
#SBATCH --mail-user=csrkang@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL

./post_fastQC.sh data/metagenome/trimm_results data/metagenome/post_trim_fastQC
