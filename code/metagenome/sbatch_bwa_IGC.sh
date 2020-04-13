#!/bin/bash

#SBATCH --job-name=bwa_IGC
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G
#SBATCH --time=24:00:00
#SBATCH --account=class
#SBATCH --mail-user=csrkang@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --qos=nolimit

code/bwa_mem_IGC.sh \
    data/metagenome/bwa_GRCh38_results \
    data/metagenome/bwa_IGC_results \
    data/metagenome/bwa_DB/IGC/IGC
