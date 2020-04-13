#!/bin/bash

#SBATCH --job-name=bwa_GRCh38
#SBATCH --cpus-per-task=1
#SBATCH --mem=39G
#SBATCH --time=24:00:00
#SBATCH --account=class
#SBATCH --mail-user=csrkang@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --qos=nolimit

code/bwa_mem_GRCh38.sh \
  data/metagenome/trimm_results \
  data/metagenome/bwa_GRCh38_results \
  data/metagenome/bwa_DB/GRCh38/GRCh38
