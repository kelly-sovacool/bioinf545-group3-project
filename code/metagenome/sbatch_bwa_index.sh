#!/bin/bash

#SBATCH --job-name=bwa_index_IGC
#SBATCH --cpus-per-task=16
#SBATCH --mem=39G
#SBATCH --time=24:00:00
#SBATCH --account=class
#SBATCH --mail-user=csrkang@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --qos=nolimit

bwa index -p data/metagenome/bwa_DB/IGC/IGC -a bwtsw data/metagenome/bwa_DB/IGC/IGC.fa 
