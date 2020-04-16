#!/bin/bash

#SBATCH --job-name=makecontigs
#SBATCH --cpus-per-task=8
#SBATCH --qos=nolimit
#SBATCH --mem=55G
#SBATCH --time=60:00:00
#SBATCH --output=slurm-%j_%x.out
#SBATCH --account=class
#SBATCH --mail-user=britbrow@umich.edu
#SBATCH --mail-type=END

module load mothur
mothur contigs.batch 
