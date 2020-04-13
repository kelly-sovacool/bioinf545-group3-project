#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=48:00:00 
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --cpus-per-task=2
#SBATCH --mem=60G
#SBATCH --qos=nolimit
#SBATCH --mail-user=csrkang@umich.edu
#SBATCH --job-name=bwa_IGC
#########################################################################################

function Fq2Bam {

  BASENAME=$1

  bwa mem -t 16 data/metagenome/bwa_DB/IGC/IGC ${BASENAME}_1_unmapped.fastq.gz ${BASENAME}_2_unmapped.fastq.gz | samtools sort -o ${BASENAME}_sorted.bam

}
#########################################################################################

export -f Fq2Bam

ls data/metagenome/bwa_GRC_sbatch/*_1_unmapped.fastq.gz | awk -F "_1_unmapped.fastq.gz" '{print $1}' | parallel -j 2 "Fq2Bam {}"
