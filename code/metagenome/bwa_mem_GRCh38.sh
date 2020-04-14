#!/bin/bash
indir=$1
outdir=$2
indexpath=$3
for infile in $indir/*1_paired.fastq.gz
do
  base=$(basename ${infile} _1_paired.fastq.gz)
  bwa mem -t 16 ${indexpath} \
          ${indir}/${base}_1_paired.fastq.gz ${indir}/${base}_2_paired.fastq.gz |
  samtools view -b - > ${outdir}/${base}_GRCh38.bam
  samtools flagstat ${outdir}/${base}_GRCh38.bam > ${outdir}/${base}_flagstat.txt
done
 


