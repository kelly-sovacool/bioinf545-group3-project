#!/bin/bash
indir=$1
outdir=$2
indexpath=$3
for infile in $indir/*1_paired.fastq.gz
do
  base=$(basename ${infile} _1_paired.fastq.gz)
  bwa mem ${indexpath} \
          ${indir}/${base}_1_paired.fastq.gz ${indir}/${base}_2_paired.fastq.gz |
  samtools view -b - | samtools collate - ${outdir}/${base}_col
  samtools flagstat ${outdir}/${base}_col.bam > ${outdir}/${base}_flagstat.txt
  samtools view -h -f 4 ${outdir}/${base}_col.bam | \
  samtools fastq -1 ${outdir}/${base}_1_unmapped.fastq.gz -2 ${outdir}/${base}_2_unmapped.fastq.gz -
done
 


