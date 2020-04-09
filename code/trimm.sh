#!/bin/bash
indir=$1
outdir=$2
adapterPath=$3
for infile in $indir/*_1.fastq.gz
do
  base=$(basename ${infile} _1.fastq.gz)
  trimmomatic PE -phred33 -threads 4 \
               ${infile} ${indir}/${base}_2.fastq.gz \
               ${outdir}/${base}_1_paired.fastq.gz ${outdir}/${base}_1_unpaired.fastq.gz \
               ${outdir}/${base}_2_paired.fastq.gz ${outdir}/${base}_2_unpaired.fastq.gz \
               ILLUMINACLIP:${adapterPath}:2:40:15 \
               LEADING:3 TRAILING:3 MINLEN:24 SLIDINGWINDOW:4:15 \
  &> ${outdir}/${base}_logX.txt
  grep -vwE '(Skipping|Using)' ${outdir}/${base}_logX.txt > ${outdir}/${base}_log.txt
  rm ${outdir}/${base}_logX.txt
done
