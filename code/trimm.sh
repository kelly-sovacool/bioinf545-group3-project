#!/bin/bash
indir=$1
outdir=$2
adapterPath=$3
for infile in $indir/*_1.fastq.gz
do
  base=$(basename ${infile} _1.fastq.gz)
  trimmomatic PE -phred33 -threads 4 \
               ${infile} ${indir}/${base}_2.fastq.gz \
               ${outdir}/${base}_paired_1.fastq.gz ${outdir}/${base}_unpaired_1.fastq.gz \
               ${outdir}/${base}_paired_2.fastq.gz ${outdir}/${base}_unpaired_2.fastq.gz \
               ILLUMINACLIP:${adapterPath}:2:40:15 \
               LEADING:3 TRAILING:3 MINLEN:24 SLIDINGWINDOW:4:15 \
  &> ${outdir}/${base}_logX.txt
  grep -vwE '(Skipping|Using)' ${outdir}/${base}_logX.txt > ${outdir}/${base}_log.txt
  rm ${outdir}/${base}_logX.txt
done
