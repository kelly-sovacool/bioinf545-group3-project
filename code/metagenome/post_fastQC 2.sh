#!/bin/bash
indir=$1
outdir=$2
for infile in $indir/*_paired.fastq.gz
do
  base=$(basename ${infile} _paired.fastq.gz)
  fastqc -o ${outdir} ${infile}
done
