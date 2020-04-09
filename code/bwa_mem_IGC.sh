#!/bin/bash
indir=$1
outdir=$2
indexpath=$3
for infile in $indir/SRR56650*1_unmapped.fastq.gz
do
  base=$(basename ${infile} _1_unmapped.fastq.gz)
  bwa mem ${indexpath} \
          ${indir}/${base}_1_unmapped.fastq.gz ${indir}/${base}_2_unmapped.fastq.gz |
  samtools view -h -f 2 - > ${outdir}/${base}_IGC.sam
done

