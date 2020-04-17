#!/bin/bash
indir=$1
for input in $indir/*_GRCh38_mapped.bam
do
  base=$(basename ${input} _GRCh38_mapped.bam)
  samtools view -bh -f 5 ${input} >  ${base}_GRCh38_unmapped.bam
done

