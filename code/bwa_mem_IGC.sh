#!/bin/bash
indir=$1
outdir=$2
indexpath=$3
for infile in $indir/*unmapped_1.fastq.gz
do
  base=$(basename ${infile} unmapped_1.fastq.gz)
  bwa mem ${indexpath} \
          ${indir}/${base}_unmapped_1.fastq.gz ${indir}/${base}_unmapped_2.fastq.gz |
  samtools view -b -f 2 - > ${outdir}/${base}_IGC.bam
done

