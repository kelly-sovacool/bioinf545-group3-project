#!/bin/bash
indir=$1
outdir=$2
for infile in ${indir}*.bam
do
  base=$(basename ${infile} _col.bam)
  samtools fasta ${infile} | cat | metaphlan2.py --input_type fasta --bowtie2out ${outdir}/${base}.bowtie2out.bz2 > ${outdir}/${base}_mtphln2.txt
done

