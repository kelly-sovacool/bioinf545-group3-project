#!/bin/bash
indir=$1
for infile in $indir/*1_unpaired.fastq.gz
do
  base=$(basename ${infile} _1_unpaired.fastq.gz)
  mv ${infile} ${indir}/${base}_unpaired_1.fastq.gz
done
for infile in $indir/*2_unpaired.fastq.gz
do
  base=$(basename ${infile} _2_unpaired.fastq.gz)
  mv ${infile} ${indir}/${base}_unpaired_2.fastq.gz
done
