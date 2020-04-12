#!/bin/bash
list_file=$1
outdir=$2
for SRA in $(cat $list_file)
do
    prefetch $SRA
    fasterq-dump --split-files $SRA -O $outdir #--gzip 
#Brittany: I commented out --gzip because fasterq doesn't recognize gzip, only works for fastq
done
