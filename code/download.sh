#!/bin/bash
list_file=$1
outdir=$2
for SRA in $(cat $list_file)
do
    prefetch $SRA
    fasterq-dump --split-files $SRA -O $outdir --gzip 

done
