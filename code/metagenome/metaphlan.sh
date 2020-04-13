#!/bin/bash
indir=$1
outdir=$2
for infile in ${indir}*.bam
do
  base=$(basename ${infile} _GRCh38.bam)
  samtools fasta ${infile} | cat |
           metaphlan2.py --input_type fasta --nproc 4 --bowtie2out ${outdir}/${base}.bowtie2out.bz2 \
           > ${outdir}/${base}_mtphln2.txt
done

merge_metaphlan_tables.py ${outdir}/*.txt > ${outdir}/merged.txt
grep -E "(s__)|(^clade_name)" merged.txt | grep -v "t__" | sed 's/^.*s__//g' > ${outdir}/merged_species.txt
grep -E "(g__)|(^clade_name)" merged.txt | grep -v "t__" | sed 's/^.*g__//g' | sed '/|s_/d' > ${outdir}/merged_genus.txt
grep -E "(f__)|(^clade_name)" merged.txt | grep -v "t__" | sed 's/^.*f__//g' | sed '/|g_/d' > ${outdir}/merged_family.txt
grep -E "(p__)|(^clade_name)" merged.txt | grep -v "t__" | sed 's/^.*p__//g' | sed '/|c_/d' > ${outdir}/merged_phylum.txt

