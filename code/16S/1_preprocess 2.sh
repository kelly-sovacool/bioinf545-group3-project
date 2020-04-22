#!/bin/bash

workdir=$1
procs=$2
silva=data/references/silva.v4.align
reffasta=data/references/trainset14_032015.pds.fasta
reftax=data/references/trainset14_032015.pds.tax
mock=data/references/HMP_MOCK.fasta

set.dir(input=data/raw/, output=data/16S/processed);
make.contigs(file=data/16S/crc.files, inputDir=data/raw/, processors=${procs});
summary.seqs(fasta=current);
screen.seqs(fasta=current, group=current, summary=current, maxambig=0, maxlength=275);
summary.seqs(fasta=current);
unique.seqs(fasta=current);
summary.seqs(fasta=current, name=current);
count.seqs(name=current, group=current);
summary.seqs(fasta=current, count=current);
align.seqs(fasta=current, reference=${silva});
summary.seqs(fasta=current, count=current);
screen.seqs(fasta=current, count=current, summary=current, start=1968, end=11550, maxhomop=8);
summary.seqs(fasta=current,count=current);
filter.seqs(fasta=current, vertical=T, trump=.);
unique.seqs(fasta=current, count=current);
pre.cluster(fasta=current, count=current, diffs=2);
chimera.uchime(fasta=current, count=current, dereplicate=t);
remove.seqs(fasta=current, accnos=current);
summary.seqs(fasta=current,count=current);
classify.seqs(fasta=current, count=current, reference=${reffasta}, taxonomy=${reftax}, cutoff=80);
remove.lineage(fasta=current, count=current, taxonomy=current, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota);
get.groups(fasta=current, count=current, groups=mock1-mock2);
seq.error(fasta=current, count=current, reference=${mock}, aligned=F)