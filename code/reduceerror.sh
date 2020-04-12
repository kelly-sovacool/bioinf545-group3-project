#Commands were run  after contigs are made
#Commands were run interactively in mothur, so imperfect script
#Follows guidelines from Miseq_SOP

mothur

summary.seqs(fasta=stability.trim.contigs.fasta)

screen.seqs(fasta=stability.trim.contigs.fasta, group=stability.contigs.groups, maxambig=0, maxlength=275)

#Do another summary check of the screened data
summary.seqs(fasta=stability.trim.contigs.good.fasta)

#Eliminating repeating sequences but does track number of repeats
unique.seqs(fasta=stability.trim.contigs.good.fasta)

#Shows number of a particular sequence in each sample
count.seqs(name=stability.trim.contigs.good.names, group=stability.contigs.good.groups)

#silva file is specific to V4 region of 16S
#next run align.sh file, followed by improveseqs.sh
