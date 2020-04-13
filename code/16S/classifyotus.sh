#Removes mock communities and classifies OTUs
#Performed interactively in mothur

mothur

get.groups(count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, groups=mock1)
 
remove.groups(count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.rdp1.wang.pick.taxonomy, groups=mock1)

get.groups(count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.pick.count_table, fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.fasta, groups=mock2)

remove.groups(count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.pick.count_table, fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.fasta, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.rdp1.wang.pick.pick.taxonomy, groups=mock2)

dist.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.pick.fasta, cutoff=0.03)

cluster(column=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.pick.dist, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.pick.pick.count_table, cutoff=0.03)
 
make.shared(list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.pick.opti_mcc.list, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.pick.pick.count_table, label=0.03)
 
classify.otu(list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.pick.opti_mcc.list, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.pick.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.rdp1.wang.pick.pick.pick.taxonomy, label=0.03)

