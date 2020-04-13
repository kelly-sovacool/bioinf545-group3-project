#Run interactively in mothur
mothur

#may need to add normalization depending on results?
#default is to normalize by sample with smallest sequences
summary.single(shared=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.pick.opti_mcc.shared)
 
rarefaction.single(shared=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.pick.opti_mcc.shared)
