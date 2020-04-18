#Ran mothur interactively
mothur

#performed with default normalization
dist.shared(shared=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.pick.opti_mcc.shared, calc=braycurtis-jclass)
 
pcoa(phylip=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.pick.opti_mcc.braycurtis.0.03.lt.dist)
 
