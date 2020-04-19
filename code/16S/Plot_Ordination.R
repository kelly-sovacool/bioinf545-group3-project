#Plot Ordination
OTU_ord <- ordinate(OTU, "NMDS", "bray")
p1 = plot_ordination(OTU, OTU_ord, type="taxa", color="Phylum", title="taxa")
print(p1)

GP.ord <- ordinate(GP1, "NMDS", "bray")
p1 = plot_ordination(GP1, GP.ord, type="taxa", color="Phylum", title="taxa")
print(p1)
