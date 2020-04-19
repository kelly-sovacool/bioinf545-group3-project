library(phyloseq)
library(ggplot2)
library(gridExtra)
library(vegan)
library(tidyverse)
library(here)

#load data
#alpha_Summary <- read.delim("~/Desktop/Projects/bioinf545-group3-project/data/16S/mothur_output/alpha_diversity/groups.summary")
#make the taxa the first column
#alpha_Summary <- alpha_Summary[,c(2:18)]

#make an OTU table
# get OTU table from github
otumat = read.delim("~/Desktop/Projects/bioinf545-group3-project/data/16S/mothur_output/stability.opti_mcc.shared", header = TRUE)

#remove first column and transpose data to match OTU format
otumat <- t(otumat[,4:90])

rownames(otumat) <- paste0("OTU", 1:nrow(otumat))
colnames(otumat) <- paste0("Sample", 1:ncol(otumat))

OTU = otu_table(otumat, taxa_are_rows = TRUE)

OTU_ord <- ordinate(OTU, "NMDS", "bray")
p1 = plot_ordination(OTU, OTU_ord, type="taxa", color="Phylum", title="taxa")
print(p1)

GP.ord <- ordinate(GP1, "NMDS", "bray")
p1 = plot_ordination(GP1, GP.ord, type="taxa", color="Phylum", title="taxa")
print(p1)


#set up fake taxonomy table
taxmat = matrix(sample(letters, 70, replace = TRUE), nrow = nrow(otumat), ncol = 7)
rownames(taxmat) <- rownames(otumat)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxmat

TAX = tax_table(taxmat)
#set up real taxonomy table

taxonomy <- read_tsv(file=here("data","16S","mothur_output","Hannigan_mcc.0.03.cons.taxonomy")) %>%
  rename_all(tolower) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="\\(\\d*\\)", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern=";$", replacement="")) %>%
  separate(taxonomy, into=c("kingdom", "phylum", "class", "order", "family", "genus"), sep=";") %>%
  select(-size)


rownames(taxonomy) <- taxonomy$otu

taxonomy <- tax_table(taxonomy %>% select(-otu))

physeq <- phyloseq(OTU, taxonomy)

plot_bar(physeq, fill = "family")








# convert the otu_table() within a phyloseq object to a vegan compatible data object
psotu2veg <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

#set ggplot theme
theme_set(theme_bw())
pal = "Set1"
scale_colour_discrete <-  function(palname=pal, ...){
  scale_colour_brewer(palette=palname, ...)
}
scale_fill_discrete <-  function(palname=pal, ...){
  scale_fill_brewer(palette=palname, ...)
}

#prepare data, prune OTUs not presenset in any samples.
#this is an example of how
as <- prune_taxa(taxa_sums(alpha_Summary) > 0, alpha_Summary)

#plots, hre is a plot of richness
plot_richness(alpha_Summary)
plot_richness(gp, measures=c("Choao1", "Shannon"))

sample_data(gp)$human <- get_variable(gp, "SampleType") %in% c("Feces", "Mock", "Skin", "Tongue")

#specify sample variable on group/organize samples in the horizontal x axis
plot_richness(gp, x="human", color="SampleType", measures=c("Choao1", "Shannon"))




