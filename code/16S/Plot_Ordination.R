# use vegan
library("phyloseq")
library("ggplot2")
library("plyr")
library("here")

theme_set(theme_bw())

otumat <- read.delim(file = here("data", "16S", "mothur_output", "stability.opti_mcc.shared"), header = TRUE) %>%
  select(-label, -numOtus) %>%
  pivot_longer(-Group, names_to = "otu", values_to = "abundance")

# set up taxonomy table

taxonomy <- read_tsv(file = here("data", "16S", "mothur_output", "Hannigan_mcc.0.03.cons.taxonomy")) %>%
  rename_all(tolower) %>%
  mutate(taxonomy = str_replace_all(string = taxonomy, pattern = "\\(\\d*\\)", replacement = "")) %>%
  mutate(taxonomy = str_replace_all(string = taxonomy, pattern = ";$", replacement = "")) %>%
  separate(taxonomy, into = c("kingdom", "phylum", "class", "order", "family", "genus"), sep = ";") %>%
  select(-size) %>%
  pivot_longer(-otu, names_to = "taxon_levels", values_to = "taxa")

# merge taxonomy and otu together
abun <- inner_join(otumat, taxonomy) %>%
  group_by(Group, taxa) %>%
  mutate(abunsum = sum(abundance))



# prune_taxa(tax_table(otumat[, "abundance"]))

# Plot Ordination
OTU_ord <- ordinate(otumat, "NMDS", "bray")
p1 <- plot_ordination(otumat, OTU_ord, type = "taxa", color = "Phylum", title = "taxa")
print(p1)

GP.ord <- ordinate(GP1, "NMDS", "bray")
p1 <- plot_ordination(GP1, GP.ord, type = "taxa", color = "Phylum", title = "taxa")
print(p1)
