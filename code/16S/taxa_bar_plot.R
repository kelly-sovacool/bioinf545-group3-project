library(phyloseq)
library(ggplot2)
library(gridExtra)
library(vegan)
library(tidyverse)
library(here)

otumat <- read.delim("~/Desktop/Projects/bioinf545-group3-project/data/16S/mothur_output/stability.opti_mcc.shared", header = TRUE) %>%
    select(-label, -numOtus) %>%
    pivot_longer(-Group, names_to = "otu", values_to = "abundance")

#set up taxonomy table

taxonomy <- read_tsv(file=here("data","16S","mothur_output","Hannigan_mcc.0.03.cons.taxonomy")) %>%
    rename_all(tolower) %>%
    mutate(taxonomy=str_replace_all(string=taxonomy, pattern="\\(\\d*\\)", replacement="")) %>%
    mutate(taxonomy=str_replace_all(string=taxonomy, pattern=";$", replacement="")) %>%
    separate(taxonomy, into=c("kingdom", "phylum", "class", "order", "family", "genus"), sep=";") %>%
    select(-size) %>%
    pivot_longer(-otu, names_to = "taxon_levels", values_to = "taxa")

#merge taxonomy and otu together
abun <- inner_join(otumat, taxonomy)



