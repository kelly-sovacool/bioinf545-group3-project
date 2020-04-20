library(tidyverse)
library(here)

otumat <- read.delim(file=here("data","16S","mothur_output","stability.opti_mcc.shared"), header = TRUE) %>%
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
abun <- inner_join(otumat, taxonomy) %>%
    group_by(Group, taxa) %>%
    mutate(abunsum = sum(abundance))
#plot abun function

plot_abundance <- function(taxon_level) {
    taxa_barplot <- abun %>% filter(taxon_levels == taxon_level) %>%
    ggplot(aes(x = Group, y = abunsum, fill = taxa)) +
    geom_bar(stat = "identity", position = "stack") +
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    labs(x = "sample", y = "abundance", title = paste0(taxon_level, " Abundance"))
ggsave(filename = here("figures",paste0("taxa_barplot_", taxon_level, ".png")),
       width = 10, height = 7,  plot = taxa_barplot)
}

map(c("kingdom", "phylum", "class"), plot_abundance)
