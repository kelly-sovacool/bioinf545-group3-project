
library(tidyverse)
library(here)

mapSampleToColor <- function(annotations) {
  # Assign color in heatmap.2() based on disease class
  # INPUT: metadata (SraRunTable.txt)
  # Healthy = green
  # Adenoma = blue
  # Cancer = red
  colorsVector <- ifelse(annotations["DiseaseClass"] == "Healthy",
    "green", ifelse(annotations["DiseaseClass"] == "Adenoma",
      "blue", "red"
    )
  )
  return(colorsVector)
}



# Metadata with sample information
metadata <- read.table(here("data", "SraRunTable.txt"), header = T, sep = ",", stringsAsFactors = F) %>%
  filter(samp_mat_process == "WholeMetagenome") %>% # Filter out virome samples
  arrange(factor(DiseaseClass, levels = c("Negative", "Healthy", "Adenoma", "Cancer"))) %>%
  select(subjectid, Run, DiseaseClass)


# Heat map construction
plot_heatmap <- function(taxonlevel){
  taxheatmap <- read.table(here("data", "metagenome","metaphlan2_results", paste0("merged_",taxonlevel,".txt")), header = T, sep = "\t", stringsAsFactors = F) %>%
    pivot_longer(cols = starts_with("SRR"), names_to = "Run", values_to = "abun") %>%
    mutate(Run = gsub("_mtphln2", "", Run))  %>%
    inner_join(metadata) %>%
  ggplot(aes(x = subjectid, y = clade_name, fill = abun)) +
  geom_tile() +
  scale_fill_distiller(type = "div", na.value = "white", direction = -1, palette = "RdYlBu") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
  }
plotlist<- lapply(c("phylum", "family", "genus", "species"), plot_heatmap)

ggsave(filename = here("figures","metag_taxa_heatmap_phylum.png"),
       width = 10, height = 7,  plot = plotlist[[1]])

ggsave(filename = here("figures","metag_taxa_heatmap_family.png"),
       width = 10, height = 7,  plot = plotlist[[2]])

ggsave(filename = here("figures","metag_taxa_heatmap_genus.png"),
       width = 10, height = 7,  plot = plotlist[[3]])


