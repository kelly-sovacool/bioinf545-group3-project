
library(tidyverse)
library(here)
library(cowplot)
library(grid)

mapSampleToColor <- function(annotations) {
  # Assign color in heatmap.2() based on disease class
  # INPUT: metadata (SraRunTable.txt)
  # Healthy = green
  # Adenoma = blue
  # Cancer = red
  colorsVector <- ifelse(annotations["DiseaseClass"] == "Healthy",
                         "green", ifelse(annotations["DiseaseClass"] == "Adenoma",
                                         "blue", "red"
                                         ))
  return(colorsVector)
}

# Metadata with sample information
metadata <- read.table(here("data", "SraRunTable.txt"), header = T, sep = ",", stringsAsFactors = F) %>%
  filter(samp_mat_process == "WholeMetagenome") %>% # Filter out virome samples
  arrange(factor(DiseaseClass, levels = c("Negative", "Healthy", "Adenoma", "Cancer"))) %>%
  select(subjectid, Run, DiseaseClass) %>%
  filter(Run != "SRR5665024" & # remove negative control and samples with extremely low read count
           Run != "SRR5665060" &
           Run != "SRR5665081" &
           Run != "SRR5665121" &
           Run != "SRR5665138" &
           Run != "SRR5665117" &
           Run != "SRR5665160")
metadata$subjectid <- factor(metadata$subjectid, levels = metadata$subjectid)

# Heat map construction
plot_heatmap <- function(taxonlevel){
  taxheatmap <- read.table(here("data", "metagenome","metaphlan2_results",
                                paste0("merged_",taxonlevel,".txt")), header = T, sep = "\t", stringsAsFactors = F) %>%
    pivot_longer(cols = starts_with("SRR"), names_to = "Run", values_to = "abun") %>%
    mutate(Run = gsub("_mtphln2", "", Run))  %>%
    inner_join(metadata) %>%
    ggplot(aes(x = subjectid, y = clade_name, fill = abun)) +
    geom_tile() +
    scale_fill_distiller(type = "div", na.value = "white", direction = -1, palette = "RdYlBu") +
    theme_classic() +
    scale_x_discrete(limits = levels(metadata$subjectid)) +
    labs(x = "Subject ID", y = taxonlevel, fill = "Abundance (%)") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 2),
          axis.line = element_blank(),
          axis.ticks = element_blank())
}
plotlist<- lapply(c("Phylum", "Family", "Genus", "Species"), plot_heatmap)

plot_phy_fam <- plot_grid(plot_grid(NULL, plotlist[[1]], labels = c("(a)",""), ncol = 2, rel_widths = c(1,8)),
         plot_grid(plotlist[[2]], NULL, labels = c("(b)",""), ncol = 2, rel_widths = c(1,0)),
         ncol = 1, rel_heights = c(2,7)
         )

rect <- rectGrob(
  x = 0.9,
  y = 1,
  width = unit(2, "in"),
  height = unit(2, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fill = "white", col = "white")
)

plot_phy_fam_paneled <- ggdraw(plot_phy_fam) + draw_grob(rect)

ggsave(filename = here("figures", "metag_taxa_heatmap_phy_fam.png"),
       width = 12, height = 7, plot = plot_phy_fam_paneled)

#ggsave(filename = here("figures","metag_taxa_heatmap_phylum.png"),
#       width = 14, height = 2,  plot = plotlist[[1]])

#ggsave(filename = here("figures","metag_taxa_heatmap_family.png"),
#       width = 14, height = 7,  plot = plotlist[[2]])

#ggsave(filename = here("figures","metag_taxa_heatmap_genus.png"),
#       width = 10, height = 7,  plot = plotlist[[3]])
