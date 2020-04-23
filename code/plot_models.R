library(here)
library(plotROC)
library(tidyverse)
otu <- read_tsv(here("data", "model", "rf_model_bacteria.tsv")) %>%
    mutate(seqtype = "16S Bacteria")
ovu <- read_tsv(here("data", "model", "rf_model_virus.tsv")) %>%
    mutate(seqtype = "Virus")
models <- rbind(otu, ovu)
plot_roc <- ggplot(models, aes(d = obs,
                   m = Healthy,
                   color = seqtype)) +
    geom_roc(n.cuts=0) +
    style_roc(theme = theme_gray()) +
    coord_equal() +
    #theme_classic() +
    scale_color_manual(values = c("#1F78B4", "#B2DF8A")) +
    theme(
        legend.position = c(0.75, 0.25),
        legend.title=element_blank()
    ) +
    ggtitle("Random Forest Model Performance")
ggsave(here("figures", "auroc.png"), plot_roc, height = 3.75, width = 3.75)
