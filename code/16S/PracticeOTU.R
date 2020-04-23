library("phyloseq")
library("ggplot2")

# load data
data("GlobalPatterns")

# set ggplot theme
theme_set(theme_bw())
pal <- "Set1"
scale_colour_discrete <- function(palname = pal, ...) {
  scale_colour_brewer(palette = palname, ...)
}
scale_fill_discrete <- function(palname = pal, ...) {
  scale_fill_brewer(palette = palname, ...)
}

# prepare data, prune OTUs not presenset in any samples.
# this is an example of how
gp <- prune_taxa(taxa_sums(GlobalPatterns) > 0, GlobalPatterns)

# plots, hre is a plot of richness
plot_richness(gp)
plot_richness(gp, measures = c("Choao1", "Shannon"))

sample_data(gp)$human <- get_variable(gp, "SampleType") %in% c("Feces", "Mock", "Skin", "Tongue")

# specify sample variable on group/organize samples in the horizontal x axis
plot_richness(gp, x = "human", color = "SampleType", measures = c("Choao1", "Shannon"))
