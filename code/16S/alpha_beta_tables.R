library(here)
library(tidyverse)
library(vegan)
# input is the otu table
# map is the metadata table

input <- read_tsv(file = here("data", "16S", "mothur_output", "stability.opti_mcc.shared"))
input$Group <- gsub("Adenoma", "A", input$Group) %>%
  gsub("Ademona", "A", .) %>%
  gsub("ADenoma", "A", .) %>%
  gsub("Cancer", "C", .) %>%
  gsub("Healthy", "H", .) %>%
  gsub("(.*)_\\w+", "\\1", .)


# input<- pivot_longer(input, col = Otu0001:Otu3904)
input <- input %>% select(-label, -numOtus)

map <- read.table(here::here("data", "SraRunTable.txt"),
  header = T,
  sep = ",", stringsAsFactors = F
) %>%
  filter(samp_mat_process == "WholeMetagenome") %>% # Filter out virome samples
  arrange(subjectid) %>% # rearrange by patient type
  arrange(factor(DiseaseClass,
    levels = c("Negative", "Healthy", "Adenoma", "Cancer")
  ))

## alpha diversity function
AlphaDiversity <- function(intable, mapping) {
  ttable <- as.matrix((input)[, -1])
  rownames(ttable) <- input$Group

  class(ttable) <- "numeric"
  divtable <- diversity(ttable, index = c("shannon"))
  divtable <- data.frame(divtable)
  colnames(divtable) <- "shannon"
  divtable$subjectid <- row.names(divtable)

  # Merge in the mapping file
  mapsub <- map[, c("subjectid", "DiseaseClass")]
  merged <- merge(divtable, mapsub, by.x = "subjectid", by.y = "subjectid")
  plot <- ggplot(merged, aes(x = DiseaseClass, y = shannon, colour = DiseaseClass)) +
    theme_classic() +
    theme(
      axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
      axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid")
    ) +
    geom_jitter() +
    scale_colour_brewer(palette = "Set1") +
    geom_boxplot(outlier.colour = NA, alpha = 0) +
    ggtitle("16S Shannon Diversity")
  return(plot)
}

BetaDiversity <- function(intable, mapping) {
  # Fix the names in the headers
  # names(intable) <- gsub("_R2", "", names(intable))
  # Get rid of the missing values
  # tabledrop <- intable[, -grep("X", colnames(intable))]
  # ttable <- as.matrix(t(tabledrop)[-1,])
  # class(ttable) <- "numeric"

  ttable <- as.matrix((input)[, -1])
  rownames(ttable) <- input$Group
  class(ttable) <- "numeric"

  # Create dist matrix
  gothedistance <- vegdist(ttable, method = "bray")
  ordnmds <- metaMDS(gothedistance, k = 2)
  ordnmdsfit <- data.frame(MDS1 = ordnmds$points[, 1], MDS2 = ordnmds$points[, 2])
  ordnmdsfit$ID <- rownames(ordnmdsfit)
  # Merge in the mapping file
  mapsub <- map[, c("subjectid", "DiseaseClass")]
  merged <- merge(ordnmdsfit, mapsub, by.x = "ID", by.y = "subjectid")
  plot <- ggplot(merged, aes(x = MDS1, y = MDS2, colour = DiseaseClass)) +
    theme_classic() +
    theme(
      axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
      axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid")
    ) +
    geom_point(size = 1.5) +
    scale_colour_brewer(palette = "Set1") +
    ggtitle("16S Bray-Curtis")
  return(plot)
}


############
# Run Alpha Diversity Function #
############
AlphaDiversity(input, map)
BetaDiversity(input, map)
## where they save it - change to our project
pdf(file = here("figures", paste0("alphadiversity.pdf"), width = 4, height = 4))
AlphaDiversity(input, map)

pdf(file = here("figures", paste0("betadiversity.pdf"), width = 5, height = 4))
BetaDiversity(input, map)
dev.off()
