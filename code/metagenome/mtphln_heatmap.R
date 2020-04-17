install.packages("ggplot2")
install.packages("dplyr")
install.packages("tidyr")

library(ggplot2)
library(dplyr)
library(tidyr)

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

setwd("/data/metagenome/metaphlan2_results")

# Metadata with sample information
SraRun <- read.table("../../SraRunTable.txt", header = T, sep = ",", stringsAsFactors = F)
SraRun <- SraRun %>%
  filter(samp_mat_process == "WholeMetagenome") %>% # Filter out virome samples
  arrange(subjectid) %>% # rearrange by patient type
  arrange(factor(DiseaseClass, levels = c("Negative", "Healthy", "Adenoma", "Cancer")))

# Abundance (%) data (could be any level, species/genus/etc., switch out 'merged_species.txt')
taxAbundance <- read.table("merged_species.txt", header = T, sep = "\t", stringsAsFactors = F)
rownames(taxAbundance) <- taxAbundance$clade_name # Have organism name as row name
taxAbundance <- taxAbundance %>%
  select(-NCBI_tax_id, -clade_name) # Collect abundance values
colnames(taxAbundance) <- gsub("_mtphln2", "", colnames(taxAbundance)) # Truncate sample names

# Reorder columns of abundance data based on order in meta data
colOrder <- as.data.frame(SraRun$Run)
msCols <- as.data.frame(colnames(taxAbundance))
colnames(colOrder) <- "Run"
colnames(msCols) <- "Run"
colOrder <- inner_join(colOrder, msCols, by = "Run")
rm(msCols)
taxAbundance <- taxAbundance %>% select(paste(colOrder$Run)) # Reordering
taxAbundance <- as.matrix(taxAbundance) # Need matrix for heatmap.2()

# Assign sample colors to meta data
sampleColors <- mapSampleToColor(SraRun)
SraRun$sampleColors <- sampleColors
SraRunCol <- SraRun
SraRunCol <- inner_join(colOrder, SraRunCol, by = "Run")
sampleColors <- SraRunCol$sampleColors
rm(colOrder)
rm(SraRunCol)

# Construct heatmap
heatmap.2(taxAbundance,
  ColSideColors = sampleColors,
  trace = "none",
  scale = "row",
  col = brewer.pal(11, "RdBu"),
  dendrogram = "col"
)
