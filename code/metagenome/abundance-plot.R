# EdgeR code for the creation of KEGG abundance plots.
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("edgeR")
library("edgeR")
library("dplyr")
library("purrr")
library("tidyr")


# Import metadata to order samples by name.
SraRun <- read.table(here::here("data", "SraRunTable.txt"),
  header = T,
  sep = ",", stringsAsFactors = F
)

SraRun <- SraRun %>%
  filter(samp_mat_process == "WholeMetagenome") %>% # Filter out virome samples
  arrange(subjectid) %>% # rearrange by patient type
  arrange(factor(DiseaseClass,
    levels = c("Negative", "Healthy", "Adenoma", "Cancer")
  ))

colOrder <- as.character(SraRun$Run) # used to order the dataframe of sample keggCounts.

# Import keggCount data of all samples and return dataframe with all sample info.
fileNames <- Sys.glob(here::here(
  "data", "metagenome", "gene_abundance_results",
  "*keggCount.txt"
))

keggCounts <- lapply(fileNames, function(i) {
  read.csv(i, header = F, stringsAsFactors = F)
})

keggCounts <- keggCounts %>% reduce(full_join, by = "V1") # full join on kegg ID.

colnames(keggCounts) <-
  c("KeggNo", gsub("_keggCount.txt", "", fileNames)) # rename columns to sample name.

keggCounts <- keggCounts %>%
  select(paste(colOrder)) # reorder by sample type.

keggCounts[is.na(keggCounts)] <- 0 # replace NA with 0

write.table(keggCounts, file = here::here("data", "metagenome", "all_kegg_counts.csv"), sep = ",", row.names = TRUE)

rm(fileNames)
rm(colOrder)

# Create DGEList object with groups for the treatments C and T
group <- c(
  rep("N", sum(SraRun$DiseaseClass == "Negative")),
  rep("H", sum(SraRun$DiseaseClass == "Healthy")),
  rep("A", sum(SraRun$DiseaseClass == "Adenoma")),
  rep("C", sum(SraRun$DiseaseClass == "Cancer"))
)
cds <- DGEList(keggCounts, group = group)

# Filter out genes with low counts, keeping those rows where the count
# per million (cpm) is at least 1 in at least 6 samples:
keep <- rowSums(cpm(cds) > 1) >= 6
cds <- cds[keep, ]

# Normalize
cds <- calcNormFactors(cds)
cds <- estimateCommonDisp(cds)
cds <- estimateTagwiseDisp(cds, prior.df = 10)

# Draw the MDS plot
plotMDS(cds, main = "MDS Plot for Count Data", labels = colnames(cds$counts))

#### Did not update following code >>> Need to do pairwise?

# Find Differentially Expressed genes
DEgenes <- exactTest(cds, pair = c("C", "T"))
summary(decideTestsDGE(DEgenes, p.value = 0.05))
DEgene.table <- topTags(de.tgw, n = nrow(DEgenes$table))$table
write.table(DEgene.table,
  file = here::here("data", "metagenome", "DEgenes.csv"),
  sep = ",", row.names = TRUE
)
