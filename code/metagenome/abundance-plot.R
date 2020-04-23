# EdgeR code for the creation of KEGG abundance plots.
library("edgeR")
library("dplyr")
library("purrr")
library("tidyr")
library("ggplot2")
# library("vegan")

# Import metadata to order samples by name.
SraRun <- read.table(here::here("data", "SraRunTable.txt"),
                     header = T,
                     sep = ",", stringsAsFactors = F
) %>%
  filter(samp_mat_process == "WholeMetagenome") %>% # Filter out virome samples
  filter(Run != "SRR5665024" & # remove negative control and samples with extremely low read count
           Run != "SRR5665060" &
           Run != "SRR5665081" &
           Run != "SRR5665121" &
           Run != "SRR5665138" &
           Run != "SRR5665117" &
           Run != "SRR5665160") %>%
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

keggCounts <- keggCounts %>% reduce(full_join, by = "V1") # full join on KEGG No.

colnames(keggCounts) <-
  c("KeggNo", gsub(".*(SRR.*)_keggCount.txt", "\\1", fileNames)) # rename columns to sample name.

keggCounts <- keggCounts %>% filter(KeggNo != "unknown")

rownames(keggCounts) <- keggCounts$KeggNo

keggCounts <- keggCounts %>%
  select(paste(colOrder)) # reorder by sample type.

colnames(keggCounts) <- c(SraRun$subjectid) # add the subject ids
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

groupColors <-
  c(
    rep("black", sum(SraRun$DiseaseClass == "Negative")),
    rep("#66C2A5", sum(SraRun$DiseaseClass == "Healthy")),
    rep("#FC8D62", sum(SraRun$DiseaseClass == "Adenoma")),
    rep("#8DA0CB", sum(SraRun$DiseaseClass == "Cancer"))
  )

cds <- keggCounts %>%
  DGEList(group = group) # analyze just the count numbers
# Filter out genes with low counts, keeping those rows where the count
# per million (cpm) is at least 1 in at least 6 samples:
keep <- rowSums(cpm(cds) > 1) >= 6
cds <- cds[keep, ]

# Normalize
cds <- calcNormFactors(cds)
cds <- estimateCommonDisp(cds)
cds <- estimateTagwiseDisp(cds, prior.df = 10)

# Draw the MDS plot
pdf(file = here::here("figures", "gene_abundance_MDS.pdf"), width = 5, height = 5)
plotMDS(cds, pch = 19, col = groupColors)
legend("bottomright",fill=c("#66C2A5","#FC8D62", "#8DA0CB"),legend=unique(factor(SraRun$DiseaseClass)))
dev.off()

# Find differentially abundant genes in healthy vs cancer/adenoma
DEgenes_HA <- exactTest(cds, pair = c("H", "A"))
summary(decideTestsDGE(DEgenes_HA, p.value = 0.05))
DEgene_table_HA <- topTags(DEgenes_HA, n = nrow(DEgenes_HA$table), sort.by = "logFC")$table
DEgene_table_HA$KeggNo <- rownames(DEgene_table_HA)
DEgene_table_HA <- DEgene_table_HA %>% filter(PValue <= 0.05)

DEgenes_HC <- exactTest(cds, pair = c("H", "C"))
summary(decideTestsDGE(DEgenes_HC, p.value = 0.05))
DEgene_table_HC <- topTags(DEgenes_HC, n = nrow(DEgenes_HC$table), sort.by = "logFC")$table
DEgene_table_HC$KeggNo <- rownames(DEgene_table_HC)
DEgene_table_HC <- DEgene_table_HC %>% filter(PValue <= 0.05)

DEgenes <- inner_join(DEgene_table_HA, DEgene_table_HC, by = "KeggNo", suffix = c("_HA", "_HC"))

DEgenes_pos <- DEgenes %>%
  filter(logFC_HA > 1 & logFC_HC > 1) %>%
  select(KeggNo) %>%
  head(500)
DEgenes_neg <- DEgenes %>%
  filter(logFC_HA < -1 & logFC_HC < -1) %>%
  select(KeggNo) %>%
  head(500)

write.table(DEgenes_pos,
  file = here::here("data", "metagenome", "gene_abundance_results", "DEgenes_pos.txt"),
  sep = "\t", row.names = T, quote = F, col.names = F
)

write.table(DEgenes_neg,
  file = here::here("data", "metagenome", "gene_abundance_results", "DEgenes_neg.txt"),
  sep = "\t", row.names = T, quote = F, col.names = F
)
