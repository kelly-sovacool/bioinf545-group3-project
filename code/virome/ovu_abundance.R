library(here)
library(tidyverse)
bins <- read_csv(here("data","virome","concoct","clustering_merged.csv"), col_types = 'cc') %>%
    rename(contig = contig_id)
coverage <- read_tsv(here("data","virome","contigs","coverage_table.tsv")) %>%
    rename_at(vars(contains("SRR")), funs(str_replace(., ".*(SRR.*)_mapped.sorted", "\\1"))) %>%
    mutate(contig = str_replace(contig, "([0-9]*)\\..*", "\\1")) %>%
    pivot_longer(-contig, names_to = "Run", values_to = "abun")
metadata <- read_csv(here('data', 'SraRunTable.txt')) %>%
    filter(samp_mat_process == "Virome") %>%
    select(subjectid, Run, DiseaseClass)
ovu_abun <- coverage %>%
    inner_join(metadata, by = "Run") %>%
    inner_join(bins, by = "contig") %>%
    select(subjectid, DiseaseClass, contig, abun, cluster_id) %>%
    mutate(cluster_id = str_replace(cluster_id, "(.*)", "OVU\\1")) %>%
    group_by(subjectid, cluster_id, DiseaseClass) %>%
    summarize(abunsum = sum(abun)) %>%
    pivot_wider(names_from = cluster_id, values_from = abunsum)
write_tsv(ovu_abun, here("data", "virome", "ovu_abundance.tsv"))
