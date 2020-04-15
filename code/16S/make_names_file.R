library(dplyr)

metadata <- readr::read_csv(here::here("data", "SraRunTable.txt"))
wgs <- metadata %>% filter(samp_mat_process == "WholeMetagenome") %>%
    select(subjectid, Run) %>%
    mutate(r1 = paste0(Run, '_1.fastq.gz'),
           r2 = paste0(Run, "_2.fastq.gz")) %>%
    select(-Run)
readr::write_tsv(wgs,
                 path = here::here("data", "16S", "crc.files"),
                 col_names = FALSE)
