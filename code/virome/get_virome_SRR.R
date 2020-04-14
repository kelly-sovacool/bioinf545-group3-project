library(dplyr)
metadata <- readr::read_csv(here::here("data", "SraRunTable.txt"))
viromes <- metadata %>% filter(samp_mat_process == "Virome") %>% pull(Run)
readr::write_lines(viromes, here::here("data", "virome", "SRR_Acc_List_virome.txt"))
