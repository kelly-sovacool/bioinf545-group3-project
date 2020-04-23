library(here)
library(tidyverse)
source(here('code', 'machine_learning.R'))

ovu_abun <- read_tsv(here("data", "virome", "ovu_abundance.tsv")) %>%
    filter(DiseaseClass %in% c("Cancer", "Healthy")) %>%
    select(-subjectid)

predict_rf(ovu_abun,
           here("data", "model", "rf_model_virus.tsv"),
           here("data", "model", "conf_mat_virus.rds"),
           here("data", "model", "feat_imp_virus.rds")
)
