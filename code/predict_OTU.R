library(here)
library(tidyverse)
source(here('code', 'machine_learning.R'))

otu_abun <-
    read_tsv(here("data", '16S', 'mothur_output',
                  'stability.opti_mcc.shared')
             ) %>%
    mutate(
        subject_id = str_to_title(Group) %>%
            str_replace("([A-Z])[a-z]*([0-9]{1,2})_[_0-9]*", "\\1\\2"),
        DiseaseClass = case_when(
            startsWith(Group, 'A') ~ 'Adenoma',
            startsWith(Group, 'C') ~ 'Cancer',
            startsWith(Group, 'H') ~ 'Healthy'
        )
    ) %>%
    filter(DiseaseClass %in% c("Cancer", "Healthy")) %>%
    select(-label, -Group, -numOtus, -subject_id) %>%
    select(DiseaseClass, everything())

predict_rf(otu_abun,
           here("data", "model", "rf_model_bacteria.tsv"),
           here("data", "model", "conf_mat_bacteria.rds"),
           here("data", "model", "feat_imp_bacteria.rds")
           )
