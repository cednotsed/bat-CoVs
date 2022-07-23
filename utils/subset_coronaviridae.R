rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)

meta <- fread("data/metadata/coronaviridae_taxid11118_complete_exclude_provirus_040722/coronaviridae_metadata.csv") %>%
  filter(Length > 22000)

meta %>%
  group_by(Species) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  View()

meta_filt <- meta %>% 
  filter(Species != "Severe acute respiratory syndrome-related coronavirus") %>%
  filter(!grepl("Infectious bronchitis", GenBank_Title)) %>%
  filter(!grepl("Porcine epidemic", GenBank_Title)) %>%
  filter(Species != "Middle East respiratory syndrome-related coronavirus")

set.seed(66)

ibv <- meta %>% 
  filter(grepl("Infectious bronchitis", GenBank_Title)) %>%
  sample_n(50)
  
pork <- meta %>% 
  filter(grepl("Porcine epidemic", GenBank_Title)) %>%
  sample_n(50)

cov2 <- meta %>% 
  filter(Species == "Severe acute respiratory syndrome-related coronavirus",
         grepl("Severe acute respiratory syndrome coronavirus 2", GenBank_Title)) %>%
  sample_n(50)

cov1_like <- meta %>% 
  filter(Species == "Severe acute respiratory syndrome-related coronavirus",
         !grepl("Severe acute respiratory syndrome coronavirus 2", GenBank_Title))

other_MERS <- meta %>%
  filter(Species == "Middle East respiratory syndrome-related coronavirus") %>%
  filter(!grepl("Camelus|Homo", Host))

MERS_human <- meta %>%
  filter(Species == "Middle East respiratory syndrome-related coronavirus",
         grepl("Homo", Host)) %>%
  sample_n(25)

MERS_camel <- meta %>%
  filter(Species == "Middle East respiratory syndrome-related coronavirus",
         grepl("Camelus", Host)) %>%
  sample_n(25)



meta_final <- bind_rows(meta_filt, ibv, pork, 
                        cov2, cov1_like, other_MERS, 
                        MERS_human, MERS_camel)
nrow(meta_final)

meta_final %>%
  group_by(Species) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

fwrite(meta_final, str_glue("data/metadata/coronaviridae_taxid11118_complete_exclude_provirus_040722/coronaviridae_metadata.n{nrow(meta_final)}_subset.csv"))
fwrite(meta_final %>% select(Accession), str_glue("data/metadata/coronaviridae_taxid11118_complete_exclude_provirus_040722/coronaviridae_metadata.n{nrow(meta_final)}_subset.accessions_only.txt"))
