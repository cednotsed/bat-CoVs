rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)

meta <- fread("data/metadata/Beta_CoVs_taxid694002.050422.csv")
nrow(meta)
head(meta)
table(meta$Species)

set.seed(66)

covid <- meta %>% 
  filter(grepl("SARS-CoV-2|Severe acute respiratory syndrome coronavirus 2", GenBank_Title)) %>%
  sample_n(100)

non_covid_mers <- meta %>% 
  filter(!grepl("SARS-CoV-2|Severe acute respiratory syndrome coronavirus 2|MERS|Middle East", GenBank_Title))

# Subset MERS
mers <- meta %>% filter(grepl("MERS|Middle East", GenBank_Title))
camel <- mers %>% filter(grepl("Camelus", Host)) %>% sample_n(50)
human_mers <- mers %>% filter(grepl("Homo", Host)) %>% sample_n(50)
other_mers <- mers %>% filter(!grepl("Homo|Camelus", Host))

beta <- bind_rows(covid, non_covid_mers, camel, human_mers, other_mers)

fwrite(beta, "data/metadata/Beta_CoVs_taxid694002.050422.filtered.csv")
fwrite(beta %>% select(Accession), "data/metadata/Beta_CoVs_taxid694002.050422.filtered.accessions_only.csv")
