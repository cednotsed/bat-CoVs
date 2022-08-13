rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)

meta <- fread("data/metadata/GISAID_050422.tsv")

parsed_meta <- meta %>%
  filter(Host != "Human") %>%
  rename_all(~ tolower(gsub(" ", "_", .))) %>%
  separate(location, into = c("loc1", "loc2", "loc3"), sep = " / ") %>%
  mutate(location = paste0(loc1, " / ", loc2),
         collection_date = as.Date(collection_date)) %>%
  select(-type, -gender, -`gc-content`, 
         -`n-content`, -`is_reference?`, -pangolin_version,
         -patient_age, -additional_location_information)

parsed_filt <- parsed_meta  %>%
  filter(sequence_length >= 22000, 
         grepl("Rhino|Manis", host))

fwrite(parsed_filt, "data/metadata/GISAID_050422.filtered.tsv")
fwrite(parsed_filt %>% select(accession_id), "data/metadata/GISAID_050422.filtered.accessions_only.tsv")
