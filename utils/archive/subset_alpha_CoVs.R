rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)

meta <- fread("data/metadata/Alpha_CoVs_taxid693996.050422.csv")
nrow(meta)
head(meta)
table(meta$Species)
meta %>% filter(Accession == "MZ218060.1")
set.seed(66)

pork <- meta %>% 
  filter(grepl("Porcine epidemic", GenBank_Title)) %>%
  sample_n(100)

no_pork <- meta %>% 
  filter(!grepl("Porcine epidemic", GenBank_Title))

alpha <- bind_rows(no_pork, pork)

fwrite(alpha, "data/metadata/Alpha_CoVs_taxid693996.050422.filtered.csv")
fwrite(alpha %>% select(Accession), "data/metadata/Alpha_CoVs_taxid693996.050422.filtered.accessions_only.csv")
