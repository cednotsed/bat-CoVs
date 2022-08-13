setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)

# Identify monophyletic sarbecovirus clade in dendroscope and copy tip labels
df <- fread("data/metadata/sarbecovirus.raw.txt", sep = "\t")
sarbe <- str_split(colnames(df), pattern = "\\|", simplify = T)[, 1]
sarbe <- sarbe[sarbe != "current"]
sarbe

fwrite(tibble(sarbe), "data/metadata/sarbecovirus.accessions_only.csv", col.names = F)

# Get metadata
meta <- fread("data/metadata/beta_merged.050422.csv")
meta_filt <- meta %>% filter(accession %in% sarbe)

fwrite(meta_filt, "data/metadata/sarbecovirus.050422.csv")
