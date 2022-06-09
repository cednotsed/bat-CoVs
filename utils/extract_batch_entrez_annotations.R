rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)
require(ape)

meta_dir <- "data/metadata"
fasta_path <- "data/proteins/batch_entrez_proteins/human_infecting_CoVs.faa"
out_path <- "data/proteins/extracted_proteins/human_infecting_CoVs_spike.faa"
fasta <- read.FASTA(fasta_path, type = "AA")

names(fasta)[grepl("NC_045512", names(fasta))]
spike <- fasta[grepl("gene=S|protein=spike", names(fasta))]
isolates <- tibble(acc = names(spike)) %>%
  separate(acc, into = c(NA, "acc"), sep = "\\|") %>%
  separate(acc, into = c("accession"), sep = "\\.")

names(spike) <- isolates$accession
write.FASTA(spike, out_path)

isolates
## Rename batch entrez spike proteins ##


