rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

pangolin <- readDNAStringSet("data/genomes/gisaid/Manis_GISAID_080822.fasta")
names(pangolin) <- str_split(names(pangolin), "\\|", simplify = T)[, 2]

bats <- readDNAStringSet("data/genomes/gisaid/Rhinolophus_GISAID_080822.fasta")
names(bats) <- str_split(names(bats), "\\|", simplify = T)[, 2]

final <- c(pangolin, bats)

writeXStringSet(final, "data/genomes/gisaid/all_GISAID_080822.fna")
