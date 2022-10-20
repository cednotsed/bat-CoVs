rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)
require(ape)

# fasta <- read.FASTA("data/genomes/coronaviridae_taxid11118_complete_exclude_provirus_040722/coronaviridae.n2089_subset.fna")
# out_dir <- "data/genomes/coronaviridae_taxid11118_complete_exclude_provirus_040722/coronaviridae.n2089_subset"

fasta <- read.FASTA("data/genomes/human_infecting_CoVs/all_references/human_infecting_CoVs.fna")
out_dir <- "data/genomes/human_infecting_CoVs"

for (i in seq(length(fasta))) {
  fna <- fasta[i]
  acc <- str_split(names(fna), " ")[[1]][1]
  acc <- str_split(acc, "\\.")[[1]][1]
  
  names(fna) <- acc
  write.FASTA(fna, str_glue("{out_dir}/{acc}.fna"))
}
