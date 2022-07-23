rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)
require(ape)

fasta <- read.FASTA("data/genomes/coronaviridae_taxid11118_complete_exclude_provirus_040722/coronaviridae.n2089_subset.fna")

for (i in seq(length(fasta))) {
  fna <- fasta[i]
  acc <- str_split(names(fna), " ")[[1]][1]
  acc <- str_split(acc, "\\.")[[1]][1]
  
  names(fna) <- acc
  write.FASTA(fna, str_glue("data/genomes/coronaviridae_taxid11118_complete_exclude_provirus_040722/coronaviridae.n2089_subset/{acc}.fna"))
}
