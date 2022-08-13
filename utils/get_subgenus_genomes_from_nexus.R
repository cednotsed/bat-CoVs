rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(ape)

all_covs <- read.FASTA("data/genomes/coronaviridae_n2118_novel_n9.080822.fna")
trees <- list.files("data/trees", ".nexus", full.names = T)

for (tree_file in trees) {
  tree_tmp <- read.nexus(tree_file)
  
  # Parse accessions
  accessions <- str_split(tree_tmp$tip.label, "\\|", simplify = T)[, 1]
  accessions <- as.vector(accessions)
  accessions <- gsub("'", "", accessions)
  
  # Get genomes
  genomes <- all_covs[names(all_covs) %in% accessions]
  fasta_file <- str_split(tree_file, "/")[[1]][3]
  fasta_file <- gsub(".nexus", str_glue(".n{length(genomes)}.fna"), fasta_file)
  write.FASTA(genomes, str_glue("data/genomes/{fasta_file}"))
}
