setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)
require(ape)

prefix <- "Beta"
multi_fastas <- list()
multi_dfs <- list()

for (prefix in c("Alpha", "Beta")) {
  file_list <- list.files(str_glue("data/genomes/raw_genomes/{prefix}"), full.names = T)
  
  multi_fasta <- foreach (file = file_list, .combine = "c") %do% {
    fasta <- read.FASTA(file)
  }
  
  fwrite("")
  if(length(multi_fastas) == 0) {
    multi_fastas <- multi_fasta
  } else {
    multi_fastas <- c(multi_fastas, multi_fasta)
  }
  
  
  morsels <- foreach (file = file_list) %do% {
    fasta <- read.FASTA(file)
    prop_n <- base.freq(fasta, all = T)[["n"]]
    seq_length <- length(as.character(fasta)[[1]])
    tibble(id = names(fasta), proportion_of_Ns = prop_n, contig_length = seq_length, genus = prefix)
  }
  
  multi_df <- bind_rows(morsels)
  multi_dfs <- bind_rows(multi_dfs, multi_df)
}

multi_dfs
multi_fastas

