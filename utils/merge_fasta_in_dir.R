setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)
require(ape)

prefix <- "Beta"
multi_fastas <- list()
multi_dfs <- list()

for (prefix in c("MERS_like", "sarbecovirus_like")) {
  file_list <- list.files(str_glue("data/genomes/{prefix}"), full.names = T)
  
  multi_fasta <- foreach (file = file_list, .combine = "c") %do% {
    fasta <- read.FASTA(file)
  }
  
  write.FASTA(multi_fasta, str_glue("data/genomes/{prefix}.fasta"))
}



