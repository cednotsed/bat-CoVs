rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)
require(ape)

meta_dir <- "data/metadata"
prokka_dir <- "data/proteins/prokka_proteins"
out_dir <- "data/genomes/extracted_spikes"
fastas <- list.files(prokka_dir, ".ffn", recursive = T, full.names = T)
fastas

extracted <- c()
accs <- c("1-GH087", "2-30B", "2-GH106", 
          "4-126A", "5-129B", "MZ937003",
          "Sample-18", "Sample-25", "Sample-30",
          "Sample-37")

group_name <- "all_novel_spikes"
acc_regex <- paste0(accs, collapse = "|")
acc_paths <- fastas[grepl(acc_regex, fastas)]
acc_paths
  
spikes <- foreach(faa_path = acc_paths, .combine = "c") %do% {
  isolate_name <- str_split(faa_path, "/")[[1]][4]
  fasta <- read.FASTA(faa_path, type = "DNA")
  spike <- fasta[grepl("Spike", names(fasta))]
  names(spike) <- isolate_name
  return(spike)
}
  
write.FASTA(spikes, str_glue("{out_dir}/{group_name}.fna"))

