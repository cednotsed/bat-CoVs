setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)
require(ape)

prefixes <- list.dirs("results/zoonotic_assesment/prokka_annotations", full.names = T)
prefixes <- prefixes[!grepl(".fna", prefixes)]
prefixes <- prefixes[2:length(prefixes)]
prefixes

spikes <- foreach(prefix = prefixes, .combine = "c") %do% { 
  fna_path <- list.files(prefix, ".ffn", full.names = T)
  
  cds <- read.FASTA(fna_path)
  spike <- cds[grepl("Spike|glyco", names(cds), ignore.case = T)]
  
  prefix_name <- rev(str_split(prefix, "/")[[1]])[1]
  
  names(spike) <- str_glue("{prefix_name} spike")
  spike
}

write.FASTA(spikes, "results/in_vitro_out/novel_spike_cds.n9.fna")

spikes <- foreach(prefix = prefixes, .combine = "c") %do% { 
  fna_path <- list.files(prefix, ".faa", full.names = T)
  
  cds <- read.FASTA(fna_path, type = "AA")
  spike <- cds[grepl("Spike|glyco", names(cds), ignore.case = T)]
  
  prefix_name <- rev(str_split(prefix, "/")[[1]])[1]
  
  names(spike) <- str_glue("{prefix_name} spike")
  spike
}

write.FASTA(spikes, "results/in_vitro_out/novel_spike_proteins.n9.faa")
names(spikes)
table(as.character(spikes["5-129B spike"]))
