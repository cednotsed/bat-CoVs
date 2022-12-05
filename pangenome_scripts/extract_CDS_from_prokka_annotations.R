rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)
require(ape)

meta_dir <- "data/metadata"
prokka_dir <- "results/prokka_out"
out_dir <- "results/prokka_out/extracted_spikes"
fastas <- list.files(prokka_dir, ".ffn", recursive = T, full.names = T)
fastas

extracted <- c()
accs <- c("1-GH087_assembly", "2-30B_assembly", "2-GH106_assembly", 
          "Sample-18_assembly", "MW719567")

group_name <- "sarbecovirus_spikes"
acc_regex <- paste0(accs, collapse = "|")
acc_paths <- fastas[grepl(acc_regex, fastas)]
acc_paths
  
spike_cdses <- foreach(isolate_name = accs, .combine = "c") %do% {
  fna_path <- str_glue("{prokka_dir}/{isolate_name}.ffn")
  fna <- readDNAStringSet(fna_path)
  spike_cds <- fna[grepl("Spike", names(fna))]
  names(spike_cds) <- isolate_name
  return(spike_cds)
}
 
spike_aas <- foreach(isolate_name = accs, .combine = "c") %do% {
  faa_path <- str_glue("{prokka_dir}/{isolate_name}.faa")
  faa <- readAAStringSet(faa_path)
  spike_aa <- faa[grepl("Spike", names(faa))]
  names(spike_aa) <- isolate_name
  return(spike_aa)
} 

writeXStringSet(spike_cdses, str_glue("{out_dir}/{group_name}.fna"))
writeXStringSet(spike_aas, str_glue("{out_dir}/{group_name}.faa"))

