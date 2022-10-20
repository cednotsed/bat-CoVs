rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)
require(ape)

meta_dir <- "data/metadata"
prokka_dir <- "data/proteins/prokka_proteins"
out_dir <- "data/proteins/extracted_spikes"
fastas <- list.files(prokka_dir, ".faa", recursive = T, full.names = T)
# fastas <- fastas[!grepl("0813", fastas)]

faa_groups <- list(sarbecovirus = c("NC_045512", "NC_004718", "MZ937003", "MN996532",
                                    "1-GH087", "2-GH106", "2-30B", "Sample-18"),
                   mers = c("NC_019843", "5-129B"))
extracted <- c()
missing <- c()

for(i in seq(length(faa_groups))) {
  group_name <- names(faa_groups)[i]
  acc_regex <- paste0(faa_groups[[i]], collapse = "|")
  acc_paths <- fastas[grepl(acc_regex, fastas)]
  acc_paths
  
  spikes <- foreach(faa_path = acc_paths, .combine = "c") %do% {
    isolate_name <- str_split(faa_path, "/")[[1]][4]
    fasta <- read.FASTA(faa_path, type = "AA")
    spike <- fasta[grepl("Spike", names(fasta))]
    names(spike) <- isolate_name
    return(spike)
  }
  
  write.FASTA(spikes, str_glue("{out_dir}/{group_name}_spikes.faa"))
}

# 
# 
# 
# 
# for (fasta_path in fastas) {
#   fasta <- read.FASTA(fasta_path, type = "AA")
#   isolate <- str_split(fasta_path, "/")[[1]][4]
#   spike <- fasta[grepl("gene=S", names(fasta))]
#   
#   if(length(spike) == 0) {
#     print(str_glue("{isolate} has no spike!!"))
#     missing <- c(missing, isolate)
#   } else {
#     print(str_glue("No. of spikes for {isolate}: {length(spike)}"))
#     names(spike) <- isolate
#     
#     if (length(extracted) == 0) {
#       extracted <- spike
#     } else {
#       extracted <- c(extracted, spike)
#     }
#   }
# }
# 
# write.FASTA(extracted, out_path)


# ## Run extract_batch_entrez_annotations.R first before doing this ##
# entrez_path <- "data/proteins/extracted_proteins/human_infecting_CoVs_spike.faa"
# final_path <- "data/proteins/extracted_proteins/all_spikes.faa"
# sarbecovirus_path <- "data/proteins/extracted_proteins/sarbecovirus_spikes.faa"
# MERS_path <- "data/proteins/extracted_proteins/MERS_spikes.faa"
# batch_entrez <- read.FASTA(entrez_path, type = "AA")
# all_spikes <- c(extracted, batch_entrez)
# sarbecovirus <- all_spikes[grepl("R|MZ|NC_004718|NC_045512", names(all_spikes))]
# MERS <- all_spikes[grepl("PAGB|NC_019843", names(all_spikes))]
# write.FASTA(all_spikes, final_path)
# write.FASTA(sarbecovirus, sarbecovirus_path)
# write.FASTA(MERS, MERS_path)
  