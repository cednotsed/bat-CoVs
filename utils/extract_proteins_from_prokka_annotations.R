rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)
require(ape)

meta_dir <- "data/metadata"
prokka_dir <- "data/proteins/prokka_proteins"
out_path <- "data/proteins/extracted_proteins/novel_spikes.faa"
fastas <- list.files(prokka_dir, ".faa", recursive = T, full.names = T)
fastas <- fastas[!grepl(".tmp.", fastas)]

extracted <- c()
missing <- c()

for (fasta_path in fastas) {
  fasta <- read.FASTA(fasta_path, type = "AA")
  isolate <- str_split(fasta_path, "/")[[1]][4]
  spike <- fasta[grepl("gene=S", names(fasta))]
  
  if(length(spike) == 0) {
    print(str_glue("{isolate} has no spike!!"))
    missing <- c(missing, isolate)
  } else {
    print(str_glue("No. of spikes for {isolate}: {length(spike)}"))
    names(spike) <- isolate
    
    if (length(extracted) == 0) {
      extracted <- spike
    } else {
      extracted <- c(extracted, spike)
    }
  }
}

write.FASTA(extracted, out_path)


## Run extract_batch_entrez_annotations.R first before doing this ##
entrez_path <- "data/proteins/extracted_proteins/human_infecting_CoVs_spike.faa"
final_path <- "data/proteins/extracted_proteins/all_spikes.faa"
sarbecovirus_path <- "data/proteins/extracted_proteins/sarbecovirus_spikes.faa"
MERS_path <- "data/proteins/extracted_proteins/MERS_spikes.faa"
batch_entrez <- read.FASTA(entrez_path, type = "AA")
all_spikes <- c(extracted, batch_entrez)
sarbecovirus <- all_spikes[grepl("R|MZ|NC_004718|NC_045512", names(all_spikes))]
MERS <- all_spikes[grepl("PAGB|NC_019843", names(all_spikes))]
write.FASTA(all_spikes, final_path)
write.FASTA(sarbecovirus, sarbecovirus_path)
write.FASTA(MERS, MERS_path)
  