rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

wkdir <- "data/genomes/UK_bats/UK_bats_dedup"
fna_list <- list.files(wkdir, "fna")
fna_list

bat_seqs <- foreach(fna_path = fna_list, .combine = "c") %do% {
  print(fna_path)
  
  dna <- readDNAStringSet(str_glue("{wkdir}/{fna_path}"))
  
  if (!grepl("fna.gz", fna_path)) {
    # Append species name 
    species_name <- gsub(".fna", "", fna_path)
    nuc_names <- names(dna) 
    nuc_names <- str_split(nuc_names, " ", simplify = T)[, 1]
    nuc_names <- str_split(nuc_names, "\\.", simplify = T)[, 1]
    nuc_names <- gsub("lcl\\|", "", nuc_names)
    nuc_names <- gsub("LCL\\|", "", nuc_names)
    nuc_names <- paste(species_name, nuc_names, sep = "_")
    names(dna) <- nuc_names
    
    # Remove duplicate names
    dna <- dna[which(!duplicated(names(dna)))]
    
    return(dna)
  } else {
    species_name <- gsub(".dedup.fna.gz", "", fna_path)
    nuc_names <- names(dna) 
    nuc_names <- str_split(nuc_names, " ", simplify = T)[, 1]
    nuc_names <- paste(species_name, nuc_names, sep = "_")
    names(dna) <- nuc_names
    return(dna)
  }
}

writeXStringSet(bat_seqs, "data/genomes/UK_bats/all_UK_bat_sequences.fna")
