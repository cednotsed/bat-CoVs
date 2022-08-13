rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(ape)
require(foreach)

aln_paths <- list.files("data/alignments", ".aln", full.names = T)
aln_paths <- aln_paths[!grepl("log|insertions|trimmed", aln_paths)]

for (aln_path in aln_paths) {
  print(aln_path)
  aln <- read.dna(aln_path, 
                  format = "fasta",
                  as.matrix = T)
  n_aln <- nrow(aln)
  
  pos_to_keep <- foreach(i = seq(ncol(aln)), .combine = "c") %do% {
    small_gap <- sum(as.character(aln[, i]) == "-") / n_aln <= 0.2
    
    if(small_gap) {
        return(i)
    } else {
        return(NA)
      }
  }
  
  pos_to_keep <- pos_to_keep[!is.na(pos_to_keep)]
  length(pos_to_keep)
  trimmed <- aln[, pos_to_keep]
  save_name <- gsub(".aln", str_glue(".trim_to_{ncol(trimmed)}pos.aln"), aln_path)
  write.FASTA(trimmed, save_name)
}
