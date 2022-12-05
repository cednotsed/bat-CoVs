rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

base_dir <- "results/assembly_qced/coronaspades_out/scaffolds/"
file_list <- list.files(base_dir)

pooled_contigs <- foreach(file_name = file_list, .combine = "c") %do% {
  prefix <- gsub("_scafolds.fasta", "", file_name)
  contigs <- readDNAStringSet(str_glue("{base_dir}/{file_name}"))
  contig_names <- names(contigs)
  contig_names <- paste0(prefix, ".", contig_names)
  names(contigs) <- contig_names
  return(contigs)
}

# Remove contigs <= 500nt
pooled_filt <- pooled_contigs[width(pooled_contigs) > 500]

writeXStringSet(pooled_filt, "results/assembly_qced/pooled_contigs_gt500.fna")
writeXStringSet(pooled_contigs, "results/assembly_qced/pooled_contigs.fna")
