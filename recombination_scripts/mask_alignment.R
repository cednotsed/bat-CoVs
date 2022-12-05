rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(deepredeff)

# Functions
aaset_to_matrix <- function(aaset) {
  dat <- aasset_to_df(aaset) %>%
    column_to_rownames("name") %>%
    separate(seq, into = as.character(seq(0, width(aaset)[1])), sep = "") %>%
    select(-"0")
  return(dat)
}

matrix_to_aaset <- function(mat) {
  temp_mat <- tibble(seq = apply(mat, 1, paste0, collapse = ""))
  rownames(temp_mat) <- rownames(mat)
  aaset <- AAStringSet(x = temp_mat$seq, use.names = T)
  return(aaset)
}


aln_paths <- list.files("data/alignments", ".aln", full.names = T)
aln_paths <- aln_paths[!grepl("mask|log|insertions|trim", aln_paths)]
aln_path <- aln_paths[grepl("Sarbe", aln_paths)]
aln_path
aln <- readDNAStringSet(aln_path)
n_aln <- length(aln)
n_pos <- unique(width(aln))

aln_df <- aaset_to_matrix(aln)

# Remove gaps to reference
ref <- "NC_025217"
gap_pos <- as.character(aln_df[ref, ]) == "-"
aln_no_gaps_df <- aln_df[, !gap_pos]

for (pos in seq(n_pos)) {
  temp_pos <- aln_no_gaps_df[, pos]
  prop_gaps <- sum(temp_pos == "-") / n_aln
  
  # Replace gappy positions with Ns
  if (prop_gaps > 0.2) {
    aln_no_gaps_df[, pos] <- "N"
  }
}
  
aln_no_gaps <- matrix_to_aaset(aln_no_gaps_df)


save_name <- gsub(".aln", str_glue(".masked.aln"), aln_path)
writeXStringSet(aln_no_gaps, save_name)

