setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

iter_list <- list(sarbecovirus_like = list(file = "sarbecovirus_like.aln",
                                           spike_coords = c(21563, 25384)),
                  MERS_like = list(file = "MERS_like.aln",
                                   spike_coords = c(21456, 25517)))

morsels <- foreach(i = seq(length(iter_list))) %do% {
  file <- iter_list[[i]][["file"]]
  spike_coords <- iter_list[[i]][["spike_coords"]]
  
  seqs <- readDNAStringSet(str_glue("data/alignments/{file}"), 
                           format = "fasta")
  
  spike_seqs <- subseq(seqs, spike_coords[1], spike_coords[2])
  
  # Remove gaps
  for (k in seq(length(spike_seqs))) {
    spike_seqs[k] <- DNAStringSet(gsub("-", "", spike_seqs[k]))
  }

  spike_AA <- translate(spike_seqs, if.fuzzy.codon = "X")
  
  crumbs <- foreach (j = seq(2, length(spike_AA))) %do% {
    overall_ident <- pid(pairwiseAlignment(seqs[1], seqs[j]), type = "PID2")
    spike_nt_ident <- pid(pairwiseAlignment(spike_seqs[1], spike_seqs[j], type = "PID2"))
    spike_aa_ident <- pid(pairwiseAlignment(spike_AA[1], spike_AA[j], type = "PID3"))
    accession <- names(seqs)[j]
    
    tibble(accession = accession,
           reference = names(seqs)[1],
           group = names(iter_list)[i],
           overall_nt_homology = overall_ident,
           spike_nt_homology = spike_nt_ident, 
           spike_AA_homology = spike_aa_ident)
  }
  
  bind_rows(crumbs)
  
}

final <- bind_rows(morsels)

View(final)
