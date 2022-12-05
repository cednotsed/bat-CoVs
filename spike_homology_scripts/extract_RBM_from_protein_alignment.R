setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)
require(ape)
require(deepredeff)
require(Biostrings)

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

sarbe <- readAAStringSet("data/alignments/spike_protein_alignments/sarbecovirus_spikes.faa.aln", format = "fasta")
sarbe_df <- aaset_to_matrix(sarbe)

# Remove gaps
ref <- "NC_045512.faa"
gap_pos <- as.character(sarbe_df[ref, ]) == "-"
sarbe_no_gaps_df <- sarbe_df[, !gap_pos]
sarbe_no_gaps <- matrix_to_aaset(sarbe_no_gaps_df)
writeXStringSet(sarbe_no_gaps, "data/alignments/spike_protein_alignments/spikes_aligned_to_SARS-CoV-2.faa.aln")

# Choose RBD location
RBM <- subseq(sarbe_no_gaps, 438, 506)
# writeXStringSet(RBM, "data/alignments/sarbecovirus_RBM.faa.aln")

# Choose SARS-CoV-2 contact residues
cov2_contact_res <- c(417, 446, 449,
         453, 455, 456, 
         475, 486, 487,
         489, 493, 496,
         498, 500, 501,
         502, 505)
cov2_contact_df <- sarbe_no_gaps_df[, cov2_contact_res]
cov2_contact <- matrix_to_aaset(cov2_contact_df)
# writeXStringSet(cov2_contact, "data/alignments/spike_protein_alignments/SARS-Cov-2_contact_residues.faa.aln")

cov1_contact_res <- c(402, 426, 436, 
                      440, 442, 472,
                      473, 475, 479, 
                      484, 486, 487, 
                      488, 491)
ref <- "NC_004718.faa"
gap_pos <- as.character(sarbe_df[ref, ]) == "-"
sarbe_no_gaps_cov1_df <- sarbe_df[, !gap_pos]
sarbe_no_gaps_cov1 <- matrix_to_aaset(sarbe_no_gaps_cov1_df)
writeXStringSet(sarbe_no_gaps_cov1, "data/alignments/spike_protein_alignments/spikes_aligned_to_SARS-CoV-1.faa.aln")

cov1_contact_df <- sarbe_no_gaps_cov1_df[, cov1_contact_res]
cov1_contact <- matrix_to_aaset(cov1_contact_df)
# writeXStringSet(cov1_contact, "data/alignments/spike_protein_alignments/SARS-CoV-1_contact_residues.faa.aln")

# ## MERS ##
# MERS <- readAAStringSet("data/alignments/spike_protein_alignments/mers_spikes.faa.aln", format = "fasta")
# MERS_df <- aaset_to_matrix(MERS)
# 
# # Remove gaps
# ref <- "NC_019843"
# gap_pos <- as.character(MERS_df[ref, ]) == "-"
# MERS_no_gaps_df <- MERS_df[, !gap_pos]
# 
# MERS_contact_res <- c(499, 501, 502,
#                       506, 510, 511, 
#                       513, 537, 538,
#                       539, 540, 542,
#                       553, 555)
# MERS_contact_df <- MERS_no_gaps_df[, MERS_contact_res]
# MERS_contact <- matrix_to_aaset(MERS_contact_df)
# MERS_contact
# 
# writeXStringSet(MERS_contact, "data/alignments/spike_protein_alignments/mers_contact_residues.faa.aln")
# 
# 
# test <- str_split("RVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFHHHHHH",
#                   pattern = "")[[1]]
# test[417]
# test[1, 417:418]
