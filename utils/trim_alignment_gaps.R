setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(ape)
require(foreach)

aln <- read.dna("data/alignments/sarbecovirus.050422.aln", 
                format = "fasta",
                as.matrix = T)

positions <- foreach (pos = seq(ncol(aln)), .combine = "c") %do% {
  base_freq <- base.freq(aln[, pos], all = T)
  ifelse(any((base_freq > 0.2)[5:length(base_freq)]), F, T)
}

trimmed_aln <- aln[, positions]
write.FASTA(trimmed_aln, "data/alignments/sarbecovirus.050422.trimmed.aln")
