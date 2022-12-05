rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)
require(DECIPHER)

primer_list <- list(
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7120317/
  vijgen_F = "ACWCARHTVAAYYTNAARTAYGC", vijgen_R = "TCRCAYTTDGGRTARTCCCA",
  # https://www.sciencedirect.com/science/article/pii/S1386653220301335
  xiu_F = "CCAARTTYTAYGGHGGNTGG", xiu_R = "TGTTGNGARCARAAYTCATGNGG",
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3298303/
  watanabe_F = "GTTGGGACTATCCTAAGTGTGA", watanabe_R = "CCATCATCAGATAGAATCATCATA",
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8067199/
  holbrook_F = "GGTGGGAYTAYCCHAARTGYGA", holbrook_R = "CCRTCATCAGAHARWATCAT")
  
# Disambiguate primers
primers <- DNAStringSet(unlist(primer_list))
primer_names <- names(primers)
disamb <- Disambiguate(primers)

disamb_primers <- foreach(i = seq(length(disamb)), .combine = "c") %do% {
  primer_temp <- disamb[[i]]
  
  temp_names <- paste0(rep(primer_names[i], length(primer_temp)),
                       ".",
                       seq(length(primer_temp)))
  names(primer_temp) <- temp_names
  primer_temp
}

dsz_names <- c("dsz_F.1", "dsz_F.2", "dsz_R.1", "dsz_R.2", "dsz_R.3")
dsz_primers <- DNAStringSet(c("TTATGGGTTGGGATTATC",
                              "TGATGGGATGGGACTATC",
                              "TCATCACTCAGAATCATCA",
                              "TCATCAGAAAGAATCATCA",
                              "TCGTCGGACAAGATCATCA"
                            ))
names(dsz_primers) <- dsz_names

final_primers <- c(dsz_primers, disamb_primers) 

final_primers

writeXStringSet(final_primers, "data/genomes/all_primers.fna")


