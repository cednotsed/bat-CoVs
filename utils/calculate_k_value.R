rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)
require(stringr)

dna <- readDNAStringSet("data/genomes/coronaviridae_n2118_novel_n9.080822.fna")

length(unique(names(dna))) == length(names(dna))

# Average length of genome
l <- mean(width(dna))

# Probability of observing a k-mer by chance
q <- 0.005

# no. of letters in alphabet
n_alphabet <- 4

# Calculate k value for mash
k <- ceiling(log(l * (1 - q)/ q, base = n_alphabet))
print(k)



