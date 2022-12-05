rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(ape)
require(Biostrings)


rhgb01_acc <- fread("data/metadata/for_recombination_analysis/rhgb01-like.txt", header = F)$V1
sars_acc <- fread("data/metadata/for_recombination_analysis/sars-like.txt", header = F)$V1
other_acc <- fread("data/metadata/for_recombination_analysis/other_sarbecoviruses.txt", header = F)$V1


tree <- read.tree("data/trees/for_recombination_analysis/Sarbecovirus.trim_to_28261pos.n218.contree")

all_acc <- tree$tip.label

rhgb01_tree <- drop.tip(tree, all_acc[!(all_acc %in% rhgb01_acc)])
sars_tree <- drop.tip(tree, all_acc[!(all_acc %in% sars_acc)])
other_tree <- drop.tip(tree, all_acc[!(all_acc %in% other_acc)])

rhgb01_tree <- root(rhgb01_tree, "NC_025217", resolve.root = F)
sars_tree <- root(sars_tree, "NC_025217", resolve.root = F)
other_tree <- root(other_tree, "NC_025217", resolve.root = F)

write.tree(rhgb01_tree, str_glue("data/trees/for_recombination_analysis/rhgb01-like.subtree.n{Ntip(rhgb01_tree)}.tree"))
write.tree(sars_tree, str_glue("data/trees/for_recombination_analysis/sars-like.subtree.n{Ntip(sars_tree)}.tree"))
write.tree(other_tree, str_glue("data/trees/for_recombination_analysis/other_sarbecoviruses.subtree.n{Ntip(other_tree)}.tree"))


# Get alignments
all_alns <- readDNAStringSet("data/alignments/for_recombination_analysis/Sarbecovirus.trim_to_28261pos.n218.aln", format = "fasta")

rhgb01_aln <- all_alns[names(all_alns) %in% rhgb01_tree$tip.label]
sars_aln <- all_alns[names(all_alns) %in% sars_tree$tip.label]
other_aln <- all_alns[names(all_alns) %in% other_tree$tip.label]

writeXStringSet(rhgb01_aln, str_glue("data/alignments/for_recombination_analysis/rhgb01-like.subtree.n{Ntip(rhgb01_tree)}.aln"))
writeXStringSet(sars_aln, str_glue("data/alignments/for_recombination_analysis/sars-like.subtree.n{Ntip(sars_tree)}.aln"))
writeXStringSet(other_aln, str_glue("data/alignments/for_recombination_analysis/other_sarbecoviruses.subtree.n{Ntip(other_tree)}.aln"))

