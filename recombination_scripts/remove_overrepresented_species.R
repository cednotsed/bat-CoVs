rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)
require(ape)

tree <- read.tree("data/trees/Sarbecovirus_subtree.n534.trim_to_28261pos.tree.contree")
tree

meta <- fread("data/metadata/all_meta.n2127.080822.csv")
             
meta_filt <- tibble(accession = tree$tip.label) %>% 
  left_join(meta)

SARS <- deframe(meta_filt %>% 
  filter(grepl("SARS coronavirus", genbank_title) &
           !grepl("Bat SARS coronavirus", genbank_title, ignore.case = T)) %>% 
  select(accession))

SARS2 <- deframe(meta_filt %>% 
  filter(grepl("SARS-CoV-2|Severe acute respiratory syndrome coronavirus 2", 
               genbank_title, 
               ignore.case = T)) %>% 
           select(accession))

to_remove <- c(SARS, SARS2)
to_remove <- to_remove[!(to_remove %in% c("NC_004718", "MW206198"))]

tree_filt <- drop.tip(tree, to_remove)
tree_filt <- root(tree_filt, "NC_025217", resolve.root = T)
# write.tree(tree_filt, "data/trees/for_recombination_analysis/Sarbecovirus.trim_to_28261pos.n218.contree")

aln <- readDNAStringSet("data/alignments/Sarbecovirus_subtree.n534.masked.aln")

aln_filt <- aln[tree_filt$tip.label]

for(acc in names(aln_filt)) {
  aln_filt[[acc]] <- gsub("M|R|W|S|Y|K|V|H|D|B", "N", aln_filt[[acc]])
}

all(names(aln_filt) == tree_filt$tip.label)
writeXStringSet(aln_filt, 
          "data/alignments/for_recombination_analysis/Sarbecovirus.n218.masked.aln")


mash_tree <- read.tree("data/trees/coronaviridae_NJ_mash.tree")
to_keep <- names(aln_filt)
to_drop <- mash_tree$tip.label[!(mash_tree$tip.label %in% to_keep)]

mash_tree_filt <- drop.tip(mash_tree, to_drop)
mash_tree_filt <- root(mash_tree_filt, "NC_025217", resolve.root = T)
# write.tree(mash_tree_filt, "results/recombination_out/Sarbecovirus.trim_to_28261pos.n218.mash.tree")
