rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(ape)

tree_paths <- list.files("data/trees/", ".treefile", full.names = T)
tree_paths <- tree_paths[!grepl("ANNOT", tree_paths)]
meta <- fread("data/metadata/all_meta.n2127.080822.csv")

# Annotate tree tip labels
for(tree_path in tree_paths) {
  tree <- read.tree(tree_path)

    # Match metadata
  tips <- tree$tip.label
  meta_parsed <- meta[match(tips, meta$accession), ]
  tip_label <- meta_parsed %>% mutate(tip_label = str_glue("{accession}|{host}|{genbank_title}"))
  tree$tip.label <- tip_label$tip_label
  save_path <- gsub(".treefile", ".ANNOT.treefile", tree_path)
  write.tree(tree, save_path)
  # tip annotation
}
