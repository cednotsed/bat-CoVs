rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(ape)
require(foreach)

prokka_dir <- "data/proteins/prokka_proteins/coronaviridae.n2089_subset/"
df <- fread("results/pangenome/gene_presence_absence.csv")
core_genes <- c("rep", "S", "M", "N")
core_df <- df %>% 
  filter(Gene %in% core_genes)
to_remove <- colnames(core_df)[2:14]
core_parsed <- core_df %>%
  select(!all_of(to_remove)) %>%
  pivot_longer(!Gene, names_to = "accession_id", values_to = "cds")
core_parsed
acc_list <- deframe(core_parsed %>% distinct(accession_id))
acc_list

# Build core genomes
foreach(acc = acc_list) %do% {
  
  # Get CDS ids
  cds_df <- core_parsed %>%
    filter(accession_id == acc)
  cds_df <- tibble(Gene = core_genes) %>%
    left_join(cds_list)
  cds_list <- cds_df$cds
  
  # Load CDS fna
  ffn <- read.FASTA(str_glue("{prokka_dir}/{acc}.ffn"))
  
  # Build single genome
  foreach(cds = cds_list) %do% {
    gene <- ffn[grepl(cds, names(ffn))]
  }
  
}

View(df)
