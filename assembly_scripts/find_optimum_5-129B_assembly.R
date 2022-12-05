rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

file_dir <- "results/prokka_out/"
files <- list.files(file_dir, "5-129B", full.names = T)
files <- files[grepl(".tsv", files)]
files <- files[!grepl("assembly", files)]
files <- files[!grepl("3prime", files)]

morsels <- foreach(file = files) %do% {
  fread(file) %>%
    mutate(candidate = gsub(file_dir, "", file),
           presence = T, .before = 1) %>%
    mutate(gene = ifelse(gene == "", str_glue("hyp_{length_bp}"), gene))
}

bind_rows(morsels) %>% 
  ggplot(aes(x = gene, y = candidate, fill = presence)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Gene annotation", y = "Candidate genome", fill = "Gene presence")

ggsave("results/assembly/constructed_genomes/5-129B_candidates/optimal_5-129B_candidate.png",
       dpi = 600, width = 5, height = 5)
