rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)

blastout <- fread("results/blast_out/blast_to_UK_bats.tsv")
meta <- fread("data/metadata/sequencing_metadata_260622.csv") %>%
  mutate(sample_id = ifelse(grepl("-", sample_id), sample_id, paste0("Sample-", sample_id)))
colnames(blastout) <- c("qseqid","sseqid","pident","aln_length","mismatch",
                        "gapopen", "qstart", "qend", "sstart", "send",
                        "evalue", "bitscore")

parsed <- blastout %>%
  separate(qseqid, into = c("sample_id", NA, "contig_id"), sep = "\\.") %>%
  mutate(sample_id = gsub("_scaffolds", "", sample_id)) %>%
  separate(sseqid, into = c("Accession"), sep = "\\.") %>%
  separate(contig_id, into = c(NA, NA, NA, "contig_length"), sep = "_", remove = F) %>%
  mutate(contig_length = as.numeric(contig_length),
         bitscore = as.numeric(bitscore)) %>%
  separate(Accession, into = c("genus", "species"), sep = "_") %>%
  mutate(best_species_hit = paste(genus, species))

best_hits <- parsed %>%
  filter(pident > 90) %>%
  group_by(best_species_hit, sample_id) %>%
  summarise(max_id = max(pident))

meta %>%
  left_join(best_hits) %>%
  View()

fread("data/metadata/sample_collection_metadata.parsed.csv") %>% 
  filter(collection_date != "Unknown") %>%
  mutate(collection_date = as.Date(collection_date, "%d/%m/%Y")) %>%
  summarise(min = min(collection_date),
            max = max(collection_date))

fread('data/metadata/sequencing_metadata_260622.csv') %>%
  filter(pcr_positive == "y") %>%
  distinct(host_species)
