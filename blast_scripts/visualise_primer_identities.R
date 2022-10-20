rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)

meta <- fread("data/metadata/sequencing_metadata_260622.csv") %>%
  mutate(Accession = ifelse(grepl("-", sample_id), 
                            sample_id, 
                            paste0("Sample-", sample_id))) %>%
  select(-sample_id)

# Parse output from blasting primers against novel genomes
blastout <- fread("results/blast_out/blast_primers_against_novel_genomes.tsv")
colnames(blastout) <- c("qseqid","sseqid","pident","aln_length","mismatch",
                        "gapopen", "qstart", "qend", "sstart", "send",
                        "evalue", "bitscore")
parsed <- blastout %>%
  separate(sseqid, into = c("Accession"), sep = "\\.") %>%
  right_join(meta, "Accession")


genomes_assembled <- c("1-GH087", "2-GH106", "Sample-25", 
                       "Sample-37", "4-126A", "Sample-30",
                       "2-30B", "5-129B", "Sample-18")
parsed %>% 
  filter(Accession %in% genomes_assembled) %>%
  arrange(Accession) %>%
  select(Accession, qseqid, pident, pcr_positive) %>% 
  ggplot(aes(x = Accession, y = pident, fill = qseqid)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(cols = vars(pcr_positive), scales = "free_x") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

min(meta$collection_date)
max(meta$collection_date)
meta
