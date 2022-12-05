rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)

parsed <- fread("results/blast_out/blast_contigs.parsed.csv") %>%
  filter(contig_length > 500)
# parsed <- fread("results/blast_out/blast_qc_contigs.parsed.csv")

# Get best hits for each contig
contig_list <- deframe(parsed %>% distinct(contig_id))

best_morsels <- foreach(contig = contig_list) %do% {
  parsed %>% 
    filter(contig_id == contig) %>%
    arrange(desc(bitscore)) %>%
    head(1)
}

best_hits <- bind_rows(best_morsels) %>%
  filter(!grepl("MHV|Murine", genbank_title)) %>%
  arrange(desc(contig_length))

best_hits %>%
  group_by(sample_id, accession) %>%
  summarise(max = max(contig_length)) %>%
  View()

fwrite(best_hits, "results/blast_out/best_contig_hits_to_reference.csv")
best_hits
