rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

primers <- readDNAStringSet("data/genomes/all_primers.fna")
blastout <- fread("results/blast_out/blast_primers_against_novel_genomes.tsv")
colnames(blastout) <- c("qseqid","sseqid","pident","aln_length","mismatch",
                        "gapopen", "qstart", "qend", "sstart", "send",
                        "evalue", "bitscore")

parsed <- blastout %>% 
  # filter(sseqid %in% c("Sample-25", "Sample-30")) %>%
  mutate(n_matches = aln_length - mismatch,
         primer_len = ifelse(grepl("PC2S2", qseqid), 18, 19)) %>%
  mutate(prop_matches = n_matches / primer_len) %>%
  separate(qseqid, into = c("study"), sep = "_", remove = F) %>%
  mutate(primer_type = ifelse(grepl("_F", qseqid), "F", "R"))

primer_lengths <- tibble(qseqid = names(primers), primer_length = width(primers))

parsed %>% 
  left_join(primer_lengths) %>%
  mutate(prop_match = (aln_length - mismatch) / primer_length) %>%
  group_by(study, primer_type, sseqid) %>%
  summarise(max_prop = max(prop_match)) %>%
  mutate(primer = str_glue("{study}_{primer_type}")) %>%
  ggplot(aes(x = primer, y = sseqid, fill = max_prop)) +
  geom_tile() +
  geom_text(aes(label = round(max_prop, 2))) +
    labs(x = "Primer", y = "Sample", fill = "Prop. matches") +
    scale_fill_gradient(low = "slateblue", high = "indianred")

parsed
distinct_pairs <- parsed %>% 
  distinct(qseqid, sseqid)

# Get best hits
morsels <- foreach(i = seq(nrow(distinct_pairs))) %do% {
  qseq <- distinct_pairs[i, ]$qseqid
  sseq <- distinct_pairs[i, ]$sseqid
  
  parsed %>% 
    filter(sseqid == sseq, qseqid == qseq) %>%
    arrange(desc(prop_matches)) %>%
    head(1)
}

meta <- fread("data/metadata/sequencing_metadata_260622.csv") %>%
  mutate(sseqid = ifelse(grepl("-", sample_id), sample_id, paste0("Sample-", sample_id))) %>%
  select(sseqid, host_species, genbank_title)

# bind_rows(morsels) %>%
#   left_join(meta) %>%
#   mutate(primer_type = ifelse(grepl("PC2As", qseqid), "Reverse", "Forward")) %>%
#   ggplot(aes(x = qseqid, y = sseqid, fill = prop_matches)) +
#   geom_tile() +
#   geom_text(aes(label = n_matches)) +
#   facet_grid(rows = vars(host_species), cols = vars(primer_type), scales = "free") +


ggsave("results/surveillance_out/primer_match.png", dpi = 600, height = 8, width = 8)


  
  
