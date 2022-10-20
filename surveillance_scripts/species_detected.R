rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)

meta <- fread("data/metadata/all_meta.n2127.080822.csv")
meta2 <- fread("data/metadata/sequencing_metadata_260622.csv") %>%
  mutate(sample_id = ifelse(grepl("-", sample_id), sample_id, paste0("Sample-", sample_id)))
blastout <- fread("results/blast_out/best_contig_hits_to_reference.csv")
sum_df <- blastout %>% 
  filter(contig_length > 500) %>% 
  mutate(Species = ifelse(grepl("Bat coronavirus|Bat alpha|Myotis bat", GenBank_Title), 
                          str_glue("{Host} coronavirus"), Species)) %>% 
  group_by(sample_id, Species) %>% 
  summarise(total_length = sum(contig_length), n_contigs = n_distinct(contig_id)) %>%
  right_join(meta2) %>%
  arrange(Species)

sum_df %>% 
  mutate(sample_id = factor(sample_id, unique(sum_df$sample_id)),
         Species = factor(Species, rev(unique(sum_df$Species)))) %>% 
  ggplot(aes(x = Species, y = sample_id)) +
  geom_tile() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
