rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

meta2 <- fread("data/metadata/sequencing_metadata_260622.csv") %>%
  mutate(sample_id = ifelse(grepl("-", sample_id), sample_id, paste0("Sample-", sample_id)))
ass_results <- fread("results/assembly/assembly_results.csv")
best_df <- fread("results/blast_out/best_contig_hits_to_reference.csv") %>%
  select(-host_species)

# Get best match
best_morsels <- foreach(sample_name = unique(best_df$sample_id)) %do% {
  best_df %>% 
    filter(sample_id == sample_name) %>%
    arrange(desc(bitscore)) %>%
    head(1)
}

plot_df <- bind_rows(best_morsels) %>%
  left_join(meta2 %>% select(sample_id, host_species)) %>%
  arrange(desc(host_species))

plot_df %>%
  mutate(sample_id = factor(sample_id, unique(plot_df$sample_id))) %>%
  ggplot(aes(x = genbank_title, y = sample_id, fill = host_species)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save genbank genomes that are similar for alignment and pident analysis
mers_acc <- plot_df %>% 
  filter(grepl("Middle", genbank_title)) %>%
  distinct(accession)

bat_acc <- plot_df %>% 
  filter(grepl("Bat coronavirus isolate", genbank_title)) %>%
  distinct(accession)

all_fnas <- readDNAStringSet("data/genomes/coronaviridae_n2118.080822.with_novel.fna")

mers <- all_fnas[mers_acc$Accession]
bat_covs <- all_fnas[bat_acc$Accession]

writeXStringSet(mers, "data/genomes/similar_genbank_matches/MERS-like_hits.fna")
writeXStringSet(bat_covs, "data/genomes/similar_genbank_matches/alpha-like_hits.fna")

# Merge genbank titles
plot_df <- bind_rows(best_morsels) %>%
  left_join(meta2 %>% select(sample_id, host_species)) %>% 
  mutate(final_species = case_when(grepl("M.dau", genbank_title) ~ "M. daubentonii bat coronavirus",
                                   grepl("P.pyg", genbank_title) ~ "P. pygmaeus bat coronavirus",
                                   grepl("RhGB01", genbank_title) ~ "Sarbecovirus RhGB01",
                                   grepl("Middle", genbank_title) ~ "MERS-related bat coronavirus")) %>% 
  arrange(desc(host_species)) 

plot_df
plot_df %>%
  left_join(ass_results, "sample_id") %>%
  mutate(sample_id = factor(sample_id, unique(plot_df$sample_id))) %>%
  ggplot(aes(y = final_species, x = sample_id, fill = host_species)) +
  geom_tile(color = "grey") +
  geom_text(aes(label = genome_name)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        legend.title = element_blank()) +
  guides(fill=guide_legend(nrow=3,byrow=TRUE)) +
  labs(x = "Sample name", y = "Best identity for assembled scaffold") +
  coord_flip()

ggsave("results/surveillance_out/viral_sharing.pdf", dpi = 600, width = 5, height = 7.5)

