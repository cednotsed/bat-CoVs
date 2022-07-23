rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)

best_df <- fread("results/blast_out/best_contig_hits_to_reference.csv") %>%
  as_tibble()

genomes_in_samples <- best_df %>% 
  filter(!grepl("MHV|Murine hepatitis", GenBank_Title, ignore.case = T)) %>%
  group_by(sample_id, Accession) %>%
  summarise(total_length = sum(contig_length)) %>%
  filter(total_length > 5000) %>%
  select(-total_length)

plt_morsels <- foreach(i = seq(nrow(genomes_in_samples))) %do% {
  s_id <- genomes_in_samples[i, ]$sample_id
  acc <- genomes_in_samples[i, ]$Accession
  contig_df <- best_df %>% 
    filter(sample_id == s_id,
           Accession == acc)
  contig_df %>%
    mutate(sample_contig = str_glue("{sample_id}-->{GenBank_Title}")) %>%
      select(sample_contig, contig_id, Accession, pident, sstart, send) %>%
      pivot_longer(c(sstart, send), names_to = "start_end", values_to = "coords") %>%
      ggplot(aes(x = coords, y = contig_id, color = pident)) +
      facet_grid(rows = vars(sample_contig)) +
      geom_point() +
      geom_line()
  
  ggsave(str_glue("results/assembly/contig_on_reference_plots/{s_id}_{acc}.png"), 
         dpi = 300,
         height = 5, width = 8)
}
  
best_df %>% 
  filter(sample_id == "Sample-38") %>%
  View()
  