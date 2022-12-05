rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)
require(Hmisc)

df <- fread("results/classification_out/abundance_matrix/abundance_matrix.F.tsv") %>%
  mutate(sample_id = gsub(".tsv", "", sample), .before = 1) %>%
  select(-sample, -unclassified) %>%
  rename_all(tolower) %>%
  rename_all(function(x) gsub(" ", "_", x)) %>%
  as_tibble()

meta <- fread("results/surveillance_out/sample_locations.csv") %>%
  select(sample_id = "Sequencing ID", host_species) %>%
  filter(sample_id != "") %>%
  mutate(sample_id = ifelse(grepl("-", sample_id),
                            sample_id, 
                            paste0("Sample-", sample_id)))

# Convert to relative abundance
otu_to_RA <- function(df) {
  mat <- as.matrix(df[, colnames(df) != "sample_id"])
  RA_df <- as.data.frame(mat / rowSums(mat))
  RA_df <- add_column(RA_df, df$sample_id, .before = 1)
  colnames(RA_df)[1] <- "sample_id"
  
  return(RA_df)
}

RA_df <- otu_to_RA(df)

# Convert to presence absence
RA_to_PA <- function(RA_df, PA_threshold) {
  prev_RA <- RA_df %>% column_to_rownames("sample_id")
  prev_RA[prev_RA <= PA_threshold] <- 0
  prev_RA[prev_RA > PA_threshold] <- 1
  prev_RA <- prev_RA %>% rownames_to_column("sample_id")
  return(prev_RA)
}

PA_df <- RA_to_PA(RA_df, 0.001)

# Set absent taxa to zero
RA_df_zeroed <- RA_df %>% 
  select(all_of(colnames(PA_df))) %>%
  column_to_rownames("sample_id")

prev_bool_df <- PA_df %>%
  column_to_rownames("sample_id")

for(i in seq(ncol(RA_df_zeroed))) {
  RA_df_zeroed[!prev_bool_df[, i], i] <- 0
}

RA_df_zeroed <- RA_df_zeroed %>% 
  rownames_to_column("sample_id")

# Remove absent taxa
filter_taxa_by_presence <- function(prev_df, presence_t) {
  prev_temp <- prev_df %>% column_to_rownames("sample_id")
  taxa_counts <- apply(prev_temp, 2, sum)
  to_keep <- names(taxa_counts)[taxa_counts > presence_t]
  return(prev_df %>% select(all_of(c("sample_id", to_keep))))
}

PA_filt <- filter_taxa_by_presence(PA_df, 0)

# Visualise results
dsrna <- c("Amalgaviridae", "Birnaviridae", "Chrysoviridae",
           "Cystoviridae", "Endornaviridae", "Hypoviridae",
           "Megabirnaviridae", "Partitiviridae", "Picobirnaviridae",
           "Reoviridae", "Totiviridae", "Quadriviridae")

plusrna <- c("Arteriviridae", "Coronaviridae", "Mesoniviridae",
             "Roniviridae", "Dicistroviridae", "Iflaviridae",
             "Marnaviridae", "Picornaviridae", "Secoviridae",
             "Alphaflexiviridae", "Betaflexiviridae", "Gammaflexiviridae",
             "Tymoviridae", "Alphatetraviridae", "Alvernaviridae",
             "Astroviridae", "Barnaviridae", "Benyviridae", 
             "Botourmiaviridae", "Bromoviridae", "Caliciviridae",
             "Carmotetraviridae", "Closteroviridae", "Flaviviridae",
             "Fusariviridae", "Hepeviridae", "Hypoviridae",
             "Leviviridae", "Luteoviridae", "Polycipiviridae",
             "Narnaviridae", "Nodaviridae", "Permutotetraviridae",
             "Potyviridae", "Sarthroviridae", "Statovirus",
             "Togaviridae", "Tombusviridae", "Virgaviridae")
negrna <- c("Qinviridae", "Aspiviridae", "Chuviridae", 
            "Bornaviridae", "Filoviridae", "Mymonaviridae",
            "Nyamiviridae", "Paramyxoviridae", "Pneumoviridae",
            "Rhabdoviridae", "Sunviridae", "Yueviridae", 
            "Arenaviridae", "Cruliviridae", "Feraviridae", 
            "Fimoviridae", "Hantaviridae", "Jonviridae",
            "Nairoviridae", "Peribunyaviridae", "Phasmaviridae",
            "Phenuiviridae", "Tospoviridae", "Tilapineviridae",
            "Amnoonviridae", "Orthomyxoviridae")

plot_df <- RA_df_zeroed %>%
  select(all_of(colnames(PA_filt))) %>%
  pivot_longer(!sample_id, names_to = "taxa", values_to = "rel_a") %>%
  left_join(meta) %>%
  group_by(taxa, host_species) %>%
  summarise(freq = sum(rel_a > 0)) %>%
  mutate(freq = ifelse(freq == 0, NA, freq),
         taxa = capitalize(taxa)) %>%
  mutate(viral_type = case_when(taxa %in% dsrna ~ "dsRNA",
                                taxa %in% plusrna ~ "+ssRNA",
                                taxa %in% negrna ~ "-ssRNA")) %>%
  mutate(viral_type = as.factor(viral_type)) %>%
  filter(!is.na(viral_type))

plot_df %>% 
  ggplot(aes(y = host_species, x = taxa, fill = freq)) +
  geom_tile(color = "black") +
  facet_grid(cols = vars(viral_type),
             scales = "free_x",
             space = "free_x") +
  scale_fill_gradient(low = "blue", 
                       high = "red",
                       na.value = "white") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Family", y = "Host species", fill = "Freq.")

ggsave("results/classification_out/RNA_virus_heatmap.311022.pdf", dpi = 600, height = 5, width = 8.5)
