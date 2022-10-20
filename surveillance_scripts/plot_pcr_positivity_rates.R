rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)
require(randomcoloR)

meta <- fread("data/metadata/all_meta.n2127.080822.csv")

df <- fread("results/surveillance_out/pcr_surveillance_results.csv")

plot_df <- df %>%
  mutate(prev = n_positive / n_collected * 100,
         perc_seq = n_sequenced / n_collected * 100) %>%
  separate(species, into = c("genus", "species"), sep = " ", remove = F) %>%
  mutate(genus_short = substr(genus, 1, 1)) %>%
  mutate(species = str_glue("{genus_short}. {species}")) %>%
  arrange(desc(prev))

pal <- c("darkcyan", "seagreen4", "darkgoldenrod3", 
             "#2e4057", "#d1495b", "cornflowerblue", 
             "maroon4", "darkorchid4", "slateblue4",
             "black", "darkgrey")
pal <- col_pal[1:n_distinct(plot_df$species)]

plot_df %>%
  mutate(species = factor(species, unique(plot_df$species)),
         genus = factor(genus, unique(plot_df$genus))) %>%
  ggplot(aes(x = species, y = prev, fill = genus)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = pal) +
  # geom_text(aes(x = species,
  #               y = -5,
  #               label = str_glue("n = {n_collected}"))) +
  geom_text(aes(x = species,
                y = prev,
                label = n_positive),
            color = "orchid4",
            vjust = 0) +
  theme_bw() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 30, 
                                   hjust = 1,
                                   face = "italic"),
        legend.text = element_text(face = "italic")) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = "UK bat species", 
       y = "RT-qPCR positivity (%)",
       fill = "Host genus") +
  ylim(0, 120)

ggsave("results/surveillance_out/pcr_positivity_rate.pdf", 
       dpi = 600,
       width = 11,
       height = 3)
  


