rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)
require(randomcoloR)

df <- fread("results/surveillance_out/sample_locations.csv", header = T) %>% 
  rename_all(~tolower(gsub(" ", "_", .x))) %>%
  mutate(collection_date = as.Date(collection_date, format = "%d/%m/%Y")) %>%
  filter(sequencing_id != "")

df_counts <- df %>%
  group_by(host_species) %>%
  summarise(n = n())

plot_df <- df %>%
  left_join(df_counts) %>%
  arrange(desc(n)) %>%
  mutate(host_species_lab = str_glue("{host_species} (n={n})"))

# Color palette
pal <- distinctColorPalette(n_distinct(df$host_species))
pal

# Plot temporal distribution
plot_df %>%
  filter(sequencing_id != "") %>%
  mutate(host_species_lab = factor(host_species_lab, unique(plot_df$host_species_lab))) %>%
  ggplot(aes(x = collection_date, fill = host_species_lab)) +
  geom_density() +
  geom_point(aes(y = -0.005), color = "black", pch = 21) +
  expand_limits(y=-0.01) +
  facet_grid(rows = vars(host_species_lab),
             switch = "y",
             scale = "free") +
  scale_fill_manual(values = pal) +
  theme_bw() + 
  theme(legend.position = "none",
        panel.spacing.x = unit(0, "lines") , 
        panel.spacing.y = unit(0.1, "lines"),
        strip.text.y.left = element_text(angle = 0, face = "italic"),
        strip.background = element_rect(fill = "white"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank()) +
  labs(x = "Sampling date")

  ggsave("results/surveillance_out/temporal_distribution.pdf", dpi = 600, width = 5.5, height = 4.5)

  plot_df %>% summarise(range(collection_date))

