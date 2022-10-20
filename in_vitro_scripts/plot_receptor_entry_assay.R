rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)

# Species receptor entry assay
df <- fread("results/in_vitro_out/species_cell_entry_assay.csv")

plot_df <- df %>%
  pivot_longer(!c(exp_no, receptor), names_to = "spike", values_to = "RLU") %>%
  group_by(receptor, spike) %>%
  summarise(mean_RLU = mean(RLU),
            sd_RLU = sd(RLU)) %>%
  mutate(type = case_when(grepl("GB", spike) ~ "Novel genomes", 
                          grepl("MLV", spike) ~ "Positive control",
                          grepl("Blank", spike) ~ "Negative control",
                          TRUE ~ "Reference CoVs"))

col_pal <- c("darkcyan", "seagreen4", "darkgoldenrod3", 
             "#2e4057", "#d1495b", "cornflowerblue", 
             "maroon4", "darkorchid4", "slateblue4",
             "black", "darkgrey")

col_pal <- col_pal[1:n_distinct(plot_df$spike)]

plot_df %>%
  mutate(spike = factor(spike, c("MLV-A", "MERS-CoV", 
                                 "SARS-CoV-2", "PAGB01", "RHGB02", "Blank")),
         type = factor(type, c("Negative control", "Positive control", 
                               "Reference CoVs", "Novel genomes"))) %>%
  ggplot(aes(x = spike, y = mean_RLU, fill = spike)) +
  geom_bar(stat = "identity") +
  facet_grid(rows = vars(receptor),
             scale = "free_x",
             space = "free_x") +
  geom_errorbar(aes(x = spike,
                    ymin = mean_RLU - sd_RLU,
                    ymax = mean_RLU + sd_RLU),
                width = 0.1,
                size = 0.2,
                colour = "black") +
  scale_fill_manual(values = col_pal) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("results/in_vitro_out/species_entry_pseudovirus_assays.pdf", dpi = 600, width = 4, height = 7)
