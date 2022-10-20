rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)

# Species receptor entry assay
df <- fread("results/in_vitro_out/human_cell_assay.csv")
df2 <- fread("results/in_vitro_out/species_cell_entry_assay.csv")
plot_df <- df %>%
  pivot_longer(!c(spike, exp_no), names_to = "cell_type", values_to = "RLU") %>%
  group_by(cell_type, spike) %>%
  summarise(mean_RLU = mean(RLU),
            sd_RLU = sd(RLU)) %>%
  mutate(type = case_when(grepl("GB", spike) ~ "Novel genomes", 
                          grepl("MLV", spike) ~ "Positive control",
                          grepl("Blank", spike) ~ "Negative control",
                          TRUE ~ "Reference CoVs"),
         assay = "human_entry")

plot_df2 <- df2 %>%
  dplyr::rename(cell_type = receptor) %>%
  pivot_longer(!c(exp_no, cell_type), names_to = "spike", values_to = "RLU") %>%
  group_by(cell_type, spike) %>%
  summarise(mean_RLU = mean(RLU),
            sd_RLU = sd(RLU)) %>%
  mutate(type = case_when(grepl("GB", spike) ~ "Novel genomes", 
                          grepl("MLV", spike) ~ "Positive control",
                          grepl("Blank", spike) ~ "Negative control",
                          TRUE ~ "Reference CoVs"),
         assay = "species_entry") %>%
  bind_rows(plot_df)


col_pal <- c("darkcyan", "seagreen4", "darkgoldenrod3", 
             "#2e4057", "#d1495b", "cornflowerblue", 
             "maroon4", "darkorchid4", "slateblue4",
             "black", "darkgrey")

col_pal <- col_pal[1:n_distinct(plot_df$spike)]

plot_df2 %>%
  filter(assay == "human_entry") %>%
  mutate(spike = factor(spike, c("MLV-A", "MERS-CoV", 
                                 "SARS-CoV-2", "PAGB01", "RHGB02", "Blank")),
         type = factor(type, c("Negative control", "Positive control", 
                               "Reference CoVs", "Novel genomes")),
         cell_type = factor(cell_type, c("293T-hACE2", "Calu-2", "Caco-2",
                                         "Fruitbat ACE2", "Horseshoe bat ACE2", 
                                         "Little brown bat ACE2", "Human ACE2", "Human DPP4",
                                         "Blank"))) %>%
  filter(!(cell_type %in% c("Fruitbat ACE2", "Little brown bat ACE2"))) %>%
  ggplot(aes(x = spike, y = mean_RLU, fill = spike)) +
  geom_bar(stat = "identity") + 
  facet_grid(cols = vars(cell_type),
             scale = "free",
             space = "free_y") +
  geom_errorbar(aes(x = spike, 
                    ymin = mean_RLU - sd_RLU, 
                    ymax = mean_RLU + sd_RLU), 
                width = 0.1, 
                size = 0.2,
                colour = "black") +
  scale_fill_manual(values = col_pal) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank()) +
  labs(x = "Spike", y = "Mean RLU") +
  coord_flip()


ggsave("results/in_vitro_out/human_cell_entry.pdf", dpi = 600, width = 5, height = 2)


plot_df2 %>%
  filter(assay == "species_entry") %>%
  mutate(spike = factor(spike, c("MLV-A", "MERS-CoV", 
                                 "SARS-CoV-2", "PAGB01", "RHGB02", "Blank")),
         type = factor(type, c("Negative control", "Positive control", 
                               "Reference CoVs", "Novel genomes")),
         cell_type = factor(cell_type, c("293T-hACE2", "Calu-2", "Caco-2",
                                         "Fruitbat ACE2", "Horseshoe bat ACE2", 
                                         "Little brown bat ACE2", "Human ACE2", "Human DPP4",
                                         "Blank"))) %>%
  filter(!(cell_type %in% c("Fruitbat ACE2", "Little brown bat ACE2"))) %>%
  ggplot(aes(x = spike, y = mean_RLU, fill = spike)) +
  geom_bar(stat = "identity") + 
  facet_grid(cols = vars(cell_type),
             scale = "free",
             space = "free_y") +
  geom_errorbar(aes(x = spike, 
                    ymin = mean_RLU - sd_RLU, 
                    ymax = mean_RLU + sd_RLU), 
                width = 0.1, 
                size = 0.2,
                colour = "black") +
  scale_fill_manual(values = col_pal) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank()) +
  labs(x = "Spike", y = "Mean RLU") +
  coord_flip()


ggsave("results/in_vitro_out/species_receptor_entry.pdf", dpi = 600, width = 5, height = 2)


plot_df2 %>%
  filter(spike %in% c("RHGB02", "Blank"),
         cell_type == "Human ACE2")
test <- df2 %>%
  filter(receptor %in% c("Human DPP4", "Blank")) %>%
  select(RHGB02, receptor)

anova(lm(test$RHGB02 ~ test$receptor))
df2