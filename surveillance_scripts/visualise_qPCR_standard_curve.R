rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

df <- fread("results/surveillance_out/standard_curve_trial_4_Superscript_091122.csv", header = T) %>%
  rename_with(~ paste0("col_", .x))
df
plot_df <- df %>%
  select(all_of(c("col_V1", "col_V2", paste0("col_", seq(7))))) %>% 
  pivot_longer(!c(col_V1, col_V2), names_to = "standard", values_to = "CT") %>%
  filter(col_V1 %in% c("A", "B", "C"), col_V2 == "Cq") %>%
  mutate(standard = case_when(standard == "col_1" ~ 6.83 / 1, 
                              standard == "col_2" ~ 6.83 / 2,
                              standard == "col_3" ~ 6.83 / 4,
                              standard == "col_4" ~ 6.83 / 8,
                              standard == "col_5" ~ 6.83 / 16,
                              standard == "col_7" ~ 6.83 / 1,
                              standard == "col_6" ~ 0)) %>%
  # filter(standard != 0) %>%
  mutate(CT = 2 ^ as.numeric(CT))


plot_df
cor(as.numeric(plot_df$standard), as.numeric(plot_df$CT))
plot_df %>%
  filter(standard != 0) %>%
ggplot(aes(x = standard, y = as.numeric(CT), color = col_V1)) +
  geom_point()
    
