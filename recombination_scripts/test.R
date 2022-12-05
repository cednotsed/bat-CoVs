setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)

cf <- fread("results/recombination_out/clonalframeml_per_branch/Sarbecovirus.trim_to_28261pos.n218.em.txt") %>%
  rename_all(~tolower(gsub(" ", "_", .x))) %>%
  filter(branch != "Mean")

cf_parsed <- cf %>%
  select(parameter, branch, posterior_mean) %>%
  pivot_wider(names_from = "parameter", values_from = "posterior_mean") %>%
  rename_all(~gsub("\\/", "_", .x)) %>%
  rename(one_delta = "1_delta") %>%
  mutate(r_m_cf = R_theta * (1/one_delta) * nu) %>%
  rename(node = branch) %>%
  select(node, r_m_cf)

gub <- fread("results/recombination_out/gubbins_out/Sarbecovirus.trim_to_28261pos.n218.per_branch_statistics.tsv") %>%
  rename_all(~tolower(gsub(" ", "_", .x))) %>%
  rename_all(~gsub("\\/", "_", .x))

gub_parsed <- gub %>% select(node, r_m_gub = r_m)

plot_df <- cf_parsed %>% 
  left_join(gub_parsed) %>%
  filter(!grepl("NODE", node)) %>%
  filter(!is.na(r_m_gub))

median(cf_parsed$r_m_cf)
median(gub_parsed$r_m_gub)

plot_df %>%
pivot_longer(!node, names_to = "method", values_to = "r_m") %>%
  ggplot(aes(x = r_m, y = method, color = method)) +
    geom_density()

plot_df