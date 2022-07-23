rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

dat <- fread("results/mapping_out/4-126A_CoV2.coverage.txt", sep = "\t")


dat %>% 
  right_join(tibble(V2 = seq(29891))) %>%
  mutate(V3 = replace_na(V3, 0)) %>%
  ggplot(aes(x = V2, y = log(V3, base = 10))) +
    geom_line() +
    labs(title = "4-126A (M. daubetonii) map to WIV04 reference (SARS-CoV-2)",
         y = "log10(coverage)", 
         x = "Pos. on WIV04")
  