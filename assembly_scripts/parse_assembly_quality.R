setwd("../Desktop/git_repos/bat-CoVs/")
require(data.table)
require(tidyverse)
require(ape)
require(foreach)

dir_list <- list.files("results/assembly/checkv_out/")
dir_list <- dir_list[!(dir_list %in% c("Sample-15", "Sample-39"))]

morsels <- foreach(prefix = dir_list) %do% {
  df <- fread(str_glue("results/assembly/checkv_out/{prefix}/completeness.tsv"))
  qual <- fread(str_glue("results/assembly/checkv_out/{prefix}/quality_summary.tsv"))
  df %>%
    left_join(qual) %>%
    add_column(sample_id = prefix, .before = 1)
}

final <- bind_rows(morsels)

final %>% 
  group_by(sample_id, checkv_quality) %>%
  summarise(n = n())

final_filt <- final %>%
  filter(aai_expected_length > 20000 & aai_expected_length < 35000) %>%
  filter(aai_top_hit != "GCA_000862345.1") %>%
  select(sample_id, contig_length, completeness, 
         gene_count, aai_top_hit, aai_expected_length) %>%
  arrange(sample_id, aai_top_hit)



