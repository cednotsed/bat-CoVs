rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)
require(randomcoloR)
meta <- fread("data/metadata/sequencing_metadata_260622.csv") %>%
  mutate(sample_id = ifelse(grepl("-", sample_id), sample_id, paste0("Sample-", sample_id)))
dir <- "results/mapping_out/coverage_files"
files <- list.files(dir, full.names = T)

file <- files[1]

morsels <- foreach(file = files) %do% {
  tmp <- fread(file)
  colnames(tmp) <- c("sample_id", "pos", "depth")
  return(tmp)
}

plot_df <- bind_rows(morsels) %>% 
  left_join(meta)

pal <- distinctColorPalette(length(unique(plot_df$genbank_title)))

plot_df %>%
  ggplot(aes(x = pos, y = log(depth, 10), color = genbank_title)) +
  facet_grid(rows = vars(genbank_title)) +
  geom_line() +
  scale_color_manual(values = pal) +
  labs(x = "Position on assembly", y = "log10(read coverage)") +
  theme_bw() +
  theme(legend.position = "none")
  
ggsave("results/mapping_out/assembly_coverage.pdf", dpi = 600, height = 6, width = 6)
  
