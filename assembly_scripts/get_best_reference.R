rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

parsed <- fread("results/blast_out/blast_contigs.parsed.csv")

length_df <- parsed %>% 
  filter(!grepl("MHV|Murine", GenBank_Title)) %>%
  group_by(sample_id, Accession, GenBank_Title) %>%
  summarise(total_length = sum(contig_length), sum_bit = sum(bitscore)) %>%
  ungroup()
# View(length_df %>% filter(sample_id == "1-GH087"))
match_list <- length_df %>% distinct(sample_id)

# Get sample_id-->Accession match with largest combined contig length
best_morsels <- foreach(i = seq(nrow(match_list))) %do% {
  sample_name <- match_list[i, ]$sample_id

  length_temp <- length_df %>% 
    filter(sample_id == sample_name)
  max_total <- max(length_temp$total_length)
  max_sum_bit <- max(length_temp$sum_bit)
  length_temp %>% 
    filter(total_length == max_total,
           sum_bit == max_sum_bit)
}

best_df <- bind_rows(best_morsels)
View(best_df)
  # select(sample_id, Accession, GenBank_Title, t)

fwrite(best_df, "results/assembly/best_references_per_sample.csv")

# Get relevant references for genome construction
all_fna <- Biostrings::readDNAStringSet("data/genomes/coronaviridae_taxid11118_complete_exclude_provirus_040722/coronaviridae.n2089_subset.fna")
acc_list <- deframe(best_df %>% distinct(Accession))

# Do mapping of contigs to each sample
relevant_refs <- foreach(acc = acc_list, .combine = "c") %do% {
  all_fna[grepl(acc, names(all_fna))]
}

# Parse ref names
ref_names <- names(relevant_refs)
ref_names <- deframe(tibble(x = ref_names) %>%
                       separate(x, into = c("ref_names"), sep = "\\."))
names(relevant_refs) <- ref_names
Biostrings::writeXStringSet(relevant_refs, 
                            str_glue("data/genomes/coronaviridae_taxid11118_complete_exclude_provirus_040722/relevent_references.n{length(relevant_refs)}.fna"))

## SANITY CHECK ##
# Get all murine/rat accessions for alignment #
# View(length_df %>% distinct(GenBank_Title))
# 
# mice_acc <- deframe(length_df %>% 
#   filter(grepl("Murine|Rat|MHV", GenBank_Title)) %>%
#   distinct(Accession))

