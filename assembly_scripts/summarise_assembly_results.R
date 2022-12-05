rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

# Merge table (Rerun after checkv and BLASTn analysis)
meta <- fread("data/metadata/sequencing_metadata_260622.csv") %>%
  mutate(sample_id = ifelse(grepl("-", sample_id), sample_id, paste0("Sample-", sample_id))) %>%
  select(pcr_positive)

ass_results <- fread("results/assembly/assembly_results.csv")

checkv <- fread("results/assembly/checkv_out/all_novel_genomes.n9.fna/completeness.tsv") %>%
  left_join(fread("results/assembly/checkv_out/all_novel_genomes.n9.fna/quality_summary.tsv"), "contig_id") %>%
  select(sample_id = contig_id, aai_completeness, checkv_quality, warnings) %>%
  mutate(aai_completeness = round(aai_completeness))

blast_meta <- fread("results/blast_out/best_contig_hits_to_reference.csv") %>%
  rename_all(~tolower(gsub(" ", "_", .x))) %>%
  mutate(qcov = qend - qstart + 1,
         scov = abs(send - sstart + 1),
         subject_length = NA) %>%
  select(sample_id, genbank_title, pident, accession, bitscore, qcov, scov, subject_length)


all_genomes <- readDNAStringSet("data/genomes/coronaviridae_n2118_novel_n9.080822.fna")


blast_morsels <- foreach(sample_name = unique(blast_meta$sample_id)) %do% {
  temp <- blast_meta %>% 
    filter(sample_id == sample_name) %>%
    arrange(desc(bitscore)) %>%
    head(1)
  
  # Add genome length of subject
  acc <- temp$accession
  acc_length <- width(all_genomes[acc])
  
  temp %>% 
    mutate(subject_length = ifelse(accession == acc, acc_length, subject_length))
}

blast_parsed <- bind_rows(blast_morsels) %>%
  select(-bitscore)

# Include results
complete_results <- ass_results %>%
  left_join(checkv)
  
final <- blast_parsed %>% 
  left_join(complete_results) %>% 
    mutate(pident = round(pident, 1),
           perc_query_aligned = round(qcov / genome_length * 100, 1),
           perc_subj_aligned = round(scov / subject_length * 100, 1)) %>% 
  select(sample_id, genome_name, closest_accession = accession, 
         closest_genome_name = genbank_title,
         is_stitched,
         novel_genome_length = genome_length, 
         pident, perc_N, 
         perc_query_aligned, perc_subj_aligned,
         aai_completeness, checkv_quality)

final %>% View()
# incomplete_parsed <- all_contigs %>% 
#   filter(!(contig_id %in% to_remove)) %>%
#   group_by(sample_id) %>%
#   summarise(genome_length = sum(contig_length)) %>%
#   mutate(stitched = T)
# 
# final_genome_meta <- complete %>%
#   select(sample_id, genome_length = contig_length) %>%
#   mutate(stitched = F) %>%
#   bind_rows(incomplete_parsed) %>%
#   left_join(meta) %>%
#   left_join(checkv) %>%
#   left_join(blast_parsed) %>%
#   arrange(host_species, sample_id) %>%
#   mutate(pident = round(pident, 1),
#          perc_query_aligned = round(qcov / genome_length * 100, 1)) %>%
#   select(-qcov) %>%
#   relocate(genome_length, .after = 12)
# 
# final_genome_meta
fwrite(final, "results/assembly/final_assembly_results.csv")
