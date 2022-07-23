rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

best_hits <- fread("results/blast_out/best_contig_hits_to_reference.csv") %>% 
  as_tibble()
parsed <- fread("results/blast_out/blast_contigs.parsed.csv")

# Get best scaffolds for each sample
match_list <- best_hits %>% 
  distinct(sample_id)

best_sample_match <- foreach(i = seq(nrow(match_list))) %do% {
  best_hits %>%
    filter(sample_id == match_list[i, ]$sample_id) %>%
    arrange(desc(contig_length)) %>%
    head(1)
}

best_match_df <- bind_rows(best_sample_match)

# Get complete genomes
complete <- best_match_df %>% filter(contig_length > 28000)

complete_samples <- unique(complete$sample_id)
for(sample_name in complete_samples) {
  contigs <- readDNAStringSet(str_glue("results/assembly/coronaspades_out/scaffolds/{sample_name}_scaffolds.fasta"))
  contig_name <- complete %>% filter(sample_id == sample_name)
  contig_name <- contig_name$contig_id
  complete_genome <- contigs[contig_name]
  names(complete_genome) <- sample_name
  writeXStringSet(complete_genome, str_glue("results/assembly/constructed_genomes/{sample_name}_assembly.fna"))
}

# Get best candidate scaffold > 10000 as primary
incomplete <- best_match_df %>% filter(contig_length < 28000)
primary_scaffold <- incomplete %>% filter(contig_length > 2000)
to_stitch <- primary_scaffold$sample_id
to_stitch_filt <- to_stitch[to_stitch != "3-32A"]

# Plot scaffolds on reference
plot_morsels <- foreach (sample_name = to_stitch) %do% {
  plot_df <- best_hits %>% 
    filter(sample_id == sample_name) %>% 
    select(contig_id, sstart, send) %>%
    arrange(sstart)
  
  plot_df %>%
    mutate(contig_id = factor(contig_id, unique(plot_df$contig_id))) %>%
    pivot_longer(!contig_id, names_to = "start_end", values_to = "pos") %>%
    ggplot(aes(x = pos, y = contig_id, color = contig_id)) +
    geom_point() +
    geom_line() +
    labs(x = "Aligned position on reference", y = "Scaffold", title = sample_name) +
    theme(legend.position = "none")
}

ggpubr::ggarrange(plotlist = plot_morsels, nrow = length(to_stitch))
ggsave("results/assembly/contig_on_reference_plots/to_stitch.pdf", width = 8, height = 10)

# Contigs to remove
to_remove <- c("NODE_648_length_260_cluster_677_candidate_1_domains_1", 
               "NODE_52_length_233_cluster_81_candidate_1_domains_1", 
               "NODE_42_length_262_cluster_49_candidate_1_domains_1",
               "NODE_27_length_226_cluster_34_candidate_1_domains_1")

# Remove contigs
pool_morsels <- foreach(sample_name = to_stitch_filt) %do% {
  best_hits %>% 
    filter(sample_id == sample_name) %>% 
    select(sample_id, Accession, contig_id, sstart, send, contig_length) %>%
    arrange(sstart)
}

all_contigs <- bind_rows(pool_morsels) %>% 
  filter(!(contig_id %in% to_remove))

# Get contig pools for genomes
for(sample_name in to_stitch_filt) {
  sample_contigs_df <- all_contigs %>%
    filter(sample_id == sample_name) %>%
    mutate(plus = ifelse(sstart < send, T, F)) %>%
    mutate(actual_start = ifelse(plus, sstart, send)) %>%
    arrange(actual_start)
  
  ref_name <- primary_scaffold %>%
    filter(sample_id == sample_name)
  
  ref_name <- ref_name$Accession
  
  ref <- readDNAStringSet(str_glue("data/genomes/coronaviridae_taxid11118_complete_exclude_provirus_040722/coronaviridae.n2089_subset/{ref_name}.fna"))
  contigs <- readDNAStringSet(str_glue("results/assembly/coronaspades_out/scaffolds/{sample_name}_scaffolds.fasta"))
  
  # Forward sequences
  contigs_F <- contigs[names(contigs) %in% (sample_contigs_df %>% filter(plus))$contig_id]
  contigs_R <- contigs[names(contigs) %in% (sample_contigs_df %>% filter(!plus))$contig_id]
  contigs_R_rc <- reverseComplement(contigs_R)
  
  justified_contigs <- c(contigs_F, contigs_R_rc)
  
  # Stitch genomes
  contig_order <- sample_contigs_df$contig_id
  seq_string <- foreach(contig_name = contig_order, .combine = "c") %do% {
    justified_contigs[[contig_name]]
  }

  final_seq <- DNAStringSet(seq_string)
  names(final_seq) <- sample_name
  final_seq
  ref_seq <- readDNAStringSet(str_glue("results/assembly/constructed_genomes/references/{ref_name}.fna"))
  writeXStringSet(c(ref_seq, final_seq), str_glue("results/assembly/constructed_genomes/for_inspection/{sample_name}-{ref_name}.fna"))
  writeXStringSet(final_seq, str_glue("results/assembly/constructed_genomes/{sample_name}_assembly.fna"))
  # writeXStringSet(justified_contigs, str_glue("results/assembly/constructed_genomes/contig_pools/{sample_name}_contigs_{ref_name}.fna"))
  # writeXStringSet(ref, str_glue("results/assembly/constructed_genomes/references/{ref_name}.fna"))
}

# Save assemblies into single file
ass_names <- list.files("results/assembly/constructed_genomes", ".fna", full.names = T)
novel_genomes <- foreach(ass_name = ass_names, .combine = "c") %do% {
  readDNAStringSet(ass_name)
}

writeXStringSet(novel_genomes, "results/assembly/all_novel_genomes.fna")

# Merge table (Rerun after checkv and BLASTn analysis)
meta <- fread("data/metadata/sequencing_metadata_260622.csv") %>%
  mutate(sample_id = ifelse(grepl("-", sample_id), sample_id, paste0("Sample-", sample_id)))

checkv <- fread("results/assembly/checkv_out/all_novel/completeness.tsv") %>%
  left_join(fread("results/assembly/checkv_out/all_novel/quality_summary.tsv"), "contig_id") %>%
  select(sample_id = contig_id, aai_completeness, checkv_quality, warnings) %>%
  mutate(aai_completeness = round(aai_completeness))

blast_meta <- fread("results/blast_out/blast_all_novel_genomes.parsed.csv") %>% 
  select(sample_id, GenBank_Title, pident, Accession, aln_length, bitscore)

blast_morsels <- foreach(sample_name = unique(blast_meta$sample_id)) %do% {
  blast_meta %>% 
    filter(sample_id == sample_name) %>%
    arrange(desc(bitscore)) %>%
    head(1)
}

blast_parsed <- bind_rows(blast_morsels) %>%
 select(-bitscore)

incomplete_parsed <- all_contigs %>% 
  filter(!(contig_id %in% to_remove)) %>%
  group_by(sample_id) %>%
  summarise(genome_length = sum(contig_length)) %>%
  mutate(stitched = T)

final_genome_meta <- complete %>%
  select(sample_id, genome_length = contig_length) %>%
  mutate(stitched = F) %>%
  bind_rows(incomplete_parsed) %>%
  left_join(meta) %>%
  left_join(checkv) %>%
  left_join(blast_parsed) %>%
  arrange(host_species, sample_id) %>%
  mutate(pident = round(pident, 1))

fwrite(final_genome_meta, "results/assembly/assembly_results.csv")
