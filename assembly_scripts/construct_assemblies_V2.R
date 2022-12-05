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

complete_morsels <- foreach(sample_name = complete_samples) %do% {
  ref_df <- complete %>%
    filter(sample_id == sample_name)
  
  ref_name <- ref_df$accession
  ref_genbank <- ref_df$genbank_title
  ref_seq <- readDNAStringSet(str_glue("data/genomes/coronaviridae_taxid11118_complete_exclude_provirus_040722/coronaviridae.n2089_subset/{ref_name}.fna"))

  # Rename genome
  contigs <- readDNAStringSet(str_glue("results/assembly/coronaspades_out/scaffolds/{sample_name}_scaffolds.fasta"))
  contig_name <- complete %>% filter(sample_id == sample_name)
  contig_name <- contig_name$contig_id
  complete_genome <- contigs[contig_name]
  names(complete_genome) <- sample_name
  
  writeXStringSet(complete_genome, str_glue("results/assembly/constructed_genomes/{sample_name}_assembly.fna"))
  writeXStringSet(c(complete_genome, ref_seq), str_glue("results/assembly/constructed_genomes/for_inspection/{sample_name}_{ref_name}.fna"))
  
  # Align concatenated genome to reference (local alignment)
  tibble(sample_id = sample_name, 
         accession = ref_name,
         genbank_title = ref_genbank,
         perc_N = 0, 
         genome_length = width(complete_genome))
}

# Genome metrics
complete_results <- bind_rows(complete_morsels) %>%
  mutate(is_stitched = F)

## For incomplete genomes #####################################################
# Get best candidate scaffold > 10000 as primary
incomplete <- best_match_df %>% filter(contig_length < 28000)
primary_scaffold <- incomplete %>% filter(contig_length > 2000)
to_stitch <- primary_scaffold$sample_id
to_stitch_filt <- to_stitch[to_stitch != "3-32A"]

## PLOT SCAFFOLDS ON CLOSEST REFERENCE ##
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

## REMOVE SUBOPTIMAL CONTIGS BASED ON PLOT ##
to_remove <- c("NODE_648_length_260_cluster_677_candidate_1_domains_1", 
               "NODE_52_length_233_cluster_81_candidate_1_domains_1", 
               "NODE_42_length_262_cluster_49_candidate_1_domains_1",
               "NODE_27_length_226_cluster_34_candidate_1_domains_1")

pool_morsels <- foreach(sample_name = to_stitch_filt) %do% {
  best_hits %>% 
    filter(sample_id == sample_name) %>% 
    select(sample_id, accession, genbank_title, contig_id, sstart, send, qstart, qend, aln_length, contig_length) %>%
    arrange(sstart)
}

all_contigs <- bind_rows(pool_morsels) %>% 
  filter(!(contig_id %in% to_remove))

## STITCH GENOMES ##
assembly_morsels <- foreach(sample_name = to_stitch_filt) %do% {
  sample_contigs_df <- all_contigs %>%
    filter(sample_id == sample_name) %>%
    # Label complemented contigs
    mutate(plus = ifelse(sstart < send, T, F)) %>%
    # Rectify alignment coordinates for ordering contigs
    mutate(rect_start = ifelse(plus, sstart, send),
           rect_end = ifelse(plus, send, sstart)) %>%
    arrange(rect_start)
  
  # Get reference sequence
  ref_df <- primary_scaffold %>%
    filter(sample_id == sample_name)
  
  ref_name <- ref_df$accession
  ref_genbank <- ref_df$genbank_title
  
  ref_seq <- readDNAStringSet(str_glue("data/genomes/coronaviridae_taxid11118_complete_exclude_provirus_040722/coronaviridae.n2089_subset/{ref_name}.fna"))
  
  # Rectify contig sequences
  contigs <- readDNAStringSet(str_glue("results/assembly/coronaspades_out/scaffolds/{sample_name}_scaffolds.fasta"))
  contigs_F <- contigs[names(contigs) %in% (sample_contigs_df %>% filter(plus))$contig_id]
  contigs_R <- contigs[names(contigs) %in% (sample_contigs_df %>% filter(!plus))$contig_id]
  contigs_R_rc <- reverseComplement(contigs_R)
  
  rectified_contigs <- c(contigs_F, contigs_R_rc)
  
  # Concatenate contigs
  contig_order <- sample_contigs_df$contig_id
  seq_string <- foreach(contig_name = contig_order, .combine = "c") %do% {
    rectified_contigs[[contig_name]]
  }
  
  # Align concatenated genome to reference (local alignment)
  palign <- pairwiseAlignment(seq_string, ref_seq, type = "local")
  draft_seq <- alignedPattern(palign)
  
  # Replace gaps with Ns 
  draft_seq_vec <- str_split(draft_seq, "")[[1]]
  final_seq_vec <- gsub("-", "N", draft_seq_vec)
  final_seq <- DNAStringSet(paste0(final_seq_vec, collapse = ""))
  names(final_seq) <- sample_name
  
  # Calculate percentage Ns and percentage identity
  n_gaps <- sum(grepl("-", draft_seq_vec))
  perc_N <- n_gaps / width(final_seq) * 100
  
  ####### FOR 5-129B ONLY ######################################################
  if (sample_name == "5-129B") {
    require(gtools)
    all_perms <- permutations(3, 3, c(0, 1, 2), repeats.allowed = T)
    foreach(i = seq(nrow(all_perms))) %do% {
      rectified_contigs <- rectified_contigs[contig_order]
      perm <- all_perms[i, ]
      gen_129B <- c(rectified_contigs[[1]], DNAString(paste0(rep("N", perm[1]), collapse = "")),
                    rectified_contigs[[2]], DNAString(paste0(rep("N", perm[2]), collapse = "")),
                    rectified_contigs[[3]], DNAString(paste0(rep("N", perm[3]), collapse = "")),
                    rectified_contigs[[4]])
      
      gen_129B <- DNAStringSet(gen_129B)
      
      perm_name <- paste0(perm, collapse = ".")
      writeXStringSet(gen_129B, str_glue("results/assembly/constructed_genomes/5-129B_candidates/5-129B_{perm_name}.fna"))
    }
  }
  ###############################################################################
  
  # Save sequences
  writeXStringSet(c(ref_seq, final_seq), str_glue("results/assembly/constructed_genomes/for_inspection/{sample_name}-{ref_name}.fna"))
  writeXStringSet(final_seq, str_glue("results/assembly/constructed_genomes/{sample_name}_assembly.fna"))
  
  tibble(sample_id = sample_name, 
         accession = ref_name,
         genbank_title = ref_genbank,
         perc_N = perc_N, 
         genome_length = width(final_seq))
}

# Summarise final assembly metrics
assembly_df <- bind_rows(assembly_morsels) %>%
  mutate(is_stitched = T) %>%
  bind_rows(complete_results)

ass_names <- fread("results/assembly/novel_genomes_names.csv") %>%
  dplyr::rename(sample_id = accession, genome_name = genbank_title)

assembly_parsed <- assembly_df %>%
  left_join(ass_names) %>%
  relocate(genome_name, .before = 2)
  
fwrite(assembly_parsed, "results/assembly/assembly_results.csv")

