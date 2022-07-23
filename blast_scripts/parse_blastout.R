rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)

meta <- fread("data/metadata/coronaviridae_taxid11118_complete_exclude_provirus_040722/coronaviridae_metadata.n2089_subset.csv")
isolate_meta <- fread("data/metadata/sequencing_metadata_260622.csv")
blastout <- fread("results/blast_out/blast_contigs.tsv")
blastout <- fread("results/blast_out/blast_all_novel_genomes.tsv")
colnames(blastout) <- c("qseqid","sseqid","pident","aln_length","mismatch",
                        "gapopen", "qstart", "qend", "sstart", "send",
                        "evalue", "bitscore")
parsed <- blastout %>%
  separate(qseqid, into = c("sample_id", NA, "contig_id"), sep = "\\.") %>%
  mutate(sample_id = gsub("_scaffolds", "", sample_id)) %>%
  separate(sseqid, into = c("Accession"), sep = "\\.") %>%
  separate(contig_id, into = c(NA, NA, NA, "contig_length"), sep = "_", remove = F) %>%
  mutate(contig_length = as.numeric(contig_length),
         bitscore = as.numeric(bitscore)) %>%
  left_join(meta) %>%
  relocate(GenBank_Title, .before = 4)

# fwrite(parsed, "results/blast_out/blast_contigs.parsed.csv")
fwrite(parsed, "results/blast_out/blast_all_novel_genomes.parsed.csv")
