rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)

meta <- fread("data/metadata/all_meta.n2127.080822.csv") %>%
  select(-sample_id)
# isolate_meta <- fread("data/metadata/sequencing_metadata_260622.csv")
blastout <- fread("results/blast_out/blast_contigs.tsv")
# blastout <- fread("results/blast_out/blast_all_novel_genomes.tsv")
colnames(blastout) <- c("qseqid","sseqid","pident","aln_length","mismatch",
                        "gapopen", "qstart", "qend", "sstart", "send",
                        "evalue", "bitscore")
parsed <- blastout %>%
  separate(qseqid, into = c("sample_id", NA, "contig_id"), sep = "\\.") %>%
  mutate(sample_id = gsub("_scaffolds", "", sample_id)) %>%
  mutate(accession = gsub("\\.1|\\.2|\\.3|\\.4", "", sseqid)) %>%
  separate(contig_id, into = c(NA, NA, NA, "contig_length"), sep = "_", remove = F) %>%
  mutate(contig_length = as.numeric(contig_length),
         bitscore = as.numeric(bitscore)) %>% 
  left_join(meta, "accession") %>%
  relocate(genbank_title, .before = 4)

fwrite(parsed, "results/blast_out/blast_contigs.parsed.csv")
# fwrite(parsed, "results/blast_out/blast_all_novel_genomes.parsed.csv")

