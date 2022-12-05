rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)
require(doParallel)
require(ape)
require(aplot)
require(ggtree)
require(ggnewscale)

primers <- readDNAStringSet("data/genomes/all_primers.fna")
blastout <- fread("results/blast_out/blast_primers_against_all_genomes.tsv") %>%
  as_tibble()

colnames(blastout) <- c("qseqid","sseqid","pident","aln_length","mismatch",
                        "gapopen", "qstart", "qend", "sstart", "send",
                        "evalue", "bitscore")

parsed <- blastout %>% 
  # filter(sseqid %in% c("Sample-25", "Sample-30")) %>%
  mutate(n_matches = aln_length - mismatch,
         primer_len = ifelse(grepl("PC2S2", qseqid), 18, 19)) %>%
  mutate(prop_matches = n_matches / primer_len) %>%
  separate(qseqid, into = c("study"), sep = "_", remove = F) %>%
  mutate(primer_type = ifelse(grepl("_F", qseqid), "F", "R"))

primer_lengths <- tibble(qseqid = names(primers), primer_length = width(primers))

# annotate tree
tree <- read.tree("data/trees/coronaviridae_NJ_mash.rooted.tree")
novel <- c("1-GH087", "2-30B",
           "2-GH106", "4-126A",
           "5-129B", "Sample-18",
           "Sample-25", "Sample-30",
           "Sample-37")

# Retrieve novel genomes
to_drop <- tree$tip.label[!(tree$tip.label %in% novel)]
tree_filt <- drop.tip(tree, to_drop)

meta <- fread("data/metadata/all_meta.n2127.080822.csv")
meta2 <- fread("results/assembly/novel_genomes_names.csv")
meta2
meta.match <- tibble(accession = tree_filt$tip.label) %>%
  left_join(meta) %>%
  select(-genbank_title) %>%
  left_join(meta2) %>%
  mutate(accession = genbank_title) %>%
  mutate(subgenus = case_when(grepl("Rhino", host) ~"Sarbecovirus",
                              grepl("Plecot", host) ~"Merbecovirus",
                              grepl("Pipi|Myotis", host) ~"Pedacovirus"),
         genus = case_when(grepl("Rhino", host) ~"Betacoronavirus",
                           grepl("Plecot", host) ~"Betacoronavirus",
                           grepl("Pipi|Myotis", host) ~"Alphacoronavirus")) %>%
  select(accession, host, genus, subgenus)

dd <- data.frame(Accession = meta.match$accession,
                 host = meta.match$host,
                 genus = meta.match$genus,
                 subgenus = meta.match$subgenus) 

# Rename tip labels
tree_filt$tip.label <- meta.match$accession

# Plot tree
col_pal <- c("darkcyan", "seagreen4", "darkgoldenrod3", 
             "#2e4057", "#d1495b", "cornflowerblue", 
             "maroon4", "darkorchid4", "slateblue4",
             "black", "darkgrey")

plot_df <- parsed %>% 
  left_join(primer_lengths) %>%
  mutate(primer = str_glue("{study}_{primer_type}")) %>%
  mutate(primer = factor(primer)) %>%
  filter(sseqid %in% novel) %>%
  group_by(study, primer_type, sseqid, primer) %>%
  summarise(min_mismatch = min(mismatch)) %>%
  ungroup() %>%
  mutate(accession = sseqid) %>%
  left_join(meta2) %>%
  mutate(accession = factor(genbank_title, levels = meta.match$accession))

hm <- plot_df %>%
  ggplot(aes(x = primer, y = accession, fill = as.character(min_mismatch))) +
  geom_tile(color = "black") +
  scale_x_discrete(drop = F) +
  scale_y_discrete(drop = F) +
  scale_fill_manual(values = c("skyblue4", "orchid", "indianred")) +
  labs(x = "Primer", y = "Sample", fill = "Prop. matches") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        # axis.text.y = element_text(size = 10),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())
hm

p <- ggtree(tree_filt, color = "grey") %<+% dd +
  geom_tippoint(aes(fill = subgenus), color = "black", pch = 21) +  
  scale_fill_manual(values = col_pal[c(1, 3, 8)]) +
  geom_tiplab() +
  new_scale_fill() + 
  geom_fruit(geom = geom_tile, 
             aes(fill = genus),
             color = "black",
             offset = 0.5,
             width = 0.01)

pdf("results/surveillance_out/primer_match/primer_mismatch_against_novel.pdf", width = 11, height = 7.5)
hm %>% insert_left(p)
dev.off()

