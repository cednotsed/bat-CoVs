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

# # Fix GISAID accessions
# gisaid <- blastout %>% 
#   filter(grepl("EPI_", sseqid)) %>%
#   separate(sseqid, into  = c(NA, "sseqid"), sep = "\\|")
# 
# not_gisaid <- blastout %>% 
#   filter(!grepl("EPI_", sseqid))

# merged <- gisaid %>%
#   bind_rows(not_gisaid)

parsed <- blastout %>% 
  # filter(sseqid %in% c("Sample-25", "Sample-30")) %>%
  mutate(n_matches = aln_length - mismatch,
         primer_len = ifelse(grepl("PC2S2", qseqid), 18, 19)) %>%
  mutate(prop_matches = n_matches / primer_len) %>%
  separate(qseqid, into = c("study"), sep = "_", remove = F) %>%
  mutate(primer_type = ifelse(grepl("_F", qseqid), "F", "R"))

primer_lengths <- tibble(qseqid = names(primers), primer_length = width(primers))

# distinct_pairs <- parsed %>% 
#   distinct(qseqid, sseqid)
# 
# # Get best hits
# no_cores <- detectCores() - 3 
# cl <- makeCluster(no_cores)  
# 
# morsels <- foreach(i = seq(nrow(distinct_pairs)), 
#                    .packages = c("tidyverse")) %dopar% {
#   qseq <- distinct_pairs[i, ]$qseqid
#   sseq <- distinct_pairs[i, ]$sseqid
#   
#   parsed %>% 
#     filter(sseqid == sseq, qseqid == qseq) %>%
#     arrange(desc(prop_matches)) %>%
#     head(1)
# }
# stopCluster(cl)


# annotate tree
tree <- read.tree("data/trees/coronaviridae_NJ_mash.rooted.tree")
# novel <- c("RFGB01", "RHGB01", "RHGB02", "RFGB02", "MDGB01", "PAGB01", "RHGB04", "PPGB02", "MDGB02", "MDGB03")

# tree_filt <- drop.tip(tree, novel)

meta <- fread("data/metadata/all_meta.n2127.080822.csv")

# Match metadata to tips
genera_counts <- table(meta$host_genus)
genera_filt <- names(genera_counts[genera_counts > 5])
genera_filt <- genera_filt[genera_filt != ""]
genera_filt <- genera_filt[genera_filt != "Chiroptera"]

meta.match <- meta[match(tree$tip.label, meta$accession), ] %>%
  separate(accession, into = c("isolate_name"), sep = "_", remove = F) %>%
  mutate(annotation = case_when(grepl("Rhino", host) ~ "Rhinolophus",
                                grepl("Camelus", host) ~ "Camelus",
                                TRUE ~ as.character(NA))) %>%
  mutate(genus = ifelse(genus == "", NA, genus),
         new_isolate = ifelse(species == "novel" |accession == "MW719567.1", 
                              isolate_name, NA),
         is_new_isolate = ifelse(species == "novel", T, NA),
         host_genus_filt = ifelse(host_genus %in% genera_filt, 
                                  host_genus, NA),
         common_name = ifelse(common_name == "", NA, common_name),
         tip_labels = ifelse(species == "novel", isolate_name, genbank_title),
         accession_labels = ifelse(species == "novel", "novel", accession))

all(tree$tip.label == meta.match$accession)

dd <- data.frame(Accession = meta.match$accession,
                 new_isolate = meta.match$new_isolate,
                 is_new_isolate = meta.match$is_new_isolate,
                 host = meta.match$host,
                 host_genus_filt = meta.match$host_genus_filt,
                 common_name = meta.match$common_name,
                 annotation = meta.match$annotation,
                 genus = meta.match$genus,
                 species = meta.match$genbank_title,
                 tip_labels = meta.match$tip_labels,
                 accession_labels = meta.match$accession_labels)

# Plot tree
col_pal <- c("darkcyan", "seagreen4", "darkgoldenrod3", 
             "#2e4057", "#d1495b", "cornflowerblue", 
             "maroon4", "darkorchid4", "slateblue4",
             "black", "darkgrey")

set.seed(66)
random_pal <- distinctColorPalette(length(unique(dd$host_genus_filt)))

plot_df <- parsed %>% 
  # filter(!(sseqid %in% novel)) %>%
  left_join(primer_lengths) %>%
  group_by(study, primer_type, sseqid) %>%
  summarise(min_mismatch = min(mismatch)) %>%
  mutate(primer = str_glue("{study}_{primer_type}")) %>%
  ungroup() %>%
  mutate(accession = gsub("\\.1|\\.2|\\.3|\\.4", "", sseqid)) %>% 
  mutate(accession = factor(accession, levels = dd$Accession))

hm <- plot_df %>%
  ggplot(aes(x = primer, y = accession, fill = as.character(min_mismatch))) +
    geom_tile() +
    scale_y_discrete(drop = F) +
    scale_fill_manual(values = c("skyblue4", "orchid", "indianred")) +
    labs(x = "Primer", y = "Sample", fill = "Prop. matches") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        # axis.text.y = element_text(size = 1),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())

p <- ggtree(tree, color = "grey") %<+% dd +
  geom_tippoint(aes(color = genus)) +  
  scale_color_discrete(na.translate = F) +
  new_scale_color() +
  geom_tippoint(aes(color = is_new_isolate), alpha = 1, size = 2) +
  scale_color_manual(values = c("black"), na.translate = F, guide = "none")

pdf("results/surveillance_out/primer_match/primer_mismatch_against_all_coronaviruses.pdf", width = 11, height = 7.5)
hm %>% insert_left(p)
dev.off()

studies <- unique(plot_df$study)
foreach(study_name = studies) %do% {
  primer_F <- parsed %>%
    filter(study == study_name,
           primer_type == "F") %>% 
    distinct(sseqid)
  
  primer_R <- parsed %>%
    filter(study == study_name,
           primer_type == "R") %>% 
    distinct(sseqid)
  prop_detectable <- length(intersect(primer_F$sseqid, primer_R$sseqid)) / 2127
  str_glue("{study_name}: {prop_detectable * 100}")
  }

# meta <- fread("data/metadata/sequencing_metadata_260622.csv") %>%
#   mutate(sseqid = ifelse(grepl("-", sample_id), sample_id, paste0("Sample-", sample_id))) %>%
#   select(sseqid, host_species, genbank_title)
# 
# bind_rows(morsels) %>%
#   left_join(meta) %>%
#   mutate(primer_type = ifelse(grepl("PC2As", qseqid), "Reverse", "Forward")) %>%
#   ggplot(aes(x = qseqid, y = sseqid, fill = prop_matches)) +
#   geom_tile() +
#   geom_text(aes(label = n_matches)) +
#   facet_grid(rows = vars(host_species), cols = vars(primer_type), scales = "free")


# ggsave("results/surveillance_out/primer_match.all_coronaviruses.png", dpi = 600, height = 8, width = 8)




