rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(ape)
require(ggtree)
require(ggutils)
require(ggtreeExtra)

prefix <- "beta_merged.050422"

cm <- read.csv(str_glue("results/mash_out/{prefix}.tsv"), sep = "\t", 
               header = T,
               row.names = 1,
               stringsAsFactors = F)

cm <- data.matrix(cm)
tree <- nj(cm)
outgrp <- fread("data/metadata/Delta_CoVs_refseq_taxid1159901.csv")$Accession
write.tree(tree, str_glue("data/trees/{prefix}.tree"))
rooted <- root(tree, outgroup = outgrp, resolve.root = T)

# Parse GISAID tip labels
gisaid_labs <- rooted$tip.label[grepl("EPI", rooted$tip.label)]
gisaid_labs <- str_split(gisaid_labs, pattern = "\\|", simplify = T)[, 2]
rooted$tip.label[grepl("EPI", rooted$tip.label)] <- gisaid_labs

# Get tree metadata
meta <- read.csv(str_glue("data/metadata/{prefix}.csv"), check.names = F, stringsAsFactors = F)

# No. of genera for annotation
genera_counts <- table(meta$host_genus)
genera_filt <- names(genera_counts[genera_counts > 5])
genera_filt <- genera_filt[genera_filt != ""]
genera_filt <- genera_filt[genera_filt != "Chiroptera"]

# Match metadata to tips
meta.match <- meta[match(rooted$tip.label, meta$accession), ] %>%
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
         tip_labels = ifelse(species == "novel", accession, genbank_title),
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
             "maroon4", "darkorchid4", "slateblue4")

set.seed(66)
random_pal <- distinctColorPalette(length(unique(dd$host_genus_filt)))
random_pal[length(random_pal) - 1] <- "darkolivegreen4"

p_linear <- ggtree(rooted, 
                   size = 0.001,
                   branch.length = "none",
                   color = "darkslategrey",
                   options(ignore.negative.edge = TRUE)) %<+% dd +
  geom_tippoint(aes(hjust = 0.5, color = genus), alpha = 1, size = 1) +
  scale_color_discrete(na.translate = F) +
  geom_tiplab(aes(label = tip_labels), size = 0.5) +
  labs(color = "Viral genus", fill = "Host genus") +
  new_scale_color() +
  geom_tippoint(aes(color = is_new_isolate), alpha = 1, size = 1) +
  scale_color_manual(values = c("black"), na.translate = F, guide = "none") +
  geom_fruit(geom = geom_tile, 
             aes(fill = host_genus_filt),
             color = "grey",
             offset = 0.15,
             width = 10) +
  scale_fill_manual(values = random_pal, 
                    na.translate = F) +
  new_scale_fill() +
  geom_fruit(geom = geom_tile, 
             aes(fill = common_name),
             color = "grey",
             offset = 0.25,
             width = 10) +
  scale_fill_manual(values = col_pal, 
                    na.translate = F) +
  labs(fill = "Host")
  

p_circle <- ggtree(rooted, 
                   layout= "fan",
                   size = 0.001,
                   branch.length = "none",
                   color = "darkslategrey",
                   options(ignore.negative.edge = TRUE)) %<+% dd +
  geom_tippoint(aes(hjust = 0.5, color = genus), alpha = 1, size = 1) +
  scale_color_discrete(na.translate = F) +
  labs(color = "Viral genus", fill = "Host genus") +
  new_scale_color() +
  geom_tippoint(aes(color = is_new_isolate), alpha = 1, size = 1) +
  scale_color_manual(values = c("black"), na.translate = F, guide = "none") +
  geom_fruit(geom = geom_tile, 
             aes(fill = host_genus_filt),
             color = "grey",
             offset = 0.15,
             width = 10) +
  scale_fill_manual(values = random_pal, 
                    na.translate = F) +
  new_scale_fill() +
  geom_fruit(geom = geom_tile, 
             aes(fill = common_name),
             color = "grey",
             offset = 0.25,
             width = 10) +
  scale_fill_manual(values = col_pal, 
                    na.translate = F) +
  labs(fill = "Host")

ggsave("results/beta_circle.png",
       plot = p_circle,
       dpi = 300,
       width = 10,
       height = 10)

ggsave("results/beta_tree.linear.pdf",
       plot = p_linear,
       dpi = 300,
       width = 20,
       height = 20)

# Save annotated tree
annot_tree <- rooted
annot_tree$tip.label <- paste0(dd$accession_labels, "|", dd$tip_labels)
annot_tree$edge.length[annot_tree$edge.length < 0] <- 0
write.tree(annot_tree, "data/trees/beta_merged.050422.rooted.annotated.tree")

meta.match %>% filter(grepl("EPI_ISL", accession))
meta.match
rooted$tip.label[grepl("EPI", rooted$tip.label)]
# meta %>% 
#   filter(host_genus == "Homo") %>%
#   distinct(species)
# 
# meta %>% 
#   filter(host_genus == "Homo", species == "Betacoronavirus 1") %>%
#   select(genbank_title)
# 
# meta %>% filter(common_name == "")
# meta %>% filter(host_genus == "Chiroptera")
#   table(meta$host_genus)[table(meta$host_genus) > 10]
  
