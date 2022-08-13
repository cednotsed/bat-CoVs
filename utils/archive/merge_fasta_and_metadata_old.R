rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)
require(lubridate)
require(ape)

meta_dir <- "data/metadata"
genome_dir <- "data/genomes"

# Merge fasta
alpha_gen <- read.FASTA(str_glue("{genome_dir}/raw_genomes/Alpha_CoVs_taxid693996.050422.filtered.fasta"))
alpha_closest_gen <- read.FASTA(str_glue("{genome_dir}/raw_genomes/Alpha_closest_hits.fasta"))
beta_gen <- read.FASTA(str_glue("{genome_dir}/raw_genomes/Beta_CoVs_taxid694002.050422.filtered.fasta"))
delta_gen <- read.FASTA(str_glue("{genome_dir}/raw_genomes/Delta_CoVs_refseq_taxid1159901.fasta"))
gisaid_gen <- read.FASTA(str_glue("{genome_dir}/raw_genomes/GISAID_050422.filtered.fasta"))
alpha_files <- list.files(str_glue("{genome_dir}/raw_genomes/Alpha"), full.names = T)
alpha_new <- foreach(i = alpha_files, .combine = "c") %do% {
  read.FASTA(i)
}

beta_files <- list.files(str_glue("{genome_dir}/raw_genomes/Beta"), full.names = T)
beta_new <- foreach(i = beta_files, .combine = "c") %do% {
  read.FASTA(i)
}

alpha_gen_merged <- c(alpha_new, alpha_gen, alpha_closest_gen, delta_gen)
beta_gen_merged <- c(beta_new, beta_gen, delta_gen, gisaid_gen)

write.FASTA(alpha_gen_merged, str_glue("{genome_dir}/alpha_merged.050422.fasta"))
write.FASTA(beta_gen_merged, str_glue("{genome_dir}/beta_merged.050422.fasta"))

# Merge metadata
alpha <- fread(str_glue("{meta_dir}/Alpha_CoVs_taxid693996.050422.filtered.csv")) %>%
  rename_all(~ tolower(.)) %>%
  bind_rows(tibble(accession = names(alpha_new), species = "novel"))
alpha_closest <- fread(str_glue("{meta_dir}/Alpha_closest_hits.csv")) %>%
  rename_all(~ tolower(.)) %>%
  select(all_of(colnames(alpha))) %>%
  mutate(genus = "Alphacoronavirus")

beta <- fread(str_glue("{meta_dir}/Beta_CoVs_taxid694002.050422.filtered.csv")) %>%
  rename_all(~ tolower(.)) %>%
  bind_rows(tibble(accession = names(beta_new), species = "novel"))
delta <- fread(str_glue("{meta_dir}/Delta_CoVs_refseq_taxid1159901.csv")) %>%
  rename_all(~ tolower(.))

gisaid <- fread(str_glue("{meta_dir}/GISAID_050422.filtered.tsv")) %>%
  mutate(species = "GISAID CoV",
         genus = "Betacoronavirus",
         collection_date = as_datetime(collection_date, tz = "UTC"),
         submission_date = as_datetime(submission_date, tz = "UTC")) %>%
  select(accession = accession_id, 
         release_date = submission_date,
         genbank_title = virus_name,
         country = loc2, 
         length = sequence_length,
         host,
         genus,
         species)

# Merge and annotate
bats <- c("Rhinolophus", "Myotis", "Scotophilus", 
          "Chiroptera", "Triaenops", "Pipistrellus", 
          "Nyctalus", "Eptesicus", "Miniopterus",
          "Chaerephon", "Hipposideros", "Scotophilus",
          "Cynopterus", "Tylonycteris", "Macronycteris",
          "Desmodus", "Rousettus", "Pteropus",
          "Eonycteris", "Hypsugo", "Aselliscus",
          "Vespertilio", "Laephotis", "Murina")

birds <- c("Zosteropidae", "Passeridae", "Pycnonotus",
           "Lonchura", "Muscicapidae", "Ardeidae",
           "Mareca", "Gallinula", "Turdus")

rodents <- c("Apodemus", "Rattus", "Mus", 
             "Berylmys", "Eothenomys")
mustelids <- c("Mustela", "Neovison")
camelids <- c("Camelus", "Lama", "Vicugna")
civets <- c("Viverridae", "Paradoxurus", "Paguma")
pangolins <- c("Manis", "Pholidota")

alpha_merged <- bind_rows(alpha, alpha_closest, delta) %>%
  separate(host, into = c("host_genus", "host_species"), sep = " ", remove = F) %>%
  mutate(host_genus = ifelse(host_genus == "Feliformia", "Felis", host_genus)) %>%
  mutate(common_name = case_when(host_genus %in% bats ~ "Bats",
                                 grepl("bat|Bat", genbank_title) ~ "Bats",
                                 grepl(paste0(bats, collapse = "|"), genbank_title) ~ "Bats",
                                 host_genus %in% birds ~ "Birds",
                                 host_genus %in% rodents ~ "Rodents",
                                 grepl("Murine", genbank_title) ~ "Rodents",
                                 host_genus %in% mustelids ~ "Mustelids",
                                 host_genus %in% camelids ~ "Camelids",
                                 host_genus %in% civets ~ "Civets",
                                 host_genus %in% pangolins ~ "Pangolins",
                                 host_genus == "Canis" ~ "Dogs",
                                 host_genus == "Felis" ~ "Cats",
                                 host_genus == "Sus" ~ "Pigs",
                                 host_genus == "Homo" ~ "Humans",
                                 TRUE ~ as.character(NA)))

beta_merged <- bind_rows(beta, delta, gisaid) %>%
  separate(host, into = c("host_genus", "host_species"), sep = " ", remove = F) %>%
  mutate(common_name = case_when(host_genus %in% bats ~ "Bats",
                                 grepl("bat|Bat", genbank_title) ~ "Bats",
                                 grepl(paste0(bats, collapse = "|"), genbank_title) ~ "Bats",
                                 host_genus %in% birds ~ "Birds",
                                 host_genus %in% rodents ~ "Rodents",
                                 grepl("Murine", genbank_title) ~ "Rodents",
                                 host_genus %in% mustelids ~ "Mustelids",
                                 host_genus %in% camelids ~ "Camelids",
                                 host_genus %in% civets ~ "Civets",
                                 host_genus %in% pangolins ~ "Pangolins",
                                 host_genus == "Canis" ~ "Dogs",
                                 host_genus == "Felis" ~ "Cats",
                                 host_genus == "Sus" ~ "Pigs",
                                 host_genus == "Homo" ~ "Humans",
                                 TRUE ~ as.character(NA)))

fwrite(alpha_merged, str_glue("{meta_dir}/alpha_merged.050422.csv"))
fwrite(beta_merged, str_glue("{meta_dir}/beta_merged.050422.csv"))

fwrite(alpha_merged %>% select(accession), 
       str_glue("{meta_dir}/alpha_merged.050422.accessions_only.csv"))
fwrite(beta_merged %>% select(accession), 
       str_glue("{meta_dir}/beta_merged.050422.accessions_only.csv"))


total <- bind_rows(alpha_merged, beta_merged) %>% distinct()
unique(total$host_genus)[!(unique(total$host_genus) %in% c(bats, birds, rodents, 
                                                           mustelids, camelids, civets))]
table(alpha_merged$common_name)

table(alpha_merged$genus)
table(beta_merged$genus)
sum(alpha_merged$species == "novel")
sum(beta_merged$species == "novel")

# sum(grepl("EPI_", beta_merged$accession))
