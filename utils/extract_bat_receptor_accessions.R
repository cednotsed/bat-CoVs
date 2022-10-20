rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

uk_bats <- fread("results/UK_breeding_bats_list.csv")
uk_bats %>% separate(species, into = c("genus"), sep = " ") %>%
  distinct(genus)
# Extract accessions for batch entrez
dpp4 <- fread("data/metadata/receptor_genes/DPP4_list_230822.tsv")
ace2 <- fread("data/metadata/receptor_genes/ACE2_list_230822.tsv")
apn <- fread("data/metadata/receptor_genes/APN_list_290822.tsv")

bats <- c("Kerivoula", "Phoniscus", "Eudiscopus",
          "Myotis", "Submyotodon", "Harpiocephalus",
          "Harpiola", "Murina", "Antrozous",
          "Bauerus", "Rhogeessa", "Arielulus",
          "Eptesicus", "Glauconycteris", "Hesperoptenus",
          "Histiotus", "Ia", "Lasionycteris",
          "Scoteanax", "Scotomanes", "Scotorepens",
          "Thainycteris", "Aeorestes" , "Dasypterus", 
          "Lasiurus", "Nycticeius", "Parastrellus",
          "Perimyotis", "Glischropus", "Nyctalus",
          "Pipistrellus", "Scotoecus", "Scotozous", 
          "Vansonia", "Barbastella", "Corynorhinus",
          "Euderma", "Idionycteris", "Otonycteris",
          "Plecotus", "Scotophilus", "Afronycteris",
          "Afronycteris", "Cassistrellus", "Chalinolobus",
          "Falsistrellus", "Hypsugo", "Mirostrellus",
          "Neoromicia", "Nycticeinops", "Nyctophilus",
          "Pharotis", "Philetor", "Tylonycteris", 
          "Vespadelus", "Vespertilio", "Rhinolophus")

bat_regex <- paste0(bats, collapse = "|")

bat_dpp4 <- dpp4 %>% 
  filter(grepl(bat_regex, Org_name))

fwrite(bat_dpp4 %>% select(Org_name, genomic_nucleotide_accession.version),
       "data/metadata/receptor_genes/NCBI_bats_with_DPP4.csv", col.names = F)
fwrite(bat_dpp4 %>% select(genomic_nucleotide_accession.version),
       "data/metadata/receptor_genes/NCBI_bats_with_DPP4.accessions_only.txt", col.names = F)

bat_ace2 <- ace2 %>% 
  filter(grepl(bat_regex, Org_name))

fwrite(bat_ace2 %>% select(Org_name, genomic_nucleotide_accession.version),
       "data/metadata/receptor_genes/NCBI_bats_with_ACE2.csv", col.names = F)
fwrite(bat_ace2 %>% select(genomic_nucleotide_accession.version),
       "data/metadata/receptor_genes/NCBI_bats_with_ACE2.accessions_only.txt", col.names = F)

bat_apn <- apn %>% 
  filter(grepl(bat_regex, Org_name))

fwrite(bat_apn %>% select(Org_name, genomic_nucleotide_accession.version),
       "data/metadata/receptor_genes/NCBI_bats_with_APN.csv", col.names = F)
fwrite(bat_apn %>% select(genomic_nucleotide_accession.version),
       "data/metadata/receptor_genes/NCBI_bats_with_APN.accessions_only.txt", col.names = F)



# # Extract CDS from whole genomes
# ace2_fna <- readDNAStringSet("data/genomes/UK_bats/receptor_genes/bat_genomes_with_ACE2.fna")
# ace2_seqs <- ace2_fna[grepl("ACE2", names(ace2_fna), ignore.case = T)]
# writeXStringSet(ace2_seqs, "data/genomes/UK_bats/receptor_genes/bat_ACE2.fna")
# 
# dpp4_fna <- readDNAStringSet("data/genomes/UK_bats/receptor_genes/bat_genomes_with_DPP4.fna")
# dpp4_seqs <- dpp4_fna[grepl("dpp4", names(dpp4_fna), ignore.case = T)]
# writeXStringSet(dpp4_seqs, "data/genomes/UK_bats/receptor_genes/bat_DPP4.fna")

# test <- names(dpp4_seqs)
# test <- str_split(test, "\\.", simplify = T)[, 1]
# test <- gsub("lcl\\|", "", test)
# bat_dpp4 %>% 
#   separate(genomic_nucleotide_accession.version, into = c("acc"), sep = "\\.") %>%
#   select(Org_name, acc) %>%
  
