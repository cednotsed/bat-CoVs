rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)
require(lubridate)
require(ape)

meta_dir <- "data/metadata"
genome_dir <- "data/genomes"
novel_dir <- "results/assembly/constructed_genomes"

## Merge metadata ##
cov_meta <- fread("data/metadata/coronaviridae_taxid11118_complete_exclude_provirus_040722/coronaviridae_metadata.n2089_subset.csv") %>%
  rename_all(~ tolower(gsub(" ", "_", .)))
gisaid_meta <- fread("data/metadata/gisaid/Manis_GISAID_080822.tsv") %>%
  bind_rows(fread("data/metadata/gisaid/Rhinolophus_GISAID_080822.tsv")) %>%
  rename_all(~ tolower(gsub(" ", "_", .)))

# Parse gisaid metadata
gisaid_parsed <- gisaid_meta %>%
  separate(virus_name, into = c(NA, NA, NA, "genome_id"), sep = "/", remove = F) %>%
  mutate(genome_id = case_when(grepl("VM22000138_HKUVOC0589P2", virus_name) ~ "VM22000138_HKUVOC0589P2",
                               grepl("Guangdong/1/2019", virus_name) ~ "PCoV-1",
                               TRUE ~ genome_id))

# Get only gisaid genomes not deposited on NCBI
to_keep <- foreach(gen = gisaid_parsed$genome_id, .combine = "c") %do% {
  bool_var <- sum(grepl(gen, cov_meta$GenBank_Title))
  if (bool_var == 0) {
    acc <- deframe(gisaid_parsed %>% 
      filter(genome_id == gen) %>% 
      select(accession_id))
    return(acc)
  } else{
    return(NA)
  }
}

to_keep <- to_keep[!is.na(to_keep)]
gisaid_filt <- gisaid_parsed %>%
  filter(accession_id %in% to_keep) %>%
  select(accession = accession_id, geo_location = location, 
         host, collection_date, genbank_title = virus_name) %>%
  mutate(species = "GISAID CoV", genus = "Betacoronavirus")

# Create novel metadata
novel_list <- list.files(novel_dir, ".fna")
novel_list
novel <- foreach(i = novel_list, .combine = "c") %do% {
  read.FASTA(str_glue("{novel_dir}/{i}"))
}
novel_meta <- tibble(accession = names(novel))
sample_meta <- fread("data/metadata/sequencing_metadata_260622.csv") %>%
  rename(accession = sample_id, host = host_species) %>%
  mutate(accession = ifelse(grepl("-", accession), accession, paste0("Sample-", accession)))

novel_parsed_meta <- novel_meta %>%
  left_join(sample_meta) %>%
  select(-pcr_positive) %>%
  mutate(genus = "novel", species = "novel",
         geo_location = "United Kingdom",
         sample_id = accession) %>%
  mutate(accession = genbank_title,
         genbank_title = str_glue("{sample_id}/{accession}"))

# Parse novel genome names
novel_names <- deframe(tibble(sample_id = names(novel)) %>%
  left_join(novel_parsed_meta) %>%
  select(accession))
names(novel) <- novel_names

# Parse metadata
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

all_meta <- gisaid_filt %>% 
  bind_rows(novel_parsed_meta, cov_meta) %>%
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

fwrite(all_meta, str_glue("{meta_dir}/all_meta.n{nrow(all_meta)}.080822.csv"))

fwrite(all_meta %>% select(accession), 
       str_glue("{meta_dir}/all_meta.n{nrow(all_meta)}.080822.accessions_only.csv"))

## Merge genomes ##
cov_genomes <- read.FASTA(str_glue("{genome_dir}/coronaviridae_taxid11118_complete_exclude_provirus_040722/coronaviridae.n2089_subset.fna"))
cov_names <- names(cov_genomes)
cov_names <- str_split(cov_names, pattern = "\\.", simplify = T)[, 1]
names(cov_genomes) <- cov_names

gisaid_genomes <- c(read.FASTA(str_glue("{genome_dir}/gisaid/Manis_GISAID_080822.fasta")),
  read.FASTA(str_glue("{genome_dir}/gisaid/Rhinolophus_GISAID_080822.fasta")))
acc_names <- names(gisaid_genomes)
acc_names <- str_split(acc_names, pattern = "\\|", simplify = T)[, 2]
names(gisaid_genomes) <- acc_names
gisaid_genomes <- gisaid_genomes[names(gisaid_genomes) %in% gisaid_filt$accession]
  
all_genomes <- c(novel, gisaid_genomes, cov_genomes)

# Check if all genomes/metadata are parsed
all(names(all_genomes) %in% all_meta$accession)

write.FASTA(all_genomes, str_glue("{genome_dir}/coronaviridae_n2118_novel_n9.080822.fna"))
length(all_genomes)
