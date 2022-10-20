rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)
require(sf)
require(ggmaps)
require(rnrfa)

df <- fread("data/metadata/sample_collection_metadata.csv", header = T) %>%
  mutate(host_species = gsub("pipustrellus", "pipistrellus", host_species)) %>%
  mutate(host_species = gsub("Pipistellus", "Pipistrellus", host_species)) %>%
  mutate(host_species = gsub("Pipistellus", "Pipistrellus", host_species)) %>%
  mutate(host_species = gsub("Pipistrellys", "Pipistrellus", host_species)) %>%
  mutate(host_species = ifelse(grepl("pygmaeus", host_species), 
                               "Pipistrellus pygmaeus",
                               host_species)) %>%
  mutate(host_species = ifelse(grepl("Nathusius", host_species), 
                               "Pipistrellus nathusii",
                               host_species)) %>%
  mutate(host_species = ifelse(grepl("daubenton", host_species, ignore.case = T), 
                               "Myotis daubentonii",
                               host_species)) %>%
  mutate(host_species = ifelse(grepl("auritus", host_species), 
                               "Plecotus auritus",
                               host_species)) %>%
  mutate(host_species = ifelse(grepl("serotine", host_species), 
                               "Eptesicus serotinus",
                               host_species)) %>%
  mutate(host_species = ifelse(grepl("natteri", host_species), 
                               "Myotis nattereri",
                               host_species)) %>%
  mutate(host_species = ifelse(grepl("sp", host_species), 
                               "Unindentified species",
                               host_species)) %>%
  mutate(host_species = ifelse(grepl("presumed|possible|pipistrellus OR Plecotus|cryptic group", host_species, ignore.case = T), 
                               "Unindentified species",
                               host_species)) %>%
  filter(!grepl("kuhlii", host_species))


fwrite(df, "data/metadata/sample_collection_metadata.parsed.csv")

df %>% 
  group_by(host_species) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

# df_filt <- df %>%
#   filter(!grepl("Beer|back|Colyton|region|area|unknown|Burghfield|Bradley|Chudleigh|Fingringhoe|Sevenoaks|DS16|Buckfastleigh", coord, ignore.case = T)) %>%
#   filter(!(coord %in% c("46608 14387", "48217, 239539", "047539578, 51.51708757")))
#   
# latlong <- df_filt %>% 
#   filter(grepl(",", coord)) %>%
#   separate(coord, into = c("lat", "lon"), sep = ", ")
# 
# others <- df_filt %>% 
  
# Identify erroneous coord
# for (coord in others$coord) {
#   print(coord)
#   print(osg_parse(coord, coord_system = "WGS84"))
# }
# 
# wgs <- osg_parse(others$coord, coord_system = "WGS84")
# wgs_parsed <- tibble(lat = wgs$lat, lon = wgs$lon)
# 
# View(wgs_parsed)
fread("data/metadata/sample_collection_metadata.csv", header = T) %>%
  filter(grepl("noctula", host_species))
