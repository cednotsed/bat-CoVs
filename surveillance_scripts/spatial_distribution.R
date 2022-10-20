rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)
require(randomcoloR)
require(sf)
require(raster)
require(ggrepel)

# Get GADM UK regions
uk <- st_as_sf(getData("GADM", country="GBR", level=1)) %>%
  st_set_crs(4326)

# Get sf from city information
# From https://www.latlong.net/category/cities-235-15.html
cities <- data.frame(
  city = c("London", "Brighton", "Plymouth",
           "Birmingham", "Leeds", "Glasgow", 
           "Sheffield", "Bradford", "Manchester", 
           "Edinburgh", "Liverpool", "Cardiff",
           "Belfast"),
  lat = c(51.509865, 50.827778, 50.376289,
          52.489471, 53.801277, 55.860916,
          53.383331, 53.799999, 53.483959,
          55.953251, 53.400002, 51.481583,
          54.607868),
  lon = c(-0.118092, -0.152778, -4.143841,
           -1.898575, -1.548567, -4.251433,
           -1.466667, -1.750000, -2.244644,
           -3.188267, -2.983333, -3.179090,
           -5.926437))  

cities <- st_as_sf(cities, coords = c("lon", "lat"),
         crs = 4326)

cities <- cities %>%
  bind_cols(st_coordinates(cities))

# Get bat collection locations
loc_dat <- fread("results/surveillance_out/sample_locations.csv", header = T) %>%
  separate(coords, into = c("lat", "lon"), sep = ", ") %>%
  mutate(lat = as.numeric(lat),
         long = as.numeric(lon)) %>% 
  filter(!is.na(lat)) %>%
  select(sample_id, lat, lon)

loc_dat <- st_as_sf(loc_dat, coords = c("lon", "lat"),
                   crs = 4326)

loc_dat <- loc_dat %>%
  bind_cols(st_coordinates(loc_dat))

ggplot() + 
  geom_sf(data = uk, aes(fill = NAME_1)) +
  geom_sf(data = cities) +
  geom_sf(data = loc_dat, 
          pch = 24, 
          fill = "red", 
          color = "black", 
          size = 1) +
  geom_text_repel(data = cities, 
                  aes(x = X, y = Y, label = city),
                  size = 3) +
  coord_sf(xlim = c(-8, 2), 
           ylim = c(50, 60.5), 
           crs = 4326) +
  theme_bw() +
  scale_fill_manual(values = c("paleturquoise2", "darksalmon", "darkseagreen4", "tan")) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank()) +
  labs(fill = "Region")

ggsave("results/surveillance_out/spatial_distribution.pdf", 
       dpi = 600,
       width = 4,
       height = 5)

