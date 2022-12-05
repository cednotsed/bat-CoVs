setwd("../Desktop/git_repos/bat-CoVs/")
require(raster)
require(rasterVis)
require(tidyverse)

# Load probability maps for all species
species <- c("sp1","sp2","sp3","sp4","sp5","sp6",
             "sp7","sp8","sp9","sp10","sp11","sp12",
             "sp13","sp14","sp15","sp16",
             "sp17")
species <- paste0("results/SDM_out/raw_rasters/", species, "_final.img")
species <- setNames(species,
                    c("Barbestella_barbastellus","Myotis_alcathoe","Myotis_bechsteinii",
                      "Myotis_brandtii","Plecotus_auritus","Pipistrellus_pipistrellus",
                      "Myotis_daubentonii","Rhinolophus_ferrumequinum","Plecotus_austriacus",
                      "Nyctalus_leisleri","Rhinolophus_hipposideros","Pipistrellus_nathusii",
                      "Moyits_nattereri","Nyctalus_noctula","Eptesicus_serotinus",
                      "Pipistrellus_pygmaeus","Myotis_mystacinus"))


maps <- Map(raster, species[16:17])
# Make the combined maps
binarise_raster <- function(m) {
  m_temp <- m
  m_temp[m_temp < 0.8] <- 0
  m_temp[m_temp > 0] <- 1
  return(m_temp)
}

m1 <- binarise_raster(maps[[1]])
m2 <- binarise_raster(maps[[2]])
m1[m1 == 1] <- -10
m2[m2 == 1] <- 5

combined <- m1 + m2

# For testing
combined_test <- sampleRegular(combined, size = 100000, asRaster = TRUE)

presence_df <- as.data.frame(combined_test, xy = TRUE, na.rm = TRUE) %>%
  filter(x > 0 & x < 800000, y < 1100000) %>%
  mutate(layer = case_when(layer == 0 ~ "base",
                           layer == -10 ~ names(maps)[[1]],
                           layer == -5 ~ "overlap",
                           layer == 5 ~ names(maps)[[2]]))

# Calculate prop. of overlap
prop_df <- presence_df %>%
  group_by(layer) %>%
  summarise(n = n()) %>%
  filter(layer != "base") %>%
  mutate(prop_total = n / sum(n)) %>%
  filter(layer == "overlap")
prop_overlap <- round(prop_df$prop_total, 3)

p <- ggplot(presence_df,
            aes(x = x, y = y, fill = layer)) +
  geom_tile() +
  geom_text(x = 550000, y = 1000000, label = str_glue("Prop. overlap = {prop_overlap}"),
            size = 2) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_fill_manual(values = c("grey", "#0072B2", "#CC79A7", "#009E73")) 
p
ggsave("results/SDM_out/test_plot.pdf", plot = p, height = 5, width = 5, dpi = 300)





# combined <- ratify(combined)
# rat <- levels(combined)[[1]]
# rat$legend <- c(names(maps)[[1]], "overlap", "base", names(maps[[2]]))
# levels(combined) <- rat
# pdf("results/SDM_out/test_plot.pdf", width = 5, height = 8)
#   levelplot(combined), col = c("grey", "#0072B2", "#CC79A7", "#009E73"))
# dev.off()
# # Make the combined maps
# binarise_raster <- function(m) {
#   m[m < 0.5] <- 0
#   m[m > 0] <- 1
#   return(m)
# }
# m1 <- binarise_raster(maps[[1]])
# m2 <- binarise_raster(maps[[2]])
# 
# combined <- m1 * m2
# combined <- crop(combined, extent(0, 700000, 0, 1100000))
# raster::plot(combined, axes = F)



