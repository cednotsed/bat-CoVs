setwd("../Desktop/git_repos/bat-CoVs/")
require(raster)
require(rasterVis)
require(tidyverse)

# Loads probability maps for species we need
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


maps <- Map(raster, species[c(12, 13, 14, 15, 16)])
# Make the combined maps
binarise_raster <- function(m) {
  m_temp <- m
  m_temp[m_temp < 0.5] <- 0
  m_temp[m_temp > 0] <- 1
  return(m_temp)
}

for (i in seq(length(maps))) {
  maps[[i]] <- binarise_raster(maps[[i]])
}

# Get intersect of all maps
combined <- maps[[1]]

for (i in seq(length(maps))) {
  combined <- combined * maps[[i]]
}

# For testing
combined_test <- sampleRegular(combined, size = 100000, asRaster = TRUE)

presence_df <- as.data.frame(combined_test, xy = TRUE, na.rm = TRUE) %>%
  filter(x > 0 & x < 800000, y < 1100000) %>%
  mutate(layer = ifelse(layer == 1, T, F))

p <- ggplot(presence_df,
            aes(x = x, y = y, fill = layer)) +
  geom_tile() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_fill_manual(values = c("grey", "#CC79A7"))
p
ggsave("results/SDM_out/test_all_overlap.pdf", plot = p, height = 5, width = 5, dpi = 300)





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



