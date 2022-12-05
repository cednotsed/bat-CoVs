setwd("../Desktop/git_repos/bat-CoVs/")
require(raster)

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


maps <- Map(raster, species[16:17])
maps
plot(maps[[1]])
# Make the combined maps
m1 <- maps[[1]]
m2 <- maps[[2]]
m2[m2 < 0.5] <- 0
m2[m2 > 0] <- 1

plot(m2)


combined <- maps[["Myotis_alcathoe"]] * maps[["Myotis_alcathoe"]]
combined

plot(maps[["Myotis_bechsteinii"]])
# Load the HPD data
# - This is explicitly projected onto the British National Grid
# - It is at 1000m resolution, so needs to be upsampled.
hpd <- raster('human-pop/data/UK_residential_population_2011_1_km.asc')
proj4string(hpd)

# You could _interpolate_ population here (method  = 'bilinear') but
# I've currently just duplicated the 1000m cell value into 10x10 100 metre cells

hpd100 <- disaggregate(hpd, fact = 10, method = '')

# The two maps have different extents, so crop combined to the smaller 
# extent of hpd100 
extent(hpd100)
extent(combined)
combined_crop <- crop(combined, extent(hpd100))

# Force set the CRS of combined crop
crs(combined_crop) <- crs(hpd100)

# Get the product
bat_times_hpd <- combined_crop * hpd100

plot(bat_times_hpd)

# Save the files
writeRaster(combined_crop, "bats_infectious_raw.img", overwrite=TRUE)
writeRaster(combined_crop/5, "bats_infectious_rescaled.img", overwrite=TRUE)


writeRaster(bat_times_hpd, "~/bats_infectious_human.img", overwrite=TRUE)
writeRaster(bat_times_hpd / 5, "~/bats_infectious_human.img", overwrite=TRUE)

