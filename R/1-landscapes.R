library(terra)
library(gen3sis)

biome_rast_dir <-   file.path(
  "../gen3sis_tests/007-speciation_spatial_variation",
  "out/rast/no_antarctica"
)
r <- rast(file.path(biome_rast_dir, "Time_0_biome_array.tiff"))
mask <- !is.na(r)

if (!dir.exists("./out")) {
  dir.create("./out", )
}
writeRaster(mask, "./out/world_mask.tiff")

# Edit mask in QGIS
mask <- rast("./out/america_mask.tiff")

rast_out_dir <- "./out/biome_rast"
if (!dir.exists(rast_out_dir)) {
  dir.create(rast_out_dir, recursive = TRUE)
}

for (t in 0:5000) {
  biome_file_name <- paste0("Time_", t, "_biome_array.tiff")
  r <- rast(file.path(biome_rast_dir, biome_file_name))
  r_masked <- mask(r, mask, maskvalues = 0)
  r_cropped <- crop(r_masked, c(-180, -22, -90, 90))
  writeRaster(
    r_cropped,
    file.path(rast_out_dir, paste0("Time_", t, "_biome_rast_america.tiff"))
  )
}

## Cost function
dist_km <- function(source, habitable_src, dest, habitable_dest) {
  if (!all(habitable_src, habitable_dest)) {
    2 / 1000
  } else {
    1 / 1000
  }
}

# There is no tectonic plate movement in this simulation, so distances
# do not change between timesteps.  Create distance matrices for one
# timestep only.
rast_list <- list(
  biome = file.path(rast_out_dir, "Time_0_biome_rast_america.tiff")
)

ldsc_outdir <- "./landscapes/america"

create_input_landscape(
  landscapes = rast_list,
  cost_function = dist_km,
  directions = 8,
  output_directory = ldsc_outdir,
  timesteps = "0kya",
  crs = "EPSG:4326",
  calculate_full_distance_matrices = TRUE,
  verbose = TRUE,
  overwrite = TRUE
)

first_timestep <- 5000
last_timestep <- 0
timestep_interval <- 1

timestep_names <- paste0(
  seq(last_timestep, first_timestep, timestep_interval), "kya"
)

full_dist_dir <- file.path(ldsc_outdir, "distances_full")
link_names <- file.path(
  full_dist_dir,
  paste0("distances_full_", 1:(length(timestep_names) - 1), ".rds")
)
if (!all(file.symlink("distances_full_0.rds", link_names))) {
  warning("Some symbolic links could not be created")
}

# Manually create the "landscapes.rds" object (normally this is not
# needed)
rast_list <- list(
  # Sort files by timestep
  biome = gtools::mixedsort(
    list.files(rast_out_dir, "_biome_", full.names = TRUE)
  )
)

# NOTE: Run this only if all biome rasters have been created beforehand
if (length(rast_list$biome) == 5001) {
  timestep_indices <- seq(last_timestep, first_timestep, timestep_interval)
  rast_list$biome <- rast_list$biome[timestep_indices + 1]
}

# Cell coordinates
ldsc_biome <- crds(rast(rast_list$biome[1]), df = TRUE, na.rm = FALSE)
# Important: setting the row names explicitly to assign an ID to each
# cell is needed for gen3sis to work properly.  When the data.frame is
# converted to a matrix the row names are only kept as dimnames when
# they are explicitly set.
row.names(ldsc_biome) <- seq_len(nrow(ldsc_biome))

# The index would be useful in this loop if there were more
# environmental variables
for (i in seq_along(rast_list$biome)) {
  biome_val <- as.data.frame(rast(rast_list$biome[i]), na.rm = FALSE)
  ldsc_biome <- cbind(ldsc_biome, biome_val)
}

names(ldsc_biome)[-(1:2)] <- timestep_names

ldsc <- list(biome = ldsc_biome)
saveRDS(ldsc, file.path(ldsc_outdir, "landscapes.rds"))
