# this script reads in or downloads raw data, clips, buffers, resamples, and
# reprojects as needed, then writes out tidy .shp or .tif files.

# set new_dl == TRUE if you wish to download new data (where available).
# otherwise, it will read in existing data from the ./raw_data/ directory.

# load packages ####
library(tidyverse)
library(sf)
library(elevatr)
library(terra)

# set options ####
# boost terra memory usage and display progress bar by default
terra::terraOptions(memfrac = 0.8, progress = 1)

# load data ####
## CBI (composite burn index) ####
cbi_full_extent <- rast("./raw_data/burn_severity/ravg_2022_cbi4.tif")
# save spatial reference for use with other data
proj <- crs(cbi_full_extent) # EPSG9822 NAD83 CONUS Albers, unit: meters

## fire perimeters ####
hpcc <- read_sf("./raw_data/fire_perimeters/Perimeters.shp") %>%
  # subset to HPCC fires
  subset(.$poly_Incid %in% c("Calf Canyon", "Hermits Peak")) %>%
  # ensure valid geometries
  st_make_valid %>%
  # reproject
  st_transform(crs = proj)

# now that we have polygons, use them to clip the cbi raster
cbi <- crop(cbi_fullextent, hpcc)

## vegetation treatments ####
veg_treatments <- read_sf("./raw_data/NMVeg.gdb") %>%
  # ensure valid geometries
  st_make_valid %>%
  # reproject
  st_transform(crs = proj) %>%
  # subset by ROI. treatment centroid must be within HPCC burn boundary.
  filter(apply(st_intersects(st_centroid(.), hpcc, sparse = FALSE), 1, any))

# buffer the vegetation treatments by 10m
veg_rast <- st_buffer(veg_treatments, 10) %>%
  # rasterize, using the CBI raster as a template
  rasterize(cbi)

# use veg_rast to create masked CBI raster
# set inverse = TRUE to mask out the vegetation treatment areas *if* you didn't
# already invert the veg_rast raster
masked_cbi <- terra::mask(cbi, veg_rast, inverse = TRUE)

## roads ####
# read in road segments, combine, reproject, and crop to HPCC bounds
roads <- lapply(c("./raw_data/transportation/Trans_RoadSegment_0.shp",
                  "./raw_data/transportation/Trans_RoadSegment_1.shp"),
                st_read) %>%
  # combine road segments
  do.call(rbind, .) %>%
  # reproject to match CBI
  st_transform(crs = proj) %>%
  # ensure valid geometries
  st_make_valid %>%
  # crop to bounds of CBI raster
  st_crop(st_bbox(cbi))

# distance to road raster
roads_distance_rast <- rasterize(roads, cbi) %>%
  distance

## topography ####
# can download from source
new_dl <- FALSE # set to TRUE to download new elevation data
if (new_dl == TRUE) {
  elev <- get_elev_raster(locations = cbi, z = 12) %>%
    # convert RasterLayer to SpatRaster
    rast %>%
    # project to match CBI.
    project(cbi)
} else {
  elev <- rast("./raw_data/HPCC_elevations.tif") %>%
    project(cbi)
}

# calculate terrain metrics
terrain_vars <- c("slope", "aspect", "TRI", "TPI")
terrain_metrics <- terrain(elev, terrain_vars)

## environmental site potential ####
# Nate: do we want to write in an option to download this like we did for the
# DEM? Could use
# https://cran.r-project.org/web//packages//rlandfire/vignettes/rlandfire.html
lf_site_potential <-
  rast("./raw_data/landfire_environmental_site_potential/Tif/us_140esp.tif") %>%
  project(cbi, method = "near")

# TODO check downstream code to make sure factor levels are correct

## climate ####
clim_vars <- c("tmin", "vpdmax", "ppt")
for (i in seq_along(clim_vars)) {
  # read in each climate variable. note: ppt remains in version M4, others at
  # version M5. this is why the filename changes between vars.
  if (i < 3) {
    clim_var <- rast(paste0("./raw_data/prism_climate/PRISM_", clim_vars[i],
                            "_30yr_normal_800mM5_annual_asc.asc")) %>%
      # project to match CBI
      project(cbi)
  } else {
    clim_var <- rast(paste0("./raw_data/prism_climate/PRISM_", clim_vars[i],
                            "_30yr_normal_800mM4_annual_asc.asc")) %>%
      project(cbi)
  }
  # assign to global environment
  assign(clim_vars[i], clim_var)
}

# Write out ####
## CBI ####
# cropped CBI raster
writeRaster(cbi,
            filename = "./processed_data/clipped_burn_raster.tif",
            overwrite = TRUE)
# cropped and masked CBI raster
writeRaster(masked_cbi, "./processed_data/masked_raster.tif", overwrite = TRUE)

## fire perimeters ####
write_sf(hpcc, dsn = "./processed_data/burn_perimeter.shp")

## vegetation treatments ####
st_write(veg_treatments,
         dsn = "./processed_data/vegetation_treatments_hpcc_new.shp",
         append = FALSE)

## roads ####
# cropped and reprojected road segments
st_write(roads, "./processed_data/burn_scar_roads.shp", append = TRUE)
# distance to road raster
writeRaster(roads_distance_rast,
            "./processed_data/distance_to_road.tif",
            overwrite = TRUE)

## topography ####
# elevation
writeRaster(elev, "./processed_data/elev_down.tif", overwrite = TRUE)
map(terrain_vars, function(x) {
  # write each terrain metric to a file
  writeRaster(terrain_metrics[[x]],
              filename = paste0("./processed_data/", x, "_down.tif"),
              overwrite = TRUE)
})

## environmental site potential ####
writeRaster(lf_site_potential,
            "./processed_data/lf_site_potential_new.tif",
            overwrite = TRUE)

## climate ####
map(clim_vars, function(x) {
  # write each climate variable to a file
  writeRaster(get(x),
              filename = paste0("./processed_data/", x, ".tif"),
              overwrite = TRUE)
})
