# this script follows the methods described in Parks 2014 Mapping day-of-burning
# with coarse-resolution satellite fire detection data.

# load packages ####
library(tidyverse)
library(sf)
library(terra)
library(FNN)
library(conflicted)
# set options ####

# load data ####
# HPCC burn perimeter
burn_perimeter <- read_sf("./processed_data/burn_perimeter.shp")

# composite burn index (CBI) raster
cbi <- rast("./processed_data/clipped_burn_raster.tif")

#  fire detections
d <- read_sf(
  "./raw_data/satelite_fire_detection/fire_archive_M-C61_605340.shp"
) %>%
  st_transform(st_crs(burn_perimeter)) %>%
  st_crop(burn_perimeter) %>%
  mutate(ACQ_DATE = ymd(ACQ_DATE), jdate = as.numeric(format(ACQ_DATE, "%j")))

# make interpolated raster ####
## precompute coords and distance matrix ####
# create sparse DOB raster, using CBI raster as a template
dob <- rasterize(d, cbi, field = "jdate")
# coords of each cell
coords <- xyFromCell(dob, 1:ncell(dob))
# coords of all NA cells (no hotspot detected)
na_coords <- coords[is.na(values(dob)), ]
# coords of all non-NA cells
detection_coords <- coords[!is.na(values(dob)), ]
# values of all non-NA cells (julian day of earliest hotspot detection)
dob_values <- values(dob)[!is.na(values(dob))]
# build KNN index
knn_result <- get.knnx(detection_coords, na_coords, k = 5)

## define interpolation function ####
weighted_dob <- function(target_idx) {
  # use the indices to pull the distances of the nearest neighbors
  nn_dist <- knn_result[[2]][target_idx, ]
  # use the indices to pull the values of the nearest neighbors
  days_of_burn <- dob_values[knn_result[[1]][target_idx, ]]
  # calculate "weights" following Parks 2014.
  w_i <- (1 / (abs(days_of_burn - mean(days_of_burn)) + 1)) * nn_dist
  # use weights to calculate the interpolated value.
  # Nate: fyi style guide says to use implicit returns
  round(sum((w_i * days_of_burn)) / sum(w_i), 0)
}

## map function over all NA cells ####
# Nate: this works, but we could code-golf it further with purrr::modify_if().
# I'm not very familiar with the modify() family, so I'm going to stop here.
dob[is.na(values(dob))] <- map(seq_len(nrow(na_coords)),
                               \(x) weighted_dob(x)) %>%
  unlist

## mask result ####
dob_masked <- mask(dob, burn_perimeter)

plot(dob_masked, main = "Ordinal date of burn")
# write out ####
writeRaster(dob_masked,
            filename = "./processed_data/julian_day_of_burn.tiff",
            overwrite = TRUE)
