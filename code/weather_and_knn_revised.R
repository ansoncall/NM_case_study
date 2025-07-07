# this script prepares weather data and identifies candidate control plots for
# the burn severity analysis. it (optionally) downloads gridmet data, processes
# it, and extracts relevant features for each day of the burn. It also creates a
# grid of control plots based on the burn perimeter and extracts environmental
# features for these plots. this script should take less than an hour to run on
# most laptops - ~32 minutes on a 2022ish laptop with 32 GB ram.

# load packages ####
library(tidyverse)
library(downloader)
library(sf)
library(terra)
library(exactextractr)

# set options ####
# boost terra memory usage and never display progress bar
terraOptions(memfrac = 0.8, progress = 0)

# load data ####
## composite burn index (CBI) ####
# Nate: you were using the full cbi raster (covering all of north america) here?
cbi <- rast("./processed_data/clipped_burn_raster.tif") %>%
  # 9 is the code for "unmappable". set to NA for downstream analysis.
  subst(9, NA) %>%
  # 0 == outside the fire perimeter and 1 == inside, but unburned. here, we
  # consider these to be the same.
  subst(0, 1)

## day of burn (DOB) ####
dob <- rast("./processed_data/julian_day_of_burn.tiff")
# get unique dob values, sort them, and convert to Date objects
dates <- unique(dob) %>%
  pull(1) %>%
  sort

## HPCC burn perimeter ####
burn_perimeter <- read_sf("./processed_data/burn_perimeter.shp")

## weather ####
# define vector of weather variables
weather_vars <- c("rmax", "vs", "th", "fm100", "fm1000", "tmmx")
# optional: download new weather data to disk. set new_dl == TRUE to download.
new_dl <- FALSE
if (new_dl == TRUE) {
  map(weather_vars,
      \(x) {
        download(
          url = paste0("http://www.northwestknowledge.net/metdata/data/",
                       x,
                       "_2022.nc"),
          destfile = paste0("raw_data/", x, "_2022.nc"),
          mode = "wb"
        )
      })
}
# load the weather data from disk
weather_rasts <- map(
  weather_vars,
  \(x) {
    # read in rasters, subset to get the bands for the relevant days only
    rast(paste0("./raw_data/", x, "_2022.nc"), lyrs = dates)
  }
) %>%
  # stack all bands: 6 vars * 71 days = 426 bands
  rast %>%
  project(cbi)
# tidy layer names and set time value to dob
names(weather_rasts) <- rep(weather_vars, each = length(dates))
time(weather_rasts) <- rep(dates, 6)
# map over the rasters and mask them by the appropriate date subset of the dob
# raster
weather_masked <- map(
  dates,
  # i == the ordinal day of burn
  \(i) {
    # create dob mask
    msk <- ifel(dob == i, 1, NA)
    # select the raster layers with the target i
    weather_rasts_subset <- weather_rasts[[time(weather_rasts) == i]]
    # mask them
    masked_weather_rasts <- mask(weather_rasts[[time(weather_rasts) == i]], msk)
  },
  .progress = "Masking weather rasters..."
) %>%
  # stack result into single multiband raster
  do.call(c, .)
# composite the partial rasters for each output variable
weather_composites <- map(
  weather_vars,
  \(x) {
    # select all the rasters for the variable x
    var_subset <- weather_masked[[grep(x, names(weather_masked))]]
    # use any function to summarize layers (there is only one non-NA value per
    # pixel and layer). median() is reasonably fast.
    out <- app(var_subset, median, na.rm = TRUE, wopt = list(names = x))
  },
  .progress = "Building composites..."
) %>%
  do.call(c, .)
# view composite rasters. Nate: check
weather_composites <-  rast("./processed_data/weather_composite.tif")
plot(weather_composites, main = weather_vars)

# control plots ####
## create geometries ####

# Create a grid of control plots. Control plots will be 10 acre circles to match
# the validation plots. They will be built on a square grid.The area of square
# grid cell that that contains a 10 acre is 12.732 acres (confirmed. -Anson)

# total burned area in acres
# Nate: shouldn't we be buffering the burn perimeter so we don't end up with
# circular plots hanging over the edge of the burn perimeter? I didn't see you
# doing this in this script. Probably more efficient to buffer up front so you
# don't have to run an extra geoprocessing operation on thousands of polygons
# later on.
total_burn_area <- sum(st_area(st_buffer(burn_perimeter,
                                         dist = -113.49694))) / 4046.68

# the area in 12.732 acre units--the maximum possible number of control plots,
# not accounting for square-packing in non-square areas.
how_may_cells <- as.numeric(total_burn_area) / 12.732

# sample points in the burn perimeter to create a regular grid of control plot
# centers. size parameter is used to specify an *approximate* number of points
# to return.
gridded_points <- st_sample(st_buffer(burn_perimeter, dist = -113.49694),
                            size = round(how_may_cells, 0),
                            type = "regular")
# TODO maybe get rid of plots that fall within treatments at this stage to speed
# up downstream stuff

# buffer points by 113.49... meters to create circular 10 acre polygons.
# (pi*113.49694^2)/4046.86 == 10 acres. (checked. -Anson)
gridded_plots <- st_buffer(gridded_points, dist = 113.49694) %>% st_as_sf

# Nate: check
# mapview::mapview(burn_perimeter) +
#   mapview::mapview(gridded_plots, color = "red", alpha = 0.5) # nolint

## load additional data ####
# topography
topo_vars <- c("elev", "aspect", "TRI", "TPI", "slope")
topo_rasts <- map(
  topo_vars,
  \(x) {
    # load rasters, project to cbi crs, and rename layers
    rast(paste0("./processed_data/", x, "_down.tif")) %>%
      project(cbi) %>%
      setNames(x)
  }
) %>%
  do.call(c, .)

# roads, LandFire site potential
roads_distance <- rast("./processed_data/distance_to_road.tif")
site_potential <- rast("./processed_data/lf_site_potential_new.tif")
# set active category of LF site potential to 2, which is the fine-scale
# zone*esp*esplf categorization
activeCat(site_potential) <- 2

# compile full multiband raster
# all layers, excluding the categorical site_potential layer
all_layers <- c(
  topo_rasts, roads_distance, weather_composites, cbi
) %>%
  # set names for the layers
  setNames(c(topo_vars, "road", weather_vars, "cbi"))

# Nate: check
plot(all_layers, main = names(all_layers))
plot(site_potential, main = "esp")

## extract values ####

# Nate: previously, esp values were extracted in the same way as the other
# continuous values, i.e. using fun = "mean". This seems inappropriate for
# factor values. Here, we take the mode to extract the most common cover type
# within each control plot. This loses some information about the secondary
# cover types and the proportions of each cover type, but it maintains
# compatibility with the downstream analysis and makes more sense than averaging
# the numerical values, which are actually codes for qualitative cover type
# descriptions.

# # aside: it seems like normalizeWeights is unnecessary; see examples in help.
# # Also, it was breaking the code and idk why. Says normalizeWeights argument
# # was unused even though I was specifying the raster namespace. See:
# nolint start
# test <- all_layers["cbi"]
# # broken:
# raster::extract(test, st_as_sf(gridded_plots), fun = mean, na.rm = TRUE,
#                 weights = TRUE, exact = TRUE, normalizeWeights = TRUE,
#                 small = TRUE)
# # works (took me 40 minutes to run):
# tictoc::tic()
# testout <- raster::extract(
#   test, st_as_sf(gridded_plots), fun = mean, na.rm = TRUE, weights = TRUE,
#   exact = TRUE, small = TRUE
# )
# tictoc::toc()
# nolint end
# aside 2: switching to exactextractr::exact_extract for performance. I spot
# checked with CBI and this will return almost exactly the same result you had
# previously.

## extract modal esp ####
# TODO: remove esp stuff from gridded plots.
# esp is not used in propensity score matching, so we dont even need to do this.
gridded_plots$esp <- exact_extract(site_potential,
                                   gridded_plots,
                                   fun = "mode",
                                   weights = "area",
                                   progress = TRUE)

## extract mean aspect ####
# Nate: aspect is in degrees. This is going to create averaging problems when we
# extract to the control plots. We should convert to "northness" and "eastness"
# before averaging, then convert back to degrees after averaging.
plot(all_layers["aspect"], main = "aspect")
### convert aspect to uv component vectors ####
# define function to convert aspect to uv component vectors
aspect_to_uv <- function(aspect_degrees) {
  # Convert to radians
  aspect_radians <- aspect_degrees * pi / 180
  # Compute x and y components
  u <- cos(aspect_radians)
  v <- sin(aspect_radians)
  c(u, v)
}
# apply uv function to aspect layer
aspect_uv <- all_layers["aspect"] %>%
  app(aspect_to_uv, wopt = list(names = c("aspect_u", "aspect_v")))
plot(aspect_uv, main = names(aspect_uv))

# calculate means of u and v
uv <- exact_extract(aspect_uv,
                    gridded_plots,
                    fun = "mean",
                    weights = "area",
                    stack_apply = TRUE,
                    force_df = TRUE,
                    progress = TRUE)
### convert mean uv back to bearings ####
# radians, counterclockwise from x axis (east)
angles_rad <- atan2(uv$mean.aspect_v, uv$mean.aspect_u)
# convert to degrees, clockwise from north. assign to gridded_plots.
gridded_plots$aspect <- (450 - angles_rad * 180 / pi) %% 360

# extract means for other variables ####
var_means <- exact_extract(all_layers[[-2]], # exclude aspect
                           gridded_plots,
                           fun = "mean",
                           weights = "area",
                           stack_apply = TRUE,
                           colname_fun = \(values, ...) values,
                           progress = TRUE)
gridded_plots <- cbind(gridded_plots, var_means)

# write out ####
## weather rasters ####
writeRaster(weather_composites,
            filename = "./processed_data/weather_composite.tif",
            overwrite = TRUE)
## gridded control plots ####
write_sf(gridded_plots, dsn = "./processed_data/gridded_candidate_plots.shp")
tictoc::toc()
