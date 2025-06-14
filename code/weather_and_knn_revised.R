# this script prepares weather data and identifies candidate control plots for
# the burn severity analysis. it (optionally) downloads gridmet data, processes
# it, and extracts relevant features for each day of the burn. It also creates a
# grid of control plots based on the burn perimeter and extracts environmental
# features for these plots.

# load packages ####
library(tidyverse)
library(downloader)
library(sf)
library(terra)
library(conflicted)
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
plot(cbi)

## day of burn (DOB) ####
dob <- rast("./processed_data/julian_day_of_burn.tiff")
plot(dob)
# get unique dob values, sort them, and convert to Date objects
dates <- unique(dob) %>%
  pull(1) %>%
  sort

## HPCC burn perimeter ####
burn_perimeter<-read_sf("./processed_data/burn_perimeter.shp")

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
          mode = 'wb'
        )
      })
}
# load the weather data from disk
weather_rasts <- map(
  weather_vars,
  \(x) {
    # read in rasters, subset to get the bands for the relevant days only
    rast(paste0("./raw_data/", x, "_2022.nc"), lyrs = dates)
  }) %>%
  # stack all bands: 6 vars * 71 days = 426 bands
  rast %>%
  project(cbi)

# tidy layer names and set time value to dob
names(weather_rasts) <- rep(weather_vars, each = length(dates))
time(weather_rasts) <- rep(dates, 6)

# map over the rasters and mask them by the appropriate date subset of the dob
# raster.
weather_masked <- map(
  dates,
  # i == the ordinal day of burn
  \(i) {
    # i = 104
    # create dob mask
    msk <- ifel(dob == i, 1, NA)
    # plot(msk)
    # select the raster layers with the target i
    weather_rasts_subset <- weather_rasts[[time(weather_rasts) == i]]
    # plot(weather_rasts_subset,
    #      main = paste("Weather data for day of burn:", i))
    # mask them
    masked_weather_rasts <- mask(weather_rasts[[time(weather_rasts) == i]], msk)
    # plot(masked_weather_rasts,
    #      main = paste("masked:", i))
  },
  .progress = "Masking weather rasters..."
) %>%
  do.call(c, .)

# composite the partial rasters for each output variable
weather_composites <-  map(
  weather_vars,
  \(x) {
    # select all the rasters for the variable x
    var_subset <- weather_masked[[grep(x, names(weather_masked))]]
    # plot(var_subset[[1:5]],
    #      main = paste("var:", x, " day:", dates[1:5]))
    # get first non-na layer value (there should only be one)
    out <- app(var_subset, mean, na.rm = TRUE, wopt = list(names = x))
  },
  .progress = "Building composites..."
) %>%
  do.call(c, .)

plot(weather_composites)


# create control plots ####

# Create a grid of control plots. Control plots will be 10 acre circles to match
# the validation plots. They will be built on a grid.The area of square grid
# cell that that contains a 10 acre is 12.732 acres

burn_perimiter<-read_sf("./processed_data/burn_perimiter.shp")

total_burn_area<-sum(st_area(burn_perimiter))/4046.68 ## total burn acres

how_may_cells<-as.numeric(total_burn_area)/12.732

gridded_points<-st_sample(burn_perimiter,size=round(how_may_cells,0),type="regular")

gridded_plots<-st_buffer(st_as_sf(gridded_points),dist=113.49694)

mapview(gridded_plots)

## grid plots are made now need to extract all of the data for each gridded plot

elev_down<-raster("./processed_data/elev_down.tif")
aspect_down<-raster("./processed_data/aspect_down.tif")
TRI_down<-raster("./processed_data/TRI_down.tif")
TPI_down<-raster("./processed_data/TPI_down.tif")
slope_down<-raster("./processed_data/slope_down.tif")
roads_distance<-raster("./processed_data/distance_to_road.tif")
site_potential<-raster("./processed_data/lf_site_potential_new.tif")
ppt<-raster("./processed_data/ppt.tif")
tmin<-raster("./processed_data/tmin.tif")
tmmx<-raster("./processed_data/tmmx.tif")
th<-raster("./processed_data/th.tif")
vpdmax<-raster("./processed_data/vpdmax.tif")
rmax<-raster("./processed_data/rmax.tif")
vs<-raster("./processed_data/vs.tif")
fm100<-raster("./processed_data/fm100.tif")
fm1000<-raster("./processed_data/fm1000.tif")


masked_cbi<-raster("./processed_data/masked_raster.tif")
masked_cbi[masked_cbi==9]<-NA
masked_cbi[masked_cbi==0]<-1
# Nate: terra::extract is supposed to be faster. Maybe a good place to use mclapply as well.
gridded_plots$elevation<-raster::extract(elev_down,st_as_sf(gridded_plots),fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)
gridded_plots$aspect<-raster::extract(aspect_down,st_as_sf(gridded_plots),fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)
gridded_plots$tri<-raster::extract(TRI_down,st_as_sf(gridded_plots),fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)
gridded_plots$tpi<-raster::extract(TPI_down,st_as_sf(gridded_plots),fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)
gridded_plots$slope<-raster::extract(slope_down,st_as_sf(gridded_plots),fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)
gridded_plots$road<-raster::extract(roads_distance,st_as_sf(gridded_plots),fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)
gridded_plots$esp<-raster::extract(site_potential,st_as_sf(gridded_plots),fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)
gridded_plots$cbi<-raster::extract(masked_cbi,st_as_sf(gridded_plots),fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)
gridded_plots$fm1000<-raster::extract(fm1000,st_as_sf(gridded_plots),fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)
gridded_plots$fm100<-raster::extract(fm100,st_as_sf(gridded_plots),fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)
gridded_plots$ppt<-raster::extract(ppt,st_as_sf(gridded_plots),fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)
gridded_plots$tmin<-raster::extract(tmin,st_as_sf(gridded_plots),fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)
gridded_plots$tmmx<-raster::extract(tmmx,st_as_sf(gridded_plots),fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)
gridded_plots$th<-raster::extract(th,st_as_sf(gridded_plots),fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)
gridded_plots$vpdmax<-raster::extract(vpdmax,st_as_sf(gridded_plots),fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)
gridded_plots$rmax<-raster::extract(rmax,st_as_sf(gridded_plots),fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)
gridded_plots$vs<-raster::extract(vs,st_as_sf(gridded_plots),fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)
gridded_plots$fm100<-raster::extract(fm100,st_as_sf(gridded_plots),fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)
gridded_plots$fm100<-raster::extract(fm100,st_as_sf(gridded_plots),fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)
gridded_plots$fm1000<-raster::extract(fm1000,st_as_sf(gridded_plots),fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)


write_sf(gridded_plots,dsn="gridded_candidate_plots.shp")
