## this script works through the burn severity, topography, and climate data to
## transform and load the data and do any initial clipping that is required.

############## Load packages

library(ggplot2)
library(dplyr)
library(sf)
library(raster)
library(stars)
library(elevatr)
library(mapview)

################ Load data

# load CBI data
cbi <- raster("./raw_data/burn_severity/ravg_2022_cbi4.tif")
# save spatial reference for use with other data
proj <- crs(cbi)

# load fire perimeters
hpcc <- read_sf("./raw_data/fire_perimeters/Perimeters.shp") %>% 
  # subset to HPCC fires
  subset(.$poly_Incid %in% c("Calf Canyon","Hermits Peak")) %>% 
  # ensure valid geometries
  st_make_valid %>% 
  # reproject
  st_transform(crs = proj)

# load vegetation treatment data
veg_treatments <- read_sf("./raw_data/NMVeg.gdb") %>%
  # ensure valid geometries
  st_make_valid %>%
  # reproject
  st_transform(crs = proj) %>%
  # subset by ROI. treatment centroid must be within HPCC burn boundary.
  filter(apply(st_intersects(st_centroid(.), hpcc, sparse = FALSE), 1, any))

# crop CBI data to ROI
cbi_crop <- crop(cbi, hpcc)

# write out projected perimeter data
write_sf(hpcc, dsn="./processed_data/burn_perimeter.shp")
# write out cropped CBI raster
writeRaster(cbi, 
            filename="./processed_data/clipped_burn_raster.tif",
            overwrite=TRUE)
# write out cropped vegetation treatment data
st_write(veg_treatments, 
         dsn="./processed_data/vegetation_treatments_hpcc_new.shp",
         append=FALSE)

######################### Raster of Vegetation Treatments with a 10m buffer
# define raster dimensions
ycell<-cbi_crop@nrows
xcell<-cbi_crop@ncols
bounds<-st_bbox(cbi_crop)

# confirm measurements are in meters
if (st_crs(veg_treatments)$units != "m") {
  stop("Vegetation treatments are not in meters\n. 
       Please reproject to a CRS with meters as units.")
} else {
  message("Vegetation treatments are in meters.")
}

# buffer the vegetation treatments by 10m
veg_treatments_buffer <- st_buffer(veg_treatments, 10) %>% 
  st_union %>% 
  st_sf

# TODO not sure this is all that necessary, could just use cbi raster to
# define dimensions
empty_raster <- raster(nrows = ycell, ncols = xcell, xmn = bounds[1], 
                       xmx = bounds[3], ymn = bounds[2], ymx = bounds[4], 
                       crs = st_crs(veg_treatments_buffer))

library(fasterize)
rast_veg_full <- fasterize(veg_treatments_buffer, cbi_crop)

plot(rast_veg_full)
mapview(rast_veg_full)
## using this as a mask requires that the areas we want to remove (the treatments) be NA values
# create mask by defining raster with NA values for treatment areas
rast_veg_full[is.na(rast_veg_full[])] <- (-1) 
rast_veg_full[rast_veg_full[]>0]<-NA

terra::terraOptions(memfrac = 0.8, progress = 1)
# slow
masked_cbi <- terra::mask(cbi_crop, veg_treatments_buffer, inverse = TRUE)

plot(masked_cbi)

## save the masked raster for model training
writeRaster(masked_cbi, "./processed_data/masked_raster.tif", overwrite=TRUE)



###################### Roads ################################ 
# list shapefiles to read
roadlist <- c("./raw_data/transportation/Trans_RoadSegment_0.shp", 
              "./raw_data/transportation/Trans_RoadSegment_1.shp")
roads <- lapply(roadlist, st_read) %>% 
  # combine road segments
  do.call(rbind, .) %>% 
  # reproject to match CBI
  st_transform(crs = proj) %>% 
  # ensure valid geometries
  st_make_valid %>% 
  # crop to bounds of CBI raster
  st_crop(bounds) 

  # buffer by 10m
  # Nate: previously, you had created a buffered roads layer, but exported the
  # non-buffered one. Is this correct? If you're not using the buffered roads,
  # do you need to create the buffer in the first place?
  # st_buffer()


# Nate: the RoadSegment_0 shapefile did not download correctly, so I could not
# test this. It looks okay though. 
st_write(roads, "./processed_data/burn_scar_roads.shp", append=TRUE) 

# distance to road raster
roads_raster_distance <-raster::rasterize(roads, empty_raster) %>% 
  terra::distance

writeRaster(roads_raster_distance, 
            "./processed_data/distance_to_road.tif", 
            overwrite=TRUE)


################# Topographic data
# Nate: could not download this one either. I just took the 3DEP 1/3 arc  second
# DEM, I assume that's what you started with
elev_0 <- raster("./raw_data/USGS_13_n36w106_20250311.tif")
elev_1 <- raster("./raw_data/USGS_13_n37w106_20250311.tif")
# combine the two DEMs
elev <- merge(elev_0, elev_1) 
mapview(elev)

# todo should be terra::terrain?
slopes <- terrain(elev, "slope")
aspect <- terrain(elev, "aspect")
TRI <- terrain(elev, "TRI")
TPI <- terrain(elev, "TPI")


# define function to aggregate and reproject topographic rasters to match CBI
transform_raster <- function(from_raster, to_raster = cbi) {
  # define aggregation factors
  facts <- list(round(from_raster@nrows/to_raster@nrows), 
            round(from_raster@ncols/to_raster@nrows))
  # aggregate raster
  transformed <- terra::aggregate(from_raster, fact = facts, fun = "mean") %>% 
    # reproject to target CRS
    projectRaster(crs=st_crs(to_raster))
  return(transformed)
}

yfact <- round(elev@nrows/cbi@nrows)
xfact <- round(elev@ncols/cbi@ncols)

####### manipulate elevation data to calculate aspect, topogr
elev_down<-terra::aggregate(elev,fact=c(yfact,xfact),fun="mean")
elev_down <- transform_raster(elev, cbi)
aspect_down<-terra::aggregate(aspect,fact=c(yfact,xfact),fun="mean")
TRI_down<-terra::aggregate(TRI,fact=c(yfact,xfact),fun="mean")
TPI_down<-terra::aggregate(TPI,fact=c(yfact,xfact),fun="mean")
slope_down<-terra::aggregate(slopes,fact=c(yfact,xfact),fun="mean")

### transform everything to correct CRS
#crs(masked_cbi)<-CRS('+init=EPSG:4269')

elev_down<-projectRaster(elev_down,masked_cbi)
aspect_down<-projectRaster(aspect_down,masked_cbi)
TRI_down<-projectRaster(TRI_down,masked_cbi)
TPI_down<-projectRaster(TPI_down,masked_cbi)
slope_down<-projectRaster(slope_down,masked_cbi)


writeRaster(elev_down,"./processed_data/elev_down.tif",overwrite=TRUE)
writeRaster(aspect_down,"./processed_data/aspect_down.tif",overwrite=TRUE)
writeRaster(TRI_down,"./processed_data/TRI_down.tif",overwrite=TRUE)
writeRaster(TPI_down,"./processed_data/TPI_down.tif",overwrite=TRUE)
writeRaster(slope_down,"./processed_data/slope_down_degrees.tif",overwrite=TRUE)

rm(elev,slopes,aspect,TRI,TPI)

################## Environmental site potential 
#crs(masked_cbi)<-CRS('+init=EPSG:4269')


lf_site_potential<-raster("./raw_data/landfire_environmental_site_potential/Tif/us_140esp.tif")
plot(lf_site_potential)

#projected_cbi<-projectRaster(cbi_raster,crs=crs(lf_site_potential))                         


#plot(lf_site_potential_cropped)

## this takes a very long time to run
lf_site_potential_projected<-projectRaster(lf_site_potential,crs=crs(cbi),method="ngb")

lf_site_potential_cropped<-crop(lf_site_potential_projected,cbi)

plot(lf_site_potential_cropped)

lf_site_potential_cropped@data<-as.factor(lf_site_potential_cropped)

is.factor(lf_site_potential_cropped)

lf_site_potential_cropped_down<-terra::resample(lf_site_potential_cropped,elev_down,method="ngb") ## ngb is nearest neighbor



writeRaster(lf_site_potential_cropped_down,"./processed_data/lf_site_potential_new.tif",overwrite=TRUE)

rm(lf_site_potential,lf_site_potential_cropped)

#### climate variables

tmin<-raster("./raw_data/prism_climate/PRISM_tmin_30yr_normal_800mM5_annual_asc.asc")
tmin<-projectRaster(tmin,cbi)

vpdmax<-raster("./raw_data/prism_climate/PRISM_vpdmax_30yr_normal_800mM5_annual_asc.asc")
vpdmax<-projectRaster(vpdmax,cbi)

ppt<-raster("./raw_data/prism_climate/PRISM_ppt_30yr_normal_800mM4_annual_asc.asc")
ppt<-projectRaster(ppt,cbi)

writeRaster(tmin,"./processed_data/tmin.tif")
writeRaster(vpdmax,"./processed_data/vpdmax.tif")
writeRaster(ppt,"./processed_data/ppt.tif")
