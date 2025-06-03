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
hpcc <- read_sf("raw_data/fire_perimeters/Perimeters.shp") %>% 
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
  # intersect with the burn perimeter to remove treatments outside the burn scar
  # st_intersection(hpcc) %>% # not needed, see below 
  # remove treatments with centroids that are not in the burn scar Nate: you
  # only need to run the intersection once - if the treatment centroid isn't in
  # the burn scar, the full polygon won't intersect either.
  st_intersection(st_centroid(.), hpcc)

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
ycell<-cbi@nrows
xcell<-cbi@ncols
bounds<-st_bbox(cbi)

# confirm measurements are in meters
st_crs(veg_treatments)$units

# buffer the vegetation treatments by 10m
veg_treatments_buffer <- st_buffer(veg_treatments_cropped, 10)

# TODO not sure this is all that necessary, could just use cbi raster to
# define dimensions
empty_raster <- raster(nrows = ycell, ncols = xcell, xmn = bounds[1], 
                       xmx = bounds[3], ymn = bounds[2], ymx = bounds[4], 
                       crs = st_crs(veg_treatments_buffer))

rast_veg_full <- rasterize(veg_treatments_buffer, empty_raster)

## using this as a mask requires that the areas we want to remove (the treatments) be NA values
# create mask by defining raster with NA values for treatment areas
rast_veg_full[is.na(rast_veg_full[])] <- (-1) 
rast_veg_full[rast_veg_full[]>0]<-NA

plot(rast_veg_full)

masked_cbi<-raster::mask(cbi,rast_veg_full)

plot(masked_cbi)

## save the masked raster for model training
writeRaster(masked_cbi, "./processed_data/masked_raster.tif", overwrite=TRUE)



###################### Roads ################################ 

# level 0 road segments
roads_0 <- st_read("./raw_data/transportation/Trans_RoadSegment_0.shp") %>% 
  st_transform(crs = proj) %>% 
  st_crop(bounds)

# level 1 road segments
roads_1 <- st_read("./raw_data/transportation/Trans_RoadSegment_1.shp") %>% 
  st_transform(crs = proj) %>% 
  st_crop(bounds)

# combine road segments
roads <- rbind(roads_0, roads_1) %>% 
  # buffer by 10m
  # Nate: previously, you had created a buffered roads layer, but exported the
  # non-buffered one. Is this correct? If you're not using the buffered roads,
  # do you need to create the buffer in the first place?
  st_buffer()


# Note, this looks ok. has some forest service roads but not every single one.
st_write(roads, "./processed_data/burn_scar_roads.shp", append=TRUE) 

# distance to road raster
roads_raster_distance <-raster::rasterize(roads, empty_raster) %>% 
  terra::distance

writeRaster(roads_raster_distance, 
            "./processed_data/distance_to_road.tif", 
            overwrite=TRUE)


################# Topographic data

elev <- raster("./raw_data/HPCC_elevations.tif")
# todo should be terra::terrain?
slopes <- terrain(elev, "slope") %>% resample()
aspect <- terrain(elev, "aspect")
TRI <- terrain(elev, "TRI")
TPI <- terrain(elev, "TPI")


# define function to aggregate and reproject topographic rasters to match CBI
transform_raster <- function(from_raster, to_raster = cbi) {
  # define aggregation factors
  fact <- c(round(from_raster@nrows/to_raster@nrows), 
            round(from_raster@ncols/to_raster@nrows))
  # aggregate raster
  transformed <- terra::aggregate(from_raster, fact = fact, fun = "mean") %>% 
    # reproject to target CRS
    projectRaster(crs=st_crs(to_raster))
  return(transformed)
}




####### manipulate elevation data to calculate aspect, topogr
elev_down<-terra::aggregate(elev,fact=c(yfact,xfact),fun="mean")
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
