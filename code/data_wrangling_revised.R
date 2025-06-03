## this script works through the burn severity, topography, and climate data to transform and load the data and do any initial clipping that is required. 

############## Load packages

library(ggplot2)
library(dplyr)
library(sf)
library(raster)
library(stars)
library(elevatr)
library(mapview)

################ Load data

# Fire perimiters
bounds<-read_sf("raw_data/fire_perimiters/Perimeters.shp") 

# subset the hermits peak and calf canyon fires
hpcc <- subset(bounds, bounds$poly_Incid %in% c("Calf Canyon","Hermits Peak"))
ggplot(hpcc)+geom_sf()
hpcc<-st_make_valid(hpcc)

# open the vegetation treatment data
veg_treatments<-read_sf("./raw_data/NMVeg.gdb")
veg_treatments<-st_make_valid(veg_treatments)

# open CBI data, this is what I reference all other coordinate reference systems too
CBI<-raster("./raw_data/burn_severity/ravg_2022_cbi4.tif")

#transform the boundaries for HPCC boundaries
hpcc_projected<-st_transform(hpcc,crs(CBI))

write_sf(hpcc_projected,dsn="./processed_data/burn_perimiter.shp")

CBI<-crop(CBI,hpcc_projected)
raster::writeRaster(CBI,filename="./processed_data/clipped_burn_raster.tif",overwrite=TRUE)

veg_treatments<-st_transform(veg_treatments,crs(CBI))

##########################

# crop the NM vegetation treatment database to the burn scar
veg_treatments_cropped<-st_intersection(veg_treatments,hpcc_projected)


centroids_in_scar<-st_intersection(st_centroid(veg_treatments),hpcc_projected)

veg_treatments_cropped<-veg_treatments[veg_treatments$GlobalID %in% centroids_in_scar$GlobalID,]

mapview(hpcc_projected)+mapview(veg_treatments_cropped)

st_write(veg_treatments_cropped,dsn="./processed_data/vegetation_treatments_hpcc_new.shp",append=FALSE)
#veg_treatments_cropped_2<-st_transform(veg_treatments_cropped,crs=4269)

ycell<-CBI@nrows
xcell<-CBI@ncols
bounds<-st_bbox(CBI)
######################### Raster of Vegetation Treatments with a 10m buffer
# confirm measurements are in meters
st_crs(veg_treatments_cropped)$units

veg_treatments_buffer<-st_buffer(veg_treatments_cropped,10)

empty_raster<-raster(nrows=ycell,ncols=xcell,xmn=bounds[1],xmx=bounds[3],ymn=bounds[2],ymx=bounds[4],crs=st_crs(veg_treatments_buffer))

rast_veg_full<-raster::rasterize(veg_treatments_buffer,empty_raster)

plot(rast_veg_full)

## using this as a mask requires that the areas we want to remove (the treatments) be NA values

rast_veg_full[is.na(rast_veg_full[])] <- (-1) 
rast_veg_full[rast_veg_full[]>0]<-NA

plot(rast_veg_full)

#####################

masked_cbi<-raster::mask(CBI,rast_veg_full)

plot(masked_cbi)

## save the masked raster for model training
writeRaster(masked_cbi,"./processed_data/masked_raster.tif",overwrite=TRUE)



###################### Roads ################################ 


roads1<-st_read("./raw_data/transportation/Trans_RoadSegment_0.shp")
roads1<-st_transform(roads1,st_crs(veg_treatments_cropped))
roads1_cropped<-st_crop(roads1,st_bbox(CBI))


roads2<-st_read("./raw_data/transportation/Trans_RoadSegment_1.shp")
roads2<-st_transform(roads2,st_crs(veg_treatments_cropped))
roads2_cropped<-st_crop(roads2,st_bbox(CBI))

rm(roads1,roads2)
roads<-rbind(roads1_cropped,roads2_cropped)
roads_buffer<-st_buffer(roads,10)## 10 meter buffer around roads

st_write(roads,"./processed_data/burn_scar_roads.shp",append=TRUE) ### looks OK. has some forest service roads but not every single one.

#roads_raster<-st_rasterize(roads%>% dplyr::select(geometry), nx = 700, ny = 700)

roads_raster<-raster::rasterize(roads,empty_raster)

plot(roads_raster)

## making a raster of distance to road

roads_raster_distance<-terra::distance(roads_raster)

plot(roads_raster_distance)


writeRaster(roads_raster_distance,"./processed_data/distance_to_road.tif",overwrite=TRUE)


################# Topographic data

elev<-raster("./raw_data/HPCC_elevations.tif")

slopes<-terrain(elev,v="slope")
aspect<-terrain(elev,"aspect")
TRI<-terrain(elev,"TRI")
TPI<-terrain(elev,"TPI")


### figuring out how much to aggregate these data
xfact<-round(elev@nrows/xcell)
yfact<-round(elev@ncols/ycell)



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

#projected_cbi<-projectRaster(CBI_raster,crs=crs(lf_site_potential))                         


#plot(lf_site_potential_cropped)

## this takes a very long time to run
lf_site_potential_projected<-projectRaster(lf_site_potential,crs=crs(CBI),method="ngb")

lf_site_potential_cropped<-crop(lf_site_potential_projected,CBI)

plot(lf_site_potential_cropped)

lf_site_potential_cropped@data<-as.factor(lf_site_potential_cropped)

is.factor(lf_site_potential_cropped)

lf_site_potential_cropped_down<-terra::resample(lf_site_potential_cropped,elev_down,method="ngb") ## ngb is nearest neighbor



writeRaster(lf_site_potential_cropped_down,"./processed_data/lf_site_potential_new.tif",overwrite=TRUE)

rm(lf_site_potential,lf_site_potential_cropped)

#### climate variables

tmin<-raster("./raw_data/prism_climate/PRISM_tmin_30yr_normal_800mM5_annual_asc.asc")
tmin<-projectRaster(tmin,CBI)

vpdmax<-raster("./raw_data/prism_climate/PRISM_vpdmax_30yr_normal_800mM5_annual_asc.asc")
vpdmax<-projectRaster(vpdmax,CBI)

ppt<-raster("./raw_data/prism_climate/PRISM_ppt_30yr_normal_800mM4_annual_asc.asc")
ppt<-projectRaster(ppt,CBI)

writeRaster(tmin,"./processed_data/tmin.tif")
writeRaster(vpdmax,"./processed_data/vpdmax.tif")
writeRaster(ppt,"./processed_data/ppt.tif")
