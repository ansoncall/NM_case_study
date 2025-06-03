###### this code follows the methods described in Parks 2014 Mapping day-of-burning with coarse-resolution satellite fire-detection data.
######
library(sf)
library(raster)
library(lubridate)
library(mapview)
library(terra)
library(FNN)
library(terra)

########


d<-read_sf("./raw_data/satelite_fire_detection/fire_archive_M-C61_605340.shp")

mapview(d)

## trim data to burn scar and transform

burn_perimiter<-read_sf("./processed_data/burn_perimiter.shp")

d<-st_transform(d,st_crs(burn_perimiter))

d<-st_crop(d,burn_perimiter)

d$ACQ_DATE<-ymd(d$ACQ_DATE)

d$jdate<-as.numeric(format(d$ACQ_DATE,"%j"))

## Empty raster
CBI<-raster("./processed_data/clipped_burn_raster.tif")
CBI[CBI==9]<-NA
CBI[CBI==0]<-1


ycell<-CBI@nrows
xcell<-CBI@ncols
bounds<-st_bbox(CBI)

empty_raster<-raster(nrows=ycell,ncols=xcell,xmn=bounds[1],xmx=bounds[3],ymn=bounds[2],ymx=bounds[4],crs=st_crs(CBI)$wkt)

##

days_raster<-raster::rasterize(d,empty_raster,field="jdate")

##########

t_days<-rast(days_raster)


d_raster<-distance(t_days)

####################
coords<-crds(t_days,na.rm=TRUE)
jday_vals <- values(t_days)[!is.na(values(t_days))]
d_vals <- values(d_raster)[!is.na(values(t_days))]

knn_model <- get.knnx(coords, coords, k = 5)


filled_r <- t_days

# Loop through NA cells and compute weighted value
na_cells <- which(is.na(values(t_days)))
all_coords <- crds(t_days,na.rm=FALSE)

remaining_na_cells<-na_cells[na_cells>i]
  
for (i in remaining_na_cells) {
  target_coord <- all_coords[i, , drop = FALSE]
  
  # Find 5 nearest non-NA neighbors
  nn <- get.knnx(coords, target_coord, k = 5)
  idx <- nn$nn.index[1, ]
  
  nn.coords<-all_coords[idx,]
  
  d_i<-distance(target_coord,nn.coords)
  
  jdays <- jday_vals[idx]
  d_i <- values(d_raster)[i]  # Use distance value from d raster
  
  if (!is.na(d_i)) {
    mean_jday <- mean(jdays)
    w_i <- (1 / (abs(jdays - mean_jday) + 1)) * d_i
    values(filled_r)[i] <- round(sum((w_i*jdays))/sum(w_i),0)
  }
}

mapview(filled_r)

raster::writeRaster(filled_r,filename="./processed_data/julian_day_of_burn.tiff",overwrite=TRUE)

raster::writeRaster(mask(filled_r,burn_perimiter),filename="./processed_data/julian_day_of_burn_masked.tiff",overwrite=TRUE)

mapview(burn_perimiter)+mapview(mask(filled_r,burn_perimiter))


