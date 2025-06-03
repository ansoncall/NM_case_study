############## Load packages

library(ggplot2)
library(dplyr)
library(sf)
library(raster)
library(stars)
library(elevatr)
library(mapview)


##### Download the weather data from gridmet

downloader::download(
  url = "http://www.northwestknowledge.net/metdata/data/rmax_2022.nc",
  destfile = "raw_data/rmax_2022.nc",
  mode = 'wb'
)

downloader::download(
  url = "http://www.northwestknowledge.net/metdata/data/vs_2022.nc",
  destfile = "raw_data/vs_2022.nc",
  mode = 'wb'
)

downloader::download(
  url = "http://www.northwestknowledge.net/metdata/data/th_2022.nc",
  destfile = "raw_data/th_2022.nc",
  mode = 'wb'
)

downloader::download(
  url = "http://www.northwestknowledge.net/metdata/data/fm100_2022.nc",
  destfile = "raw_data/fm100_2022.nc",
  mode = 'wb'
)

downloader::download(
  url = "http://www.northwestknowledge.net/metdata/data/fm1000_2022.nc",
  destfile = "raw_data/fm1000_2022.nc",
  mode = 'wb'
)

downloader::download(
  url = "http://www.northwestknowledge.net/metdata/data/tmmx_2022.nc",
  destfile = "raw_data/tmmx_2022.nc",
  mode = 'wb'
)


# load CBI raster
CBI<-raster("./raw_data/burn_severity/ravg_2022_cbi4.tif")
CBI[CBI==9]<-NA # 9 is the code for unmappable and I need it to not be considered as a number
CBI[CBI==0]<-1 # in the raw data a zero is outside of the fire perimiter and a 1 is a pixel that is unchanged. For this analysis we consider these to be similar. For instance, when kriging it is better to have the information from the area outside the perimiter considered as unchanged instead of unknown. 

## load progression raster
progression<-raster("./processed_data/julian_day_of_burn_masked.tiff")

dates<-sort(unique(progression$layer))
dates_converted<-as.Date(dates,origin=as.Date("2024-01-01"))


#first_day_of_fire<-progression[progression$GDB_FROM_D==dates[1],]


### load weather data

ycell<-CBI@nrows
xcell<-CBI@ncols
bounds<-st_bbox(CBI)
######################### Raster of Vegetation Treatments with a 10m buffer

empty_raster<-raster(nrows=ycell,ncols=xcell,xmn=bounds[1],xmx=bounds[3],ymn=bounds[2],ymx=bounds[4],crs=st_crs(CBI))
empty_raster$dummy=1


burn_perimiter<-read_sf("./processed_data/burn_perimiter.shp")

progression



#https://tmieno2.github.io/R-as-GIS-for-Economists/gridMET.html

fm100<-projectRaster(raster("./raw_data/fm100_2022.nc",as.numeric(format(dates_converted[1],"%j"))),to=CBI) ## 100 hour fuel moisture
fm1000<-projectRaster(raster("./raw_data/fm1000_2022.nc",as.numeric(format(dates_converted[1],"%j"))),to=CBI) ## 1000 hour fuel moisture
vs<-projectRaster(raster("./raw_data/vs_2022.nc",as.numeric(format(dates_converted[1],"%j"))),to=CBI) ## wind speed at 10 M
rmax<-projectRaster(raster("./raw_data/rmax_2022.nc",as.numeric(format(dates_converted[1],"%j"))),to=CBI) ## rel humidity maximum
th<-projectRaster(raster("./raw_data/th_2022.nc",as.numeric(format(dates_converted[1],"%j"))),to=CBI) ## high temp wind direction at 10m
tmax<-projectRaster(raster("./raw_data/tmmx_2022.nc",as.numeric(format(dates_converted[1],"%j"))),to=CBI) ## high temp wind direction at 10m

weather_raster<-data.frame(rasterToPoints(progression==dates[1]))
weather_raster<-weather_raster[weather_raster$layer==1,]


weather_raster$fm100<-raster::extract(fm100,weather_raster[,1:2])
weather_raster$fm1000<-raster::extract(fm1000,weather_raster[,1:2])
weather_raster$vs<-raster::extract(vs,weather_raster[,1:2])
weather_raster$rmax<-raster::extract(rmax,weather_raster[,1:2])
weather_raster$th<-raster::extract(th,weather_raster[,1:2])
weather_raster$tmmx<-raster::extract(th,weather_raster[1,2])

weather_raster$day_of_burn<-dates_converted[1]

weather_raster$xy<-paste(weather_raster$x,weather_raster$y,sep="_")


## looping through the burn days and extracting the weather data for the appropriate cells. 
for (i in 2:length(dates)){

  fm100<-projectRaster(raster("./raw_data/fm100_2022.nc",as.numeric(format(dates_converted[i],"%j"))),to=CBI)
  fm1000<-projectRaster(raster("./raw_data/fm1000_2022.nc",as.numeric(format(dates_converted[i],"%j"))),to=CBI)
  vs<-projectRaster(raster("./raw_data/vs_2022.nc",as.numeric(format(dates_converted[i],"%j"))),to=CBI)
  rmax<-projectRaster(raster("./raw_data/rmax_2022.nc",as.numeric(format(dates_converted[i],"%j"))),to=CBI)
  th<-projectRaster(raster("./raw_data/th_2022.nc",as.numeric(format(dates_converted[i],"%j"))),to=CBI)
  tmmx<-projectRaster(raster("./raw_data/tmmx_2022.nc",as.numeric(format(dates_converted[i],"%j"))),to=CBI)
  

  weather_raster_new<-data.frame(rasterToPoints(progression==dates[i]))
  weather_raster_new<-weather_raster_new[weather_raster_new$layer==1,]

  weather_raster_new$fm100<-raster::extract(fm100,weather_raster_new[,1:2])
  weather_raster_new$fm1000<-raster::extract(fm1000,weather_raster_new[,1:2])
  weather_raster_new$vs<-raster::extract(vs,weather_raster_new[,1:2])
  weather_raster_new$rmax<-raster::extract(rmax,weather_raster_new[,1:2])
  weather_raster_new$th<-raster::extract(th,weather_raster_new[,1:2])
  weather_raster_new$xy<-paste(weather_raster_new$x,weather_raster_new$y,sep="_")
  weather_raster_new$tmmx<-raster::extract(tmmx,weather_raster_new[,1:2])
  
  weather_raster_new$day_of_burn<-dates_converted[i]
  
  weather_raster_new<-weather_raster_new[weather_raster_new$xy %in% weather_raster$xy==FALSE,]
  
  weather_raster<-rbind(weather_raster,weather_raster_new)
  
}

weather_raster_total<-rasterFromXYZ(weather_raster,crs=crs(CBI))

#weather_raster_total<-raster::mask(weather_raster_total,st_cast(first_day_of_fire,to="POLYGON"))

mapview(weather_raster_total$fm100)


## rename for the difference columns
writeRaster(weather_raster_total$fm1000,filename = "./processed_data/fm1000.tif",overwrite=TRUE)
writeRaster(weather_raster_total$fm100,filename = "./processed_data/fm100.tif",overwrite=TRUE)
writeRaster(weather_raster_total$vs,filename = "./processed_data/vs.tif",overwrite=TRUE)
writeRaster(weather_raster_total$rmax,filename = "./processed_data/rmax.tif",overwrite=TRUE)
writeRaster(weather_raster_total$th,filename = "./processed_data/th.tif",overwrite=TRUE)
writeRaster(weather_raster_total$tmmx,filename = "./processed_data/tmmx.tif",overwrite=TRUE)



############################################

#####  Create a grid of control plots. Control plots will be 10 acre circles to match the validation plots. They will be built on a grid.The area of square grid cell that that contains a 10 acre is 12.732 acres

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
