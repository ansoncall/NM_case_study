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




## load progression boundaries 
CBI<-raster("./raw_data/burn_severity/ravg_2022_cbi4.tif")
CBI[CBI==9]<-NA
CBI[CBI==0]<-1

progression<-read_sf("./raw_data/fire_progression/PerimeterLine2022.shp") ## https://data-nifc.opendata.arcgis.com/datasets/nifc::perimeterline2022/explore?location=35.823546%2C-105.140721%2C10.40


progression<-raster("./processed_data/julian_day_of_burn_masked.tiff")

#progression<-sf::st_transform(progression,crs=crs(CBI))

#progression<-progression[progression$SourceGlob!="{C42FB181-2486-4D92-AE97-E1348E4E032F}",]

#do_not_use<-c("Cooks Peak","Grass Mountain Fire")

#progression<-progression[progression$IncidentNa %in% do_not_use==FALSE,]

#library(stringr)

#progression$Comments<-tolower(progression$Comments)

#progression<-progression %>%
  filter(str_detect(Comments, "ir")) # should include IR flights and firewatch/ firescan



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

weather_raster<-data.frame(rasterToPoints(progression==dates[40]))
weather_raster<-weather_raster[weather_raster$layer==1,]

#weather_raster<-data.frame(rasterToPoints(raster::crop(empty_raster,first_day_of_fire,snap="in")))



weather_raster$fm100<-raster::extract(fm100,weather_raster[,1:2])
weather_raster$fm1000<-raster::extract(fm1000,weather_raster[,1:2])
weather_raster$vs<-raster::extract(vs,weather_raster[,1:2])
weather_raster$rmax<-raster::extract(rmax,weather_raster[,1:2])
weather_raster$th<-raster::extract(th,weather_raster[,1:2])
weather_raster$tmmx<-raster::extract(th,weather_raster[1,2])

weather_raster$xy<-paste(weather_raster$x,weather_raster$y,sep="_")

library(concaveman)

for (i in 2:length(dates)){
#for (i in 2:5){
  
  next_day_of_fire<-st_make_valid(progression[progression$GDB_FROM_D==dates[i],])
  
  #next_day_of_fire<-next_day_of_fire[st_is_simple(next_day_of_fire),]


  next_day_of_fire<-st_combine(st_as_sf(next_day_of_fire))
  
  
  
  polys<-st_polygonize(next_day_of_fire)
  
  next_day_of_fire<-st_buffer(next_day_of_fire,30)
  
  #next_day_of_fire<-st_union(st_make_valid(next_day_of_fire),st_simplify(st_make_valid(st_buffer(polys,0))))
  next_day_of_fire<-st_union(st_make_valid(next_day_of_fire),st_buffer(polys,0))
  
  
  
  #next_day_of_fire<-st_cast(next_day_of_fire,"MULTIPOLYGON")
  
  #next_day_of_fire_closed<-concaveman(next_day_of_fire[1,])
  
  #for (j in 2:nrow(next_day_of_fire)){
    #close_geo<-concaveman(next_day_of_fire[1,])
  #  next_day_of_fire_closed[j,]<-concaveman(next_day_of_fire[j,])
    
 # }
  
  #next_day_of_fire<-next_day_of_fire_closed
  
  fm100<-projectRaster(raster("./raw_data/fm100_2022.nc",as.numeric(format(dates[i],"%j"))),to=CBI)
  fm1000<-projectRaster(raster("./raw_data/fm1000_2022.nc",as.numeric(format(dates[i],"%j"))),to=CBI)
  vs<-projectRaster(raster("./raw_data/vs_2022.nc",as.numeric(format(dates[i],"%j"))),to=CBI)
  rmax<-projectRaster(raster("./raw_data/rmax_2022.nc",as.numeric(format(dates[i],"%j"))),to=CBI)
  th<-projectRaster(raster("./raw_data/th_2022.nc",as.numeric(format(dates[i],"%j"))),to=CBI)
  tmmx<-projectRaster(raster("./raw_data/tmmx_2022.nc",as.numeric(format(dates[i],"%j"))),to=CBI)
  
  next_day_of_fire<-st_as_sf(next_day_of_fire)
  
  weather_raster_new<-raster::crop(empty_raster,next_day_of_fire)
  weather_raster_new<-raster::mask(weather_raster_new,next_day_of_fire)
  weather_raster_new<-data.frame(rasterToPoints(weather_raster_new))
  
  
  weather_raster_new$fm100<-raster::extract(fm100,weather_raster_new[,1:2])
  weather_raster_new$fm1000<-raster::extract(fm1000,weather_raster_new[,1:2])
  weather_raster_new$vs<-raster::extract(vs,weather_raster_new[,1:2])
  weather_raster_new$rmax<-raster::extract(rmax,weather_raster_new[,1:2])
  weather_raster_new$th<-raster::extract(th,weather_raster_new[,1:2])
  weather_raster_new$xy<-paste(weather_raster_new$x,weather_raster_new$y,sep="_")
  weather_raster_new$tmmx<-raster::extract(tmmx,weather_raster_new[,1:2])
  
  weather_raster_new<-weather_raster_new[weather_raster_new$xy %in% weather_raster$xy==FALSE,]
  
  weather_raster<-rbind(weather_raster,weather_raster_new)
  
}

weather_raster_total<-rasterFromXYZ(weather_raster,crs=crs(CBI))

#weather_raster_total<-raster::mask(weather_raster_total,st_cast(first_day_of_fire,to="POLYGON"))

mapview(weather_raster_total$fm100)

fm100_final <- focal(weather_raster_total$fm100, w=matrix(1, nrow=11, ncol=11), fun=mean ,na.rm=TRUE,NAonly=TRUE)
fm1000_final <- focal(weather_raster_total$fm1000, w=matrix(1, nrow=11, ncol=11), fun=mean ,na.rm=TRUE,NAonly=TRUE)
vs_final <- focal(weather_raster_total$vs, w=matrix(1, nrow=11, ncol=11), fun=mean ,na.rm=TRUE,NAonly=TRUE)
rmax_final <- focal(weather_raster_total$rmax, w=matrix(1, nrow=11, ncol=11), fun=mean ,na.rm=TRUE,NAonly=TRUE)
th_final <- focal(weather_raster_total$th, w=matrix(1, nrow=11, ncol=11), fun=mean ,na.rm=TRUE,NAonly=TRUE)
tmmx_final<-focal(weather_raster_total$tmmx, w=matrix(1, nrow=11, ncol=11), fun=mean ,na.rm=TRUE,NAonly=TRUE)


for (i in 1:200){
  fm1000_final <- focal(fm1000_final, w=matrix(1, nrow=11, ncol=11), fun=mean ,na.rm=TRUE,NAonly=TRUE)
  fm100_final <- focal(fm100_final, w=matrix(1, nrow=11, ncol=11), fun=mean ,na.rm=TRUE,NAonly=TRUE)
  vs_final <- focal(vs_final, w=matrix(1, nrow=11, ncol=11), fun=mean ,na.rm=TRUE,NAonly=TRUE)
  rmax_final <- focal(rmax_final, w=matrix(1, nrow=11, ncol=11), fun=mean ,na.rm=TRUE,NAonly=TRUE)
  th_final <- focal(th_final, w=matrix(1, nrow=11, ncol=11), fun=mean ,na.rm=TRUE,NAonly=TRUE)
  tmmx_final <- focal(tmmx_final, w=matrix(1, nrow=11, ncol=11), fun=mean ,na.rm=TRUE,NAonly=TRUE)
}

mapview(th_final)
mapview(tmmx_final)

writeRaster(fm1000_final,filename = "./processed_data/fm1000.tif")
writeRaster(fm100_final,filename = "./processed_data/fm100.tif")
writeRaster(vs_final,filename = "./processed_data/vs.tif")
writeRaster(rmax_final,filename = "./processed_data/rmax.tif")
writeRaster(th_final,filename = "./processed_data/th.tif")
writeRaster(tmmx_final,filename = "./processed_data/tmmx.tif")



############################################

#####  Create a grid of control plots. Control plots will be 10 acre circles to match the validation plots. They will be built on a grid.I need to first calculate the area of square grid cell that that contains a 10 acre circle. It is 12.732 acres

burn_perimiter<-read_sf("./processed_data/burn_perimiter.shp")

total_burn_area<-sum(st_area(burn_perimiter))/4046.68 ## total burn acres

how_may_cells<-as.numeric(total_burn_area)/12.732

gridded_points<-st_sample(burn_perimiter,size=round(how_may_cells,0),type="regular")

gridded_plots<-st_buffer(st_as_sf(gridded_points),dist=113.49694)

mapview(gridded_plots)



elev_down<-raster("./processed_data/elev_down.tif")
aspect_down<-raster("./processed_data/aspect_down.tif")
TRI_down<-raster("./processed_data/TRI_down.tif")
TPI_down<-raster("./processed_data/TPI_down.tif")
slope_down<-raster("./processed_data/slope_down.tif")
roads_distance<-raster("./processed_data/distance_to_road.tif")
site_potential<-raster("./processed_data/lf_site_potential_new.tif")

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
gridded_plots$cbi<-raster::extract(masked_cbi,st_as_sf(gridded_plots),fun=mean,na.rm=FALSE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)


gridded_plots<-gridded_plots[is.na(gridded_plots$cbi)==FALSE,]


write_sf(gridded_plots,dsn="gridded_candidate_plots.shp")
