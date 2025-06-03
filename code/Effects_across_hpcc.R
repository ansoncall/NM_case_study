# load packages
library(raster)
library(terra)
library(sf)
library(stars)
library(ggplot2)
library(spmodel)

############# Load processed data

elev_down<-raster("./processed_data/elev_down.tif")
aspect_down<-raster("./processed_data/aspect_down.tif")
TRI_down<-raster("./processed_data/TRI_down.tif")
TPI_down<-raster("./processed_data/TPI_down.tif")
slope_down<-raster("./processed_data/slope_down.tif")
roads_distance<-raster("./processed_data/distance_to_road.tif")
site_potential<-raster("./processed_data/lf_site_potential_new.tif")
vpdmax<-raster("./processed_data/vpdmax.tif")
ppt<-raster("./processed_data/ppt.tif")
tmin<-raster("./processed_data/tmin.tif")

masked_cbi<-raster("./processed_data/masked_raster.tif")

CBI<-raster("./processed_data/clipped_burn_raster.tif")

CBI[CBI==9]<-NA
CBI[CBI==0]<-1

veg_treatments<-read_sf("./processed_data/vegetation_treatments_hpcc_new.shp")
veg_treatments<-st_make_valid(veg_treatments)
veg_treatments<-st_transform(veg_treatments,crs=crs(masked_cbi))
veg_treatments_buffer<-st_buffer(veg_treatments,60)

masked_cbi<-raster::mask(masked_cbi,st_as_sf(veg_treatments_buffer),inverse=TRUE)

burn_perimiter<-read_sf("./processed_data/burn_perimiter.shp")

spatial_rf_interative<-function(one_plot){

  size_of_perimiter<-sqrt((150*4046.86+as.numeric(st_area(one_plot)))/pi)
  
  one_plot_buffer<-st_buffer(one_plot,size_of_perimiter)
  
  clipped_burn_raster<-raster::crop(masked_cbi,st_as_sf(one_plot_buffer))
  
  clipped_burn_raster[clipped_burn_raster==0]<-1
  clipped_burn_raster[clipped_burn_raster==9]<-NA
  
  training_df<-as.data.frame(rasterToPoints(clipped_burn_raster))
  
  training_df$CBI<-training_df$masked_raster
  training_df$elev<-raster::extract(elev_down,y=training_df[,1:2])
  training_df$aspect<-raster::extract(aspect_down,y=training_df[,1:2])
  training_df$TRI<-raster::extract(TRI_down,y=training_df[,1:2])
  training_df$TPI<-raster::extract(TPI_down,y=training_df[,1:2])
  training_df$slope<-raster::extract(slope_down,y=training_df[,1:2])
  training_df$road_distance<-raster::extract(roads_distance,y=training_df[,1:2])
  training_df$env_potential<-as.factor(raster::extract(site_potential,y=training_df[,1:2]))
  training_df$ppt<-raster::extract(ppt,y=training_df[,1:2])
  training_df$tmin<-raster::extract(ppt,y=training_df[,1:2])
  training_df$ppt<-raster::extract(tmin,y=training_df[,1:2])
  training_df$vpdmax<-raster::extract(vpdmax,y=training_df[,1:2])
  
  training_df<-training_df[complete.cases(training_df),]
  
  training_df<-st_as_sf(training_df,coords=c("x","y"),crs=st_crs(site_potential))
  
  training_df$lat<-as.data.frame(rasterToPoints(clipped_burn_raster))[,1]
  training_df$lon<-as.data.frame(rasterToPoints(clipped_burn_raster))[,2]
  training_df$lat2<-training_df$lat**2
  training_df$lon2<-training_df$lon**2
  
  
  try(spatial_rf_model<-splmRF(CBI~elev+aspect+TRI+TPI+slope+road_distance+env_potential+lat+lon+tmin+vpdmax+ppt,data=training_df,spcov_type = "exponential",local=c(parallel=TRUE,ncores=10),mtry=4,min.node.size=2,sample.fraction=0.89))
  
  #

  
  predict_df<-data.frame(lat=rasterToPoints(raster::crop(site_potential,y=st_as_sf(one_plot)))[,1],
                         lon=rasterToPoints(raster::crop(site_potential,y=st_as_sf(one_plot)))[,2])
  
  predict_df$elev<-raster::extract(elev_down,y=predict_df[,1:2])
  predict_df$aspect<-raster::extract(aspect_down,y=predict_df[,1:2])
  predict_df$TRI<-raster::extract(TRI_down,y=predict_df[,1:2])
  predict_df$TPI<-raster::extract(TPI_down,y=predict_df[,1:2])
  predict_df$slope<-raster::extract(slope_down,y=predict_df[,1:2])
  predict_df$road_distance<-raster::extract(roads_distance,y=predict_df[,1:2])
  predict_df$env_potential<-as.factor(raster::extract(site_potential,y=predict_df[,1:2]))
  predict_df$ppt<-raster::extract(ppt,y=predict_df[,1:2])
  predict_df$tmin<-raster::extract(ppt,y=predict_df[,1:2])
  predict_df$ppt<-raster::extract(tmin,y=predict_df[,1:2])
  predict_df$vpdmax<-raster::extract(vpdmax,y=predict_df[,1:2])
  
  
  predict_df<-predict_df[complete.cases(predict_df),]
  
  predict_df<-st_as_sf(predict_df,coords=c("lat","lon"),crs=st_crs(site_potential))
  
  predict_df$lat<-rasterToPoints(raster::crop(site_potential,y=st_as_sf(one_plot)))[,1]
  predict_df$lon<-rasterToPoints(raster::crop(site_potential,y=st_as_sf(one_plot)))[,2]
  
  #predict_df$modeled_values<-predict(ranger_rf,data=predict_df,na.rm=TRUE)$predictions
  if (sd(training_df$CBI,na.rm=TRUE)==0){predict_df$modeled_values<-mean(training_df$CBI,na.rm=TRUE)}
  
  if (sd(training_df$CBI,na.rm=TRUE)>0){predict_df$modeled_values<-predict(spatial_rf_model,newdata=predict_df,local=FALSE)}
  
  mean_cbi<-predict_df$modeled_values[1]
  for_conversion<-data.frame(lat=predict_df$lat,lon=predict_df$lon,modeled_values=as.numeric(predict_df$modeled_values))
  
  if (nrow(predict_df)>1){
  
    if (length(unique(predict_df$lat))==1 |length(unique(predict_df$lon))==1){mean_cbi=mean(predict_df$modeled_values)
    
    
    
    }else{
    

  
  predict_raster<-rasterFromXYZ(for_conversion,crs=crs(masked_cbi))
  
  mean_cbi<-raster::extract(predict_raster$modeled_values,st_as_sf(one_plot),fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)
  }}
  
  return(list(mean_cbi,for_conversion))
  
  
}

descriptions<-tolower(unique(veg_treatments$Dscrptn))

thinning<-descriptions[grepl("thin",descriptions)]
cutting<-descriptions[grepl("cut",descriptions)]
defens<-descriptions[grepl("defens",descriptions)]
burn<-descriptions[grepl("burn",descriptions)]
fuel<-descriptions[grepl("fuel",descriptions)]
fire<-descriptions[grepl("fire",descriptions)]


treatments_to_use<-unique(c(thinning,cutting,defens,burn,fuel,fire))


veg_treatments<-veg_treatments[tolower(veg_treatments$Dscrptn) %in% treatments_to_use,]


veg_treatments$control_burn_severity<-NA
map_values<-data.frame()

l<-nrow(veg_treatments)

veg_treatments<-st_cast(veg_treatments,"POLYGON")

for (i in 1:nrow(veg_treatments)){

  print(i/l*100)
  output<-spatial_rf_interative(veg_treatments[i,])
    
    
  veg_treatments$control_burn_severity[i]<-output[[1]]
  map_values<-rbind(map_values,output[2][[1]])
  
  st_write(veg_treatments,dsn = "./results/hpcc_processed_cbi_new.shp",append=FALSE)
  write.csv(map_values,"./results/mapping_predictions.csv")
}

#### adding in propensity score matching here


veg_treatments<-read_sf("./results/hpcc_processed_cbi_new.shp")

gridded_plots<-read_sf("gridded_candidate_plots.shp")


gridded_plots<-st_transform(gridded_plots,st_crs(CBI))
#3control_plots<-st_transform(control_plots,st_crs(CBI))

gridded_plots_2<-st_difference(gridded_plots,st_union(veg_treatments))


gridded_plots_3<-st_difference(gridded_plots_2,st_union(st_buffer(validation_plots,60)))

max_area<-max(st_area(gridded_plots_3))

gridded_plots_4<-gridded_plots_3[as.numeric(st_area(gridded_plots_3))>as.numeric(max_area)-50,]

library(caret)
#### normalizing variables first

gridded_plots_4$elevation_norm<-(gridded_plots_4$elevatn-mean(gridded_plots_4$elevatn))/sd(gridded_plots_4$elevatn)
gridded_plots_4$ppt_norm<-(gridded_plots_4$nrml_pp-mean(gridded_plots_4$nrml_pp))/sd(gridded_plots_4$nrml_pp)
gridded_plots_4$vs_norm<-(gridded_plots_4$vs-mean(gridded_plots_4$vs))/sd(gridded_plots_4$vs)


nearest_model<-knnreg(cbi~elevation_norm+ppt_norm+vs_norm,data=gridded_plots_4,k=1)
summary(nearest_model)

ppt<-raster("./processed_data/ppt.tif")
elev_down<-raster("./processed_data/elev_down.tif")
vs<-raster("./processed_data/vs.tif")

predict_df<-st_as_sf(veg_treatments)

predict_df$elev<-raster::extract(elev_down,y=predict_df,fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)
predict_df$vs<-raster::extract(vs,y=st_as_sf(predict_df),fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)
predict_df$ppt<-raster::extract(ppt,y=st_as_sf(predict_df),fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)

predict_df$elevation_norm<-(predict_df$elev-mean(gridded_plots_4$elevatn))/sd(gridded_plots_4$elevatn)
predict_df$ppt_norm<-(predict_df$ppt-mean(gridded_plots_4$nrml_pp))/sd(gridded_plots_4$nrml_pp)
predict_df$vs_norm<-(predict_df$vs-mean(gridded_plots_4$vs))/sd(gridded_plots_4$vs)



predict_df$knn_1<-predict(nearest_model,as.data.frame(predict_df))


predict_df<-as.data.frame(predict_df)

veg_treatments$knn_pred<-predict_df$knn_1


st_write(veg_treatments,dsn = "./results/hpcc_processed_cbi_new.shp",append=FALSE)



