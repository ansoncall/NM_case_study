library(ggplot2)
library(dplyr)
library(sf)
library(raster)
library(stars)
library(elevatr)


validation_plots<-read_sf("./processed_data/val_points_revised.shp",fid_column_name="ID")
#write_sf(validation_plots,dsn="./processed_data/branch_validation_points.shp")
pairing<-read.csv("./processed_data/plot_matching.csv")

gridded_plots<-read_sf("gridded_candidate_plots.shp")

masked_cbi<-raster("./processed_data/masked_raster.tif")

masked_cbi[masked_cbi==9]<-NA
masked_cbi[masked_cbi==0]<-1

CBI[CBI==9]<-NA
CBI[CBI==0]<-1

plot_method_cbi<-as.data.frame(raster::extract(CBI,st_as_sf(validation_plots),fun=mean,na.rm=FALSE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE))

plot_method_cbi<-data.frame(FID=validation_plots$ID,cbi=plot_method_cbi$V1)

validation_cbi<-plot_method_cbi[plot_method_cbi$FID>352,]


control_plots<-merge(plot_method_cbi,pairing,by.x="FID",by.y="FID_1")

names(control_plots)[names(control_plots) == "cbi"] <- "origional_plot_cbi"

paired_plots<-merge(control_plots,validation_cbi,by.x="FID.of.New.plot",by.y="FID")


################# How well did the handmade controls do?



compare_df<-data.frame(plot_id=paired_plots$FID,origional_plot_cbi=paired_plots$origional_plot_cbi,control_plot_cbi=paired_plots$cbi,skip=paired_plots$skip)

compare_df<-compare_df[is.na(compare_df$origional_plot_cbi)==FALSE,]
compare_df<-compare_df[is.na(compare_df$control_plot_cbi)==FALSE,]


plot(compare_df$control_plot_cbi,compare_df$origional_plot_cbi)



## RMSE
sqrt(mean((compare_df$control_plot_cbi-compare_df$origional_plot_cbi)**2,na.rm=TRUE))

origional_plots<-validation_plots[validation_plots$ID %in% compare_df$plot_id,]


############# Simple perimiters
size_of_perimiter<-sqrt((33*4046.86)/pi)-sqrt((10*4046.86)/pi)


perimiters<-st_buffer(origional_plots,size_of_perimiter,allow_holes=TRUE)
st_area(perimiters)/4046.86
og_plots<-st_buffer(origional_plots,60)
st_area(og_plots)/4046.86


buffer_sizes<-function(plot_buffer_size){


compare_df$buffer_method<-NA

for (i in 1:50){
  one_plot<-origional_plots[origional_plots$FID_1==compare_df$plot_id[i],]
  
  area_of_plot<-st_area(one_plot)/4046.86
  
  buffer_plot<-st_buffer(one_plot,plot_buffer_size)
  area_of_plot_with_buffer<-st_area(buffer_plot)/4046.86
  
  size_of_perimiter<-sqrt(((area_of_plot+area_of_plot_with_buffer)*4046.86)/pi)-sqrt((area_of_plot)*4046.86/pi)
  
  perimiters<-st_buffer(one_plot,size_of_perimiter,allow_holes=TRUE)
  area_perimiter<-st_area(perimiters)/4046.86
  
  perimiter_extract<-raster::extract(masked_cbi,st_as_sf(perimiters),fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)
  origional_extract<-raster::extract(masked_cbi,st_as_sf(buffer_plot),fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)
  
  perimiter_only<-as.numeric(((perimiter_extract[1]*area_perimiter)-(origional_extract[1]*area_of_plot_with_buffer))/(area_perimiter-area_of_plot_with_buffer))
  
  compare_df[compare_df$plot_id==one_plot$FID_1,"buffer_method"]<-perimiter_only
  
}

rmse<-sqrt(mean((compare_df$origional_plot_cbi-compare_df$buffer_method)**2,na.rm=TRUE))
bias<-mean(compare_df$origional_plot_cbi-compare_df$buffer_method)

return(data.frame(rmse=rmse,bias=bias))
}


buffer_sizes(60)

sizes_to_test<-c(10,30,60,90,150,250,350,500)

comparison<-data.frame()

for (i in 1:length(sizes_to_test)){
  one_size<-buffer_sizes(sizes_to_test[i])
  one_size$buffer_size=sizes_to_test[i]
  comparison<-rbind(comparison,one_size)
}


ggplot(comparison,aes(x=buffer_size,y=rmse))+geom_point(size=2)+theme_classic()+xlab("Size of buffer excluded (m)")+ylab("RMSE")+theme(text=element_text(size=20))

ggplot(comparison,aes(x=buffer_size,y=bias))+geom_point(size=2)+theme_classic()+xlab("Size of buffer excluded (m)")+ylab("Bias")+theme(text=element_text(size=20))


################# What about the spatial rf models?


##### spatial rf try two

size_of_pannel<-200
size_of_buffer<-60

library(spmodel)

sensitivity_analyses<-function(size_of_pannel,size_of_buffer){

local_spatial_rf<-c()

for (i in 1:50){
  #for (i in 1:6){
  one_plot<-origional_plots[origional_plots$FID_1==compare_df$plot_id[i],]
  
  size_of_perimiter<-sqrt((size_of_pannel*4046.86)/pi)
  
  one_plot_buffer<-st_buffer(one_plot,size_of_perimiter)
  
  clipped_burn_raster<-raster::crop(CBI,st_as_sf(one_plot_buffer))
  clipped_burn_raster[clipped_burn_raster==0]<-1
  clipped_burn_raster[clipped_burn_raster==9]<-NA
  
  clipped_burn_raster<-raster::mask(clipped_burn_raster,st_buffer(one_plot,size_of_buffer),inverse=TRUE)
  
  training_df<-as.data.frame(rasterToPoints(clipped_burn_raster))
  
  training_df$CBI<-training_df$masked_raster
  training_df$elev<-raster::extract(elev_down,y=training_df[,1:2])
  training_df$aspect<-raster::extract(aspect_down,y=training_df[,1:2])
  training_df$TRI<-raster::extract(TRI_down,y=training_df[,1:2])
  training_df$TPI<-raster::extract(TPI_down,y=training_df[,1:2])
  training_df$slope<-raster::extract(slope_down,y=training_df[,1:2])
  training_df$road_distance<-raster::extract(roads_distance,y=training_df[,1:2])
  training_df$env_potential<-as.factor(raster::extract(site_potential,y=training_df[,1:2]))
  
  training_df<-training_df[complete.cases(training_df),]
  
  training_df<-st_as_sf(training_df,coords=c("x","y"),crs=st_crs(site_potential))
  
  training_df$lat<-as.data.frame(rasterToPoints(clipped_burn_raster))[,1]
  training_df$lon<-as.data.frame(rasterToPoints(clipped_burn_raster))[,2]
  
  
  spatial_rf_model<-splmRF(clipped_burn_raster~elev+aspect+TRI+TPI+slope+road_distance+env_potential+lat+lon,data=training_df,spcov_type = "exponential",local=c(parallel=TRUE,ncores=12),mtry=4,min.node.size=2,sample.fraction=0.89)
  
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
  
  predict_df<-predict_df[complete.cases(predict_df),]
  
  predict_df<-st_as_sf(predict_df,coords=c("lat","lon"),crs=st_crs(site_potential))
  
  predict_df$lat<-rasterToPoints(raster::crop(site_potential,y=st_as_sf(one_plot)))[,1]
  predict_df$lon<-rasterToPoints(raster::crop(site_potential,y=st_as_sf(one_plot)))[,2]
  
  #predict_df$modeled_values<-predict(ranger_rf,data=predict_df,na.rm=TRUE)$predictions
  predict_df$modeled_values<-predict(spatial_rf_model,newdata=predict_df,local=FALSE)
  
  for_conversion<-data.frame(lat=predict_df$lat,lon=predict_df$lon,modeled_values=as.numeric(predict_df$modeled_values))
  
  predict_raster<-rasterFromXYZ(for_conversion,crs=crs(masked_cbi))
  
  mean_cbi<-raster::extract(predict_raster$modeled_values,st_as_sf(one_plot),fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)
  
  local_spatial_rf[i]<-mean_cbi
  
  
}


rmse_output<-sqrt(mean((local_spatial_rf-compare_df$origional_plot_cbi[1:length(local_spatial_rf)])**2,na.rm=TRUE))
rmse_output

output<-data.frame(pannel_size=(ncell(clipped_burn_raster)*30*30)/4046,buffer_size=size_of_buffer,rmse=rmse_output)

return(output)
}




conditions_to_test<-data.frame(pannel=rep(c(200,1200,2200),times=length(sizes_to_test)),buffer=rep(sizes_to_test,each=3))

conditions_to_test_2<-conditions_to_test[conditions_to_test$buffer<(conditions_to_test$pannel-20),]


spatial_rf_sensitivity<-data.frame()

for (i in 1:nrow(conditions_to_test_2)){
  output<-sensitivity_analyses(conditions_to_test_2$pannel[i],conditions_to_test_2$buffer[i])
  
  spatial_rf_sensitivity<-rbind(spatial_rf_sensitivity,output)
  
  print((i/nrow(conditions_to_test_2))*100)
  
}


comparison_2<-data.frame(pannel_size="buffer",buffer_size=comparison$buffer_size,rmse=comparison$rmse)


spatial_rf_sensitivity$pannel_size<-round(spatial_rf_sensitivity$pannel_size,0)
plotting_data<-rbind(comparison_2,spatial_rf_sensitivity)

cbPalette <- c( "#E69F00", "#56B4E9", "#009E73", "#D55E00", "#CC79A7")

sensitivity<-ggplot(plotting_data,aes(x=buffer_size,y=rmse,color=pannel_size))+geom_jitter(size=2)+theme_classic()+xlab("Size of buffer around plots (m)")+ylab("RMSE")+theme(text=element_text(size=20))+ theme(legend.position = c(0.8, 0.2))+scale_color_manual(values=cbPalette,name="Training size (acres)")
sensitivity

tiff(filename=("./figures/sensitivity_analysis.tif"),units='in',compression='lzw',width=8,height=8,res=300)
sensitivity
dev.off()

write.csv(plotting_data,"./processed_data/comparing_buffer_sizes.csv")
