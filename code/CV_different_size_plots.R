library(raster)
library(terra)
library(sf)
library(stars)
library(ggplot2)

############# Load processed data

elev_down<-raster("./processed_data/elev_down.tif")
aspect_down<-raster("./processed_data/aspect_down.tif")
TRI_down<-raster("./processed_data/TRI_down.tif")
TPI_down<-raster("./processed_data/TPI_down.tif")
slope_down<-raster("./processed_data/slope_down.tif")
roads_distance<-raster("./processed_data/distance_to_road.tif")
site_potential<-raster("./processed_data/lf_site_potential_new.tif")

masked_cbi<-raster("./processed_data/masked_raster.tif")

veg_treatments<-read_sf("./processed_data/vegetation_treatments_hpcc.shp")
veg_treatments<-st_make_valid(veg_treatments)

burn_perimiter<-read_sf("./processed_data/burn_perimiter.shp")

gridded_plots<-read_sf("gridded_candidate_plots.shp")

####### extracting data from the rasters to a dataframe. Make sure these match how the rasters were created in the data_warning_revised_script
ycell<-masked_cbi@nrows
xcell<-masked_cbi@ncols


extract_df<-data.frame(lat=rep(seq(masked_cbi@extent[1],masked_cbi@extent[2],length.out=ycell),each=xcell),lon=rep(seq(masked_cbi@extent[3],masked_cbi@extent[4],length.out=ycell),times=xcell))



extract_df$CBI<-raster::extract(masked_cbi,y=as.matrix(extract_df[,1:2]))
extract_df$elev<-raster::extract(elev_down,y=as.matrix(extract_df[,1:2]))
extract_df$aspect<-raster::extract(aspect_down,y=as.matrix(extract_df[,1:2]))
extract_df$TRI<-raster::extract(TRI_down,y=as.matrix(extract_df[,1:2]))
extract_df$TPI<-raster::extract(TPI_down,y=as.matrix(extract_df[,1:2]))
extract_df$slope<-raster::extract(slope_down,y=as.matrix(extract_df[,1:2]))
extract_df$road_distance<-raster::extract(roads_distance,y=as.matrix(extract_df[,1:2]))
extract_df$env_potential<-raster::extract(site_potential,y=as.matrix(extract_df[,1:2]))


extract_df<-extract_df[extract_df$CBI!=0,]
extract_df<-extract_df[extract_df$CBI!=9,]

### site potentials are numerical codes but should be considered as factors
extract_df$env_potential<-as.factor(paste("env_",extract_df$env_potential,sep=""))

library(ranger)

predict_df<-extract_df[is.na(extract_df$CBI)==TRUE,]
extract_df_full<-extract_df

extract_df<-extract_df[complete.cases(extract_df),]

## what sizes are treatment areas

# 4046.86 m2 in an acre

plot_sizes<-data.frame(acres=c(2,5,15,50,100,150,250))
plot_sizes$radius_m<-sqrt((plot_sizes$acres*4046.86)/pi)

veg_treatments$actual_size<-st_area(veg_treatments)/4046.86
hist(veg_treatments$actual_size)
quantile(veg_treatments$Acre_US,c(0.05,0.25,0.50,0.75,0.95))

#CBI<-st_make_valid(CBI)
#CBI<-st_union(CBI)

set.seed(2076)


random_point_in_burn<-st_sample(burn_perimiter,size=nrow(plot_sizes)*60)

random_point_in_burn<-st_buffer(random_point_in_burn,dist=rep(plot_sizes$radius_m))

random_point_in_burn<-sf::st_transform(random_point_in_burn,crs=st_crs(masked_cbi))




size_cv_results<-data.frame()

##### Setup rf models
masked_cbi[masked_cbi==9]<-NA
masked_cbi[masked_cbi==0]<-1

######

library(spmodel)
library(snapKrig)

size_cv<-function (index_number){
  ## this function needs to mask the validation area, train the RF model, predict on the validation area, then compare to the actual values in the validation area.
  validation_plot<-random_point_in_burn[index_number]
  #validation_plot<-index_number
  ### actual CBI in validation area
  
  actual_burn<-raster::extract(x=masked_cbi,y=st_as_sf(validation_plot),na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)
  
  #output=data.frame(actual_burn_severity=NA,modeled_burn=NA,size_acres=as.numeric(st_area(validation_plot)/4046.86))
  
  if (sum(is.na(actual_burn[[1]]))==0){
  
  #actual_burn_severity<-mean(actual_burn[[1]])

    
    validation_mask_cbi<-raster::mask(masked_cbi,st_as_sf(st_buffer(validation_plot,60)),inverse=TRUE)
    
    size_of_perimiter<-sqrt((150*4046.86+as.numeric(st_area(validation_plot)))/pi)
    
    one_plot_buffer<-st_buffer(validation_plot,size_of_perimiter)
    
    clipped_burn_raster<-raster::crop(validation_mask_cbi,st_as_sf(one_plot_buffer))
    
    xcell<-clipped_burn_raster@ncols
    ycell<-clipped_burn_raster@nrows
    
    extract_df<-data.frame(lat=rep(seq(clipped_burn_raster@extent[1],clipped_burn_raster@extent[2],length.out=ycell),each=xcell),lon=rep(seq(clipped_burn_raster@extent[3],clipped_burn_raster@extent[4],length.out=ycell),times=xcell))
    
    
    
    extract_df$CBI<-raster::extract(validation_mask_cbi,y=as.matrix(extract_df[,1:2]))
    extract_df$elev<-raster::extract(elev_down,y=as.matrix(extract_df[,1:2]))
    extract_df$aspect<-raster::extract(aspect_down,y=as.matrix(extract_df[,1:2]))
    extract_df$TRI<-raster::extract(TRI_down,y=as.matrix(extract_df[,1:2]))
    extract_df$TPI<-raster::extract(TPI_down,y=as.matrix(extract_df[,1:2]))
    extract_df$slope<-raster::extract(slope_down,y=as.matrix(extract_df[,1:2]))
    extract_df$road_distance<-raster::extract(roads_distance,y=as.matrix(extract_df[,1:2]))
    extract_df$env_potential<-as.factor(raster::extract(site_potential,y=as.matrix(extract_df[,1:2])))
    
    predict_df<-extract_df[is.na(extract_df$CBI)==TRUE,]
    extract_df_full<-extract_df
    
    
    extract_df<-extract_df[extract_df$CBI!=0,]
    extract_df<-extract_df[extract_df$CBI!=9,]
    
    extract_df<-extract_df[complete.cases(extract_df),]

    
    
    
    #ranger_rf<-ranger(as.factor(CBI)~elev+aspect+TRI+TPI+slope+lat+lon+road_distance+env_potential,extract_df,num.trees=1000,importance="impurity",mtry=5,min.node.size=2,sample.fraction=0.89)
    
    library(spmodel)
    
    coordinates<- data.frame(lat=extract_df$lat,lon=extract_df$lon)
    
    
    extract_sf<-st_as_sf(extract_df,coords=c("lat","lon"))
    extract_sf$lat<-extract_df$lat
    extract_sf$lon<-extract_df$lon
    
    spatial_rf_model<-splmRF(CBI~elev+aspect+TRI+TPI+slope+road_distance+env_potential+lat+lon,data=extract_sf,spcov_type = "gravity",local=c(parallel=TRUE,ncores=20),mtry=4,min.node.size=2,sample.fraction=0.89)
    
    
    
  
  predict_df<-data.frame(lat=rasterToPoints(raster::crop(site_potential,y=st_as_sf(validation_plot)))[,1],
                         lon=rasterToPoints(raster::crop(site_potential,y=st_as_sf(validation_plot)))[,2])

  
  predict_df$elev<-raster::extract(elev_down,y=predict_df[,1:2])
  predict_df$aspect<-raster::extract(aspect_down,y=predict_df[,1:2])
  predict_df$TRI<-raster::extract(TRI_down,y=predict_df[,1:2])
  predict_df$TPI<-raster::extract(TPI_down,y=predict_df[,1:2])
  predict_df$slope<-raster::extract(slope_down,y=predict_df[,1:2])
  predict_df$road_distance<-raster::extract(roads_distance,y=predict_df[,1:2])
  predict_df$env_potential<-as.factor(raster::extract(site_potential,y=predict_df[,1:2]))
  
  predict_df$actual_burn<-raster::extract(x=masked_cbi,y=predict_df[,1:2],na.rm=TRUE,,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)
  
  predict_df<-predict_df[predict_df$actual_burn!=0,]
  predict_df<-predict_df[predict_df$actual_burn!=9,]
  
  predict_df<-predict_df[complete.cases(predict_df),]
  predict_df$lat2<-predict_df$lat
  predict_df$lon2<-predict_df$lon
    ### spatial rf
  predict_df<-st_as_sf(predict_df,coords=c("lat2","lon2"),crs=st_crs(site_potential))
  

  
  #predict_df$modeled_values<-predict(ranger_rf,data=predict_df,na.rm=TRUE)$predictions
  predict_df$modeled_values<-predict(spatial_rf_model,newdata=predict_df)
  
  for_conversion<-data.frame(lat=predict_df$lat,lon=predict_df$lon,modeled_values=as.numeric(as.character(predict_df$modeled_values)))
  
  predict_raster<-rasterFromXYZ(for_conversion,crs=crs(masked_cbi))
  
  mean_cbi_spatial_rf<-raster::extract(predict_raster$modeled_values,st_as_sf(validation_plot),fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)
  
  
  
  ### perimiter method
  area_of_plot<-st_area(validation_plot)/4046.86
  buffer_plot<-st_buffer(validation_plot,60)
  area_of_plot_with_buffer<-st_area(buffer_plot)/4046.86
  
  
  size_of_perimiter<-sqrt(((area_of_plot+area_of_plot_with_buffer)*4046.86)/pi)-sqrt((area_of_plot)*4046.86/pi)
  
  perimiters<-st_buffer(validation_plot,size_of_perimiter,allow_holes=TRUE)
  area_perimiter<-st_area(perimiters)/4046.86

  perimiter_extract<-raster::extract(masked_cbi,st_as_sf(perimiters),fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)
  origional_extract<-raster::extract(masked_cbi,st_as_sf(buffer_plot),fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)
  
  perimiter_only<-as.numeric(((perimiter_extract[1]*area_perimiter)-(origional_extract[1]*area_of_plot_with_buffer))/(area_perimiter-area_of_plot_with_buffer))
  
  ### kriging
  
  size_of_perimiter<-sqrt((150*4046.86+as.numeric(st_area(validation_plot)))/pi)
  
  one_plot_buffer<-st_buffer(validation_plot,size_of_perimiter)
  
  clipped_burn_raster<-raster::crop(validation_mask_cbi,st_as_sf(one_plot_buffer))
  clipped_burn_raster[clipped_burn_raster==0]<-1
  clipped_burn_raster[clipped_burn_raster==9]<-NA
  
  clipped_burn_raster<-raster::mask(clipped_burn_raster,st_as_sf(st_buffer(validation_plot,60)),inverse=TRUE)
  
  clipped_burn_raster_sk<-sk(clipped_burn_raster)
  
  if (is.na(clipped_burn_raster[1,1])==TRUE){clipped_burn_raster[1,1]<-1}
  if (is.na(clipped_burn_raster[nrow(clipped_burn_raster),1])==TRUE){clipped_burn_raster[nrow(clipped_burn_raster),1]<-1}
  if (is.na(clipped_burn_raster[1,ncol(clipped_burn_raster)])==TRUE){clipped_burn_raster[1,ncol(clipped_burn_raster)]<-1}
  if (is.na(clipped_burn_raster[nrow(clipped_burn_raster),ncol(clipped_burn_raster)])==TRUE){clipped_burn_raster[nrow(clipped_burn_raster),ncol(clipped_burn_raster)]<-1}
  
  kriging_fit<-sk_fit(clipped_burn_raster_sk,n_max=15000)
  
  
  krigged_final<-sk_cmean(clipped_burn_raster_sk,kriging_fit)
  plot(krigged_final)
  
  krigged_final_raster<-sk_export(g=krigged_final,template = 'raster')
  
  kriging_predictions<-extract(krigged_final_raster,st_as_sf(validation_plot),fun=mean)
  
  
  #############

  
  ##############

  
  ##### summary output
  
  actual_burn_severity<-raster::extract(x=masked_cbi,y=st_as_sf(validation_plot),fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)
  

  output<-data.frame(id_number=index_number,actual_burn_severity=actual_burn_severity,perimiter_cbi_value=perimiter_only,kriging_cbi=kriging_predictions,cbi_rf_spatial=mean_cbi_spatial_rf,size_acres=as.numeric(st_area(validation_plot)/4046.86))
  
  
  }

  return(output)
}

#cv_data<-read.csv("./results/plot_size_cross_validation.csv")

cv_data<-data.frame()
for (i in 1:420){
#for (i in 1:25){
    
  print(i)
  try(cv_data<-rbind(cv_data,size_cv(i)))
  write.csv(cv_data,file="./results/plot_size_cross_validation_new.csv")

}

#####
library(dplyr)

cv_data<-read.csv("./results/plot_size_cross_validation_new.csv")
cv_data$acre_bins<-round(cv_data$size_acres)

########## Adding clustering method #######
## remove the treated and validation plots
veg_treatments<-read_sf("./processed_data/vegetation_treatments_hpcc_new.shp")

gridded_plots<-st_transform(gridded_plots,st_crs(CBI))
#3control_plots<-st_transform(control_plots,st_crs(CBI))

gridded_plots_2<-st_difference(gridded_plots,st_union(veg_treatments))


gridded_plots_3<-st_difference(gridded_plots_2,st_union(st_buffer(validation_plots,60)))

mapview(gridded_plots_3)

max_area<-max(st_area(gridded_plots_3))

gridded_plots_4<-gridded_plots_3[as.numeric(st_area(gridded_plots_3))>as.numeric(max_area)-50,]

library(caret)
#### normalizing variables first

gridded_plots_4$elevation_norm<-(gridded_plots_4$elevatn-mean(gridded_plots_4$elevatn))/sd(gridded_plots_4$elevatn)
gridded_plots_4$ppt_norm<-(gridded_plots_4$nrml_pp-mean(gridded_plots_4$nrml_pp))/sd(gridded_plots_4$nrml_pp)
gridded_plots_4$vs_norm<-(gridded_plots_4$vs-mean(gridded_plots_4$vs))/sd(gridded_plots_4$vs)


ppt<-raster("./processed_data/ppt.tif")

nearest_model<-knnreg(cbi~elevation_norm+ppt_norm+vs_norm,data=gridded_plots_4,k=1)
summary(nearest_model)

predict_df<-st_as_sf(random_point_in_burn)

predict_df$elev<-raster::extract(elev_down,y=predict_df,fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)

predict_df$vs<-raster::extract(vs,y=st_as_sf(predict_df),fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)
predict_df$ppt<-raster::extract(ppt,y=st_as_sf(predict_df),fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)




predict_df$elevation_norm<-(predict_df$elev-mean(gridded_plots_4$elevatn))/sd(gridded_plots_4$elevatn)
predict_df$ppt_norm<-(predict_df$ppt-mean(gridded_plots_4$nrml_pp))/sd(gridded_plots_4$nrml_pp)
predict_df$vs_norm<-(predict_df$vs-mean(gridded_plots_4$vs))/sd(gridded_plots_4$vs)



predict_df$knn_1<-predict(nearest_model,as.data.frame(predict_df))


predict_df<-as.data.frame(predict_df)

predict_df$id_number<-row_number(predict_df)

############################################

cv_data_2<-merge(cv_data,as.data.frame(predict_df),by="id_number")


rmse_by_size<-cv_data_2 %>% group_by(acre_bins) %>% 
  dplyr::summarise(rmse_perimiter=sqrt(mean((actual_burn_severity-perimiter_cbi_value)**2,na.rm=TRUE)),
                   rmse_kriging=sqrt(mean((actual_burn_severity-kriging_cbi)**2,na.rm=TRUE)),
                   rmse_spatial_rf=sqrt(mean((actual_burn_severity-cbi_rf_spatial)**2,na.rm=TRUE)),
                   rmse_knn=sqrt(mean((actual_burn_severity-knn_1)**2,na.rm=TRUE)))

mean_error_size<-cv_data_2 %>% group_by(acre_bins) %>% 
  dplyr::summarise(me_perimiter=mean((actual_burn_severity-perimiter_cbi_value),na.rm=TRUE),
                   me_kriging=mean((actual_burn_severity-kriging_cbi),na.rm=TRUE),
                   me_spatial_rf=mean((actual_burn_severity-cbi_rf_spatial),na.rm=TRUE),
                   me_knn=mean((actual_burn_severity-knn_1),na.rm=TRUE))


rmse_by_size

library(tidyr)

rmse_by_size_2<-data.frame(acre_bins=rep(rmse_by_size$acre_bins,times=4),model=rep(c("perimiter","kriging","spatial_rf","Clustering similarity"),each=7),rmse=c(rmse_by_size$rmse_perimiter,rmse_by_size$rmse_kriging,rmse_by_size$rmse_spatial_rf,rmse_by_size$rmse_knn))

me_by_size_2<-data.frame(acre_bins=rep(rmse_by_size$acre_bins,times=4),model=rep(c("perimiter","kriging","spatial_rf","Clustering similarity"),each=7),mean_error=c(mean_error_size$me_perimiter,mean_error_size$me_kriging,mean_error_size$me_spatial_rf,mean_error_size$me_knn))

rmse_by_size_2 %>% group_by(model) %>% summarise(rmse=mean(rmse))

cb_key<-c("#009E73", "#0072B2","#D55E00","black")

rmse_size<-ggplot(rmse_by_size_2,aes(x=acre_bins,y=rmse,color=model))+geom_point(size=3)+geom_line(lty=2)+theme_classic()+ylab("RMSE")+xlab("Validation Plot Size (acres)")+theme(text=element_text(size=20))+scale_color_manual(values=cb_key)+theme(legend.position = "none")
rmse_size

mean_error_size<-ggplot(me_by_size_2,aes(x=acre_bins,y=mean_error,color=model))+geom_point(size=3)+geom_line(lty=2)+theme_classic()+ylab("Mean Error (actual-modeled)")+xlab("Validation Plot Size (acres)")+theme(text=element_text(size=20))+scale_color_manual(values=cb_key,labels=c("Propensity score matching","Kriging","Perimiter method","Local spatial RF no weather"),name="Method")+theme(legend.position=c(.70,0.85))
mean_error_size

library(cowplot)

tiff(filename=("./figures/size_cross_validation.tif"),units='in',compression='lzw',width=12,height=6,res=300)
plot_grid(rmse_size,mean_error_size,labels=c("a","b"))
dev.off()



cv_data$difference<-cv_data$actual_burn_severity-cv_data$modeled_burn

plot_size_cv<-ggplot(cv_data,aes(y=difference,x=as.factor(round(size_acres))))+geom_boxplot()+theme_classic()+
  geom_hline(yintercept = 0)+ylab("Error (actual - modeled)")+xlab("CV plot size (acres)")+
  theme(text=element_text(size=20))
plot_size_cv

cv_data$acre_bins<-paste(round(cv_data$size_acres,0),"acres",sep="_")

rmse_by_size<-cv_data %>% group_by(acre_bins) %>% 
  dplyr::summarise(rmse=sqrt(mean(difference**2,na.rm=TRUE)))

ggplot(rmse_by_size,aes(x=acre_bins,y=rmse))+geom_point()

tiff(filename=("./figures/plot_size_cv.tif"),units='in',compression='lzw',width=8,height=8,res=300)
plot_size_cv
dev.off()
