'''
To do for updating this analysis. 
Bring in roads to mask the burn area. Possibly use NLCD to make a mask based on more features, e.g., water and buildings.
also calculate closest distnace to road for every object
rerun models to see how the R2 of the RF increases.
explore how the OOB prediction error is effected by size of a treated area - expect that it will decline somewhat
compare my control estimates to Biggs
experiment with artificially changing the resolution

try to scale up to estiamte how much the effect of all the treatments was at the level of the whole fire. Fuel treatments in the burn scar
NET reduced the area of sever burn by XXX and the area of moderate burn by XXX
'''


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

### site potentials are numerical codes but should be considered as factors
extract_df$env_potential<-as.factor(paste("env_",extract_df$env_potential,sep=""))


extract_sf<-sf::st_as_sf(extract_df,coords=c("lat","lon"),crs=crs(veg_treatments))

########### model
#library(SpatialML)
library(ranger)

predict_df<-extract_df[is.na(extract_df$CBI)==TRUE,]
extract_df_full<-extract_df

extract_df<-extract_df[complete.cases(extract_df),]

ranger_rf<-ranger(CBI~elev+aspect+TRI+TPI+slope+lat+lon+road_distance+env_potential,extract_df,importance="impurity")

################### How does our predictive power change over different size validation plots

## what sizes are treatment areas

# 4046.86 m2 in an acre

plot_sizes<-data.frame(acres=c(2,5,15,50,100,150,1000))
plot_sizes$radius_m<-sqrt((plot_sizes$acres*4046.86)/pi)

veg_treatments$actual_size<-st_area(veg_treatments)/4046.86
hist(veg_treatments$actual_size)
quantile(veg_treatments$Acre_US,c(0.05,0.25,0.50,0.75,0.95))

#CBI<-st_make_valid(CBI)
#CBI<-st_union(CBI)

set.seed(207)


random_point_in_burn<-st_sample(burn_perimiter,size=nrow(plot_sizes)*10)
random_point_in_burn<-st_buffer(random_point_in_burn,dist=rep(plot_sizes$radius_m))

random_point_in_burn<-sf::st_transform(random_point_in_burn,crs=st_crs(masked_cbi))


validation_plot<-random_point_in_burn[1]

size_cv_results<-data.frame()

size_cv<-function (index_number){
  ## this function needs to mask the validation area, train the RF model, predict on the validation area, then compare to the actual values in the validation area.
  validation_plot<-random_point_in_burn[index_number]
  #validation_plot<-index_number
  ### actual CBI in validation area
  
  actual_burn<-raster::extract(x=masked_cbi,y=st_as_sf(validation_plot))
  
  if (sum(is.na(actual_burn[[1]]))==0){
  
  #actual_burn_severity<-mean(actual_burn[[1]])
  
  validation_mask_cbi<-raster::mask(masked_cbi,st_as_sf(validation_plot),inverse=TRUE)
  

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
  
  
  
  
  extract_df<-extract_df[complete.cases(extract_df),]
  extract_df$CBI<-as.factor(extract_df$CBI)
  
  ranger_rf<-ranger(CBI~elev+aspect+TRI+TPI+slope+lat+lon+road_distance+env_potential,extract_df,importance="impurity")
  
  predict_df<-data.frame(lat=rasterToPoints(raster::crop(site_potential,y=st_as_sf(validation_plot)))[,1],
                         lon=rasterToPoints(raster::crop(site_potential,y=st_as_sf(validation_plot)))[,2])

  
  predict_df$elev<-raster::extract(elev_down,y=predict_df[,1:2])
  predict_df$aspect<-raster::extract(aspect_down,y=predict_df[,1:2])
  predict_df$TRI<-raster::extract(TRI_down,y=predict_df[,1:2])
  predict_df$TPI<-raster::extract(TPI_down,y=predict_df[,1:2])
  predict_df$slope<-raster::extract(slope_down,y=predict_df[,1:2])
  predict_df$road_distance<-raster::extract(roads_distance,y=predict_df[,1:2])
  predict_df$env_potential<-raster::extract(site_potential,y=predict_df[,1:2])
  
  predict_df$actual_burn<-raster::extract(x=masked_cbi,y=predict_df[,1:2])
  
  predict_df$modeled_values<-predict(ranger_rf,data=predict_df)$predictions
  
  output<-data.frame(actual_burn_severity=mean(predict_df$actual_burn),modeled_burn=mean(as.numeric((as.character(predict_df$modeled_values)))),
                     size_acres=as.numeric(st_area(validation_plot)/4046.86))
  
  
  }
  else(
    output=data.frame(actual_burn_severity=NA,modeled_burn=NA,size_acres=as.numeric(st_area(validation_plot)/4046.86))
    
    )
  size_cv_results<-rbind(size_cv_results,output)
}



#cv_data<-data.frame()
for (i in 1:70){
  cv_data<-rbind(cv_data,size_cv(random_point_in_burn[i]))
  write.csv(cv_data,file="./results/plot_size_cross_validation.csv")
}

#####

cv_data$difference<-cv_data$actual_burn_severity-cv_data$modeled_burn

plot_size_cv<-ggplot(cv_data,aes(y=difference,x=as.factor(round(size_acres))))+geom_boxplot()+theme_classic()+
  geom_hline(yintercept = 0)+ylab("Error (actual - modeled)")+xlab("CV plot size (acres)")+
  theme(text=element_text(size=20))
plot_size_cv

tiff(filename=("./figures/plot_size_cv.tif"),units='in',compression='lzw',width=8,height=8,res=300)
plot_size_cv
dev.off()


##### trying to parall process this part

library(parallel)

cl<-makeCluster(16)
clusterExport(cl,c("random_point_in_burn"))
clusterEvalQ(cl, library("lme4"),library(ranger))

size_cv_results<-lapply(1:length(random_point_in_burn),size_cv)

size_cv_results<-as.data.frame(do.call(rbind, size_cv_results))

#cl<-makeCluster(16)



write.csv(size_cv_results,file="./processed_data/size_cv_results.csv")


test<-read.csv("./processed_data/size_cv_results.csv")

test$diff<-test$actual_burn_severity-test$modeled_burn

ggplot(test,aes(x=as.factor(round(size_acres)),y=diff))+geom_boxplot()

ggplot(test,aes(x=actual_burn_severity,y=modeled_burn,color=as.factor(round(size_acres))))+geom_point()+geom_abline(intercept = 0,slope=1)

summary(lm(actual_burn_severity~modeled_burn,data=test))
