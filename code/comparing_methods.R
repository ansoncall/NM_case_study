#load packages
library(ggplot2)
library(dplyr)
library(sf)
library(raster)
library(stars)
library(elevatr)
library(mapview)
library(ranger)
library(SpatialML)
library(tuneRanger)
library(spmodel)
library(forcats)
library(caret)
# load data

# validation plots are the randomly selected plots that the different methods
# will be evalauted across. It contains both the validation plots and the
# practitioner generated controls
validation_plots<-read_sf("./processed_data/val_points_revised.shp",fid_column_name="ID")

# pairing is the information about how to match the practitioner generated controls
pairing<-read.csv("./processed_data/plot_matching.csv")

# gridded_plots are for the propensity score matching
gridded_plots<-read_sf("gridded_candidate_plots.shp")

masked_cbi<-raster("./processed_data/masked_raster.tif")

masked_cbi[masked_cbi==9]<-NA
masked_cbi[masked_cbi==0]<-1

CBI[CBI==9]<-NA
CBI[CBI==0]<-1

## this is for the practitioner generated controls
plot_method_cbi<-as.data.frame(raster::extract(CBI,st_as_sf(validation_plots),fun=mean,na.rm=FALSE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE))

plot_method_cbi<-data.frame(FID=validation_plots$ID,cbi=plot_method_cbi$V1)
# Nate: double check that the subsetting is correct here. 
validation_cbi<-plot_method_cbi[plot_method_cbi$FID>352,]


control_plots<-merge(plot_method_cbi,pairing,by.x="FID",by.y="FID_1")

names(control_plots)[names(control_plots) == "cbi"] <- "origional_plot_cbi"

paired_plots<-merge(control_plots,validation_cbi,by.x="FID.of.New.plot",by.y="FID")


################# How well did the handmade controls do?

compare_df<-data.frame(plot_id=paired_plots$FID,origional_plot_cbi=paired_plots$origional_plot_cbi,control_plot_cbi=paired_plots$cbi,skip=paired_plots$skip)

compare_df<-compare_df[is.na(compare_df$origional_plot_cbi)==FALSE,]
compare_df<-compare_df[is.na(compare_df$control_plot_cbi)==FALSE,]

plot(compare_df$control_plot_cbi,compare_df$origional_plot_cbi)

origional_plots<-validation_plots[validation_plots$ID %in% compare_df$plot_id,]

############# Simple perimiters 

# "We also used an alternative method that defines the post hoc control as a
# simple 4-hectare buffer around the validation plot, which we call the “buffer
# method” (Vorster et al. 2023). As with the other methods, a 60-m wide buffer
# around the validation plot was removed to minimize edge effects, and a
# 4-hectare area around the excluded area was summarized to estimate the
# counterfactual burn severity (Figure S2). 
 
# Nate: I'm not following along 100% since I'm just reading the code, not
# running it. That said, it looks like you are taking a


# the correct math(?):
# original plot radius (in meters)
og_plot_radius = sqrt(10*4046.86/pi)
# add 60 m buffer
inner_ring = og_plot_radius + 60
# need to know area of original plot + buffer
inner_ring_area = pi*inner_ring^2
inner_ring_area/4046.86 
# okay, supplement S2 shows a 10 acre buffer area, but the text says 4 ha buffer
# area. I'll do the math for both. 

# now calculate the desired outer ring area (inner ring area + ac)
# (choose one)
total_area = inner_ring_area + 40000 # 40,000 m2 == 4 hectares 
# total_area = inner_ring_area + (10*4046.86) # 40,000 m2 == 4 hectares
# now calculate the outer ring radius
outer_ring_radius = sqrt(total_area/pi)
# check final area
outer_ring_area = pi*outer_ring_radius^2
final_area = outer_ring_area - inner_ring_area # checks out

# I think this is just leftover from double checking?
size_of_perimiter<-sqrt((33*4046.86)/pi)-sqrt((10*4046.86)/pi) # not sure where 33 comes from? possibly og 10 acres + 60 m buffer area + new 10 acres?
perimiters<-st_buffer(origional_plots,size_of_perimiter,allow_holes=TRUE) # what is allow_holes doing?
st_area(perimiters)/4046.86
og_plots<-st_buffer(origional_plots,60)
st_area(og_plots)/4046.86


test_outer_ring_area = pi*size_of_perimiter^2



perimiter_extract<-raster::extract(masked_cbi,st_as_sf(perimiters),na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)
origional_extract<-raster::extract(masked_cbi,st_as_sf(og_plots),na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)

compare_df$buffer_method<-NA

# another place to refactor. Can use terra::extract and parallelize. 
for (i in 1:nrow(compare_df)){
  one_plot<-origional_plots[origional_plots$FID_1==compare_df$plot_id[i],]
  
  area_of_plot<-st_area(one_plot)/4046.86 # dividing by 4046, so I assume the units are m2 and you want acres. This should be ==10.
  
  buffer_plot<-st_buffer(one_plot,60) # adding 60m buffer, adding 60m to the radius, buffer plot is now ~94565 m2 or ~23.36 acres in area. 
  area_of_plot_with_buffer<-st_area(buffer_plot)/4046.86 # should be ~23.36 acres. 
  # it looks like size_of_perimeter is the buffer size. to calc this, you want
  # to take the ~23 acre buffer_plot, add 10 acres in area, then take the radius
  # of that area, then subtract the original radius. That would give you the
  # correct parameter for the next st_buffer(). 
  
  # you have: 
  # radius( (area(ogplot) + area(ogplot_w60mbuff)) ) - radius( area(ogplot) )
  # looks correct!
  size_of_perimiter<-sqrt(((area_of_plot+area_of_plot_with_buffer)*4046.86)/pi)-sqrt((area_of_plot)*4046.86/pi) 
  
  perimiters<-st_buffer(one_plot,size_of_perimiter,allow_holes=TRUE) # what is allow_holes doing here?
  area_perimiter<-st_area(perimiters)/4046.86 # should be ~33 acres, unless allow_holes is doing something funky.
  # the following looks to be returning a 1-element vector with CBI mean of the
  # 33 acre plot, not the 10 acre perimeter plot. should leave fun=NULL to
  # return a list of all raster values within the geometry. weights and
  # normalizeWeights just make mean() return an area-weighted pixel
  # mean--addresses edge pixels--but afaik doesn't change the output type.
  perimiter_extract<-raster::extract(masked_cbi,st_as_sf(perimiters),fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE) 
  origional_extract<-raster::extract(masked_cbi,st_as_sf(buffer_plot),fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)
  # You can't recover the perimeter-plot-only mean CBI this way. 
  perimiter_only<-as.numeric(((perimiter_extract[1]*area_perimiter)-(origional_extract[1]*area_of_plot_with_buffer))/(area_perimiter-area_of_plot_with_buffer)) 

  compare_df[compare_df$plot_id==one_plot$FID_1,"buffer_method"]<-perimiter_only
  
}

# imagine the og plot is all 2's and the perimeter plot is all 9's. This matrix
# isn't scaled to the right dimensions, but that doesn't matter for this
# example. Here, the perimeter plot area ==8 and the og plot area ==1.

test_all <- matrix(c(9, 9, 9, 
                     9, 2, 9, 
                     9, 9, 9), 
                   nrow = 3)
# we know the mean of the perimeter plot should be 9. If you take the mean of the whole enchilada, you get ~8.22.
mean_all <- mean(test_all, na.rm=TRUE)
# doesn't matter if you have the plot areas correct, you can't get back to 9 with this kind of math. 
((mean_all * 8) - (2 * 1)) / (8 - 1) 
sum_all <- sum(test_all, na.rm=TRUE)
((sum_all * 8) - (2 * 1)) / (8 - 1) 

sqrt(mean((compare_df$origional_plot_cbi-compare_df$buffer_method)**2,na.rm=TRUE))


#####

### rf model

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

validation_mask_cbi<-raster::mask(masked_cbi,st_as_sf(og_plots),inverse=TRUE)

veg_treatments<-read_sf("./processed_data/vegetation_treatments_hpcc.shp")
veg_treatments<-st_make_valid(veg_treatments)

burn_perimiter<-read_sf("./processed_data/burn_perimiter.shp")

## make empty raster to extract data into
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
extract_df$ppt<-raster::extract(ppt,y=as.matrix(extract_df[,1:2]))
extract_df$tmin<-raster::extract(tmin,y=as.matrix(extract_df[,1:2]))
extract_df$tmmx<-raster::extract(tmmx,y=as.matrix(extract_df[,1:2]))
extract_df$th<-raster::extract(th,y=as.matrix(extract_df[,1:2]))
extract_df$vs<-raster::extract(vs,y=as.matrix(extract_df[,1:2]))
extract_df$vpdmax<-raster::extract(vpdmax,y=as.matrix(extract_df[,1:2]))
extract_df$fm100<-raster::extract(fm100,y=as.matrix(extract_df[,1:2]))
extract_df$fm1000<-raster::extract(fm1000,y=as.matrix(extract_df[,1:2]))
extract_df$rmax<-raster::extract(rmax,y=as.matrix(extract_df[,1:2]))

extract_df<-extract_df[extract_df$CBI!=0,]
extract_df<-extract_df[extract_df$CBI!=9,]

### site potentials are numerical codes but should be considered as factors
extract_df$env_potential<-as.factor(paste("env_",extract_df$env_potential,sep=""))


predict_df<-extract_df[is.na(extract_df$CBI)==TRUE,]
extract_df_full<-extract_df

extract_df<-extract_df[complete.cases(extract_df),]

#extract_df$CBI<-as.factor(extract_df$CBI)

#cbi.task = makeRegrTask(data = extract_df, target = "CBI")
#estimateTimeTuneRanger(cbi.task)

#tuning<-tuneRanger(cbi.task)
#tuning

ranger_rf<-ranger(CBI~elev+aspect+TRI+TPI+slope+lat+lon+road_distance+env_potential+ppt+tmin+tmmx+th+vs+vpdmax+fm100+fm1000+rmax,extract_df,num.trees=1000,importance="impurity",mtry=5,min.node.size=2,sample.fraction=0.89)
ranger_rf

### variable importance plot

var_imp_plot<-data.frame(importance=ranger_rf$variable.importance)
var_imp_plot$variable<-row.names(var_imp_plot)
var_imp_plot$category<-c("landscape","landscape","landscape","landscape","landscape","landscape","landscape","landscape","bioclimatic","bioclimatic","bioclimatic","weather","weather","weather","bioclimatic","weather","weather","weather")

var_imp_plot$variable<-fct_reorder(var_imp_plot$variable,var_imp_plot$importance)
cbPalette <- c( "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

var_imp_plot$variable <- fct_recode(var_imp_plot$variable,"TRI"="TRI","TPI"="TPI","aspect"="aspect","slope"="slope","elev" = "elevation", "lat" = "lattitude","lon"="longitude","road_distance"="road distance","ppt"="precipitation","tmin"="minimum temperature","tmmx"="maximum temperature","th"="wind direction","vs"="wind speed","vpdmax"="maximum vapor pressure deficit","fm100"="moisture 100-hr","fm1000"="moisture 1000-hr","rmax"="maximum relative humidity","env_potential"="environmental potential")

variable_importance_plot<-ggplot(var_imp_plot,aes(x=importance,y=variable,fill=category))+geom_bar(stat="identity")+theme_classic()+theme(text=element_text(size=15))+xlab("Importance (Gini index)")+ylab("")+scale_fill_manual(values=cbPalette)+scale_y_discrete(labels=c("TPI","TRI","slope","wind direction","aspect","maximum temperature","moisture 100-hr fuel","maximum relative humidity","moisture 1000-hr fuel","wind speed","minimum temperature","maximum vpd","road distance","lattitude","environmental potential","precipitation","elevation","longitude"))
variable_importance_plot

tiff(filename=("./figures/rf_importance.tif"),units='in',compression='lzw',width=10,height=8,res=300)
variable_importance_plot
dev.off()

###############

coordinates<- data.frame(lat=extract_df$lat,lon=extract_df$lon)

extract_sf<-st_as_sf(extract_df,coords=c("lat","lon"))
extract_sf$lat<-extract_df$lat
extract_sf$lon<-extract_df$lon

spatial_rf_model<-splmRF(CBI~elev+aspect+TRI+TPI+slope+lat+lon+road_distance+env_potential+ppt+tmin+tmmx+th+vs+vpdmax+fm100+fm1000+rmax,data=extract_sf,spcov_type = "gravity",local=c(parallel=TRUE,ncores=16),mtry=4,min.node.size=2,sample.fraction=0.89)


## predicting for one plot at a time
predict_for_each_plot<-function(location){
  
  predict_df<-data.frame(lat=rasterToPoints(raster::crop(site_potential,y=st_as_sf(location)))[,1],
                         lon=rasterToPoints(raster::crop(site_potential,y=st_as_sf(location)))[,2])
  
  predict_df$elev<-raster::extract(elev_down,y=predict_df[,1:2])
  predict_df$aspect<-raster::extract(aspect_down,y=predict_df[,1:2])
  predict_df$TRI<-raster::extract(TRI_down,y=predict_df[,1:2])
  predict_df$TPI<-raster::extract(TPI_down,y=predict_df[,1:2])
  predict_df$slope<-raster::extract(slope_down,y=predict_df[,1:2])
  predict_df$road_distance<-raster::extract(roads_distance,y=predict_df[,1:2])
  predict_df$env_potential<-as.factor(raster::extract(site_potential,y=predict_df[,1:2]))
  predict_df$ppt<-raster::extract(ppt,y=as.matrix(predict_df[,1:2]))
  predict_df$tmin<-raster::extract(tmin,y=as.matrix(predict_df[,1:2]))
  predict_df$tmmx<-raster::extract(tmmx,y=as.matrix(predict_df[,1:2]))
  predict_df$th<-raster::extract(th,y=as.matrix(predict_df[,1:2]))
  predict_df$vs<-raster::extract(vs,y=as.matrix(predict_df[,1:2]))
  predict_df$vpdmax<-raster::extract(vpdmax,y=as.matrix(predict_df[,1:2]))
  predict_df$fm100<-raster::extract(fm100,y=as.matrix(predict_df[,1:2]))
  predict_df$fm1000<-raster::extract(fm1000,y=as.matrix(predict_df[,1:2]))
  predict_df$rmax<-raster::extract(rmax,y=as.matrix(predict_df[,1:2]))
  
  
  predict_df<-predict_df[complete.cases(predict_df),]
  
  #predict_df<-st_as_sf(predict_df,coords=c("lat","lon"),crs=st_crs(site_potential))
  
  #predict_df$lat<-rasterToPoints(raster::crop(site_potential,y=st_as_sf(location)))[,1]
  #predict_df$lon=rasterToPoints(raster::crop(site_potential,y=st_as_sf(location)))[,2]
  
  predict_df$modeled_values<-predict(ranger_rf,data=predict_df,na.rm=TRUE)$predictions
  #predict_df$modeled_values<-predict(spatial_rf_model,newdata=predict_df)
  
  for_conversion<-data.frame(lat=predict_df$lat,lon=predict_df$lon,modeled_values=as.numeric(predict_df$modeled_values))
  
  predict_raster<-rasterFromXYZ(for_conversion,crs=crs(masked_cbi))
  
  mean_cbi<-raster::extract(predict_raster$modeled_values,st_as_sf(location),fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)
  

  return(mean_cbi)
}



predict_for_each_plot_spatial<-function(location){
  
  predict_df<-data.frame(lat=rasterToPoints(raster::crop(site_potential,y=st_as_sf(location)))[,1],
                         lon=rasterToPoints(raster::crop(site_potential,y=st_as_sf(location)))[,2])
  
  predict_df$elev<-raster::extract(elev_down,y=predict_df[,1:2])
  predict_df$aspect<-raster::extract(aspect_down,y=predict_df[,1:2])
  predict_df$TRI<-raster::extract(TRI_down,y=predict_df[,1:2])
  predict_df$TPI<-raster::extract(TPI_down,y=predict_df[,1:2])
  predict_df$slope<-raster::extract(slope_down,y=predict_df[,1:2])
  predict_df$road_distance<-raster::extract(roads_distance,y=predict_df[,1:2])
  predict_df$env_potential<-as.factor(raster::extract(site_potential,y=predict_df[,1:2]))
  predict_df$ppt<-raster::extract(ppt,y=as.matrix(predict_df[,1:2]))
  predict_df$tmin<-raster::extract(tmin,y=as.matrix(predict_df[,1:2]))
  predict_df$tmmx<-raster::extract(tmmx,y=as.matrix(predict_df[,1:2]))
  predict_df$th<-raster::extract(th,y=as.matrix(predict_df[,1:2]))
  predict_df$vs<-raster::extract(vs,y=as.matrix(predict_df[,1:2]))
  predict_df$vpdmax<-raster::extract(vpdmax,y=as.matrix(predict_df[,1:2]))
  predict_df$fm100<-raster::extract(fm100,y=as.matrix(predict_df[,1:2]))
  predict_df$fm1000<-raster::extract(fm1000,y=as.matrix(predict_df[,1:2]))
  predict_df$rmax<-raster::extract(rmax,y=as.matrix(predict_df[,1:2]))
  
  predict_df$lat2<-predict_df$lat
  predict_df$lon2<-predict_df$lon
  
  
  predict_df<-predict_df[complete.cases(predict_df),]
    
  predict_df<-st_as_sf(predict_df,coords=c("lat2","lon2"),crs=st_crs(site_potential))

  #predict_df$modeled_values<-predict(ranger_rf,data=predict_df,na.rm=TRUE)$predictions
  
  
  
  predict_df$modeled_values<-predict(spatial_rf_model,newdata=predict_df)
  
  for_conversion<-data.frame(lat=predict_df$lat,lon=predict_df$lon,modeled_values=as.numeric(as.character(predict_df$modeled_values)))
  
  predict_raster<-rasterFromXYZ(for_conversion,crs=crs(masked_cbi))
  
  mean_cbi<-raster::extract(predict_raster$modeled_values,st_as_sf(location),fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)
  
  
  return(mean_cbi)
}



compare_df$spatial_rf_prediction<-NA
compare_df$rf_prediction<-NA
for (i in 1:nrow(compare_df)){
  print(i/nrow(compare_df)*100)
  one_plot<-origional_plots[origional_plots$FID_1==compare_df$plot_id[i],]
  compare_df$spatial_rf_prediction[i]=predict_for_each_plot_spatial(one_plot)
  compare_df$rf_prediction[i]=predict_for_each_plot(one_plot)
  }

## RMSE
sqrt(mean((compare_df$rf_prediction-compare_df$origional_plot_cbi)**2,na.rm=TRUE))
sqrt(mean((compare_df$spatial_rf_prediction-compare_df$origional_plot_cbi)**2,na.rm=TRUE))

write.csv(compare_df,"./processed_data/larger_rf_model_predictions.csv")

##### Kriging

library(snapKrig)

compare_df$krigged_values<-NA

for (i in 1:nrow(compare_df)){
#for (i in 1:5){
    
  
one_plot<-origional_plots[origional_plots$FID_1==compare_df$plot_id[i],]

size_of_perimiter<-sqrt((200*4046.86)/pi)

one_plot_buffer<-st_buffer(one_plot,size_of_perimiter)

clipped_burn_raster<-raster::crop(masked_cbi,one_plot_buffer)

clipped_burn_raster<-raster::mask(clipped_burn_raster,st_as_sf(st_buffer(one_plot,60)),inverse=TRUE)

clipped_burn_raster[clipped_burn_raster==0]<-1

clipped_burn_raster[clipped_burn_raster==9]<-NA

### sk_fit fails when corner cells are NA. This sets it to 1 when it is NA. Because corners are the fartest from the treatments this should have minimal impact on predictions
if (is.na(clipped_burn_raster[1,1])==TRUE){clipped_burn_raster[1,1]<-1}
if (is.na(clipped_burn_raster[nrow(clipped_burn_raster),1])==TRUE){clipped_burn_raster[nrow(clipped_burn_raster),1]<-1}
if (is.na(clipped_burn_raster[1,ncol(clipped_burn_raster)])==TRUE){clipped_burn_raster[1,ncol(clipped_burn_raster)]<-1}
if (is.na(clipped_burn_raster[nrow(clipped_burn_raster),ncol(clipped_burn_raster)])==TRUE){clipped_burn_raster[nrow(clipped_burn_raster),ncol(clipped_burn_raster)]<-1}


clipped_burn_raster_sk<-sk(clipped_burn_raster)



kriging_fit<-sk_fit(clipped_burn_raster_sk,n_max=3000)


krigged_final<-sk_cmean(clipped_burn_raster_sk,kriging_fit)
plot(krigged_final)

krigged_final_raster<-sk_export(krigged_final)

values<-extract(krigged_final_raster,one_plot)

compare_df$krigged_values[i]<-mean(values$lyr.1)


}


sqrt(mean((compare_df$krigged_values-compare_df$origional_plot_cbi)**2,na.rm=TRUE))
write.csv(compare_df,"./processed_data/kriging_results.csv")

##### spatial rf 

local_spatial_rf<-c()

for (i in 1:nrow(compare_df)){

  one_plot<-origional_plots[origional_plots$FID_1==compare_df$plot_id[i],]
  
  size_of_perimiter<-sqrt((1200*4046.86)/pi)
  
  one_plot_buffer<-st_buffer(one_plot,size_of_perimiter)
  
  clipped_burn_raster<-raster::crop(validation_mask_cbi,st_as_sf(one_plot_buffer))
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
  training_df$ppt<-raster::extract(ppt,y=as.matrix(training_df[,1:2]))
  training_df$tmin<-raster::extract(tmin,y=as.matrix(training_df[,1:2]))
  training_df$tmmx<-raster::extract(tmmx,y=as.matrix(training_df[,1:2]))
  training_df$vpdmax<-raster::extract(vpdmax,y=as.matrix(training_df[,1:2]))

  training_df$lat<-training_df$x
  training_df$lon<-training_df$y
  
  training_df<-training_df[complete.cases(training_df),]
  
  training_df<-st_as_sf(training_df,coords=c("x","y"),crs=st_crs(site_potential))
  


  
  spatial_rf_model<-splmRF(CBI~elev+aspect+TRI+TPI+slope+road_distance+env_potential+lat+lon+ppt+tmin+vpdmax,data=training_df,spcov_type = "exponential",local=c(parallel=TRUE,ncores=8),mtry=4,min.node.size=2,sample.fraction=0.89)
  
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
  predict_df$ppt<-raster::extract(ppt,y=as.matrix(predict_df[,1:2]))
  predict_df$tmin<-raster::extract(tmin,y=as.matrix(predict_df[,1:2]))
  predict_df$tmmx<-raster::extract(tmmx,y=as.matrix(predict_df[,1:2]))
  predict_df$vpdmax<-raster::extract(vpdmax,y=as.matrix(predict_df[,1:2]))

  g<-predict_df
  
  predict_df$lat2<-predict_df$lat
  predict_df$lon2<-predict_df$lon
  
  predict_df<-predict_df[complete.cases(predict_df),]
  
  predict_df<-st_as_sf(predict_df,coords=c("lat2","lon2"),crs=st_crs(site_potential))
  

  #predict_df$modeled_values<-predict(ranger_rf,data=predict_df,na.rm=TRUE)$predictions
  predict_df$modeled_values<-predict(spatial_rf_model,newdata=predict_df)
  
  for_conversion<-data.frame(lat=predict_df$lat,lon=predict_df$lon,modeled_values=as.numeric(predict_df$modeled_values))
  
  predict_raster<-rasterFromXYZ(for_conversion,crs=crs(masked_cbi))
  
  mean_cbi<-raster::extract(predict_raster$modeled_values,st_as_sf(one_plot),fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)
  
  local_spatial_rf[i]<-mean_cbi
  
  print(i)
}


compare_df$local_spatial_rf_no_weather<-local_spatial_rf

write.csv(compare_df,"./processed_data/with_local_spatial_rf_no_weather.csv")

############# Adding weather to spatial RF

local_spatial_rf_weather<-data.frame()

for (i in 1:nrow(compare_df)){

  one_plot<-origional_plots[origional_plots$FID_1==compare_df$plot_id[i],]
  
  size_of_perimiter<-sqrt((1200*4046.86)/pi)
  
  one_plot_buffer<-st_buffer(one_plot,size_of_perimiter)
  
  clipped_burn_raster<-raster::crop(validation_mask_cbi,st_as_sf(one_plot_buffer))
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
  training_df$ppt<-raster::extract(ppt,y=as.matrix(training_df[,1:2]))
  training_df$tmin<-raster::extract(tmin,y=as.matrix(training_df[,1:2]))
  training_df$tmmx<-raster::extract(tmmx,y=as.matrix(training_df[,1:2]))
  training_df$th<-raster::extract(th,y=as.matrix(training_df[,1:2]))
  training_df$vs<-raster::extract(vs,y=as.matrix(training_df[,1:2]))
  training_df$vpdmax<-raster::extract(vpdmax,y=as.matrix(training_df[,1:2]))
  training_df$fm100<-raster::extract(fm100,y=as.matrix(training_df[,1:2]))
  training_df$fm1000<-raster::extract(fm1000,y=as.matrix(training_df[,1:2]))
  training_df$rmax<-raster::extract(rmax,y=as.matrix(training_df[,1:2]))
  
  
  
  
  training_df$lat<-as.data.frame(rasterToPoints(clipped_burn_raster))[,1]
  training_df$lon<-as.data.frame(rasterToPoints(clipped_burn_raster))[,2]
  
  training_df<-training_df[complete.cases(training_df),]

  training_df<-st_as_sf(training_df,coords=c("x","y"),crs=st_crs(site_potential))
  
    
  spatial_rf_model<-splmRF(CBI~elev+aspect+TRI+TPI+slope+road_distance+env_potential+lat+lon+ppt+tmin+tmmx+th+vs+vpdmax+fm100+fm1000+rmax,data=training_df,spcov_type = "exponential",local=c(parallel=TRUE,ncores=15),mtry=4,min.node.size=2,sample.fraction=0.89)
  
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
  predict_df$ppt<-raster::extract(ppt,y=as.matrix(predict_df[,1:2]))
  predict_df$tmin<-raster::extract(tmin,y=as.matrix(predict_df[,1:2]))
  predict_df$tmmx<-raster::extract(tmmx,y=as.matrix(predict_df[,1:2]))
  predict_df$th<-raster::extract(th,y=as.matrix(predict_df[,1:2]))
  predict_df$vs<-raster::extract(vs,y=as.matrix(predict_df[,1:2]))
  predict_df$vpdmax<-raster::extract(vpdmax,y=as.matrix(predict_df[,1:2]))
  predict_df$fm100<-raster::extract(fm100,y=as.matrix(predict_df[,1:2]))
  predict_df$fm1000<-raster::extract(fm1000,y=as.matrix(predict_df[,1:2]))
  predict_df$rmax<-raster::extract(rmax,y=as.matrix(predict_df[,1:2]))
  
  predict_df<-predict_df[complete.cases(predict_df),]

  predict_df$lat2<-predict_df$lat
  predict_df$lon2<-predict_df$lon
  
    
  predict_df<-st_as_sf(predict_df,coords=c("lat2","lon2"),crs=st_crs(site_potential))
  

  
  #predict_df$modeled_values<-predict(ranger_rf,data=predict_df,na.rm=TRUE)$predictions
  predict_df$modeled_values<-predict(spatial_rf_model,newdata=predict_df)
  
  for_conversion<-data.frame(lat=predict_df$lat,lon=predict_df$lon,modeled_values=as.numeric(predict_df$modeled_values))
  
  predict_raster<-rasterFromXYZ(for_conversion,crs=crs(masked_cbi))
  
  mean_cbi<-raster::extract(predict_raster$modeled_values,st_as_sf(one_plot),fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)
  
  output<-data.frame(plot_id=one_plot$FID_1,local_spatial_rf_weather=mean_cbi)
  
  local_spatial_rf_weather<-rbind(local_spatial_rf_weather,output)
  
  print(i)
  
}

compare_df<-merge(compare_df,local_spatial_rf_weather,by="plot_id")

compare_df$local_spatial_rf_with_weather<-compare_df$local_spatial_rf_weather

write.csv(compare_df,"./processed_data/both_local_spatials.csv")

####################### clustering based matching 

## remove the treated and validation plots
veg_treatments<-read_sf("./processed_data/vegetation_treatments_hpcc_new.shp")

gridded_plots<-st_transform(gridded_plots,st_crs(CBI))
#3control_plots<-st_transform(control_plots,st_crs(CBI))

gridded_plots_2<-st_difference(gridded_plots,st_union(veg_treatments))


gridded_plots_3<-st_difference(gridded_plots_2,st_union(st_buffer(validation_plots,60)))

mapview(gridded_plots_3)

max_area<-max(st_area(gridded_plots_3))

gridded_plots_4<-gridded_plots_3[as.numeric(st_area(gridded_plots_3))>as.numeric(max_area)-50,]


#### normalizing variables first

gridded_plots_4$elevation_norm<-(gridded_plots_4$elevatn-mean(gridded_plots_4$elevatn))/sd(gridded_plots_4$elevatn)
gridded_plots_4$ppt_norm<-(gridded_plots_4$nrml_pp-mean(gridded_plots_4$nrml_pp))/sd(gridded_plots_4$nrml_pp)
gridded_plots_4$vs_norm<-(gridded_plots_4$vs-mean(gridded_plots_4$vs))/sd(gridded_plots_4$vs)

nearest_model<-knnreg(cbi~elevation_norm+ppt_norm+vs_norm,data=gridded_plots_4,k=1)
summary(nearest_model)

predict_df<-origional_plots[origional_plots$FID_1 %in% compare_df$plot_id[1:247],]

predict_df$elev<-raster::extract(elev_down,y=predict_df,fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)

predict_df$vs<-raster::extract(vs,y=predict_df,fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)
predict_df$ppt<-raster::extract(ppt,y=predict_df,fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)

predict_df$elevation_norm<-(predict_df$elev-mean(gridded_plots_4$elevatn))/sd(gridded_plots_4$elevatn)
predict_df$ppt_norm<-(predict_df$ppt-mean(gridded_plots_4$nrml_pp))/sd(gridded_plots_4$nrml_pp)
predict_df$vs_norm<-(predict_df$vs-mean(gridded_plots_4$vs))/sd(gridded_plots_4$vs)

predict_df$knn_1<-predict(nearest_model,as.data.frame(predict_df))

knn_output<-data.frame(plot_id=predict_df$FID_1,knn=predict_df$knn_1)

compare_df<-merge(compare_df,knn_output,by="plot_id")

#
traditional_controls_rmse<-c()
traditional_controls_mean_diff<-c()
#
perimiter_rmse<-c()
perimiter_mean_diff<-c()
#
rf_rmse<-c()
rf_mean_diff<-c()
#
sp_rf_rmse<-c()
sp_rf_mean_diff<-c()
#
krig_rmse<-c()
krig_mean_diff<-c()
#
local_spatial_rf_rmse<-c()
local_spatial_rf_mean_diff<-c()
#
local_spatial_rf_weather_rmse<-c()
local_spatial_rf_weather_mean_diff<-c()
#
knn_rmse<-c()
knn_mean_diff<-c()

for (i in 1:5000){
  
  compare_df_sample<-compare_df[sample(nrow(compare_df), nrow(compare_df),replace = TRUE), ]
  
  traditional_controls_rmse[i]<-sqrt(mean((compare_df_sample$control_plot_cbi-compare_df_sample$origional_plot_cbi)**2,na.rm=TRUE))
  traditional_controls_mean_diff[i]<-mean(compare_df_sample$control_plot_cbi-compare_df_sample$origional_plot_cbi,na.rm=TRUE)
  
  perimiter_rmse[i]<-sqrt(mean((compare_df_sample$buffer_method-compare_df_sample$origional_plot_cbi)**2,na.rm=TRUE))
  perimiter_mean_diff[i]<-mean(compare_df_sample$buffer_method-compare_df_sample$origional_plot_cbi,na.rm=TRUE)
  
  rf_rmse[i]<-sqrt(mean((compare_df_sample$rf_prediction-compare_df_sample$origional_plot_cbi)**2,na.rm=TRUE))
  rf_mean_diff[i]<-mean(compare_df_sample$rf_prediction-compare_df_sample$origional_plot_cbi,na.rm=TRUE)
  
  sp_rf_rmse[i]<-sqrt(mean((compare_df_sample$spatial_rf_prediction-compare_df_sample$origional_plot_cbi)**2,na.rm=TRUE))
  sp_rf_mean_diff[i]<-mean(compare_df_sample$spatial_rf_prediction-compare_df_sample$origional_plot_cbi,na.rm=TRUE)
  
  krig_rmse[i]<-sqrt(mean((compare_df_sample$krigged_values-compare_df_sample$origional_plot_cbi)**2,na.rm=TRUE))
  krig_mean_diff[i]<-mean(compare_df_sample$krigged_values-compare_df_sample$origional_plot_cbi,na.rm=TRUE)
  
  local_spatial_rf_rmse[i]<-sqrt(mean((compare_df_sample$local_spatial_rf_no_weather-compare_df_sample$origional_plot_cbi)**2,na.rm=TRUE))
  local_spatial_rf_mean_diff[i]<-mean(compare_df_sample$local_spatial_rf_no_weather-compare_df_sample$origional_plot_cbi,na.rm=TRUE)
  
  local_spatial_rf_weather_rmse[i]<-sqrt(mean((compare_df_sample$local_spatial_rf_with_weather.local_spatial_rf_weather-compare_df_sample$origional_plot_cbi)**2,na.rm=TRUE))
  local_spatial_rf_weather_mean_diff[i]<-mean(compare_df_sample$local_spatial_rf_with_weather.local_spatial_rf_weather-compare_df_sample$origional_plot_cbi,na.rm=TRUE)
  
  knn_rmse[i]<-sqrt(mean((compare_df_sample$knn-compare_df_sample$origional_plot_cbi)**2,na.rm=TRUE))
  knn_mean_diff[i]<-mean(compare_df_sample$knn-compare_df_sample$origional_plot_cbi,na.rm=TRUE)
  
  
}



graphing_comparison<-data.frame(model_type=c("practitioner selected plots","perimeter method","random forest","spatial random forest","kriging","local spatial random forest without weather","local spatial random forest","propensity score matching"),mean_rmse=c(mean(traditional_controls_rmse),mean(perimiter_rmse),mean(rf_rmse,na.rm=TRUE),mean(sp_rf_rmse,na.rm=TRUE),mean(krig_rmse),mean(local_spatial_rf_rmse),mean(local_spatial_rf_weather_rmse),mean(knn_rmse)),lower_ci=c(quantile(traditional_controls_rmse,0.025),quantile(perimiter_rmse,0.025),quantile(rf_rmse,0.025,na.rm=TRUE),quantile(sp_rf_rmse,0.025,na.rm=TRUE),quantile(krig_rmse,0.025),quantile(local_spatial_rf_rmse,0.025,na.rm=TRUE),quantile(local_spatial_rf_weather_rmse,0.025,na.rm=TRUE),quantile(knn_rmse,0.025)),upper_ci=c(quantile(traditional_controls_rmse,0.975),quantile(perimiter_rmse,0.975),quantile(rf_rmse,0.975,na.rm=TRUE),quantile(sp_rf_rmse,0.975,na.rm=TRUE),quantile(krig_rmse,0.975),quantile(local_spatial_rf_rmse,0.975),quantile(local_spatial_rf_weather_rmse,0.975,na.rm=TRUE),quantile(knn_rmse,0.975)),mean_difference=c(mean(traditional_controls_mean_diff),mean(perimiter_mean_diff),mean(rf_mean_diff,na.rm=TRUE),mean(sp_rf_mean_diff,na.rm=TRUE),mean(krig_mean_diff),mean(local_spatial_rf_mean_diff),mean(local_spatial_rf_weather_mean_diff),mean(knn_mean_diff)),mean_diff_lower_ci=c(quantile(traditional_controls_mean_diff,0.025),quantile(perimiter_mean_diff,0.025),quantile(rf_mean_diff,0.025,na.rm=TRUE),quantile(sp_rf_mean_diff,0.025,na.rm=TRUE),quantile(krig_mean_diff,0.025),quantile(local_spatial_rf_mean_diff,0.025),quantile(local_spatial_rf_weather_mean_diff,0.025),quantile(knn_mean_diff,0.025)),mean_diff_upper_ci=c(quantile(traditional_controls_mean_diff,0.975),quantile(perimiter_mean_diff,0.975),quantile(rf_mean_diff,0.975,na.rm=TRUE),quantile(sp_rf_mean_diff,0.975,na.rm=TRUE),quantile(krig_mean_diff,0.975,na.rm=TRUE),quantile(local_spatial_rf_mean_diff,0.975,na.rm=TRUE),quantile(local_spatial_rf_weather_mean_diff,0.975,na.rm=TRUE),quantile(knn_mean_diff,0.975,na.rm=TRUE)))

graphing_comparison$model_type<-factor(graphing_comparison$model_type,levels=c("propensity score matching","random forest","spatial random forest","kriging","practitioner selected plots","perimeter method","local spatial random forest","local spatial random forest without weather"))

library(stringr)


rmse_plot<-ggplot(graphing_comparison,aes(x=model_type,y=mean_rmse))+geom_bar(stat="identity")+theme_classic()+geom_errorbar(aes(ymin=lower_ci,ymax=upper_ci),width=0.1)+ylab("RMSE of CBI")+xlab("")+theme(text=element_text(size=13))+scale_x_discrete(labels = function(model_type) str_wrap(model_type, width = 6))
rmse_plot


mean_difference_plot<-ggplot(graphing_comparison,aes(x=model_type,y=mean_difference))+geom_bar(stat="identity")+theme_classic()+geom_errorbar(aes(ymin=mean_diff_lower_ci,ymax=mean_diff_upper_ci),width=0.1)+ylab("Mean Error (actual-modeled)")+xlab("")+theme(text=element_text(size=13))+scale_x_discrete(labels = function(model_type) str_wrap(model_type, width = 6))
mean_difference_plot


library(cowplot)

write.csv(graphing_comparison,"./results/model_comparisons.csv")

write.csv(compare_df,"./processed_data/final_comparsion_results_new.csv")

tiff(filename=("./figures/comparison_of_methods.tif"),units='in',compression='lzw',width=7,height=10,res=300)
plot_grid(rmse_plot,mean_difference_plot,ncol = 1,labels=c("a","b"))
dev.off()


################################

gd<-read.csv("./processed_data/final_comparsion_results_new.csv")

model<-lm(origional_plot_cbi~knn,data=gd)
knn_plot<-ggplot(gd,aes(y=origional_plot_cbi,x=knn))+geom_abline(intercept=0,slope=1,size=2,color="gray45")+geom_smooth(method="lm")+geom_point()+theme_classic()+ylab("Validaiton plot burn severity (CBI)")+xlab("Propensity score matching (CBI)")+annotate("text",x=3.5,y=1.5,label=paste("Y=",round(model$coefficients[2],2),"* X +",round(model$coefficients[1],2)),size=4)+xlim(1,4)


model<-lm(origional_plot_cbi~local_spatial_rf_no_weather,data=gd)
lsrfnw_plot<-ggplot(gd,aes(y=origional_plot_cbi,x=local_spatial_rf_no_weather))+geom_abline(intercept=0,slope=1,size=2,color="gray45")+geom_smooth(method="lm")+geom_point()+theme_classic()+ylab("Validaiton plot burn severity (CBI)")+xlab("Local spatial RF without weather (CBI)")+annotate("text",x=3.5,y=1.5,label=paste("Y=",round(model$coefficients[2],2),"* X +",round(model$coefficients[1],2)),size=4)+xlim(1,4)

model<-lm(origional_plot_cbi~local_spatial_rf_with_weather.local_spatial_rf_weather,data=gd)
lsrf_plot<-ggplot(gd,aes(y=origional_plot_cbi,x=local_spatial_rf_with_weather.local_spatial_rf_weather))+geom_abline(intercept=0,slope=1,size=2,color="gray45")+geom_smooth(method="lm")+geom_point()+theme_classic()+ylab("Validaiton plot burn severity (CBI)")+xlab("Local spatial RF (CBI)")+annotate("text",x=3.5,y=1.5,label=paste("Y=",round(model$coefficients[2],2),"* X +",round(model$coefficients[1],2)),size=4)+xlim(1,4)

model<-lm(origional_plot_cbi~control_plot_cbi,data=gd)
prac_sel_plot<-ggplot(gd,aes(y=origional_plot_cbi,x=control_plot_cbi))+geom_abline(intercept=0,slope=1,size=2,color="gray45")+geom_smooth(method="lm")+geom_point()+theme_classic()+ylab("Validaiton plot burn severity (CBI)")+xlab("Practioner selected plots (CBI)")+annotate("text",x=3.5,y=1.5,label=paste("Y=",round(model$coefficients[2],2),"* X +",round(model$coefficients[1],2)),size=4)+xlim(1,4)

model<-lm(origional_plot_cbi~buffer_method,data=gd)
buffer_plot<-ggplot(gd,aes(y=origional_plot_cbi,x=buffer_method))+geom_abline(intercept=0,slope=1,size=2,color="gray45")+geom_smooth(method="lm")+geom_point()+theme_classic()+ylab("Validaiton plot burn severity (CBI)")+xlab("Perimeters (CBI)")+annotate("text",x=3.5,y=1.5,label=paste("Y=",round(model$coefficients[2],2),"* X +",round(model$coefficients[1],2)),size=4)+xlim(1,4)

model<-lm(origional_plot_cbi~krigged_values,data=gd)
krigged_values<-ggplot(gd,aes(y=origional_plot_cbi,x=krigged_values))+geom_abline(intercept=0,slope=1,size=2,color="gray45")+geom_smooth(method="lm")+geom_point()+theme_classic()+ylab("Validaiton plot burn severity (CBI)")+xlab("Kriging (CBI)")+annotate("text",x=3.5,y=1.5,label=paste("Y=",round(model$coefficients[2],2),"* X +",round(model$coefficients[1],2)),size=4)+xlim(1,4)

model<-lm(origional_plot_cbi~rf_prediction,data=gd)
rf_plot<-ggplot(gd,aes(y=origional_plot_cbi,x=rf_prediction))+geom_abline(intercept=0,slope=1,size=2,color="gray45")+geom_smooth(method="lm")+geom_point()+theme_classic()+ylab("Validaiton plot burn severity (CBI)")+xlab("RF model (CBI)")+annotate("text",x=3.5,y=1.5,label=paste("Y=",round(model$coefficients[2],2),"* X +",round(model$coefficients[1],2)),size=4)+xlim(1,4)

model<-lm(origional_plot_cbi~spatial_rf_prediction,data=gd)
sprf_plot<-ggplot(gd,aes(y=origional_plot_cbi,x=spatial_rf_prediction))+geom_abline(intercept=0,slope=1,size=2,color="gray45")+geom_smooth(method="lm")+geom_point()+theme_classic()+ylab("Validaiton plot burn severity (CBI)")+xlab("Spatial RF model (CBI)")+annotate("text",x=3.5,y=1.5,label=paste("Y=",round(model$coefficients[2],2),"* X +",round(model$coefficients[1],2)),size=4)+xlim(1,4)

tiff(filename=("./figures/scatter_plots_of_predictions.tif"),units='in',compression='lzw',width=8,height=12,res=300)
plot_grid(prac_sel_plot,buffer_plot,knn_plot,rf_plot,sprf_plot,krigged_values,lsrf_plot,lsrfnw_plot,ncol = 2,labels="auto",label_x=0.15,label_y=0.9)
dev.off()
