# this script tests the effect of different validation plot sizes on prediction
# accuracy.

# load packages ####
# library(raster)
library(terra)
library(sf)
library(stars)
library(ggplot2)
library(ranger)
library(spmodel)
library(snapKrig)
library(spmodel)
library(parallel)

# load data ####
## composite burn index (CBI) raster ####
cbi <- rast("./processed_data/masked_raster.tif")
# as usual, set "unmappable" -> NA and "outside of perimeter" -> 1 (unburned)
names(cbi) <- "cbi"
cbi[cbi == 9] <- NA
cbi[cbi == 0] <- 1

## HPCC burn perimeter ####
burn_perimeter <- read_sf("./processed_data/burn_perimeter.shp")

## treatments ####
veg_treatments <- read_sf(
  "./processed_data/vegetation_treatments_hpcc_new.shp"
) %>%
  # add rownum as first column
  mutate(rownum = row_number(), .before = everything())

## gridded candidate plots ####
gridded_plots <- read_sf("./processed_data/gridded_candidate_plots.shp") %>%
  st_transform(st_crs(cbi))

## predictor rasters ####
# LandFire ESP
site_potential <- rast("./processed_data/lf_site_potential_new.tif")
# set active category to 2, which is the fine-scale zone*esp*esplf
# categorization
activeCat(site_potential) <- 2
names(site_potential) <- "esp"

# additional variables
raster_filenames <- c(
  "elev_down", "aspect_down", "TRI_down", "TPI_down", "slope_down",
  "distance_to_road", "ppt", "tmin", "tmmx", "th",
  "vpdmax", "rmax", "vs", "fm100", "fm1000"
)
raster_varnames <- c(
  "elev", "aspect", "tri", "tpi", "slope", "roads_distance", "ppt",
  "tmin", "tmmx", "th", "vpdmax", "rmax", "vs", "fm100", "fm1000"
)
rasts <- map(raster_filenames, function(x) {
  # load each raster file
  rast(paste0("./processed_data/", x, ".tif"))
}) %>%
  # unlist
  rast
# set names for each layer
names(rasts) <- raster_varnames
# rename and combine other raster layers
rasts <- c(rasts, site_potential, cbi)

# TODO ppt was missing from gridded plots. Could fix this in weather_and_knn
# script. Just adding it here for now.
names(gridded_plots)
gridded_plots$ppt <- exact_extract(
  rasts$ppt,
  gridded_plots,
  fun = "mean",
  weights = "area",
  progress = TRUE
)


## what sizes are treatment areas

# 4046.86 m2 in an acre
# calculating buffer size for making plots of various acreages
plot_sizes<-data.frame(acres=c(2,5,15,50,100,150,250))
plot_sizes$radius_m<-sqrt((plot_sizes$acres*4046.86)/pi)

# TODO only place points in valid CBI areas and NOT in vegetation treatments

# random seed
set.seed(2076)

# select 60 random points of each size
random_point_in_burn<-st_sample(burn_perimiter,size=nrow(plot_sizes)*60)
# buffer them
random_point_in_burn<-st_buffer(random_point_in_burn,dist=rep(plot_sizes$radius_m))# check last arg
# reproject to match cbi raster
random_point_in_burn<-sf::st_transform(random_point_in_burn,crs=st_crs(masked_cbi))



# create empty dataframe for cross-validation results
size_cv_results<-data.frame()

# setup rf models ####
# masked_cbi[masked_cbi==9]<-NA # already done
# masked_cbi[masked_cbi==0]<-1

# this function masks the validation area, trains the rf model and knn models,
# predicts on the validation area, then compares to the actual values in the
# validation area.
size_cv <- function(index_number) {
  # choose random point by random point index number
  validation_plot<-random_point_in_burn[index_number]

  # calculate actual cbi in validation area
  actual_burn <- exact_extract(x = masked_cbi,
                               y = st_as_sf(validation_plot),
                               fun = "mean",
                               weights = "area")



  # if there are no NA values in the actual burn data, proceed
  if (sum(is.na(actual_burn[[1]]))==0){

    # mask out the cbi raster to avoid the validation area + 60m buffer
    validation_mask_cbi <- mask(masked_cbi,
                                st_as_sf(st_buffer(validation_plot, 60)), # should already be sf? TODO
                                inverse=TRUE)

    # calculate size of radius for the random forest neighborhood
    radius <- sqrt((150 * 4046.86 +
                      as.numeric(st_area(validation_plot))) / pi)
    # create a buffer around the validation plot
    # Nate: with a constant radius, the larger plots will have larger
    # neighborhoods. Is this desired? We could scale the radius for each plot
    # size to keep the neighborhood size constant.

    # TODO implement testing neighborhood size radius == 100 and non-testing
    # radius == 800

    one_plot_buffer <- st_buffer(validation_plot, radius)
    # crop the burn raster to the buffer
    clipped_burn_raster <- crop(validation_mask_cbi, st_as_sf(one_plot_buffer)) # need st as sf? TODO
    # TODO add mask on top of crop to get circular neighborhood

    # calc raster dims
    xcell <- clipped_burn_raster@ncols
    ycell <- clipped_burn_raster@nrows
    # build empty dataframe from raster dims TODO refactor
    extract_df <- data.frame(lat = rep(seq(clipped_burn_raster@extent[1],
                                           clipped_burn_raster@extent[2],
                                           length.out=ycell),
                                       each=xcell),
                             lon = rep(seq(clipped_burn_raster@extent[3],
                                           clipped_burn_raster@extent[4],
                                           length.out = ycell),
                                       times = xcell))

    # extract cbi values from the clipped burn raster
    # TODO should just be using cell numbers the whole time here, much more
    # efficient than extract on 8 different rasters. Could be using for CBI as
    # well.
    extract_df$CBI <- raster::extract(validation_mask_cbi, y = as.matrix(extract_df[, 1:2]))
    extract_df$elev <- raster::extract(elev_down, y = as.matrix(extract_df[, 1:2]))
    extract_df$aspect <- raster::extract(aspect_down, y = as.matrix(extract_df[, 1:2]))
    extract_df$TRI <- raster::extract(TRI_down, y = as.matrix(extract_df[, 1:2]))
    extract_df$TPI <- raster::extract(TPI_down, y = as.matrix(extract_df[, 1:2]))
    extract_df$slope <- raster::extract(slope_down, y = as.matrix(extract_df[, 1:2]))
    extract_df$road_distance <- raster::extract(roads_distance, y = as.matrix(extract_df[, 1:2]))
    extract_df$env_potential <- as.factor(raster::extract(site_potential, y = as.matrix(extract_df[, 1:2])))
    # Make new predict_df for ...?
    predict_df <- extract_df[is.na(extract_df$CBI) == TRUE, ]
    extract_df_full <- extract_df

    # remove rows with NA CBI values and incomplete cases
    extract_df <- extract_df[extract_df$CBI != 0, ]
    extract_df <- extract_df[extract_df$CBI != 9, ]
    extract_df <- extract_df[complete.cases(extract_df), ]



# list ####
    # local spatial rf
    # kriging
    # perimeter method
    # psm

    # convert the dataframe to an sf object for spatial rf
    coordinates <- data.frame(lat=extract_df$lat,lon=extract_df$lon)
    extract_sf <- st_as_sf(extract_df, coords=c("lat", "lon"))
    extract_sf$lat <- extract_df$lat
    extract_sf$lon <- extract_df$lon

    # fit spatial model
    spatial_rf_model <- splmRF(
      cbi ~ elev + aspect + TRI + TPI + slope + road_distance + env_potential + lat + lon,
      data = extract_sf,
      spcov_type = "gravity",
      local = c(parallel = TRUE, ncores = detectCores() - 4),
      mtry = 4,
      min.node.size = 2,
      sample.fraction = 0.89
    )

    # create empty dataframe to hold validation plot prediction
    predict_df <- data.frame(
      lat = rasterToPoints(raster::crop(site_potential,
                                        y = st_as_sf(validation_plot)))[, 1],
      lon = rasterToPoints(raster::crop(site_potential,
                                        y = st_as_sf(validation_plot)))[, 2]
    )

    # extract predictor variables for the validation plot # TODO fix redundancy here
    predict_df$elev <- raster::extract(elev_down, y = predict_df[, 1:2])
    predict_df$aspect <- raster::extract(aspect_down, y = predict_df[, 1:2])
    predict_df$TRI <- raster::extract(TRI_down, y = predict_df[, 1:2])
    predict_df$TPI <- raster::extract(TPI_down, y = predict_df[, 1:2])
    predict_df$slope <- raster::extract(slope_down, y = predict_df[, 1:2])
    predict_df$road_distance <- raster::extract(roads_distance, y = predict_df[, 1:2])
    predict_df$env_potential <- as.factor(raster::extract(site_potential, y = predict_df[, 1:2]))

    predict_df$actual_burn <- raster::extract( # TODO exact_extract
      x = masked_cbi, y = predict_df[, 1:2], na.rm = TRUE, weights = TRUE,
      exact = TRUE, normalizeWeights = TRUE, small = TRUE)

    predict_df <- predict_df[predict_df$actual_burn != 0, ]
    predict_df <- predict_df[predict_df$actual_burn != 9, ]

    predict_df <- predict_df[complete.cases(predict_df), ]
    predict_df$lat2 <- predict_df$lat
    predict_df$lon2 <- predict_df$lon

    # combine predictors into single df for rf predictions
    predict_df <- st_as_sf(predict_df,
                           coords = c("lat2", "lon2"),
                           crs = st_crs(site_potential))

    # make spatial rf predictions
    predict_df$modeled_values <- predict(spatial_rf_model, newdata = predict_df) # Nate: doesn't this have to be a spatial object?
    # prep conversion to raster
    for_conversion <- data.frame(
      lat = predict_df$lat,
      lon = predict_df$lon,
      modeled_values = as.numeric(
        as.character(predict_df$modeled_values))
    )
    # convert to raster
    predict_raster <- rasterFromXYZ(for_conversion, crs = crs(masked_cbi))
    # extract mean cbi value from the raster for the validation plot
    mean_cbi_spatial_rf <- raster::extract(
      predict_raster$modeled_values,
      st_as_sf(validation_plot),
      fun = mean,
      na.rm = TRUE, weights = TRUE,
      exact = TRUE, normalizeWeights = TRUE, small = TRUE
    )

    ## perimeter method
    # calculate area of the validation plot + buffer in acres
    area_of_plot <- st_area(validation_plot) / 4046.86
    buffer_plot <- st_buffer(validation_plot, 60)
    area_of_plot_with_buffer <- st_area(buffer_plot) / 4046.86

    # get radius of (plot area + buffer) - radius of (plot area)
    # Nate: are we aiming for the radius of an outer circle that encompasses the
    # plot area * 2 + buffer area? So the "neighborhood" size is the same as the
    # plot size?
    size_of_perimiter <- sqrt( # should be RADIUS in varname TODO
      ((area_of_plot + area_of_plot_with_buffer) * 4046.86) / pi
    ) - # radius of circle with 2 * plot area + buffer area
      sqrt((area_of_plot) * 4046.86 / pi)
    # create a circle around the entire validation + buffer + neighborhood area
    perimiters <- st_buffer(validation_plot, size_of_perimiter,
                            allow_holes = TRUE) # Nate: allow_holes is TRUE, but this is a circle, so no holes?
    # grab the total area for normalization
    area_perimiter <- st_area(perimiters) / 4046.86

    perimiter_extract <- raster::extract(
      masked_cbi, st_as_sf(perimiters), fun = mean, na.rm = TRUE,
      weights = TRUE, exact = TRUE, normalizeWeights = TRUE, small = TRUE
    )
    origional_extract <- raster::extract(
      masked_cbi, st_as_sf(buffer_plot), fun = mean, na.rm = TRUE,
      weights = TRUE, exact = TRUE, normalizeWeights = TRUE, small = TRUE
    )
    # calculate the perimiter cbi
    perimiter_only <- as.numeric(
      ((perimiter_extract[1] * area_perimiter) -
         (origional_extract[1] * area_of_plot_with_buffer)) /
        (area_perimiter - area_of_plot_with_buffer)
    )

    ## kriging
    # TODO reuse values from above
    size_of_perimiter <- sqrt(
      (150 * 4046.86 + as.numeric(st_area(validation_plot))) / pi
    )
    one_plot_buffer <- st_buffer(validation_plot, size_of_perimiter)
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



