# this script compares the accuracy of methods used to create counterfactual
# outcomes. it is a core part of the analysis.

# load packages ####
library(tidyverse)
library(sf)
library(terra)
library(ranger)
library(SpatialML)
library(tuneRanger)
library(spmodel)
library(caret)
library(exactextractr)
library(conflicted)

# set options ####
conflicts_prefer(dplyr::filter)

#load data ####
# validation and control plots ####
# validation plots are the randomly selected plots used to evaluate the accuracy
# of the different methods. "control" plots are the practitioner-generated
# control plots. both are included in this file.
validation_and_control_plots <- read_sf(
  "./processed_data/val_points_revised.shp",
  fid_column_name = "ID"
)

## pairing information ####
# "pairing" is a lookup table that matches validation plots with their
# respective practitioner-generated control plots, along with some notes. this
# table was manually created during the creation of the practitioner-generated
# control plots. note that some practitioner-generated control plots were used
# twice, such that the number of validation plots exceeds the number of
# practitioner-generated control plots.
pairing <- read.csv("./processed_data/plot_matching.csv")

## gridded candidate plots ####
# "gridded_plots" includes all possible 10-acre plots in a regular grid and is
# used for propensity score matching
gridded_plots <- read_sf("./processed_data/gridded_candidate_plots.shp")

## composite burn index (CBI) raster ####
cbi <- rast("./processed_data/masked_raster.tif")
# as usual, set "unmappable" -> NA and "outside of perimeter" -> 1 (unburned)
cbi[cbi == 9] <- NA
cbi[cbi == 0] <- 1

# analysis ####
## practitioner-generated controls ####
# extract CBI values for all random, non-treated validation and
# practitioner-generated control plots.
plot_method_cbi <- tibble(
  # FID is the plot's unique identifier
  "FID" = validation_and_control_plots$ID,
  # extract mean CBI. ignores NA values.
  "cbi" = exact_extract(!!cbi,
                        validation_and_control_plots,
                        fun = "mean",
                        weights = "area",
                        progress = TRUE),
  # calculate proportion of plot covered by NA values.
  "prop_na" = exact_extract(!!cbi,
                            validation_and_control_plots,
                            fun = "frac",
                            weights = "area",
                            # sets default value for NA values to 5, so that it
                            # is not ignored by "frac" calculation
                            default_value = 5,
                            progress = TRUE)[[5]]
  ) %>%
  # remove plots with >= 20% NA values (=9 plots)
  filter(prop_na < 0.2) %>%
  select(-prop_na)

# use "pairing" to pair validation plots with practitioner-generated controls in
# wide df. validation plots were created first and have low FIDs.
validations <- plot_method_cbi %>% filter(as.integer(FID) < 352)
controls <- plot_method_cbi %>% filter(as.integer(FID) >= 352)
compare <- pairing %>%
  # IDs must be of same type for use as join keys
  mutate(FID_validation = as.character(FID_1),
         FID_control = as.character(FID.of.New.plot)) %>%
  # get cbi values for validation plots
  left_join(plot_method_cbi,
            by = c("FID_validation" = "FID")) %>%
  left_join(plot_method_cbi, by = c("FID_control" = "FID")) %>%
  # rename to differentiate different CBI values
  rename(cbi_validation = cbi.x,
         cbi_control = cbi.y) %>%
  # finally, drop the extra cols
  select(5:8) %>%
  # remove rows with NA values
  filter(complete.cases(.))

# plot control vs. validation
plot(compare$cbi_control, compare$cbi_validation,
     xlab = "Practitioner-generated control plot",
     ylab = "Validation plot")

# Simple perimeters ####
# get geometry of original validation plots
validation_plots <- validation_and_control_plots %>%
  filter(ID %in% compare$FID_validation)
# demonstrate math for calculating "donut" statistics:
# we want to make a 10-acre donut, centered on the original plots.
# the inner radius of the donut is the radius of a 10-acre circle () + 60m:
one_plot <- validation_plots[1, ] # grab a sample plot
plot_area_m <- st_area(one_plot) # in m2, == 10 acres
inner_radius = as.numeric(sqrt(plot_area_m / pi)) + 60 # all in meters
# the outer radius is the radius of a circle with that gives an area equal to 10
# acres (40450.11078 m2) + the area of the inner circle.
inner_circle_area = inner_radius^2 * pi # in m2
outer_circle_area = inner_circle_area + as.numeric(plot_area_m) # +10 acres
outer_radius = sqrt(outer_circle_area / pi) # in m2
# the target donut shape is a ~33m wide ring:
outer_radius - inner_radius # in m
# generate donut geometries
plot_centers <- validation_plots %>% st_centroid()
donuts_outer <- plot_centers %>%
  st_buffer(dist = outer_radius) %>%
  st_geometry
donuts_inner <- plot_centers %>%
  st_buffer(dist = inner_radius) %>%
  st_geometry
donuts <- map2(donuts_outer, donuts_inner, st_difference) %>%
  st_sfc %>%
  st_sf(geometry = .) %>%
  # set projection
  st_set_crs(st_crs(validation_plots))
# assign FIDs
donuts$FID <- validation_plots$ID
# extract CBI for all donuts and add to compare
donuts$cbi_perimeter <- exact_extract(cbi,
                                       donuts,
                                       fun = "mean",
                                       weights = "area",
                                       progress = TRUE)

# add donut CBI values to compare df
compare <- left_join(compare,
                     as.data.frame(donuts),
                     by = c("FID_validation" = "FID"))

# plot
plot(compare$cbi_perimeter, compare$cbi_validation,
     xlab = "Perimeter control plot",
     ylab = "Validation plot")

# random forest ####
## load additional data ####
# LandFire ESP
site_potential <- rast("./processed_data/lf_site_potential_new.tif")
# set active category to 2, which is the fine-scale zone*esp*esplf categorization
activeCat(site_potential) <- 2
names(site_potential) <- "esp"
# HPCC burn perimeter
burn_perimiter <- read_sf("./processed_data/burn_perimiter.shp")
# cbi raster, masked by original plots + 60m buffer
validation_mask_cbi <- mask(cbi, st_buffer(validation_plots, 60), inverse=TRUE)
names(validation_mask_cbi) <- "masked_cbi"
# treatments
veg_treatments <- read_sf("./processed_data/vegetation_treatments_hpcc.shp") %>%
  st_make_valid
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
names(cbi) <- "cbi"
rasts <- c(rasts, site_potential, cbi, validation_mask_cbi)
plot(rasts)
plot(rasts[[17:18]])

# spatial rf ####
# make spatial rf dataframes
rf_dfs_spatial <- map(rf_dfs, \(x) {
  x %>%
    # add copy of lat/lon to make sf
    mutate(x = lon, y = lat) %>%
    # make sf
    st_as_sf(coords = c("x","y")) # dim = XYZ by default
})

rf_dat_spatial_sub <- slice_sample(rf_dfs_spatial$rf_train, n = 500) # replace after testing # nolint
## fit global model ####
# Nate: does it make sense to include lon + lat as predictors here?
spatial_rf_model <- splmRF(
  masked_cbi ~ elev + aspect + tri + tpi + slope + roads_distance + ppt + tmin +
    tmmx + th + vpdmax + rmax + vs + fm100 + fm1000 + esp + lon + lat,
  data = rf_dat_spatial_sub, # replace after testing
  spcov_type = "gravity",
  # leave 4 CPU cores free
  local = list(parallel = TRUE, ncores = detectCores() - 4),
  num.trees = 10, # remove after testing
  mtry = 2, # should be 4, replace after testing
  min.node.size = 2,
  sample.fraction = 0.89
)
spatial_rf_model

## extract global predictions ####
# predict CBI values for all points in the raster
# Nate: only took 17sec?
rf_dfs$rf_test$rf_global_predicted <- predict(spatial_rf_model,
                                              rf_dfs_spatial$rf_test)
# build raster
rf_global_preds_rast <- rf_dfs$rf_test %>%
  select(lon, lat, rf_global_predicted) %>%
  rast(crs = crs(cbi))
# Nate: check if you want.
# mapview::mapview(rf_global_preds_rast) + mapview::mapview(validation_plots)

# extract rf predicted CBI values for validation plots
validation_plots$cbi_rf_global <- exact_extract(
  rf_global_preds_rast,
  validation_plots,
  fun = "mean",
  weights = "area",
  progress = TRUE
)
# join to compare df
compare <- left_join(compare,
                     as.data.frame(validation_plots) %>%
                       select(ID, cbi_rf_global),
                     by = c("FID_validation" = "ID"))
# plot
plot(compare$cbi_rf_global, compare$cbi_validation,
     xlab = "Global spatial random forest predicted CBI",
     ylab = "Validation plot CBI")

## fit local model ####

# ## kriging ####
#
# library(snapKrig)
#
# compare_df$krigged_values<-NA
#
# for (i in 1:nrow(compare_df)){
# #for (i in 1:5){
#
#
# one_plot<-origional_plots[origional_plots$FID_1==compare_df$plot_id[i],]
#
# size_of_perimiter<-sqrt((200*4046.86)/pi)
#
# one_plot_buffer<-st_buffer(one_plot,size_of_perimiter)
#
# clipped_burn_raster<-raster::crop(cbi,one_plot_buffer)
#
# clipped_burn_raster<-raster::mask(clipped_burn_raster,st_as_sf(st_buffer(one_plot,60)),inverse=TRUE)
#
# clipped_burn_raster[clipped_burn_raster==0]<-1
#
# clipped_burn_raster[clipped_burn_raster==9]<-NA
#
# ### sk_fit fails when corner cells are NA. This sets it to 1 when it is NA.
# ### Because corners are the fartest from the treatments this should have minimal
# ### impact on predictions
# if (is.na(clipped_burn_raster[1,1])==TRUE){clipped_burn_raster[1,1]<-1}
# if (is.na(clipped_burn_raster[nrow(clipped_burn_raster),1])==TRUE){clipped_burn_raster[nrow(clipped_burn_raster),1]<-1}
# if (is.na(clipped_burn_raster[1,ncol(clipped_burn_raster)])==TRUE){clipped_burn_raster[1,ncol(clipped_burn_raster)]<-1}
# if (is.na(clipped_burn_raster[nrow(clipped_burn_raster),ncol(clipped_burn_raster)])==TRUE){clipped_burn_raster[nrow(clipped_burn_raster),ncol(clipped_burn_raster)]<-1}
#
#
# clipped_burn_raster_sk<-sk(clipped_burn_raster)
#
#
#
# kriging_fit<-sk_fit(clipped_burn_raster_sk,n_max=3000)
#
#
# krigged_final<-sk_cmean(clipped_burn_raster_sk,kriging_fit)
# plot(krigged_final)
#
# krigged_final_raster<-sk_export(krigged_final)
#
# values<-extract(krigged_final_raster,one_plot)
#
# compare_df$krigged_values[i]<-mean(values$lyr.1)
#
#
# }
#
#
# sqrt(mean((compare_df$krigged_values-compare_df$origional_plot_cbi)**2,na.rm=TRUE))
# write.csv(compare_df,"./processed_data/kriging_results.csv")

##### TODO clean up kriging above




# spatial rf with weather ####
# this is now done alongside the partial local spatial rf model above.

# cluster based matching ####

## remove the treated and validation plots
# read in veg treatments
veg_treatments<-read_sf("./processed_data/vegetation_treatments_hpcc_new.shp") %>%
  st_transform(st_crs(cbi))

# transform in grid of plots, already loaded earlier
gridded_plots<-st_transform(gridded_plots,st_crs(cbi))
# remove the treated plots from the gridded plots
gridded_plots_2<-st_difference(gridded_plots,st_union(veg_treatments))

# remove the validation plots
gridded_plots_3<-st_difference(gridded_plots_2,st_union(st_buffer(validation_and_control_plots,60)))
# Nate: validation and control plots contain (obviously) the validation plots,
# but also the practitioner-generated control plots. There's no need to make
# these areas off-limits here. In fact, this could be one of the reasons this
# method underperforms: the best gridded plots, based on PSM, are likely going
# to be near the practitioner-selected control. I know a lot has changed in this
# script, but I am pretty sure this error is present in the original and not
# introduced by me.

mapview::mapview(gridded_plots_3)

max_area<-max(st_area(gridded_plots_3))
# Nate: I assume this is just to get rid of partial plots.
gridded_plots_4<-gridded_plots_3[as.numeric(st_area(gridded_plots_3))>as.numeric(max_area)-50,]

## remove the treated and validation plots (refactored)
# define off-limits areas: treated areas and validation plots (no practitioner generated controls)
off_limits <- st_union(veg_treatments, st_buffer(validation_plots, 60))
mapview::mapview(off_limits)
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
