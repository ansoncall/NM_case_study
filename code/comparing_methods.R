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
## build dataframes ####
rasts_train <- rasts %>%
  # mask out validation plots with 60m buffer
  mask(st_buffer(validation_plots, 60), inverse = TRUE)
plot(rasts_train[[3]])
rasts_test <- rasts %>%
  # remove masked_cbi, as it is not used in testing
  tidyterra::select(-masked_cbi) %>%
  # keep only validation plots. add small buffer to ensure entire plot is
  # included in result, as mask() is seemingly not precise? (Nate: check)
  mask(st_buffer(validation_plots, 30))
# mapview::mapview(rasts_test[[3]]) + mapview::mapview(validation_plots) # nolint

# convert to dataframe
rf_dfs <- map(list(rasts_train, rasts_test), \(x) {
  x %>%
    # convert to sf points. only remove if all layers are NA.
    as.points(na.all = TRUE) %>%
    # convert to df, keeping point geoms as xy coords
    as.data.frame(geom = "XY") %>%
    # rename xy cols, even though they aren't exactly in lat + lon degrees
    rename(lat = y, lon = x) %>%
    # ensure site potential is a factor
    mutate(esp = as_factor(esp)) %>%
    # keep only complete cases
    filter(complete.cases(.))
}) %>%
  set_names(c("rf_train", "rf_test"))


rf_sub <- rf_dfs$rf_train %>% slice_sample(n = 500) # replace after testing
## fit model ####
# takes ~ 2 hrs.
ranger_rf <- ranger(
  masked_cbi ~ elev + aspect + tri + tpi + slope + roads_distance + ppt + tmin +
    tmmx + th + vpdmax + rmax + vs + fm100 + fm1000 + esp + lon + lat,
  rf_sub, # replace after testing
  # parameter tuning was previously completed with tuneRanger() Nate: check
  importance="impurity",
  num.trees = 10, # should be 1000, replace after testing
  mtry = 2, # should be 5, replace after testing
  min.node.size=2,
  sample.fraction=0.89,
  seed = 18003634377 # 1-800-DOD-GERS let's go Shohei!
)
ranger_rf

## variable importance plot ####
var_imp_data <- data.frame(importance = ranger_rf$variable.importance) %>%
  mutate(
    # save row names as variable names
    variable = fct_reorder(row.names(.), importance) %>%
      # recode to make printable
      fct_recode(TRI = "tri",
                 TPI = "tpi",
                 elevation = "elev",
                 latitude = "lat",
                 longitude = "lon",
                 `road distance`  = "roads_distance",
                 precipitation = "ppt",
                 `minimum temperature` = "tmin",
                 `maximum temperature` = "tmmx",
                 `wind direction` = "th",
                 `wind speed` = "vs",
                 `max. vapor pressure deficit` = "vpdmax",
                 `moisture, 100-hr` = "fm100",
                 `moisture, 1000-hr` = "fm1000",
                 `max. relative humidity` = "rmax",
                 `environmental site potential` = "esp"
      ),
    Category = case_when( # nolint
      row.names(.) %in% c("elev", "aspect", "tri", "tpi", "slope", "lat", "lon",
                          "roads_distance") ~ "landscape",
      row.names(.) %in% c("esp", "ppt", "tmin", "vpdmax") ~ "bioclimatic",
      row.names(.) %in% c("tmmx", "th", "vs", "fm100", "fm1000", "rmax") ~
        "weather",
      .default = "other" # should not occur
    )
  )
cbPalette <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
  )
var_imp_plot <- var_imp_data %>%
  ggplot(aes(x = importance, y = variable, fill = Category)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(text = element_text(size = 15)) +
  xlab("Importance (Gini index)") +
  ylab("") +
  scale_fill_manual(values = cbPalette)
var_imp_plot
# save figure
ggsave("./figures/rf_importance.tif", units = 'in', width=10, height=8, dpi=600)

## extract predictions ####
# predict CBI values for all points in the raster
rf_dfs$rf_test$rf_predicted <- predict(ranger_rf, rf_dfs$rf_test)$predictions
# build raster
rf_preds_rast <- rf_dfs$rf_test %>%
  select(lon, lat, rf_predicted) %>%
  rast(crs = crs(cbi))
# Nate: this is just predictions for the "testing" (=validation) plots. This is
# efficient, but if you want to make a wall-to-wall map we need to generate
# predictions over the full area. See:
# mapview::mapview(rf_preds_rast) + mapview::mapview(validation_plots)

# extract rf predicted CBI values for validation plots
validation_plots$cbi_rf <- exact_extract(
  rf_preds_rast,
  validation_plots,
  fun = "mean",
  weights = "area",
  progress = TRUE
)
# join to compare df
compare <- left_join(compare,
                     as.data.frame(validation_plots) %>%
                       select(ID, cbi_rf),
                     by = c("FID_validation" = "ID"))
# plot
plot(compare$cbi_rf, compare$cbi_validation,
     xlab = "Random forest predicted CBI",
     ylab = "Validation plot CBI")

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

# # TODO remove the stuff below when you're done with it
# predict_for_each_plot_spatial<-function(location){
#
#   predict_df<-data.frame(lat=rasterToPoints(raster::crop(site_potential,y=st_as_sf(location)))[,1],
#                          lon=rasterToPoints(raster::crop(site_potential,y=st_as_sf(location)))[,2])
#
#   predict_df$elev<-raster::extract(elev_down,y=predict_df[,1:2])
#   predict_df$aspect<-raster::extract(aspect_down,y=predict_df[,1:2])
#   predict_df$TRI<-raster::extract(TRI_down,y=predict_df[,1:2])
#   predict_df$TPI<-raster::extract(TPI_down,y=predict_df[,1:2])
#   predict_df$slope<-raster::extract(slope_down,y=predict_df[,1:2])
#   predict_df$road_distance<-raster::extract(roads_distance,y=predict_df[,1:2])
#   predict_df$env_potential<-as.factor(raster::extract(site_potential,y=predict_df[,1:2]))
#   predict_df$ppt<-raster::extract(ppt,y=as.matrix(predict_df[,1:2]))
#   predict_df$tmin<-raster::extract(tmin,y=as.matrix(predict_df[,1:2]))
#   predict_df$tmmx<-raster::extract(tmmx,y=as.matrix(predict_df[,1:2]))
#   predict_df$th<-raster::extract(th,y=as.matrix(predict_df[,1:2]))
#   predict_df$vs<-raster::extract(vs,y=as.matrix(predict_df[,1:2]))
#   predict_df$vpdmax<-raster::extract(vpdmax,y=as.matrix(predict_df[,1:2]))
#   predict_df$fm100<-raster::extract(fm100,y=as.matrix(predict_df[,1:2]))
#   predict_df$fm1000<-raster::extract(fm1000,y=as.matrix(predict_df[,1:2]))
#   predict_df$rmax<-raster::extract(rmax,y=as.matrix(predict_df[,1:2]))
#
#   predict_df$lat2<-predict_df$lat
#   predict_df$lon2<-predict_df$lon
#
#
#   predict_df<-predict_df[complete.cases(predict_df),]
#
#   predict_df<-st_as_sf(predict_df,coords=c("lat2","lon2"),crs=st_crs(site_potential))
#
#   #predict_df$modeled_values<-predict(ranger_rf,data=predict_df,na.rm=TRUE)$predictions
#
#
#
#   predict_df$modeled_values<-predict(spatial_rf_model,newdata=predict_df)
#
#   for_conversion<-data.frame(lat=predict_df$lat,lon=predict_df$lon,modeled_values=as.numeric(as.character(predict_df$modeled_values)))
#
#   predict_raster<-rasterFromXYZ(for_conversion,crs=crs(cbi))
#
#   mean_cbi<-raster::extract(predict_raster$modeled_values,st_as_sf(location),fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)
#
#
#   return(mean_cbi)
# }
#
#
#
# compare_df$spatial_rf_prediction<-NA
# compare_df$rf_prediction<-NA
# for (i in 1:nrow(compare_df)){
#   print(i/nrow(compare_df)*100)
#   one_plot<-origional_plots[origional_plots$FID_1==compare_df$plot_id[i],]
#   compare_df$spatial_rf_prediction[i]=predict_for_each_plot_spatial(one_plot)
#   compare_df$rf_prediction[i]=predict_for_each_plot(one_plot)
#   }
#
# ## RMSE
# sqrt(mean((compare_df$rf_prediction-compare_df$origional_plot_cbi)**2,na.rm=TRUE))
# sqrt(mean((compare_df$spatial_rf_prediction-compare_df$origional_plot_cbi)**2,na.rm=TRUE))
#
# write.csv(compare_df,"./processed_data/larger_rf_model_predictions.csv")
#
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

# fits a "local" spatial rf model for each validation plot, one at a time. map
# this function over all of the FID's in the compare dataframe.

# anything that uses a SpatRaster cannot be run in parallel, so we need to
# preprocess train and test data for each validation plot first.

### build train and test data ####
# makes train and test data for each plot
build_train_test <- function (id) {
  # id <- compare$FID_validation[1] # remove after testing
  # grab one plot
  one_plot <- validation_plots[validation_plots$FID_1 == id, ]

  ## the following section relates to neighborhood size. Nate: let's review
  ## this.
  # # get size of ?? why 1200*4046?? I suppose this leads to the "arbitrary" 607
  # # ha size noted in the manuscript?
  # size_of_perimeter <- sqrt((1200 * 4046.86) / pi) #1243.297
  one_plot_neighborhood <- st_buffer(one_plot, 1243) # hard-coding for speed
  one_plot_target <- st_buffer(one_plot, 60) # the validation area
  # st_area(one_plot_buffer)/10000 # this is ~578 hectares, ~1428 acres. not sure where 607 comes from

  # Create a mask with 1 inside the polygon, NA outside.
  # Nate: idk why, but rasterizing before masking is the only way I can get
  # terra to mask precisely. otherwise I end up with weird extra cells on the
  # boundary of the mask.
  neighborhood_mask_layer <- rasterize(one_plot_neighborhood, rasts, field=1)
  validation_mask_layer <- rasterize(one_plot_target, rasts, field=1)

  # # Nate: masking instead of clipping. working with validation_mask_cbi, where
  # # validation plots are already masked out. result should be a "donut" of
  # # raster cells.
  # masked_burn_raster <- mask(validation_mask_cbi, st_as_sf(one_plot_buffer))
  # # note that with this method, ALL validation plot areas are masked out, not
  # # just the focal validation plot area.
  # mapview::mapview(masked_burn_raster) + mapview::mapview(validation_plots)

  # # To "fix", could start with regular cbi raster and do masking for validation
  # # plot again. though I'm not sure it's really a problem anyway.
  # mapview::mapview(masked_burn_raster) + mapview::mapview(validation_plots)
  # # actually, makes more sense to start with multiband rast so we don't have to
  # # extract again--we can just convert to df. The radius is arbitrary anyway, so
  # # we don't care about extract's ability to area-weight partially intersecting
  # # cells.

  # make "local" raster with crop (reduces memory overhead)
  local <- rasts %>%
    # mask out the target validation area
    crop(neighborhood_mask_layer)

  masked_raster <- local %>%
    # if you don't want to drop cells from non-target validation areas, remove
    # the masked_cbi column:
    select(-masked_cbi) %>%
    # mask out the target validation area
    mask(validation_mask_layer, inverse = TRUE) %>%
    # mask out everything beyond what is considered "local" training data
    crop(neighborhood_mask_layer, mask = TRUE)
  # # remove after testing
  # mapview::mapview(masked_raster[[3]]) +
  #   mapview::mapview(one_plot_neighborhood) +
  #   mapview::mapview(one_plot_target)
  # split training and testing data
  train_df <- masked_raster %>%
    # convert to points. drop where na
    as.points() %>%
    # convert to df, keeping point geoms as xy coords
    as.data.frame(geom = "XY") %>%
    # add lat/lon from y/xcols, even though they aren't exactly in lat + lon degrees
    mutate(lat = y, lon = x) %>%
    # keep only complete cases. will drop all validation cells automatically,
    # since these cell values are NA in cbi_masked layer (unless that layer was
    # previously dropped, then non-target validation areas are retained).
    na.omit()
    # Nate: there are also some edge cells that are being dropped, because the
    # neighborhood extends over the edge of the burn scar for some plots. Not
    # sure this would have been noticed at any point. I suppose it's okay that
    # neighborhoods are different sizes for validation plots near the edge of
    # the burn scar, but might mention this in the manuscript or supplement
    # somewhere.

  # validation plot raster mask, more precise
  one_plot_mask <- rasterize(one_plot_target, rasts, touches = TRUE, field=1)
  test_df <- local %>%
    # definitely want to drop the masked_cbi layer so we retain the validation
    # plot area when creating the df.
    select(-masked_cbi) %>%
    # keep ONLY the target validation area this time. using 30m buffer to ensure
    # edge cells stay in at this point.
    mask(one_plot_mask) %>%
    # convert to points. only drop if all layers are NA.
    as.points() %>%
    # convert to df, keeping point geoms as xy coords
    as.data.frame(geom = "XY") %>%
    mutate(lat = y, lon = x)
  # we should have all complete complete cases, and i want to know if we
  # don't, so no final na.omit here.

  # # remove after testing
  # set.seed(123456) # for reproducibility
  # train_df <- train_df %>%
  #   slice_sample(n = 500)

  # convert to sf
  train_sf <- train_df %>%
    st_as_sf(coords = c("x","y"), crs = st_crs(rasts))

  test_sf <- test_df %>%
    st_as_sf(coords = c("x", "y"), crs = st_crs(rasts))

  return(list(
    validation_FID = id,
    train = train_sf,
    test = test_sf
  ))
}

# map function over all validation plot FIDs
tictoc::tic()
train_test_each_plot <- purrr::map(compare$FID_validation,
                                   build_train_test,
                                   .progress = TRUE)
tictoc::toc() # ~30 mins for subsampled data, ~1hr for full data.

# takes train and test sfs and returns prediction sf
fit_local_spatial_rf <- function(train_test_list) {

  train_sf <- train_test_list$train
  test_sf <- train_test_list$test

  # Nate: I get that part of the exercise is to see what happens when we don't
  # have to approximate the spatial structure of the data. But when/where is
  # this approximation happening in the global spatial rf model? It looks like
  # its in the splm() help, and IIRC the approximation methods kick in when n >
  # 5000 **OR** when "local=" parameters are set. Is that why we have a
  # "neighborhood" of the set size? It gives a neighborhood of ~4000 obs. If
  # this is correct, and I think it is, then we might be accidentally triggering
  # an approximation method by setting "local=" parameters here ("...If a list
  # is provided, the following arguments detail the big data approximation:..").
  # It seems that the parallelization options apply parallel processing across
  # the groups used in the approximation methods, and we don't want
  # approximation methods here, so there should be nothing to parallelize. I'm
  # setting local=FALSE here just to be double sure, though it should default to
  # that based on the sample size.
  # tictoc::tic()
  local_spatial_rf_model <- splmRF(
    cbi ~ elev + aspect + tri + tpi + slope + roads_distance + ppt +
      tmin + tmmx + th + vpdmax + rmax + vs + fm100 + fm1000 + esp + lon + lat,
    data = train_sf,
    spcov_type = "exponential", # this is the default already.
    local = FALSE, # no spatial approximation!
    mtry = 4,
    min.node.size = 2,
    sample.fraction = 0.89)
  # tictoc::toc() # 5 minutes to fit one model.

  # predict for the validation plot. for speed, just focus on the pixels in the
  # validation plot area.

  # tictoc::tic()
  predictions <- predict(local_spatial_rf_model, newdata = test_sf)
  # tictoc::toc() # 11 seconds to predict one validation plot

  # Nate: if we were to cut a corner, we could just return the mean at this
  # stage instead of converting back to rast and taking the area-weighted cell
  # means with exact_extract(). Otherwise we just return the predictions as a
  # vector here.

  # mean(predictions)
  return(list(
    validation_FID = train_test_list$validation_FID,
    predictions = predictions))

}

# map this over the list of training and testing data to fit models and extract
# predictions
tictoc::tic()
all_preds <- purrr::map(train_test_each_plot,
                        fit_local_spatial_rf,
                        .progress = TRUE)
tictoc::toc() # ~7 mins for subsampled (n=500) training data.

# # map this over the list of training and testing data to fit models and extract
# # predictions, now in parallel
# tictoc::tic()
# plan(multisession, workers = n_cores)
# all_preds <- purrr::map(train_test_each_plot,
#                         fit_local_spatial_rf,
#                         .progress = TRUE)
# tictoc::toc()


get_mean_predictions <- function(preds_vec, train_test_each_plot) {
  # get id
  id <- train_test_each_plot$validation_FID
  # get test_sf
  test_sf <- train_test_each_plot$test
  # add predictions to test_sf
  test_sf$cbi_local_spatial_rf <- preds_vec
  # convert predictions to raster and extract by exact plot polygon
  predict_rast <- test_sf %>%
    select(lon, lat, cbi_local_spatial_rf) %>%
    # convert to rast
    rast(crs = crs(cbi))
    # extract
  predict_mean <- exact_extract(predict_rast,
                                one_plot,
                                fun = "mean",
                                weights = "area")
  # final output is the mean cbi value of the validation plot
  return(list(id = id, predictions = predict_mean))
}

tictoc::tic()
local_spatial_results <- map(compare$FID_validation, fit_local_spatial_rf)
tictoc::toc() # 24 minutes with "test" subset of 500 obs, 1 core.
# this has taken multiple days



# now, parallelize
install.packages("furrr")
library(furrr)
tictoc::tic()
# get n cores
n_cores <- parallel::detectCores() - 4 # leave 4 cores free
# set up parallel processing
plan(multisession, workers = n_cores)
# fit local spatial rf model for each validation plot in parallel
local_spatial_results <- future_map(
  compare$FID_validation,
  fit_local_spatial_rf,
  .options = furrr_options(packages = c("terra",
                                         "sf",
                                         "tidyterra",
                                         "tidyverse",
                                         "exactextractr"),
                           seed = 1548966)
)
tictoc::toc() # xx minutes with "test" subset of 500 obs, 16 cores.

# test furrr and terra
testlist <- list("test1" = 1)

test <- future_map(
  testlist,
  function(x) {
    # must write rast to and pass filename to future_map
    vs <- rast("./processed_data/vs.tif")
    # this is a test function that uses terra::mask
    y <- mask(vs, st_as_sf(validation_plots[x, ]), inverse = TRUE)
    writeRaster(y,
              filename = paste0("./processed_data/test_", x, ".tif"),
              overwrite = TRUE)
  },
  .options = furrr_options(
    packages = c("terra", "sf", "tidyterra", "tidyverse", "exactextractr"),
    seed = 155)
)
library(tidyverse)
rast("./processed_data/test_1.tif") %>% plot

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

  predict_raster<-rasterFromXYZ(for_conversion,crs=crs(cbi))

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


gridded_plots_3<-st_difference(gridded_plots_2,st_union(st_buffer(validation_and_control_plots,60)))

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
