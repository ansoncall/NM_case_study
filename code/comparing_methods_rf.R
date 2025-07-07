# this script is used to generate predicted counterfactual values for the random
# forest methods. These methods are computationally expensive, and this script
# takes a long time to run.

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
library(mapview)
library(units)
library(magrittr)
# set options ####
conflicts_prefer(dplyr::filter, purrr::set_names)

# testing ####
# set testing == TRUE to run all functions on a smaller subset of the data for
# testing purposes. Otherwise, this script will take a very long time to run.
testing <- TRUE

# load data ####
## composite burn index (CBI) raster ####
cbi <- rast("./processed_data/masked_raster.tif")
# as usual, set "unmappable" -> NA and "outside of perimeter" -> 1 (unburned)
names(cbi) <- "cbi"
cbi[cbi == 9] <- NA
cbi[cbi == 0] <- 1

## validation plots ####
validation_plots <- read_sf(
  "./processed_data/val_points_revised.shp",
  fid_column_name = "ID"
) %>%
  filter(as.integer(ID) < 352) %>%
  st_transform(crs = crs(cbi))

## HPCC burn perimeter ####
# TODO unused? remove?
burn_perimiter <- read_sf("./processed_data/burn_perimiter.shp")

## treatments ####
# Nate: this should be hpcc_new.shp but was hpcc.shp in the original. Not sure
# what damage this might have done.
veg_treatments <- read_sf("./processed_data/vegetation_treatments_hpcc_new.shp")

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

# rf ####
## build dataframes ####
rasts_train <- rasts %>%
  # mask out validation plots with 60m buffer and treated areas
  mask(st_buffer(validation_plots, 60), inverse = TRUE) %>%
  mask(veg_treatments, inverse = TRUE)

rasts_test <- rasts %>%
  # remove cbi, as it is not used in testing
  tidyterra::select(-cbi) %>%
  # keep only validation plots. add small buffer to ensure entire plot is
  # included in result.
  mask(st_buffer(validation_plots, 100))

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
}, .progress = TRUE) %>%
  set_names(c("rf_train", "rf_test"))


if (testing == TRUE) {
  rf_dfs$rf_train %<>% slice_sample(n = 500) # nolint
}

## fit model ####
ranger_rf <- ranger(
  cbi ~ elev + aspect + tri + tpi + slope + roads_distance + ppt + tmin +
    tmmx + th + vpdmax + rmax + vs + fm100 + fm1000 + esp + lon + lat,
  rf_dfs$rf_train, # replace after testing
  # parameter tuning was previously completed with tuneRanger() Nate: check
  importance = "impurity",
  num.trees = 10, # should be 1000, replace after testing
  mtry = 2, # should be 5, replace after testing
  min.node.size = 2,
  sample.fraction = 0.89,
  seed = 18003634377 # 1-800-DOD-GERS let's go Shohei!
)
ranger_rf

## variable importance plot ####
var_imp_data <- data.frame(importance = ranger_rf$variable.importance) %>%
  mutate(
    # save row names as variable names
    variable = fct_reorder(row.names(.), importance) %>%
      # recode to make printable
      fct_recode(TRI = "tri", # nolint
                 TPI = "tpi", # nolint
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
                 `environmental site potential` = "esp"),
    Category = case_when( # nolint
      row.names(.) %in% c("elev", "aspect", "tri", "tpi", "slope", "lat", "lon",
                          "roads_distance") ~ "landscape",
      row.names(.) %in% c("esp", "ppt", "tmin", "vpdmax") ~ "bioclimatic",
      row.names(.) %in% c("tmmx", "th", "vs", "fm100", "fm1000", "rmax") ~
        "weather",
      .default = "other" # should not occur
    )
  )
cb_palette <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
)
var_imp_plot <- var_imp_data %>%
  ggplot(aes(x = importance, y = variable, fill = Category)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(text = element_text(size = 15)) +
  xlab("Importance (Gini index)") +
  ylab("") +
  scale_fill_manual(values = cb_palette)
var_imp_plot
# save figure
ggsave("./figures/rf_importance.tif",
       units = "in",
       width = 10,
       height = 8,
       dpi = 600)

## extract predictions ####
# predict CBI values for all points in the raster
rf_dfs$rf_test$rf_predicted <- predict(ranger_rf, rf_dfs$rf_test)$predictions
# build raster
rf_preds_rast <- rf_dfs$rf_test %>%
  select(lon, lat, rf_predicted) %>%
  rast(crs = cbi)

# # Nate: note weird raster alignment. Also note that some validation plots are
# # technically outside the burned area and do not have weather data.
# # nolint start
# mapview(rf_preds_rast) + mapview(validation_plots) + mapview(veg_treatments) +
#   mapview(rasts[[12]])
# # nolint end

# extract rf predicted CBI values for validation plots
validation_plots$cbi_rf <- exact_extract(
  rf_preds_rast,
  validation_plots,
  fun = "mean",
  weights = "area",
  progress = TRUE
)

# spatial rf ####
# make spatial rf dataframes
rf_dfs_spatial <- map(rf_dfs, \(x) {
  x %>%
    # add copy of lat/lon to make sf
    mutate(x = lon, y = lat) %>%
    # make sf
    st_as_sf(coords = c("x", "y"), crs = st_crs(cbi)) # dim = XYZ by default
})

if (testing == TRUE) {
  rf_dfs_spatial$rf_train %<>% slice_sample(n = 500) # nolint
}

## fit global model ####
spatial_rf_model <- splmRF(
  cbi ~ elev + aspect + tri + tpi + slope + roads_distance + ppt + tmin +
    tmmx + th + vpdmax + rmax + vs + fm100 + fm1000 + esp + lon + lat,
  data = rf_dfs_spatial$rf_train,
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
# predict CBI values for all points in the test raster
rf_dfs$rf_test$rf_global_predicted <- predict(spatial_rf_model,
                                              rf_dfs_spatial$rf_test)
# build raster
rf_global_preds_rast <- rf_dfs$rf_test %>%
  select(lon, lat, rf_global_predicted) %>%
  rast(crs = crs(cbi))
# # Nate: check if you want.
# mapview(rf_global_preds_rast) + mapview::mapview(validation_plots) # nolint

# extract rf predicted CBI values for validation plots
validation_plots$cbi_rf_global <- exact_extract(
  rf_global_preds_rast,
  validation_plots,
  fun = "mean",
  weights = "area",
  progress = TRUE
)

# local spatial rf ####
# this includes "local spatial" rf models with and without weather variables.

# makes train and test data for each plot
build_train_test <- function(id, test = testing) {
  # grab one plot
  one_plot <- validation_plots[validation_plots$FID_1 == id, ]

  # Nate: 1243 buffer was arbitrarily defined, here equal to roughly a based on
  # a ~1200 ac circle: sqrt((1200 acres * 4046.86 m^2/acre) / pi)

  # rethinking this, I think we should target a sample size that keeps us
  # well under splmRF's default 5000 cell max threshold (though you can change
  # the threshold if you want). We're actually getting 6420 cells with the above
  # values. see:

  # sum(neighborhood_mask[], na.rm = TRUE) # 6420 with 1243 m buffer # nolint

  # I'm going to suggest a buffer of 1000m, giving a neighborhood of ~4223
  # cells. 105 of those will be masked out since they're in the validation area.
  # We need to update the manuscript regardless, because the old "607 ha" value
  # was wrong anyway.

  # note that this buffer size REALLY impacts the memory bottlenecking and speed
  # of the splmRF fitting.

  # define masks for target and neighborhood. smaller neighborhood if testing ==
  # TRUE
  target_poly <- st_buffer(one_plot, 60)

  if (test == TRUE) {
    buffer_size <- 100
  } else {
    buffer_size <- 1000
  }

  suppressWarnings( # ignore attribute variable warning
    neighborhood_poly <- st_difference(st_buffer(one_plot, buffer_size),
                                       target_poly)
  )

  # crop rasts to minimize raster read time
  rasts_crop <- rasts %>%
    crop(st_bbox(neighborhood_poly))

  # rasterize both masks in one go: assign 1 to target, 2 to neighborhood
  masks_raster <- dplyr::bind_rows(
    target_poly %>% mutate(mask_val = 1),
    neighborhood_poly %>% mutate(mask_val = 2)
  ) %>%
    rasterize(rasts_crop, field = "mask_val")

  # extract the masks
  target_mask <- ifel(masks_raster == 1, 1, NA)
  neighborhood_mask <- ifel(masks_raster == 2, 1, NA)

  # apply masks to predictor variable rasters
  train_pts <- rasts_crop %>%
    mask(neighborhood_mask) %>%
    as.points
  # retain coords as attributes
  train_xy <- crds(train_pts)
  train_pts$lon <- train_xy[, 1]
  train_pts$lat <- train_xy[, 2]
  # convert to sf
  train <- train_pts %>%
    st_as_sf(coords = c("x", "y"), crs = st_crs(rasts), remove = FALSE)
  # repeat for test data
  test_pts <- rasts_crop %>%
    mask(target_mask) %>%  # will also include 60m buffer but this is ok.
    as.points
  test_xy <- crds(test_pts)
  test_pts$lon <- test_xy[, 1]
  test_pts$lat <- test_xy[, 2]
  test <- test_pts %>%
    st_as_sf(coords = c("x", "y"), crs = st_crs(rasts))

  return(list(
    validation_FID = id,
    train = train,
    test = test
  ))
}

# map function over all validation plot FIDs
train_test_each_plot <- purrr::map(validation_plots$FID_1,
                                   build_train_test,
                                   .progress = TRUE)
# check sizes of train and test data
walk(train_test_each_plot, function(x) {
  # check that the train and test data are not empty
  if (nrow(x$train) < 10 || nrow(x$test) < 10) {
    message(sprintf(
      "[FID #%d] Warning: %d test values for this plot",
      x$validation_FID, nrow(x$test)
    )
    )
  }
})

# Nate: these three plots are in unmappable regions of the CBI map, so we need
# to drop them. I think these were probably dropped silently in other analyses,
# but here we need to do it explicitly or the splmRF() call will error out.
train_test_each_plot[[96]]$validation_FID
train_test_each_plot[[157]]$validation_FID
train_test_each_plot[[230]]$validation_FID

# drop these two plots (note discrepancy between list index and FID values. This
# will drop FIDs 96, 156, and 229, which are the problem plots)
train_test_each_plot <- train_test_each_plot[-c(96, 157, 230)]

# takes train and test sfs and returns prediction sf
fit_local_spatial_rf <- function(train_test_list, test = testing) {
  # Nate: ran into a problem when developing the "testing" workflow that sort of
  # explains why local rf performs so well--in a way that's pretty easy to
  # understand. I think we've already realized this, but it's really obvious in
  # this instance. When the training set is small enough, there are some cases
  # where all pixels have the same CBI value. There is no variability there, and
  # thus the predictors are meaningless. This is the extreme case, but similar
  # things occur when some CBI values become rare or absent, e.g. there are zero
  # or only a few CBI values == 1 in the training set. This naturally leads to a
  # model fit that never (or rarely) has to "consider" the possibility of a 1
  # ever occurring. This natural bias in the model doesn't reflect any causal
  # mechanism, just the inherent spatial autocorrelation, which is another way
  # of saying that the local-spatial RF is doing it's job. The narrowing of the
  # possible outcome space is a benefit to prediction accuracy.

  # However, to even get this function to work for those instances with no
  # variability in the training set, we need to recognize that the model will
  # error out and catch those errors so we can still return the predictions (in
  # the case of, e.g., all training data CBI == 2, the prediction would just be
  # a field of 2's).

  local_spatial_rf_model <- try(splmRF(
    cbi ~ elev + aspect + tri + tpi + slope + roads_distance + ppt +
      tmin + tmmx + th + vpdmax + rmax + vs + fm100 + fm1000 + esp + lon + lat,
    data = train_test_list$train,
    spcov_type = "exponential", # this is the default already.
    local = FALSE, # no spatial approximation!
    mtry = 4,
    min.node.size = 2,
    sample.fraction = 0.89
  ), silent = TRUE)

  # Check if the "no variability" error has occurred. This happens when the
  # training data has no variation in CBI value. First check if an error
  # occurred, then check that the error message matches the expected one.
  if (inherits(local_spatial_rf_model, "try-error") &&
      any(grepl("The response has no variability",
                attr(local_spatial_rf_model, "condition")$message))
  ) {
    # if the model failed, return a vector of NAs
    message(
      sprintf(
        "[FID #%d] Warning: No variability in CBI.
        Returning uniform predictions.",
        train_test_list$validation_FID
      )
    )
    mean_cbi <- mean(train_test_list$train$cbi, na.rm = TRUE)
    return(list(
      validation_FID = train_test_list$validation_FID,
      predictions = rep(mean_cbi, nrow(train_test_list$test)),
      predictions_noweather = rep(mean_cbi, nrow(train_test_list$test))
    ))
  }

  # Only go on to fit the noweather model if the first model fit was successful.
  local_spatial_rf_noweather <- splmRF(
    cbi ~ elev + aspect + tri + tpi + slope + roads_distance + ppt +
      tmin + vpdmax + esp + lon + lat,
    data = train_test_list$train,
    spcov_type = "exponential", # this is the default already.
    local = FALSE, # no big data (spatial) approximation!
    mtry = 4,
    min.node.size = 2,
    sample.fraction = 0.89
  )

  # predict for the validation plot. for speed, just focus on the pixels in the
  # validation plot area.
  predictions <- predict(local_spatial_rf_model, newdata = train_test_list$test)

  predictions_noweather <- predict(local_spatial_rf_noweather,
                                   newdata = train_test_list$test)

  list(
    validation_FID = train_test_list$validation_FID,
    predictions = predictions,
    predictions_noweather = predictions_noweather
  )
}

# map this over the list of training and testing data to fit models and extract
# predictions
# Nate: tried parallelizing this, but I could never get it to work faster than
# serial execution. I think the bottleneck is memory.
all_preds <- purrr::map(train_test_each_plot,
                        fit_local_spatial_rf,
                        .progress = TRUE)

# extracts raster predictions from validation plot areas
get_mean_predictions <- function(preds_vec, train_test_one_plot) {
  # get id
  id <- train_test_one_plot$validation_FID
  # get preds vec id
  p_id <- preds_vec$validation_FID
  # check that the id matches
  if (id != p_id) {
    stop(sprintf("ID mismatch: %d != %d", id, p_id))
  }
  # get geo of validation plot
  one_plot <- validation_plots[validation_plots$FID_1 == id, ]
  # get test_sf
  test_sf <- train_test_one_plot$test
  # add predictions to test_sf
  test_sf$cbi_local_spatial_rf <- preds_vec$predictions
  test_sf$cbi_local_spatial_rf_noweather <- preds_vec$predictions_noweather

  # convert predictions to raster and extract by exact plot polygon
  predict_rast <- test_sf %>%
    vect %>%
    rasterize(rasts, field = "cbi_local_spatial_rf")
  predict_rast_noweather <- test_sf %>%
    vect %>%
    rasterize(rasts, field = "cbi_local_spatial_rf_noweather")

  # extract
  predict_mean <- exact_extract(predict_rast,
                                one_plot,
                                fun = "mean",
                                weights = "area")
  predict_mean_noweather <- exact_extract(predict_rast_noweather,
                                          one_plot,
                                          fun = "mean",
                                          weights = "area")
  # final output is the mean cbi value of the validation plot
  list(id = id,
       predictions = predict_mean,
       predictions_nw = predict_mean_noweather)
}

mean_preds <- purrr::map2(all_preds,
                          train_test_each_plot,
                          get_mean_predictions,
                          .progress = TRUE)

# bind predictions into data frame
mean_preds_df <- do.call(rbind, mean_preds) %>%
  as.data.frame() %>%
  rename(validation_FID = id)
mean_preds_df
# export
write_csv(mean_preds_df, "./processed_data/local_spatial_rf_predictions.csv")
