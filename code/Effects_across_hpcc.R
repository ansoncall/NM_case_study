# evaluate treatment effects across the HPCC burn area, using a couple of the
# different counterfactual methods.

# load packages ####
library(terra)
library(sf)
library(tidyverse)
library(exactextractr)
library(spmodel)
library(caret)

# set options ####
testing <- TRUE # subsamples to reduce run time.

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

# local spatial random forest ####
# this includes "local spatial" rf models with and without weather variables.

# makes train and test data for each plot. works with treatment polys.
build_train_test <- function(trt_rownum, test = testing) {
  # grab one plot (one treatment poly).
  # TODO fix "no visible binding for rownum"
  one_plot <- veg_treatments %>% filter(rownum == trt_rownum)

  # define buffers
  target_poly <- st_buffer(one_plot, 60)

  if (test == TRUE) {
    buffer_size <- 100
  } else {
    buffer_size <- 800
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
    trt_rownum = trt_rownum,
    train = train,
    test = test
  ))
}

# map function over all validation plot FIDs
train_test_each_plot <- purrr::map(veg_treatments$rownum,
                                   build_train_test,
                                   .progress = TRUE)

# check sizes of train and test data
small_train_test <- c()
walk(train_test_each_plot, function(x) {
  # check that the train and test data are not empty
  if (nrow(x$train) < 5 || nrow(x$test) == 0) {
    message(sprintf(
      "[rownum #%d] Warning: %d test values and %d train for this plot",
      x$trt_rownum, nrow(x$test), nrow(x$train)
    )
    )
    small_train_test <<- c(small_train_test, x$trt_rownum)
  }
})
# remove treatments with small train or test data
`%ni%` <- Negate(`%in%`)
train_test_each_plot_sub <- Filter(\(x) x$trt_rownum %ni% small_train_test,
                                   train_test_each_plot)

# takes train and test sfs and returns prediction sf
fit_local_spatial_rf <- function(train_test_list) {
  # train_test_list <- train_test_each_plot[[1]] # testing # nolint
  local_spatial_rf_model <- try(splmRF(
    # no weather included here
    cbi ~ elev + aspect + tri + tpi + slope + roads_distance + ppt +
      tmin + vpdmax + esp + lon + lat,
    data = train_test_list$train,
    spcov_type = "exponential", # this is the default already.
    local = list(parallel = TRUE, ncores = parallel::detectCores() - 4),
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
    # if the model failed due to no variability in predictor data, return a
    # prediction based on the mean (only) CBI value in the training data.
    message(
      sprintf(
        c("[FID #%d] Warning: No variability in CBI.",
          " Returning uniform predictions."),
        train_test_list$validation_FID
      )
    )
    mean_cbi <- mean(train_test_list$train$cbi, na.rm = TRUE)
    out <- list(
      trt_rownum = train_test_list$trt_rownum,
      predictions = rep(mean_cbi, nrow(train_test_list$test))
    )
  } else {
    # predict for the validation plot using the rf model.
    predictions <- predict(local_spatial_rf_model,
                           newdata = train_test_list$test)
    out <- list(
      trt_rownum = train_test_list$trt_rownum,
      predictions = predictions
    )
  }
  out
}

# map this over the list of training and testing data to fit models and extract
# predictions
all_preds <- purrr::map(train_test_each_plot_sub,
                        fit_local_spatial_rf,
                        .progress = TRUE)

# extracts raster predictions from validation plot areas.
get_mean_predictions <- function(preds_vec, train_test_one_plot) {
  # get id
  id <- train_test_one_plot$trt_rownum
  # get preds vec id
  p_id <- preds_vec$trt_rownum
  # check that the id matches
  if (id != p_id) {
    stop(sprintf("ID mismatch: %d != %d", id, p_id))
  }
  # get geo of validation plot
  one_plot <- veg_treatments[veg_treatments$rownum == id, ]
  # get test_sf
  test_sf <- train_test_one_plot$test
  # add predictions to test_sf
  test_sf$cbi_local_spatial_rf <- preds_vec$predictions

  # convert predictions to raster and extract by exact plot polygon
  predict_rast <- test_sf %>%
    vect %>%
    rasterize(rasts, field = "cbi_local_spatial_rf")

  # extract
  predict_mean <- exactextractr::exact_extract(predict_rast,
                                               one_plot,
                                               fun = "mean",
                                               weights = "area")
  # final output is the mean cbi value of the validation plot
  list(id = id,
       mean_predictions = predict_mean)
}

mean_preds <- purrr::map2(all_preds,
                          train_test_each_plot_sub,
                          get_mean_predictions,
                          .progress = TRUE)

# bind predictions into data frame
mean_preds_df <- do.call(rbind, mean_preds) %>%
  as.data.frame() %>%
  rename(rownum = id,
         cbi_local_spatial_rf = predictions) %>%
  # unlist results
  mutate(
    rownum = unlist(rownum),
    cbi_local_spatial_rf = unlist(cbi_local_spatial_rf)
  )

# export
write_csv(mean_preds_df, "./processed_data/trt_rf_preds.csv")

# cluster-based matching ####

# remove the gridded plots that touch off-limits areas
touch_mat <- st_intersects(gridded_plots, veg_treatments, sparse = FALSE)
non_touching_idx <- which(rowSums(touch_mat) == 0)

gridded_plots_clean <- gridded_plots[non_touching_idx, ] %>%
  # normalize variables
  mutate(across(c(elev, ppt, vs),
                ~ (. - mean(.)) / sd(.),
                .names = "{col}_norm"))

# find nearest neighbor
nearest_model <- knnreg(cbi ~ elev_norm + ppt_norm + vs_norm,
                        data = gridded_plots_clean,
                        k = 1)
summary(nearest_model)

# Save means and sds of predictors for all gridded plots. Required to
# z-transform treatment polygon values at the same scale.
vars <- c("elev", "vs", "ppt")
mu_sigma <- lapply(vars, function(v) {
  c(mu = mean(gridded_plots_clean[[v]]),
    sigma = sd(gridded_plots_clean[[v]]))
})
names(mu_sigma) <- vars

# extract mean values of predictors for each treatment polygon. bind columns to
# dataframe and apply z transform.
mean_predictor_vals <- map(vars, function(x) {
  # extract mean values for each predictor variable
  veg_treatments[[x]] <- exact_extract(rasts[[x]], veg_treatments,
                                       fun = "mean", weights = "area")
}) %>%
  set_names(vars) %>%
  bind_cols() %>%
  mutate(
    elev_norm = (elev - mu_sigma$elev["mu"]) / mu_sigma$elev["sigma"],
    vs_norm   = (vs   - mu_sigma$vs["mu"])   / mu_sigma$vs["sigma"],
    ppt_norm  = (ppt  - mu_sigma$ppt["mu"])  / mu_sigma$ppt["sigma"]
  )
# bind to veg_treatments
veg_treatments <- bind_cols(veg_treatments, mean_predictor_vals)
# predict from knn model
veg_treatments$knn_pred <- predict(nearest_model, veg_treatments)

# write over veg_treatments shapefile. now includes cluster-based matching
# method cbi predictions as attribute data.
st_write(veg_treatments, dsn = "./results/hpcc_processed_cbi_new.shp",
         append = FALSE)
