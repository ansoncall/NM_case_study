# this script tests the effect of different validation plot sizes on prediction
# accuracy. much of the code is borrowed from comparing_methods.R and
# comparing_methods_rf.R. see those scripts for more code comments.

# load packages ####
library(tidyverse)
library(terra)
library(sf)
library(ggplot2)
library(snapKrig)
library(spmodel)
library(parallel)
library(caret)
library(exactextractr)
library(cowplot)

# set options ####
# set testing = TRUE to reduce neighborhood size for local spatial rf models.
# will reduce runtime.
testing <- TRUE

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

# wrangle ####
# calculate buffer size for making plots of various acreages
plot_sizes <- data.frame(acres = c(2, 5, 15, 50, 100, 150, 250))
plot_sizes$radius_m <- sqrt((plot_sizes$acres * 4046.86) / pi)

# only place points in burned areas and not in vegetation treatments
# vectorize cbi raster to define burned area
cbi_poly <- as.polygons(cbi) %>%
  st_as_sf %>%
  # don't include unburned areas
  filter(cbi != 1) %>%
  st_union %>%
  # avoid the very edge of unburned areas or unmappable cbi areas
  st_buffer(-60)

# remove veg treatments from burned area poly
veg_buffer <- veg_treatments %>% st_union %>% st_buffer(60)
valid_area <- st_difference(cbi_poly, veg_buffer)

# random seed
set.seed(2076)

# select 60 random points of each size
random_plots <- st_sample(valid_area, size = nrow(plot_sizes) * 60) %>%
  # buffer and reproject
  st_buffer(dist = rep(plot_sizes$radius_m)) %>%
  st_transform(crs = st_crs(cbi))

# rf models ####
build_train_test <- function(id, test = testing) {
  one_plot <- random_plots[id] %>% st_as_sf
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
  rasts_crop <- rasts %>%
    crop(st_bbox(neighborhood_poly))
  masks_raster <- dplyr::bind_rows(
    target_poly %>% mutate(mask_val = 1),
    neighborhood_poly %>% mutate(mask_val = 2)
  ) %>%
    rasterize(rasts_crop, field = "mask_val")
  target_mask <- ifel(masks_raster == 1, 1, NA)
  neighborhood_mask <- ifel(masks_raster == 2, 1, NA)
  train_pts <- rasts_crop %>%
    mask(neighborhood_mask) %>%
    as.points
  train_xy <- crds(train_pts)
  train_pts$lon <- train_xy[, 1]
  train_pts$lat <- train_xy[, 2]
  train <- train_pts %>%
    st_as_sf(coords = c("x", "y"), crs = st_crs(rasts), remove = FALSE)
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

train_test_data <- map(seq_along(random_plots),
                       build_train_test,
                       .progress = TRUE)

# check sizes of train and test data. random points are placed in valid areas,
# but polgyons may extend outside of those areas after buffering.
walk(train_test_data, function(x) {
  if (nrow(x$train) < 10 || nrow(x$test) < 10) {
    message(sprintf(
      "[FID #%d] Warning: %d test values for this plot",
      x$validation_FID, nrow(x$test)
    )
    )
  }
})

fit_local_spatial_rf <- function(train_test_list, test = testing) {

  local_spatial_rf_model <- try(splmRF(
    # this is the "no weather" model.
    cbi ~ elev + aspect + tri + tpi + slope + roads_distance + ppt +
      tmin + vpdmax + esp + lon + lat,
    data = train_test_list$train,
    spcov_type = "exponential",
    local = list(parallel = TRUE, ncores = detectCores() - 4),
    mtry = 4,
    min.node.size = 2,
    sample.fraction = 0.89
  ), silent = TRUE)

  if (inherits(local_spatial_rf_model, "try-error") &&
      any(grepl("The response has no variability",
                attr(local_spatial_rf_model, "condition")$message))
  ) {
    message(
      sprintf(
        c("[FID #%d] Warning: No variability in CBI.",
          " Returning uniform predictions."),
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

  predictions <- predict(local_spatial_rf_model, newdata = train_test_list$test)

  list(
    validation_FID = train_test_list$validation_FID, predictions = predictions
  )
}

# apply function to fit rf models and extract predicted values.
all_preds <- purrr::map(train_test_data, fit_local_spatial_rf, .progress = TRUE)

# extracts raster predictions from validation plot areas
get_mean_predictions <- function(preds_vec, id) {
  # get id. just a number from seq_along(train_test_data)
  train_test_one_plot <- train_test_data[[id]]
  # get preds vec id. #TODO rename build_train_test output so its no longer
  # "validation_FID" since we're not working with validation plots here.
  p_id <- preds_vec$validation_FID
  # check that the id matches
  if (id != p_id) {
    stop(sprintf("ID mismatch: %d != %d", id, p_id))
  }
  # get geo of one test plot
  one_plot <- random_plots[id, ]
  # get test_sf
  test_sf <- train_test_one_plot$test
  # add predictions to test_sf
  test_sf$cbi_local_spatial_rf <- preds_vec$predictions

  # convert predictions to raster and extract by exact plot polygon
  predict_rast <- test_sf %>%
    vect %>%
    rasterize(rasts, field = "cbi_local_spatial_rf")

  # extract
  predict_mean <- exact_extract(predict_rast,
                                one_plot,
                                fun = "mean",
                                weights = "area")
  # final output is the mean cbi value of the validation plot
  list(id = id,
       predictions = predict_mean)
}

# apply function to extract mean predictions.
mean_preds <- purrr::map2(all_preds,
                          seq_along(train_test_data),
                          get_mean_predictions,
                          .progress = TRUE)

# bind predictions into data frame
mean_preds_df <- do.call(rbind, mean_preds) %>%
  as.data.frame() %>%
  # unlist results and rename
  mutate(id = unlist(id), cbi_local_spatial_rf = unlist(predictions))

# cluster-based matching ####
# define off-limits areas: treated areas and validation plots
off_limits <- rbind(st_as_sf(st_geometry(veg_treatments),
                             crs = st_crs(random_plots)),
                    st_buffer(random_plots, 60))
# remove the gridded plots that touch off-limits areas
touch_mat <- st_intersects(gridded_plots, off_limits, sparse = FALSE)
non_touching_idx <- which(rowSums(touch_mat) == 0)

gridded_plots_clean <- gridded_plots[non_touching_idx, ] %>%
  # normalize variables
  mutate(across(c(elev, ppt, vs),
                ~ (. - mean(.)) / sd(.),
                .names = "{col}_norm"))

# find nearest neighbors
nearest_model <- knnreg(cbi ~ elev_norm + ppt_norm + vs_norm,
                        data = gridded_plots_clean,
                        k = 1)

# build df to use for prediction. needs the same normalized predictors.
predict_df <- st_as_sf(random_plots) %>%
  mutate(
    # extract values to validation plots
    elev = exact_extract(rasts$elev, ., fun = "mean", weights = "area"),
    vs = exact_extract(rasts$vs, ., fun = "mean", weights = "area"),
    ppt = exact_extract(rasts$ppt, ., fun = "mean", weights = "area"),
    # normalize
    elev_norm = (elev - mean(gridded_plots_clean$elev)) /
      sd(gridded_plots_clean$elev),
    ppt_norm = (ppt - mean(gridded_plots_clean$ppt)) /
      sd(gridded_plots_clean$ppt),
    vs_norm = (vs - mean(gridded_plots_clean$vs)) /
      sd(gridded_plots_clean$vs)
  ) %>%
  mutate(
    # fit model
    cbi_knn = predict(nearest_model, .)
  ) %>%
  # keep only predictions
  select(cbi_knn) %>%
  st_drop_geometry

# maybe solution? check row order
mean_preds_df <- cbind(mean_preds_df, predict_df)

# kriging ####
# takes a validation plot, crops the cbi raster to its neighborhood and masks
# out the plot itself, and applies kriging.
krig_plot <- function(id, test = testing) {# use seq_along(random_plot_in_burn)
  # grab a plot based on id
  one_plot <- random_plots[id, ] # nolint
  if (test == TRUE) {
    buffer_size <- 100
  } else {
    buffer_size <- 800
  }

  one_plot_buffer <- st_buffer(one_plot, buffer_size)

  clipped_burn_raster <- crop(cbi, one_plot_buffer) %>%
    mask(st_as_sf(st_buffer(one_plot, 60)), inverse = TRUE)

  corner_indices <- list(
    c(1, 1),
    c(nrow(clipped_burn_raster), 1),
    c(1, ncol(clipped_burn_raster)),
    c(nrow(clipped_burn_raster), ncol(clipped_burn_raster))
  )

  for (idx in corner_indices) {
    if (is.na(clipped_burn_raster[idx[1], idx[2]])) {
      clipped_burn_raster[idx[1], idx[2]] <- 1
    }
  }

  vals <- values(clipped_burn_raster, mat = FALSE)
  vals <- vals[!is.na(vals)]
  all_same <- length(unique(vals)) == 1
  if (all_same) {
    # if all values are the same, just return that value
    message(
      sprintf(
        c("[FID #%d] Warning: No variability in CBI.",
          " Returning uniform predictions."),
        id
      )
    )
    out <- vals[[1]]
  } else {
    # apply kriging
    clipped_burn_raster_sk <- sk(clipped_burn_raster)
    kriging_fit <- sk_fit(clipped_burn_raster_sk, n_max = 3000, quiet = TRUE)
    kriged_final <- sk_cmean(clipped_burn_raster_sk, kriging_fit)
    # build raster and extract values
    kriged_final_raster <- sk_export(kriged_final)
    out <- exact_extract(kriged_final_raster, one_plot, fun = "mean")
  }

  # return
  list(id, cbi_krig = out
  )
}

# apply. t
krig_results <- map(seq_along(random_plots), krig_plot, .progress = TRUE)

mean_preds_df <- krig_results %>%
  do.call(rbind, .) %>%
  as.data.frame %>%
  select(cbi_krig) %>%
  cbind(mean_preds_df, .)

# simple perimeters ####
# TODO

# perimeter area will be equal to the plot area.

# calculate correct outer radius depending on plot size
plot_sizes_with_perimeters <- plot_sizes %>%
  mutate(
    # get plot area in m
    plot_area_m = acres * 4046.86,
    # get "inner radius" (plot radius + 60m buffer)
    inner_radius = as.numeric(sqrt(plot_area_m / pi)) + 60,
    # get area of the plot + 60m buffer
    inner_circle_area = pi * inner_radius^2,
    # get area of the "outer circle" (plot + buffer + perimeter plot)
    outer_circle_area = inner_circle_area + plot_area_m,
    # get outer circle radius
    outer_radius = sqrt(outer_circle_area / pi)
  )
# repeat rows 60 times to make length match number of created plots
ps_allplots <- do.call(rbind, replicate(60,
                                        plot_sizes_with_perimeters,
                                        simplify = FALSE)) %>%
  # cbind to plots
  cbind(random_plots %>% st_sf)

# create "donut" perimeter plots of appropriate sizes, using
# plot_sizes_as_perimeters as key to defining buffer parameters

make_donuts <- function(plot_idx) {
  # plot is one row of ps_allplots. use map()
  plot <- ps_allplots[plot_idx, ] %>% st_sf
  inner <- plot %>% st_buffer(60) %>% st_geometry
  outer <- st_centroid(plot) %>% st_buffer(plot$outer_radius) %>% st_geometry
  donut <- st_difference(outer, inner) %>% st_sf
  true_cbi <- exact_extract(cbi, plot, fun = "mean", weights = "area")
  perimeter_cbi <- exact_extract(cbi, donut, fun = "mean", weights = "area")
  list(true_cbi = true_cbi, perimeter_cbi = perimeter_cbi)
}

# apply function. suppress annoying warnings from st_ functions
perimeter_pred <- suppressWarnings(
  map(seq_len(nrow(ps_allplots)), make_donuts, .progress = TRUE)
)
# add results to mean_preds_df
mean_preds_df <- perimeter_pred %>%
  do.call(rbind, .) %>%
  as.data.frame %>%
  cbind(mean_preds_df, .)

# also add size of plot in acres as column
mean_preds_df <- mean_preds_df %>%
  cbind(., acres = ps_allplots$acres) %>%
  # fix list cols
  # TODO stop this from happening in the first place
  mutate(across(everything(), unlist))

# write out
write_csv(mean_preds_df, "./results/plot_size_cross_validation.csv")

# plots ####
# read fresh csv so you can start from here if you want to
cv_data <- read_csv("./results/plot_size_cross_validation.csv")

# build rmse and mean error summary tables, summarizing by model and size
rmse_by_size <- cv_data %>%
  group_by(acres) %>%
  summarise(
    rmse_perimeter = sqrt(mean((true_cbi - perimeter_cbi)**2, na.rm = TRUE)),
    rmse_kriging = sqrt(mean((true_cbi - cbi_krig)**2, na.rm = TRUE)),
    rmse_spatial_rf = sqrt(mean((true_cbi - cbi_local_spatial_rf)**2,
                                na.rm = TRUE)),
    rmse_knn = sqrt(mean((true_cbi - cbi_knn)**2, na.rm = TRUE))
  )

mean_error_size <- cv_data %>%
  group_by(acres) %>%
  summarise(
    me_perimeter = mean((true_cbi - perimeter_cbi), na.rm = TRUE),
    me_kriging = mean((true_cbi - cbi_krig), na.rm = TRUE),
    me_spatial_rf = mean((true_cbi - cbi_local_spatial_rf), na.rm = TRUE),
    me_knn = mean((true_cbi - cbi_knn), na.rm = TRUE)
  )

# make "long" versions
rmse_long <- rmse_by_size %>% pivot_longer(2:5, names_to = "model")
me_long <- mean_error_size %>% pivot_longer(2:5, names_to = "model")

# print summary table by model across all sizes
rmse_long %>%
  group_by(model) %>%
  summarise(rmse = mean(value))

me_long %>%
  group_by(model) %>%
  summarise(mean_error = mean(value))

# plot
rmse_size <- ggplot(rmse_long, aes(x = acres, y = value, color = model)) +
  geom_point(size = 3) +
  geom_line(lty = 2) +
  theme_classic() +
  labs(y = "RMSE", x = "Validation Plot Size (acres)") +
  theme(text = element_text(size = 20)) +
  scale_color_manual(values = c("#009E73", "#0072B2", "#D55E00", "black")) +
  theme(legend.position = "none")

mean_error_size <- ggplot(me_long, aes(x = acres, y = value, color = model)) +
  geom_point(size = 3) +
  geom_line(lty = 2) +
  theme_classic() +
  labs(y = "Mean Error (actual-modeled)", x = "Validation Plot Size (acres)") +
  theme(text = element_text(size = 20), legend.position = c(.70, 0.85)) +
  scale_color_manual(values = c("#009E73", "#0072B2", "#D55E00", "black"),
                     labels = c("Propensity score matching",
                                "Kriging",
                                "Perimeter method",
                                "Local spatial RF no weather"),
                     name = "Method")

## write out plots ####
tiff(filename = ("./figures/size_cross_validation.tif"), units = "in",
     compression = "lzw", width = 12, height = 6, res = 300)
plot_grid(rmse_size, mean_error_size, labels = c("a", "b"))
dev.off()
