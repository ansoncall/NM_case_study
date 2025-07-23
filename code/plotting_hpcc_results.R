# this script makes simple models and plots the results of Effects_across_hpcc.R

# load packages ####
library(ggspatial)
library(sf)
library(terra)
library(tidyverse)
library(exactextractr)
library(cowplot)

# load data ####
## composite burn index (CBI) raster ####
cbi <- rast("./processed_data/clipped_burn_raster.tif")
# as usual, set "unmappable" -> NA and "outside of perimeter" -> 1 (unburned)
names(cbi) <- "cbi"
cbi[cbi == 9] <- NA
cbi[cbi == 0] <- 1
## vegetation treatments ####
veg_treatments <- read_sf("./results/hpcc_processed_cbi_new.shp")
## rf predictions ####
rf_preds <- read_csv("./processed_data/trt_rf_preds.csv")
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

# wrangle ####

# categorize fuel treatments based on the description field
descriptions <- tolower(unique(veg_treatments$Dscrptn))
# extract treatment descriptions based on keywords
thinning <- descriptions[grepl("thin", descriptions)]
cutting <- descriptions[grepl("cut", descriptions)]
defens <- descriptions[grepl("defens", descriptions)]
burn <- descriptions[grepl("burn", descriptions)]
fuel <- descriptions[grepl("fuel", descriptions)]
fire <- descriptions[grepl("fire", descriptions)]
# define unique set of treatment descriptions containing any of the keywords
treatments_to_use <- unique(c(thinning, cutting, defens, burn, fuel, fire))
# define categories "burning" and "thinning"
burning <- unique(c(burn, fuel))
thin <- unique(c(thinning, cutting))
# add a new column to veg_treatments to indicate whether the treatment is a burn
# or thin
veg_treatments <- veg_treatments %>%
  mutate(thin_burn = ifelse(tolower(Dscrptn) %in% burning, "burn", "thin")) %>%
  # subset veg_treatments to only those with descriptions in treatments_to_use
  filter(tolower(Dscrptn) %in% treatments_to_use)

# extract CBI values for each treatment polygon
veg_treatments$actual_severity <- exact_extract(cbi,
                                                veg_treatments,
                                                fun = "mean",
                                                weights = "area")
# add the RF predictions to the veg_treatments data frame
veg_treatments <- veg_treatments %>%
  left_join(rf_preds, by = "rownum") %>%
  # rename the RF prediction column
  rename(cntrl__ = cbi_local_spatial_rf) %>%
  # calculate treatment effect as the difference between actual severity and
  # counterfactual prediction
  mutate(effect = actual_severity - cntrl__,
         effect_knn = actual_severity - knn_pred) %>%
  # filter out rows where Year_Cl is "needs input"
  filter(Year_Cl != "needs input")

# t tests: is the estimated effect of treatment non-zero?
t.test(veg_treatments$effect)
t.test(veg_treatments$effect_knn)

# lm: does the treatment effect size vary with expected burn severity?
summary(lm(effect ~ cntrl__, data = veg_treatments))
# with treatment size?
summary(lm(effect ~ Acre_US, data = veg_treatments))
# with treatment type?
summary(lm(effect ~ thin_burn, data = veg_treatments))

# What prop of treated area burned at lower severity than counterfactually
# predicted?
lower <- veg_treatments[veg_treatments$cntrl__ >
                          veg_treatments$actual_severity, ]
sum(st_area(lower)) / sum(st_area(veg_treatments))
# What prop of treated area burned at higher severity than counterfactually
# predicted?
high <- veg_treatments[veg_treatments$actual_severity > 3, ]
sum(st_area(high)) / sum(st_area(veg_treatments))

# plots ####
first_plot_colors <- c("#E69F00", "#56B4E9", "#009E73")
second_plot_colors <- c("#0072B2", "#D55E00", "#CC79A7")

####

effect_hist <- ggplot(veg_treatments, aes(x = effect)) +
  geom_histogram(fill = "#0072B2", alpha = 0.5) +
  geom_histogram(data = veg_treatments,
                 aes(x = effect_knn),
                 fill = "#D55E00",
                 alpha = 0.5) +
  theme_classic() +
  ylab("Number of treatments") +
  xlab(expression("Effect of Treatment (" * Delta * "CBI)")) +
  theme(text = element_text(size = 20)) +
  geom_vline(xintercept = 0, linetype = "dashed")

effect_hist

scatter_effect <- ggplot(veg_treatments, aes(x = cntrl__, y = effect)) +
  geom_point(size = 2) +
  theme_classic() +
  xlab("Untreated burn seveirty (CBI)") +
  ylab(expression("Effect of Treatment (" * Delta * "CBI)")) +
  theme(text = element_text(size = 20)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_smooth(method = "lm")

scatter_effect

tiff(filename = "./figures/hpcc_results_simple.tif", units = "in",
     compression = "lzw", width = 8, height = 12, res = 300)
plot_grid(effect_hist, scatter_effect, ncol = 1, labels = "auto",
          label_x = 0.96, label_y = 0.95, label_size = 20, align = "v")
dev.off()

########### Make a map
# TODO tidy this part

# using veg_treatments_buffer?
masked_cbi <- mask(masked_cbi, veg_treatments_buffer, inverse=TRUE)

rf_preds <- read.csv("./processed_data/trt_rf_preds.csv")

masked_df <- as.points(masked_cbi)

names(modeled_values_across_hpcc)[1] <- "x"
names(modeled_values_across_hpcc)[2] <- "y"
names(modeled_values_across_hpcc)[3] <- "masked_raster"
masked_df_full<-rbind(masked_df,modeled_values_across_hpcc)

counterfactual_raster<-rasterFromXYZ(masked_df_full)

one_point<-data.frame(lat=-832592.7,lon=1458797)
one_point<-st_as_sf(one_point,coords=c("lat","lon"),crs=crs(masked_cbi))
one_region<-st_buffer(one_point,1100)
cropped_counterfactual_raster<-raster::crop(counterfactual_raster,st_as_sf(one_region))

cropped_actual_raster<-raster::crop(CBI,st_as_sf(one_region))

cropped_veg_treatments<-st_crop(veg_treatments,one_region)

plotting_actual<-rasterToPoints(cropped_actual_raster)
plotting_counter<-rasterToPoints(cropped_counterfactual_raster)

# what to do:

# identify the treatment specified in the plot
CBI<-rast("./raw_data/burn_severity/ravg_2022_cbi4.tif")
point <- data.frame(lat=-832592.7,lon=1458797, id = 1) %>% st_as_sf(coords=c("lon", "lat"), crs=st_crs(CBI))
mapview::mapview(point) + mapview::mapview(veg_treatments)

# define bbox of mapped region. should be a little bigger than what's actually
# mapped. just need points for top left and bottom right corners.
plot_bbox <- data.frame(
  id = 1:2,
  lon = c(-105.32816, -105.2966),
  lat = c(35.82756, 35.79876)
) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_bbox %>%
  st_as_sfc %>%
  st_transform(st_crs(rasts))

# select all treatments within bbox
veg_treats_for_plot <- st_filter(veg_treatments, plot_bbox) %>%
  st_transform(st_crs(rasts))

# refit the local spatial rf for plot subset
# makes train and test data for each plot. works with treatment polys.
build_train_test <- function(trt_rownum, test = FALSE) {
  # grab one plot (one treatment poly).
  # TODO fix "no visible binding for rownum"
  one_plot <- veg_treats_for_plot %>% filter(rownum == trt_rownum)

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
    mask(target_mask) %>%   # will also include 60m buffer but this is ok.
    as.points(na.all = TRUE)

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
train_test_each_plot <- purrr::map(veg_treats_for_plot$rownum,
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
# map over the result to generate prediction rasters from the list of predicted
# values for each treatment
pred_rasts <- purrr::map2(
  train_test_each_plot_sub,
  all_preds,
  \(x, y) {
    sf_test <- x$test # spatialpoints of test data
    preds <- y$predictions # list of predictions

    sf_test <- cbind(sf_test, preds)
    rasterize(sf_test, rasts,  field = "preds")

  }
)
# composite the rasters (some of the overlap) taking the minimum value

# Nate: check. there are 17 treatments in the map pane and many overlap. each
# model is fit independently. The effects of multiple treatments are probably
# not additive, but we don't have a model for this. Minimum value composite
# seems reasonable to me.
# If rasters have different extents, mosaic first
min_composite <- app(sds(pred_rasts), "min", na.rm = TRUE)

# make the graphs so they look the same

bbox_view <- st_bbox(st_buffer(veg_treats_for_plot, 60))
counterfactual_graph <- ggplot() +
  geom_spatraster(data = crop(cbi, bbox_view), aes(fill = cbi)) +
  geom_spatraster(data = crop(mask(min_composite, veg_treats_for_plot),
                              bbox_view),
                  aes(fill = last)) +
  geom_sf(data = veg_treats_for_plot, fill = NA, color = "gray15",
          linewidth = 0.5, alpha = 0.6) +
  scale_fill_gradient(guide="none",low="steelblue",high="firebrick",
                      na.value = "transparent")+
  coord_sf(xlim = c(bbox_view["xmin"], bbox_view["xmax"]),
           ylim = c(bbox_view["ymin"], bbox_view["ymax"])) +
  theme_classic() +
  labs(x = NULL, y = NULL) +
  ggtitle("Counterfactual Burn Severity")

counterfactual_graph

actual_graph <- ggplot() +
  geom_spatraster(data = crop(cbi, bbox_view), aes(fill = cbi)) +
  geom_sf(data = veg_treats_for_plot, fill = NA, color = "gray15",
          linewidth = 0.5, alpha = 0.6) +
  scale_fill_gradient(guide="none",low="steelblue",high="firebrick")+
  coord_sf(xlim = c(bbox_view["xmin"], bbox_view["xmax"]),
           ylim = c(bbox_view["ymin"], bbox_view["ymax"])) +
  theme_classic() +
  labs(x = NULL, y = NULL) +
  ggtitle("Actual Burn Severity") +
  annotation_scale(location = "bl", width_hint = 0.4) +
  annotation_north_arrow(
    location = "bl", which_north = "true", pad_x = unit(0.3, "in"),
    pad_y = unit(0.2, "in"), style = north_arrow_fancy_orienteering
  )

actual_graph

tiff(filename = ("./figures/map_comparison.tif"), units = "in",
     compression = "lzw", width = 10, height = 5, res = 300)
plot_grid(counterfactual_graph, actual_graph, ncol = 2, labels = "auto",
          label_size = 20, align = "v")
dev.off()
