# evaluate treatment effects across the HPCC burn area, using a couple of the
# different counterfactual methods.

# load packages ####
library(terra)
library(sf)
library(tidyverse)
library(spmodel)

# set options ####
testing <- TRUE # subsamples to reduce run time.

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
burn_perimeter <- read_sf("./processed_data/burn_perimeter.shp")

## treatments ####
# Nate: this should be hpcc_new.shp but was hpcc.shp in the original. Not sure
# what damage this might have done.
veg_treatments <- read_sf(
  "./processed_data/vegetation_treatments_hpcc_new.shp"
  ) %>%
  # add rownum as first column
  mutate(rownum = row_number(), .before = everything())

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

# spatial random forest ####
# Nate: why are we redoing this? why don't we just load in the predictions we've already generated in the comparing_methods script?
# Answer: we need to do it for the treated areas, not for the validation plots.
# local spatial rf ####
# this includes "local spatial" rf models with and without weather variables.

# makes train and test data for each plot. works with treatment polys.
build_train_test <- function(trt_rownum, test = testing) {
  # grab one plot (one treatment poly).
  one_plot <- veg_treatments %>% filter(rownum == trt_rownum)

  # Nate: looks like you had a buffer size of 150 here. Much smaller than what
  # was used in the compare_methods script. Not sure what to use here.
  target_poly <- st_buffer(one_plot, 60)

  if (test == TRUE) {
    buffer_size <- 100
  } else {
    buffer_size <- 200
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
# remove treatments with small train or test data Nate: I don't know if this was
# your solution, as i didn't parse your code too carefully.  do you recall?
`%ni%` <- Negate(`%in%`)
train_test_each_plot_sub <- Filter(\(x) x$trt_rownum %ni% small_train_test, train_test_each_plot)

# takes train and test sfs and returns prediction sf
fit_local_spatial_rf <- function(train_test_list) {
  # train_test_list <- train_test_each_plot[[1]] # testing
  local_spatial_rf_model <- try(splmRF(
    # Nate: no weather included here. I think this is correct - check?
    cbi ~ elev + aspect + tri + tpi + slope + roads_distance + ppt +
      tmin + vpdmax + esp + lon + lat,
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
        c("[FID #%d] Warning: No variability in CBI.",
          " Returning uniform predictions."),
        train_test_list$validation_FID
      )
    )
    mean_cbi <- mean(train_test_list$train$cbi, na.rm = TRUE)
    return(list(
      trt_rownum = train_test_list$trt_rownum,
      predictions = rep(mean_cbi, nrow(train_test_list$test)),
      predictions_noweather = rep(mean_cbi, nrow(train_test_list$test))
    ))
  }

  # predict for the validation plot. for speed, just focus on the pixels in the
  # validation plot area.
  predictions <- predict(local_spatial_rf_model, newdata = train_test_list$test)

  list(
    trt_rownum = train_test_list$trt_rownum,
    predictions = predictions
  )
}

# map this over the list of training and testing data to fit models and extract
# predictions
all_preds <- purrr::map(train_test_each_plot_sub,
                        fit_local_spatial_rf,
                        .progress = TRUE)

# extracts raster predictions from validation plot areas
get_mean_predictions <- function(preds_vec, train_test_one_plot) {
  # preds_vec <- all_preds[[31]] # testing
  # train_test_one_plot <- train_test_each_plot_sub[[31]] # testing
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
  predict_mean <- exact_extract(predict_rast,
                                one_plot,
                                fun = "mean",
                                weights = "area")
  # final output is the mean cbi value of the validation plot
  list(id = id,
       predictions = predict_mean)
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
head(mean_preds_df)
# export
write_csv(mean_preds_df, "./processed_data/trt_rf_preds.csv")

# TODO tidy below

descriptions<-tolower(unique(veg_treatments$Dscrptn))

thinning<-descriptions[grepl("thin",descriptions)]
cutting<-descriptions[grepl("cut",descriptions)]
defens<-descriptions[grepl("defens",descriptions)]
burn<-descriptions[grepl("burn",descriptions)]
fuel<-descriptions[grepl("fuel",descriptions)]
fire<-descriptions[grepl("fire",descriptions)]


treatments_to_use<-unique(c(thinning,cutting,defens,burn,fuel,fire))


veg_treatments<-veg_treatments[tolower(veg_treatments$Dscrptn) %in% treatments_to_use,]


veg_treatments$control_burn_severity<-NA
map_values<-data.frame()

l<-nrow(veg_treatments)

veg_treatments<-st_cast(veg_treatments,"POLYGON")

for (i in 1:nrow(veg_treatments)){

  print(i/l*100)
  output<-spatial_rf_interative(veg_treatments[i,])


  veg_treatments$control_burn_severity[i]<-output[[1]]
  map_values<-rbind(map_values,output[2][[1]])

  st_write(veg_treatments,dsn = "./results/hpcc_processed_cbi_new.shp",append=FALSE)
  write.csv(map_values,"./results/mapping_predictions.csv")
}

#### adding in propensity score matching here


veg_treatments<-read_sf("./results/hpcc_processed_cbi_new.shp")

gridded_plots<-read_sf("gridded_candidate_plots.shp")


gridded_plots<-st_transform(gridded_plots,st_crs(CBI))
#3control_plots<-st_transform(control_plots,st_crs(CBI))

gridded_plots_2<-st_difference(gridded_plots,st_union(veg_treatments))


gridded_plots_3<-st_difference(gridded_plots_2,st_union(st_buffer(validation_plots,60)))

max_area<-max(st_area(gridded_plots_3))

gridded_plots_4<-gridded_plots_3[as.numeric(st_area(gridded_plots_3))>as.numeric(max_area)-50,]

library(caret)
#### normalizing variables first

gridded_plots_4$elevation_norm<-(gridded_plots_4$elevatn-mean(gridded_plots_4$elevatn))/sd(gridded_plots_4$elevatn)
gridded_plots_4$ppt_norm<-(gridded_plots_4$nrml_pp-mean(gridded_plots_4$nrml_pp))/sd(gridded_plots_4$nrml_pp)
gridded_plots_4$vs_norm<-(gridded_plots_4$vs-mean(gridded_plots_4$vs))/sd(gridded_plots_4$vs)


nearest_model<-knnreg(cbi~elevation_norm+ppt_norm+vs_norm,data=gridded_plots_4,k=1)
summary(nearest_model)

ppt<-raster("./processed_data/ppt.tif")
elev_down<-raster("./processed_data/elev_down.tif")
vs<-raster("./processed_data/vs.tif")

predict_df<-st_as_sf(veg_treatments)

predict_df$elev<-raster::extract(elev_down,y=predict_df,fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)
predict_df$vs<-raster::extract(vs,y=st_as_sf(predict_df),fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)
predict_df$ppt<-raster::extract(ppt,y=st_as_sf(predict_df),fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)

predict_df$elevation_norm<-(predict_df$elev-mean(gridded_plots_4$elevatn))/sd(gridded_plots_4$elevatn)
predict_df$ppt_norm<-(predict_df$ppt-mean(gridded_plots_4$nrml_pp))/sd(gridded_plots_4$nrml_pp)
predict_df$vs_norm<-(predict_df$vs-mean(gridded_plots_4$vs))/sd(gridded_plots_4$vs)



predict_df$knn_1<-predict(nearest_model,as.data.frame(predict_df))


predict_df<-as.data.frame(predict_df)

veg_treatments$knn_pred<-predict_df$knn_1


st_write(veg_treatments,dsn = "./results/hpcc_processed_cbi_new.shp",append=FALSE)



