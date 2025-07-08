# this script compares the accuracy of methods used to create counterfactual
# outcomes. it is a core part of the analysis. because the random forest models
# take some time to fit, they are part of a separate script.

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
library(snapKrig)
library(conflicted)

# set options ####
conflicts_prefer(dplyr::filter)
testing <- TRUE # reduces sample size for testing purposes.

#load data ####
## composite burn index (CBI) raster ####
cbi <- rast("./processed_data/masked_raster.tif")
# as usual, set "unmappable" -> NA and "outside of perimeter" -> 1 (unburned)
cbi[cbi == 9] <- NA
cbi[cbi == 0] <- 1
## HPCC burn perimeter ####
burn_perimeter <- read_sf("./processed_data/burn_perimeter.shp") %>%
  st_transform(st_crs(cbi))
## vegetation treatments ####
veg_treatments <- read_sf(
  "./processed_data/vegetation_treatments_hpcc_new.shp"
)

## validation and control plots ####
# validation plots are the randomly selected plots used to evaluate the accuracy
# of the different methods. "control" plots are the practitioner-generated
# control plots. both are included in this file.
validation_and_control_plots <- read_sf(
  "./processed_data/val_points_revised.shp",
  fid_column_name = "ID"
) %>%
  st_transform(st_crs(cbi))
# "validation" plots only
validation_plots <- validation_and_control_plots %>%
  filter(as.integer(ID) < 352)
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
gridded_plots <- read_sf("./processed_data/gridded_candidate_plots.shp") %>%
  st_transform(st_crs(cbi))
## additional raster data ####
# LandFire ESP
site_potential <- rast("./processed_data/lf_site_potential_new.tif")
# set active category to 2, which is the fine-scale zone*esp*esplf
# categorization
activeCat(site_potential) <- 2
names(site_potential) <- "esp"
# cbi raster, masked by original plots + 60m buffer
validation_mask_cbi <- mask(cbi,
                            st_buffer(validation_plots, 60),
                            inverse = TRUE)
names(validation_mask_cbi) <- "masked_cbi"
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

# practitioner-generated controls ####
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
# demonstrate math for calculating "donut" statistics:
# we want to make a 10-acre donut, centered on the original plots.
# the inner radius of the donut is the radius of a 10-acre circle () + 60m:
one_plot <- validation_plots[1, ] # grab a sample plot
plot_area_m <- st_area(one_plot) # in m2, == 10 acres
inner_radius <- as.numeric(sqrt(plot_area_m / pi)) + 60 # all in meters
# the outer radius is the radius of a circle with that gives an area equal to 10
# acres (40450.11078 m2) + the area of the inner circle.
inner_circle_area <- inner_radius^2 * pi # in m2
outer_circle_area <- inner_circle_area + as.numeric(plot_area_m) # +10 acres
outer_radius <- sqrt(outer_circle_area / pi) # in m2
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


# kriging ####
# takes a validation plot, crops the cbi raster to its neighborhood and masks
# out the plot itself, and applies kriging.
krig_plot <- function(validation_plot_id) {
  # grab a plot based on id
  # one_plot <- validation_plots %>% filter(FID_1 == 18) # for testing # nolint
  # TODO fix "no visible binding"
  one_plot <- validation_plots %>% filter(FID_1 == validation_plot_id) # nolint
  # Nate, this was the original size. Could change or leave it, just make sure
  # its documented in the paper correctly. Hardcoding for speed here.
  # size_of_neighborhood <- sqrt((200*4046.86)/pi) # nolint
  if (testing == TRUE) {
    # for testing, use a smaller neighborhood. this runs in ~10 minutes.
    buffer_size <- 200
  } else {
    buffer_size <- 507 # this runs in ~5 hrs.
  }

  one_plot_buffer <- st_buffer(one_plot, buffer_size)

  clipped_burn_raster <- crop(cbi, one_plot_buffer) %>%
    mask(st_as_sf(st_buffer(one_plot, 60)), inverse = TRUE)

  # sk_fit fails when corner cells are NA. This sets it to 1 when it is NA.
  # Because corners are the farthest from the treatments this should have
  # minimal impact on predictions

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

  # Nate: it also fails when all non-NA cells have the same value, which can
  # happen easily when the neighborhood size is small. here, we check to see if
  # that is the case, and if so we skip over kriging and just assign that
  # singular value as the prediction.
  vals <- values(clipped_burn_raster, mat = FALSE)
  vals <- vals[!is.na(vals)]
  all_same <- length(unique(vals)) == 1
  if (all_same) {
    # if all values are the same, just return that value
    message(
      sprintf(
        c("[FID #%d] Warning: No variability in CBI.",
          " Returning uniform predictions."),
        validation_plot_id
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
  list(
    id = one_plot$FID_1,
    cbi_krig = out
  )
}

# apply. t
krig_results <- map(validation_plots$FID_1, krig_plot, .progress = TRUE) %>%
  bind_rows() %>%
  rename(plot_id = id, cbi_krig = cbi_krig)

# write out (if desired)
write.csv(krig_results, "./processed_data/kriging_results.csv")

# join to compare df
compare <- left_join(compare,
                     krig_results %>% mutate(plot_id = as.character(plot_id)),
                     by = c("FID_validation" = "plot_id")) %>%
  # rename for clarity
  rename(cbi_krigged = cbi_krig)

# cluster based matching ####
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
# define off-limits areas: treated areas and validation plots
off_limits <- rbind(st_as_sf(st_geometry(veg_treatments),
                             crs = st_crs(validation_plots)),
                    st_as_sf(st_geometry(st_buffer(validation_plots, 60)),
                             crs = st_crs(validation_plots)))
# remove the gridded plots that touch off-limits areas
touch_mat <- st_intersects(gridded_plots, off_limits, sparse = FALSE)
non_touching_idx <- which(rowSums(touch_mat) == 0)


gridded_plots_clean <- gridded_plots[non_touching_idx, ] %>%
  # normalize variables
  mutate(across(c(elev, ppt, vs),
                ~ (. - mean(.)) / sd(.),
                .names = "{col}_norm"))
# Nate: validation and control plots contain (obviously) the validation plots,
# but also the practitioner-generated control plots. There's no need to make
# these areas off-limits here. In fact, this could be one of the reasons this
# method underperforms: the best gridded plots, based on PSM, are likely going
# to be near the practitioner-selected control. I know a lot has changed in this
# script, but I am pretty sure this error is present in the original and not
# introduced by me.
# mapview::mapview(gridded_plots_clean) +
#   mapview::mapview(off_limits, color = "red", alpha = 0.5) # nolint

# find nearest neighbors
nearest_model <- knnreg(cbi ~ elev_norm + ppt_norm + vs_norm,
                        data = gridded_plots_clean,
                        k = 1)
summary(nearest_model)

# build df to use for prediction. needs the same normalized predictors.
predict_df <- validation_plots %>%
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
  select(ID, cbi_knn)

compare <- left_join(compare,
                     predict_df,
                     by = c("FID_validation" = "ID")) %>%
  rename(geometry = geometry.x,
         matched_plot_geometry = geometry.y)
compare %>% names
# write out
compare %>%
  st_drop_geometry %>%
  select(FID_validation, starts_with("cbi")) %>%
  write.csv("./processed_data/compare_results.csv")
