# this script makes simple models and plots the results of Effects_across_hpcc.R

# load packages ####
library(rnaturalearth) # map data
library(rnaturalearthdata) # map data
library(ggspatial)
library(sf)
library(terra)
library(tidyverse)
library(exactextractr)
library(cowplot)

# load data ####
## composite burn index (CBI) raster ####
cbi <- rast("./processed_data/masked_raster.tif")
# as usual, set "unmappable" -> NA and "outside of perimeter" -> 1 (unburned)
names(cbi) <- "cbi"
cbi[cbi == 9] <- NA
cbi[cbi == 0] <- 1
## vegetation treatments ####
veg_treatments <- read_sf("./results/hpcc_processed_cbi_new.shp")
## rf predictions ####
rf_preds <- read_csv("./processed_data/trt_rf_preds.csv")

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
