# load packages ####
library(tidyverse)
library(cowplot)
library(ggpmisc)

# read in results ####
compare_df <- read.csv("./processed_data/compare_results.csv") %>%
  rename(validation_FID = FID_validation)
compare_rf_df <- read.csv("./processed_data/rf_preds.csv")
compare_rf_local_df <- read.csv(
  "./processed_data/local_spatial_rf_predictions.csv"
)
# join
compare_all <- compare_df %>%
  left_join(compare_rf_df, by = "validation_FID") %>%
  left_join(compare_rf_local_df, by = "validation_FID")

# rmse and mean difference functions
calc_rmse <- function(predicted, actual) {
  sqrt(mean((predicted - actual)^2, na.rm = TRUE))
}
calc_mean_diff <- function(predicted, actual) {
  mean(predicted - actual, na.rm = TRUE)
}

# Model names and columns
model_info <- list(
  pg_control = "cbi_control",
  perimeter = "cbi_perimeter",
  rf = "cbi_rf",
  sp_rf = "cbi_rf_global_spatial",
  krig = "cbi_krigged",
  local_spatial_rf = "cbi_local_spatial_rf_noweather",
  local_spatial_rf_weather = "cbi_local_spatial_rf",
  knn = "cbi_knn"
)

# bootstrap resampling and metric calculation
results_df <- map_dfr(
  1:5000,
  function(i) {
    sample_df <- compare_all %>%
      slice_sample(n = nrow(.), replace = TRUE)

    # calculate all metrics for this bootstrap sample
    model_metrics <- map_dfc(names(model_info), function(model) {
      pred_col <- model_info[[model]]
      predicted <- sample_df[[pred_col]]
      actual <- sample_df$cbi_validation

      tibble(!!paste0(model, "_rmse") := calc_rmse(predicted, actual),
             !!paste0(model, "_mean_diff") := calc_mean_diff(predicted, actual))
    }
    )

    model_metrics
  },
  .progress = TRUE
)

head(results_df)

graphing_comparison <- results_df %>%
  summarise(across(everything(), list(
    mean = mean,
    lower_ci = ~ quantile(., 0.025),
    upper_ci = ~ quantile(., 0.975)
  ), .names = "{col}__{fn}")) %>%
  # Nate: this is sort of hacky and the pivoting/string splitting could probably
  # be done in a smarter way, but it works.
  pivot_longer(everything(),
               names_to = c("model_type", ".value"),
               names_sep = "__") %>%
  # split off model_type column at "mean" or "rmse_diff". I honestly tried to
  # get a more elegant solution with separate_wider_regex but couldn't figure it
  # out.
  separate_wider_delim(
    model_type,
    names = c("model_type", "metric"),
    delim = regex("_me|_rm")
  ) %>%
  mutate(
    model_type = recode(
      model_type,
      pg_control = "practitioner selected plots",
      perimeter = "perimeter method",
      rf = "random forest",
      sp_rf = "spatial random forest",
      krig = "kriging",
      local_spatial_rf = "local spatial random forest without weather",
      local_spatial_rf_weather = "local spatial random forest",
      knn = "propensity score matching"
    ),
    metric = recode(
      metric,
      se = "rmse",
      an_diff = "mean_diff"
    )
  )

rmse_plot <- ggplot(graphing_comparison %>% filter(metric == "rmse"),
                    aes(x = model_type, y = mean)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.1) +
  ylab("RMSE of CBI") +
  xlab("") +
  theme(text = element_text(size = 13)) +
  scale_x_discrete(labels = \(x) str_wrap(x, width = 6))
rmse_plot

mean_difference_plot <- ggplot(graphing_comparison %>%
                                 filter(metric == "mean_diff"),
                               aes(x = model_type, y = mean)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci),
                width = 0.1) +
  ylab("Mean Error (actual-modeled)") +
  xlab("") +
  theme(text = element_text(size = 13)) +
  scale_x_discrete(labels = \(x) str_wrap(x, width = 6))

mean_difference_plot

write.csv(graphing_comparison, "./results/model_comparisons.csv")

write.csv(compare_all, "./processed_data/final_comparsion_results_new.csv")

tiff(filename = "./figures/comparison_of_methods.tif",
     units = "in", compression = "lzw", width = 7, height = 10, res = 300)
plot_grid(rmse_plot, mean_difference_plot, ncol = 1, labels = c("a", "b"))
dev.off()

# lm plots ####
# build a correlation plot
make_plot <- function(xvar, xname) {
  # TODO fix "no visible binding"
  ggplot(compare_all, aes(y = cbi_validation, x = !!sym(xvar))) + # nolint
    geom_abline(intercept = 0, slope = 1, size = 2, color = "gray45") +
    geom_smooth(method = "lm") +
    geom_point() +
    theme_classic() +
    ylab("Validaiton plot burn severity (CBI)") +
    xlab(xname) +
    # Nate: can avoid using separate lm() call with stat_poly_eq
    stat_poly_eq(use_label("eq"),
                 label.y = 0.15,
                 label.x = 0.85,
                 coef.digits = 3) +
    xlim(1, 4)
}

xvars <- compare_all %>% select(-c(X:cbi_validation)) %>% names
xnames <- c(
  "Practitioner-generated controls (CBI)",
  "Simple perimeter method (CBI)",
  "Kriging method (CBI)",
  "Propensity score matching (CBI)",
  "Random forest (CBI)",
  "Spatial random forest (CBI)",
  "Local spatial random forest (CBI)",
  "Local spatial random forest without weather (CBI)"
)

plots <- map2(xvars, xnames, make_plot)
plots

tiff(filename = ("./figures/scatter_plots_of_predictions.tif"), units = "in",
     compression = "lzw", width = 8, height = 12, res = 300)
plot_grid(plotlist = plots, ncol = 2, labels = "auto",
          label_x = 0.15, label_y = 0.9)
dev.off()
