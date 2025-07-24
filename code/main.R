# this script is a simple wrapper used to source the other scripts in sequence.
# environment is cleared between scripts to avoid varname conflicts.

# note: be sure to set the "testing" flag to the desired T/F value (wherever it
# occurs) before running any script. setting testing==FALSE will fully replicate
# the results with original sample sizes and spatial neighborhood sizes, but
# will take a very long time (hours to days) to run.

source("./data_wrangling.R", echo = TRUE)
rm(list = ls())
gc()

source("./day_of_burn.R", echo = TRUE)
rm(list = ls())
gc()

source("./day_of_burn.R", echo = TRUE)
rm(list = ls())
gc()

source("./weather_and_knn.R", echo = TRUE)
rm(list = ls())
gc()

source("./comparing_methods_rf.R", echo = TRUE)
rm(list = ls())
gc()

source("./comparing_methods.R", echo = TRUE)
rm(list = ls())
gc()

source("./rmse_and_plots.R", echo = TRUE)
rm(list = ls())
gc()

source("./effects_across_hpcc.R", echo = TRUE)
rm(list = ls())
gc()

source("./plotting_hpcc_results.R", echo = TRUE)
rm(list = ls())
gc()

source("./cv_different_sizes_plots.R", echo = TRUE)
rm(list = ls())
gc()
