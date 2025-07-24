# Source code for:

### Comparing methods for estimating counterfactual burn severity and evaluating treatment effectiveness in the Hermits Peak Calf Canyon Fire, NM, USA

#### Nathan Tomczyk et al. 2025

[ntomczyk\@nmhu.edu](mailto:ntomczyk@nmhu.edu){.email}

This repository contains all R code necessary to reproduce the published results. Scripts should be run in the following order:

1.  `data_wranging.R`
2.  `day_of_burn.R`
3.  `weather_and_knn.R`
4.  `comparing_methods_rf.R`
5.  `comparing_methods.R`
6.  `rmse_and_plots.R`
7.  `effects_across_hpcc.R` \# runs PSM and local spatial rf (noweather) to generate counterfactual maps.
8.  `plotting_hpcc_results.R` \# figures to plot results from 6.
9.  `cv_different_sizes_plots.R` \# this is basically the same as #5 but uses different plot sizes.
10. `buffer_and_clip_size.R` (optional)

Also provided is a `main.R` file which is a simple utility to source the files in the correct sequence.

## Testing flags

Some scripts are computationally intensive and can take hours or days to run. Most of the long computation times are the result of large sample sizes and operations that scale exponentially with spatial area. To allow faster (though still slow) unit testing, some scripts contain a `testing` variable that can be set to `TRUE` or `FALSE`. Setting `testing <- TRUE` will implement subsampling and reduce the spatial extent used in key operations, such as fitting local-spatial random forest models. After you've confirmed that the code works on your machine, set `testing <- FALSE` in all scripts and re-run the analysis to ensure the correct sample sizes and spatial extents are used.

## Overview of source files

### main.R

A utility to source all other scripts in sequence.

### day_of_burn.R

this script reads in or downloads raw data, clips, buffers, resamples, and reprojects as needed, then writes out tidy .shp or .tif files.

### weather_and_knn.R

This script prepares weather data and identifies candidate control plots for the burn severity analysis. It (optionally) downloads gridmet data, processes it, and extracts relevant features for each day of the burn. It also creates a grid of control plots based on the burn perimeter and extracts environmental features for these plots. This script should take less than an hour to run on most laptops - \~32 minutes on a ca. 2022 laptop with 32 GB ram.

### comparing_methods_rf.R

This script is used to generate predicted counterfactual values for the random forest methods. These methods are computationally expensive, and this script takes a long time to run. Set `testing <- TRUE` to ensure the script is running error-free before attempting to replicate results with the full sample sizes and spatial extents.

### comparing_methods.R

This script generates predicted counterfactual values for non-random forest methods. These predictions are used to compare the accuracy of each method and are a core part of the analysis.

### rmse_and_plots.R

This script calculates the bias and root means square error of counterfactual predictions from different methods. Summary plots and scatter plots of predicted vs. actual values are generated.

### effects_across_hpcc.R

This script is used to estimate counterfactual values of actual treated areas (rather than randomly placed validation plots), using a subset of the methods detailed in the paper.

### plotting_hpcc_results.R

This script is used to generate maps and plots from the output of the previous script.

### cv_different_sizes_plots.R

This tests the effect of varying the validation plot size. Counterfactual predictions are generated and evaluated across a range of seven different validation plot sizes.

### buffer_and_clip_size.R

This script is used to examine edge effects, justifying the choice of 60 m as a "buffer size" to create an exclusion zone around validation plots and treated areas. This script has not been styled or refactored and may be difficult to parse; please reach out for assistance if needed.
