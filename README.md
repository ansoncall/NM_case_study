# NM_case_study

Run the scripts in the following order:

1. Data_wranging_revised
2. day_of_burn
3. weather_and_knn_revised
4a. comparing_methods_rf
4b. comparing_methods
5. rmse_and_plots
6. effects_across_hpcc  # runs PSM and local spatial rf (noweather) to generate counterfactual maps.
7. plotting_hpcc_results # figures to plot results from 5. 
8. CV_different_sizes_plots # this is basically the same as #4 but uses different plot sizes. 
9. buffer_and_clip_size
