library(raster)
library(sf)
library(mapview)

CBI<-raster("./raw_data/burn_severity/ravg_2022_cbi4.tif")
burn_perimiter<-read_sf("./processed_data/burn_perimiter.shp")


set.seed(207)

CBI_sample<-st_buffer(st_sample(burn_perimiter,1),1500)

st_area(CBI_sample)/4046

CBI_cropped<-clipped_burn_raster<-raster::crop(CBI,st_as_sf(CBI_sample))
CBI_cropped[0]<-NA
CBI_cropped[9]<-NA

plot(CBI_cropped)
CBI_cropped_df<-as.data.frame(rasterToPoints(CBI_cropped))


num_points<-1500


  
  
  sampled_points<-CBI_cropped_df[sample(nrow(CBI_cropped_df), num_points,replace=TRUE), ]
  sampled_points_sf<-st_as_sf(sampled_points,coords=c("x","y"),crs=st_crs(CBI))
  plot(sampled_points_sf)
  
  sampled_points_sf$size<-runif(num_points,2,200)
  
  sampled_points_sf$radius<-sqrt((sampled_points_sf$size*4046.86)/pi)
  treatments<-st_buffer(sampled_points_sf,sampled_points_sf$radius)
  

  treatments$total_severity<-raster::extract(CBI,treatments,fun=sum)
  treatments$avg_severity<-treatments$total_severity/as.numeric(st_area(treatments))

  
treatments<-treatments[order(treatments$avg_severity, decreasing = TRUE), ]

optimal_set<-treatments[1,]

total_treatment_size<-as.numeric(st_area(CBI_sample)/4046*0.2)
total_treatment_size

for (i in 2:nrow(treatments)){
  treatment_to_consider<-treatments[i,]
  intersection<-st_intersection(optimal_set,treatments[i,])
  area_intersection<-as.numeric(st_area(intersection))

  if (length(intersection$size)==0 & (sum(optimal_set$size)+treatment_to_consider$size)<total_treatment_size){
    optimal_set<-rbind(optimal_set,treatments[i,])
    
  }
  
  
}


mapview(CBI_cropped)+mapview(optimal_set,color="red")

mapview(CBI_cropped)+mapview(treatments,color="red")

two_acre_min<-optimal_set
fourty_acre_min<-optimal_set

### no treatments
CBI_df<-as.data.frame(rasterToPoints(CBI_cropped))
mean(CBI_df$ravg_2022_cbi4)
pre<-nrow(CBI_df[CBI_df$ravg_2022_cbi4==4,])

two_acre_mask<-raster::mask(CBI_cropped,st_as_sf(two_acre_min),inverse=TRUE)
plot(two_acre_mask)
two_acre_mask[is.na(two_acre_mask[])] <- 1 
two_acre_df<-as.data.frame(rasterToPoints(two_acre_mask))
mean(two_acre_df$ravg_2022_cbi4)
nrow(two_acre_df[two_acre_df$ravg_2022_cbi4==4,])/pre

fourty_acre_mask<-raster::mask(CBI_cropped,st_as_sf(fourty_acre_min),inverse=TRUE)
plot(fourty_acre_mask)
fourty_acre_mask[is.na(fourty_acre_mask[])] <- 1 
fourty_acre_df<-as.data.frame(rasterToPoints(fourty_acre_mask))
mean(fourty_acre_df$ravg_2022_cbi4)
nrow(fourty_acre_df[fourty_acre_df$ravg_2022_cbi4==4,])/pre
