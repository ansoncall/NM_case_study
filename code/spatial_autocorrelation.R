library(gstat)
library(sf)
library(raster)

library(dplyr)
library(spdep)
library(ggplot2)


burn_perimiter<-read_sf("./processed_data/burn_perimiter.shp")

elev_down<-raster("./processed_data/elev_down.tif")
aspect_down<-raster("./processed_data/aspect_down.tif")
TRI_down<-raster("./processed_data/TRI_down.tif")
TPI_down<-raster("./processed_data/TPI_down.tif")
slope_down<-raster("./processed_data/slope_down.tif")
roads_distance<-raster("./processed_data/distance_to_road.tif")
site_potential<-raster("./processed_data/lf_site_potential_new.tif")

CBI<-raster("./processed_data/clipped_burn_raster.tif")
CBI[CBI==9]<-NA
CBI[CBI==0]<-1


radius_m<-sqrt((15000*4046.86)/pi)

set.seed(207)

random_point_in_burn<-st_sample(burn_perimiter,size=1)

random_point_in_burn<-st_buffer(random_point_in_burn,dist=radius_m)

random_point_in_burn<-sf::st_transform(random_point_in_burn,crs=st_crs(elev_down))



elev_down_small<-raster::crop(elev_down,st_as_sf(random_point_in_burn))
aspect_down_small<-raster::crop(aspect_down,st_as_sf(random_point_in_burn))
TRI_down_small<-raster::crop(TRI_down,st_as_sf(random_point_in_burn))
TPI_down_small<-raster::crop(TPI_down,st_as_sf(random_point_in_burn))
slope_down_small<-raster::crop(slope_down,st_as_sf(random_point_in_burn))
roads_small<-raster::crop(roads_distance,st_as_sf(random_point_in_burn))
esp_small<-raster::crop(site_potential,st_as_sf(random_point_in_burn))
cbi_small<-raster::crop(CBI,st_as_sf(random_point_in_burn))

#mapview(elev_down_small)

elev_down_points<-rasterToPoints(elev_down_small)
aspect_down_points<-rasterToPoints(aspect_down_small)
TRI_down_points<-rasterToPoints(TRI_down_small)
TPI_down_points<-rasterToPoints(TPI_down_small)
roads_down_points<-rasterToPoints(roads_small)
esp_down_points<-rasterToPoints(esp_small)
cbi_down_points<-rasterToPoints(cbi_small)


vario.elev<-variogram(elev_down~1,data=as.data.frame(elev_down_points),locations=~x+y)
vario.aspect<-variogram(aspect_down~1,data=as.data.frame(aspect_down_points),locations=~x+y)
vario.tpi<-variogram(TPI_down~1,data=as.data.frame(TPI_down_points),locations=~x+y)
vario.tri<-variogram(TRI_down~1,data=as.data.frame(TRI_down_points),locations=~x+y)
vario.roads<-variogram(distance_to_road~1,data=as.data.frame(roads_down_points),locations=~x+y)
vario.esp<-variogram(lf_site_potential_new~1,data=as.data.frame(esp_down_points),locations=~x+y)
vario.cbi<-variogram(clipped_burn_raster~1,data=as.data.frame(cbi_down_points),locations=~x+y)

vario.elev$var<-"elevation"
vario.aspect$var<-"aspect"
vario.tpi$var<-"TPI"
vario.tri$var<-"TRI"
vario.roads$var<-"Roads"
vario.esp$var<-"ESP"
vario.cbi$var<-"CBI"

vario.elev$gamma<-vario.elev$gamma/max(vario.elev$gamma)
vario.aspect$gamma<-vario.aspect$gamma/max(vario.aspect$gamma)
vario.tpi$gamma<-vario.tpi$gamma/max(vario.tpi$gamma)
vario.tri$gamma<-vario.tri$gamma/max(vario.tri$gamma)
vario.roads$gamma<-vario.roads$gamma/max(vario.roads$gamma)
vario.esp$gamma<-vario.esp$gamma/(max(vario.esp$gamma))
vario.cbi$gamma<-vario.cbi$gamma/max(vario.cbi$gamma)

all.variograms<-bind_rows(vario.elev,vario.aspect,vario.cbi,vario.esp,vario.roads,vario.tpi,vario.tri)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
gamma_distance<-ggplot(all.variograms,aes(x=dist,y=gamma,color=var))+geom_line(linewidth=2)+scale_color_manual(values=cbbPalette)+theme_classic()+theme(text=element_text(size=20))+xlab("Distance (m)")+ylab("Gamma / max(Gamma)")+guides(color=guide_legend(title="Variable"))


ggsave(filename="./figures/autocorrelation.tif",gamma_distance,units="in",width=8,height=8)

Moran(elev_down)
Moran(aspect_down)
Moran(TRI_down)
Moran(TPI_down)
Moran(site_potential)
Moran(roads_distance)
Moran(slope_down)
Moran(CBI)

### All moran i are higher than 0.9