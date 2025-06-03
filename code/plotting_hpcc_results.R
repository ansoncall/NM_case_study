library(rnaturalearth) # map data
library(rnaturalearthdata)# map data
library(ggspatial) 
library(sf)
library(raster)
library(ggplot2)
library(dplyr)

veg_treatments<-read_sf("./results/hpcc_processed_cbi_new.shp")
CBI<-raster("./processed_data/clipped_burn_raster.tif")
CBI[CBI==9]<-NA
CBI[CBI==0]<-1

## selecting only the fuel treatments

descriptions<-tolower(unique(veg_treatments$Dscrptn))

thinning<-descriptions[grepl("thin",descriptions)]
cutting<-descriptions[grepl("cut",descriptions)]
defens<-descriptions[grepl("defens",descriptions)]
burn<-descriptions[grepl("burn",descriptions)]
fuel<-descriptions[grepl("fuel",descriptions)]
fire<-descriptions[grepl("fire",descriptions)]

treatments_to_use<-unique(c(thinning,cutting,defens,burn,fuel,fire))

burning<-unique(c(burn,fuel))

thin<-unique(c(thinning,cutting))

veg_treatments$thin_burn<-"thin"
veg_treatments[tolower(veg_treatments$Dscrptn)%in%burning,"thin_burn"]<-"burn"

veg_treatments<-veg_treatments[tolower(veg_treatments$Dscrptn) %in% treatments_to_use,]

veg_treatments$actual_severity<-raster::extract(CBI,st_as_sf(veg_treatments),fun=mean,na.rm=TRUE,weights=TRUE,exact=TRUE,normalizeWeights=TRUE,small=TRUE)

veg_treatments$effect<-veg_treatments$actual_severity-veg_treatments$cntrl__
veg_treatments$effect_knn<-veg_treatments$actual_severity-veg_treatments$knn_pred

veg_treatments_complete<-veg_treatments[veg_treatments$Year_Cl!="needs input",]

t.test(veg_treatments_complete$effect)

t.test(veg_treatments_complete$effect_knn)

summary(lm(effect~cntrl__,data=veg_treatments_complete))

summary(lm(effect~Acre_US,data=veg_treatments_complete))

summary(lm(effect~thin_burn,data=veg_treatments_complete))

### What prop of area burned at lower severity?

lower<-veg_treatments_complete[veg_treatments_complete$cntrl__>veg_treatments_complete$actual_severity,]

sum(st_area(lower))/sum(st_area(veg_treatments_complete))

high<-veg_treatments_complete[veg_treatments_complete$actual_severity>3,]

sum(st_area(high))/sum(st_area(veg_treatments_complete))


first_plot_colors <- c("#E69F00", "#56B4E9", "#009E73")
second_plot_colors <- c("#0072B2", "#D55E00", "#CC79A7")

####

effect_hist<-ggplot(veg_treatments_complete,aes(x=effect))+geom_histogram(fill="#0072B2",alpha=0.5)+geom_histogram(data=veg_treatments_complete,aes(x=effect_knn),fill="#D55E00",alpha=0.5)+theme_classic()+ylab("Number of treatments")+xlab(expression("Effect of Treatment ("*Delta*"CBI)"))+theme(text=element_text(size=20))+geom_vline(xintercept = 0,linetype="dashed")
effect_hist

scatter_effect<-ggplot(veg_treatments_complete,aes(x=cntrl__,y=effect))+geom_point(size=2)+theme_classic()+xlab("Untreated burn seveirty (CBI)")+ylab(expression("Effect of Treatment ("*Delta*"CBI)"))+theme(text=element_text(size=20))+geom_hline(yintercept = 0,linetype="dashed")+geom_smooth(method="lm")

library(cowplot)

tiff(filename=("./figures/hpcc_results_simple.tif"),units='in',compression='lzw',width=8,height=12,res=300)
plot_grid(effect_hist,scatter_effect,ncol=1,labels="auto",label_x=0.96,label_y=0.95,label_size = 20,align='v')
dev.off()

