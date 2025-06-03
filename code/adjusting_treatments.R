veg_treatments<-read_sf("./results/hpcc_processed_cbi_new.shp")


#############
descriptions<-tolower(unique(veg_treatments$Dscrptn))

thinning<-descriptions[grepl("thin",descriptions)]
cutting<-descriptions[grepl("cut",descriptions)]
defens<-descriptions[grepl("defens",descriptions)]
burn<-descriptions[grepl("burn",descriptions)]
fuel<-descriptions[grepl("fuel",descriptions)]
fire<-descriptions[grepl("fire",descriptions)]


treatments_to_use<-unique(c(thinning,cutting,defens,burn,fuel,fire))


veg_treatments<-veg_treatments[tolower(veg_treatments$Dscrptn) %in% treatments_to_use,]




###################

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

veg_treatments_complete<-veg_treatments[tolower(veg_treatments$Dscrptn) %in% treatments_to_use,]


##############

one_treatment<-veg_treatments_complete[1,]
overlaps<-st_intersection(one_treatment,veg_treatments_complete[-1,])
difference<-st_difference(one_treatment,st_union(veg_treatments_complete[-1,]))


difference$retreat<-"no"
overlaps$retreat="yes"

if (nrow(overlaps)>0){
  for (i in 1:nrow(overlaps)){
    count<-0
    
    indexes<-c()
    
    if (overlaps$thin_burn[i]!=overlaps$thin_burn.1[i]) {overlaps$thin_burn<-"Thin+burn"
    count<-count+1}
    
    if (overlaps$Year_Cl[i]>overlaps$Year_Cl.1[i]){
      overlaps$Year_Cl[i]<-overlaps$Year_Cl.1[i]
      count<-count+1
      indexes<-c(indexes,i)
    }
    
    
  }
  if (count>0){bind_rows(difference,overlaps[i,])
    
    
    difference<-bind_rows(difference,overlaps[indexes,1:30])
  }
}

adjusted_veg_treatments<-difference








#for (j in 2:6){
for (j in 2:nrow(veg_treatments_complete)){
  
  one_treatment<-veg_treatments_complete[j,]
  overlaps<-st_intersection(one_treatment,veg_treatments_complete[-j,])
  difference<-st_difference(one_treatment,st_union(veg_treatments_complete[-j,]))
  
  
  difference$retreat<-"no"
  overlaps$retreat="yes"
  
  if (nrow(overlaps)>0){
    for (i in 1:nrow(overlaps)){
      count<-0
      
      indexes<-c()
      
      if (overlaps$thin_burn[i]!=overlaps$thin_burn.1[i]) {overlaps$thin_burn<-"Thin+burn"
      count<-count+1}
      
      if (overlaps$Year_Cl[i]=>overlaps$Year_Cl.1[i]){
        count<-count+1
        indexes<-c(indexes,i)
        }

     
    }
    if (count>0){bind_rows(difference,overlaps[i,])
    
    
    difference<-bind_rows(difference,overlaps[indexes,1:30])
    }
  }
  
  adjusted_veg_treatments<-bind_rows(adjusted_veg_treatments,difference)
  
  
}

