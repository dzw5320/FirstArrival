library(ncdf4)
ncin=nc_open("psl_Amon_GISS-E2-H_rcp85_r1i1p1_215101-220012.nc")# Get CMIP5 data. Not included in the repository due to volume. 
#Download it here- https://dataserver.nccs.nasa.gov/thredds/fileServer/CMIP5/NASA/GISS/rcp85/E2-H_rcp85_r1i1p1/psl_Amon_GISS-E2-H_rcp85_r1i1p1_205101-210012.nc
print(ncin)
#Rotate from 0 to 360 to -180 to 180
library(raster)

a <- brick("psl_Amon_GISS-E2-H_rcp85_r1i1p1_215101-220012.nc", varname="psl")
rotated_psl <- rotate(a)


#For the 50 years under consideration calculate the March SLP for 
#Reykjavik (lat 64.135666, lon -21.862675) and ponta delgada (lat 37.73333, lon -25.66667)
station_latlon<-data.frame( LONGITUDE=c(-22.73391,-25.7104 ),  LATITUDE=c(65.074,37.7493 ))
# station_latlon<-data.frame( LATITUDE=c(64.135666,37.73333 ), LONGITUDE=c(-21.862675,-25.66667 ))
coordinates(station_latlon)=~LONGITUDE +LATITUDE


Rey_slp<-data.frame(Year=seq(2151:2200), SLP=extract(rotated_psl, station_latlon[1,])[,seq(3, by = 12, to = 600)])
PD_slp<-data.frame(Year=seq(2151:2200), SLP=extract(rotated_psl, station_latlon[2,])[,seq(3, by = 12, to = 600)])
save(Rey_slp,PD_slp, file="climate_Rey_PD.Rdata" )
#######################
# Get mean and sd of Reykjavik and Ponta Delgada stations

Rey<-read.csv("SLP_Rey_March.csv")
Rey$SLP.mb.<-Rey$SLP.mb.*100
names(Rey)[4]<-c("SLP_Pa")
mean(Rey$SLP_Pa[which(Rey$Date>=1864 & Rey$Date <=1983)])
# 100611.3 mean
# 771.3932 sd


PD<-read.csv("SLP_PD_March.csv")
PD$SLP.mb.<-PD$SLP.mb.*100
names(PD)[4]<-c("SLP_Pa")

#Need to find PD data for years 1864-1894 will use the following gridded data set.

#PD_hist<-brick("slp.mean.nc", varname="slp")
#PD_hist <- rotate(PD_hist)
#hist_extr<-extract(PD_hist, station_latlon[2,])[seq(3, by = 12, to = 2643)]
load("PD_hist_extract.Rdata")
Pd_hist_df<-data.frame(Year=seq(1800,to = 2020, by = 1), SLP=hist_extr)
Pd_hist_df$SLP<-Pd_hist_df$SLP*100
PD_hist_df_early<-Pd_hist_df$SLP[Pd_hist_df$Year>=1864 & Pd_hist_df$Year<1894]
mean(c(PD_hist_df_early, PD$SLP_Pa[PD$Date>=1894 & PD$Date <=1983]),na.rm=TRUE)

sd(c(PD_hist_df_early, PD$SLP_Pa[PD$Date>=1894 & PD$Date <=1983]),na.rm=TRUE)
# 101860.5 mean
# 592.8365 sd

Rey_slp_norm<-(Rey_slp$SLP-100611.3)/771.3932
PD_slp_norm<-(PD_slp$SLP-101860.5)/592.8365

NAO_climate<-PD_slp_norm-Rey_slp_norm

save(Rey_slp, PD_slp, Rey_slp_norm, PD_slp_norm, NAO_climate, file="climate_NAO.Rdata")

##############
# tas 
# rotate .nc file
library(weathermetrics)
a <- brick("C:\\Users\\dhanu\\OneDrive\\FirstArrival\\SpatialMaxID\\Climate Data\\tas_Amon_GISS-E2-H_rcp85_r1i1p1_215101-220012.nc", varname="tas")
rotated_tas <- rotate(a)

coordinates(us_stationswithmean)=~Longitude +Latitude

station_ann_anom<-kelvin.to.celsius(extract(rotated_tas, us_stationswithmean)[,seq(3, by = 12, to = 600)], round = 2)*100-us_stationswithmean$mean

grid_anom<-matrix(NA, nrow = 2592, ncol=50)

for(i in c(1:nrow(grid_anom))){
  cat(i, "\n")
  if(length(which(us_stationswithmean$grid==i))>0){
    
  sub_1=station_ann_anom[which(us_stationswithmean$grid==i),]
  
  if(is.null(dim(sub_1))){
    
    grid_anom[i,]<-sub_1
  }else{
    
  grid_anom[i,]<-colMeans(sub_1)
  }
  
  }
}

save(grid_anom, file="grid_anom.Rdata")


# create function that gives the grid cell number given lat and loncreate

get_tempanorm<-function(data, lat, lon){
  
  lat_idx<-((which(!lat<=lat_seq)[1])-1)
  lon_idx<-((which(!lon>=lon_seq)[1])-1)
  # cat(lat_idx, lon_idx, "\n")
  return(data[lat_idx,(lon_idx)])
}


County<-County[,c(1:6, 23:25)]

County_temp_anom<-matrix(NA, nrow=nrow(County), ncol=50)


for(r in c(1:nrow(County))){
  cat(r, "\n")
  County_temp_anom[r, ]<-grid_anom[get_gridnum(grid, County$LATITUDE[r], County$LONGITUDE[r]), ]
  
  
}

County<-cbind(County, County_temp_anom)

anom_col_names<-unlist(lapply(c(1:50), function(x){paste("temp", 2150+x, sep="")}))

names(County)[10:59]<-anom_col_names



plot(NAO_climate, type="l")

save(County, NAO_climate, file="climate_run_data.Rdata")


