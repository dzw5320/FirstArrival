# Spatial interpolation using kriging
library(rgdal)
library(tmap)
library(gstat)


#Run the model with 114 counties with atleast 10 years of data to get 2000 samples (short tun). Use the last iteration 
# for kriging

#Interpolating mu gp
mu_gp<-mh_gp_mu$smp_mu_gp[2000,]#last iteration of mu_gp

df<-data.frame(mu_gp=mu_gp, X=s[,2], Y=s[,1])
coordinates(df) = ~X+Y

f.1 <- as.formula(mu_gp ~ X + Y) 
var.smpl <- variogram(f.1, data=df)

dat.fit  <- fit.variogram(var.smpl,vgm(model="Exp"))

plot(var.smpl, dat.fit)

# Perform the krige interpolation (note the use of the variogram model
# created in the earlier step)

grd=data.frame(X=County$LONGITUDE, Y=County$LATITUDE)
coordinates(grd)=~X+Y

dat.krg <- krige( f.1, locations=df, newdata=grd, model=dat.fit)


library(maptools)  # Used for conversion from SPDF to ppp
library(raster)    # Used to clip out thiessen polygons
library(ggplot2)

#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('red','blue'))

#This adds a column of color values
# based on the y values
Col <- rbPal(10)[as.numeric(cut(dat.krg$var1.pred,breaks = 10))]

plot(dat.krg,col=Col, pch=18, cex=2)


mu_gp=dat.krg$var1.pred




##################################
####lK

lK_fin=mh_gp_mu$smp_lK[2000,,]
lK_krig<-matrix(NA,nrow = 869, ncol = 8)

for(i in c(1:8)){
  
  df<-data.frame(lK=lK_fin[,i], X=s[,2], Y=s[,1])
  coordinates(df) = ~X+Y
  
  f.1 <- as.formula(lK ~ X + Y) 
  var.smpl <- variogram(f.1, data=df)
  
  dat.fit  <- fit.variogram(var.smpl,vgm(model="Exp", psill =mh_gp_mu$smp_par[2000,3] ,range = mh_gp_mu$smp_par[2000,4]))
  # dat.fit  <- fit.variogram(var.smpl,vgm(model="Exp"))
  
  plot(var.smpl, dat.fit)
  
  # Perform the krige interpolation (note the use of the variogram model
  # created in the earlier step)
  
  grd=data.frame(X=County$LONGITUDE, Y=County$LATITUDE)
  coordinates(grd)=~X+Y
  
  dat.krg <- krige( f.1, locations=df, newdata=grd, model=dat.fit)
  Col <- rbPal(10)[as.numeric(cut(dat.krg$var1.var,breaks = 10))]
  
  plot(dat.krg,col=Col, pch=18, cex=2)
  
  lK_krig[,i]<-dat.krg$var1.pred
  
  
}

save(lK_krig, mu_gp, file="C:\\Users\\dhanu\\OneDrive\\FirstArrival\\SpatialMaxID\\Results11\\krig-data.Rdata")


