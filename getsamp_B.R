sourceCpp("~//SpatialMaxID//stablemix//src//exp_stable.cpp")
load("climate_run_data.Rdata")
library(stablemix)


load("samp_data.Rdata")#from fittingmodel_climate.R


samp_lB_climate<-matrix(NA, nrow = 2218, ncol=8*50)

for(i in c(1:nrow(samp_lB_climate))){
  cat(i, "\n")
  samp_lB_climate[i,]<-rhexpstab(8*50, smp_alpha[i], smp_alpha[i],smp_theta[i])
  
  
}

save(samp_lB_climate,file="~/work/SpatialMaxID/Results13/samp_lB_climate.Rdata" )



