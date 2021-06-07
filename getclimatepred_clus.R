# Run get_climate_data.R
# Run fittingmodel_climate.R
# Run getsamp_B.R
# Run getclimatepred.R





load("climate_run_data.Rdata")#from get_climate_data.R
load("samp_data.Rdata")#from fittingmodel_climate.R
load("samp_lB_climate.Rdata")#from getsamp_B.R
load("marg_par_climate.Rdata")#from getclimatepred.R
sourceCpp("stablemix//src//exp_stable.cpp")# from stablemix sorce files
library(stablemix)

s <- as.matrix(County[,c(5,4)])
colnames(s)<-c("", "")


alpha_samps<-smp_alpha[c(900:2218)]
theta_samps<-smp_theta[c(900:2218)]
lK_samps<-smp_lK[c(900:2218),,]
lB_samps<-samp_lB_climate[c(900:2218),]
lB_samps<-log(lB_samps)
y_array<-rep(NA, 1319*50*869)
dim(y_array)<-c(1319,50, 869)
num_mcmc<-length(alpha_samps)


for(i in c(1:num_mcmc)){
  cat(i, "\n")
  lB <- matrix(lB_samps[i,], nrow = 8, ncol = 50)
  lzp <- rpostpred(50, s, alpha_samps[i], lK_samps[i,,], lB)
  y_array[i,,] <- lsm2gev(lzp, alpha_samps[i], theta_samps[i], lK_samps[i,,], t(loc_array[i,,]), t(exp(scale_array[i,,])), t(shape_array[i,,]))
  
  }

save(y_array, file="pred_climate.Rdata")

y_mean<-rep(NA, 50*869)
dim(y_mean)<-c(50, 869)
for(y in c(1:50)){
  for(c in c(1:869)){
    y_mean[y,c]<-mean(y_array[,y,c])
  }
}
y_mean_t<-t(y_mean)


y_array_min<-(y_array-102)*(-1)

for(i in c(1:1319)){
  for(j in c(1:50)){
  y_diff_array[i,j,]<-y_array_min[i,j,]-y_array_min_act[i,16,]
  }
}


summary(as.vector(y_sd))


y_sd<-rep(NA, 50*869)
dim(y_sd)<-c(50, 869)
for(y in c(1:50)){
  for(c in c(1:869)){
    y_sd[y,c]<-sd(y_diff_array[,y,c])
  }
}



y_sd_t<-t(y_sd)

County<-cbind(County, y_sd_t)
names_pred<-unlist(lapply(c(2101:2150), function(x){paste("predsd", x, sep="")}))
names(County)[110:159]<-names_pred

save(County, file="pred_climate_county.Rdata")

