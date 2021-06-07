load("climate_run_data.Rdata")
load("samp_data.Rdata")
load("samp_lB_climate.Rdata")


tempdata<-County[,c(10:59)]
tempdata<-as.matrix(tempdata)

tempdata_vec<-tempdata[,1]

for(i in c(2:50)){
  
  tempdata_vec<-c(tempdata_vec, tempdata[,i])
  
}

s <- as.matrix(County[,c(5,4)])
colnames(s)<-c("", "")

df <- data.frame( lat = rep(s[,1], 50), lon=rep(s[,2], 50), elev=rep(County$ELEVATION, 50),forestcov=rep(County$FORESTCOVER, 50),waterprop=rep(County$WATERPROP, 50), popdata=rep(County$POPDATA, 50), tempdata=tempdata_vec ,
                  NAO=rep(NAO_climate,each=869))

meantempdata<-15.45167
sdtempdata<-211.3122
meanNAO<-0.525
sdNAO<-2.140024  

df_scale<-scale(df[,c(1:6)], center=TRUE, scale=TRUE)
df_scale<-as.data.frame(df_scale)
df_scale$tempdata<-(df$tempdata-meantempdata)/sdtempdata
df_scale$NAO<-(df$NAO-meanNAO)/sdNAO

#Adjusting betas for mu gp
X=cbind(matrix(1, nrow=nrow(df_scale), ncol=1),as.matrix(df_scale))
X=X[c(1:869), c(1:7)]

mu_samps_adj<-matrix(NA, nrow = 1319, ncol = 7)
for(i in c(1:nrow(mu_samps_adj))){
  cat(i, "\n")
  mu_samps_adj[i,]<-smp_locpar[900+i-1,1:7]-solve(t(X)%*%X)%*%t(X)%*%smp_mu_gp[(900+i-1), ]
}




#Compute loc, scale and shape matrices

num_mcmc<-nrow(mu_samps_adj)
mu_samps<-cbind(mu_samps_adj, smp_locpar[900:2218,c(8,9)])
scale_samps<-smp_scalepar[900:2218,]
shape_samps<-smp_shapepar[900:2218]
mu_gp<-smp_mu_gp[900:2218,]

loc_array<-rep(NA, num_mcmc*869*50)
dim(loc_array)<-c(num_mcmc,869,50)

scale_array<-rep(NA, num_mcmc*869*50)
dim(scale_array)<-c(num_mcmc,869,50)

shape_array<-rep(NA, num_mcmc*869*50)
dim(shape_array)<-c(num_mcmc,869,50)



X_all=cbind(matrix(1, nrow=869*50, ncol=1),as.matrix(df_scale))




for(r in c(1:869)){#each county
  cat(r, "\n")
  for(y in c(1:50)){#for each year
    year_range=c((869*(y-1)+1):(869*y))
    X_all_yr=X_all[year_range, ]
  for(c in c(1:num_mcmc)){# each mcmc sample after burnout
    loc_array[c,r,y]<-sum(mu_samps[c,]*X_all_yr[r,])+mu_gp[c,r]
    scale_array[c,r,y]<-sum(scale_samps[c,]*X_all_yr[r,c(1:6,8,9)])
    shape_array[c,r,y]<-shape_samps[c]
  }
    
  }  
}

save(loc_array, scale_array, shape_array, file="marg_par_climate.Rdata")




