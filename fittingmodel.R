#Running the final model-8 basis functions, sigma(s) as a fixed linear 

load("magwarnewcovs2.Rdata")
library(stablemix)

source("fun1.R")
load("krig-data.Rdata")#Refer kriging.R on how to generate this

#Function to run mcmc
example_sim <- function(nmcmc, simset, s, y, df) {
  
  s <- s
  lK <- lK_krig
  lB <- lB
  
  
  
  
  # Run MCMC ----------------------------------------------------------------
  init <-
    list(
      alpha = simset$alpha,
      theta = simset$theta,
      gvar = simset$gvar,
      gscl = simset$gscl,
      #loc = mar_betas[["loc"]],
      loc=matrix(c(32.05185, -2.902534, -0.74165440,-1.91728676,-0.54775291,0.03412806,0.22638220,0.32990968,0.93402868 ), nrow=9, ncol=1),
      scale = matrix(c(1.81198730,0.05151185,0.04067684,0.21409325,-0.01155387,0.02463049,0.06567208,0.00306442), nrow=8, ncol=1),
      shape = -0.05
    )
  tune_var <- list(
    alpha = 0.01,
    theta = 0.01,
    gvar = c(0.1),
    gscl = c(0.1),
    loc = c(0.01, 0.01, 0.01, 0.01, 0.01, 0.01,0.01,0.01, 0.01),
    scale = c(0.01, 0.01, 0.01, 0.01,0.01,0.01,0.01, 0.01),
    shape = c(0.01),
    lB = 0.01,
    lK = 0.01,
    gps=list(gp=0.01, gvar=0.01, gscl=0.01)
  )
  run_time <- system.time(
    mh <- mcmc_lnorm_basis(
      nmcmc,
      s,
      y,
      lB,
      lK,
      init,
      tune_var,
      df,
      loc = ~ lat+lon+ elev + forestcov + waterprop+popdata+tempdata+NAO,
      scale = ~ lat+lon+ elev + forestcov + waterprop+tempdata+NAO,
      shape = ~ 1,
      mar_gp_which=c("mu"),
      mar_gp_par=list(mu = list(gvar = 2.92955e-22, gscl = 2.641677)),
      mar_gp_init=list(mu=as.matrix(mu_gp, nrow=1, ncol=869), sigma=matrix(0, nrow=1, ncol=869 ), xi=matrix(0,  nrow=1, ncol=869)),
      thin_int = 50,
      parallelize = TRUE
    )
  )
  mh$run_time <- run_time
  mh$s <- s
  return(mh)
}


#Latitude/Longitude coordinates
s <- as.matrix(County[,c(5,4)])
colnames(s)<-c("", "")

#Transforming arrival time into maxima
early_arr_max<-(early_arr*(-1))+102

#initializing some of the model parameters.
simset <- list(
  alpha = 0.4879494,               # alpha-stable index (this has to be)
  theta = 0.002396685,               # exponential tilting parameter
  gvar = 86.21692 ,                 # GP Var
  gscl = 0.5764775 ,                # GP Scale
  n = 16,                    # number of replicates (years)
  nloc = ncol(early_arr_max),                # number of spatial locations
  L = 8                    # number of random effects per year
)







#Fit all together to raw data

#get values for tempdata

tempdata<-County[,c(7:22)]
tempdata<-as.matrix(tempdata)

tempdata_vec<-tempdata[,1]

for(i in c(2:16)){
  
  tempdata_vec<-c(tempdata_vec, tempdata[,i])
  
}




#NAO values for March of each year
NAO_March<-c(1.5,   
             -3.0,   
             -1.8,   
             3.1,  
             1.6,  
             1.7,   
             -1.5,  
             0.4,   
             0.9,   
             -4.3,  
             2.2,  
             3.1,  
             1.4,  
             1.5,  
             -1.0,   
             2.6  )

#dataframe with all covariates
df <- data.frame( lat = rep(s[,1], 16), lon=rep(s[,2], 16), elev=rep(County$ELEVATION, 16),forestcov=rep(County$FORESTCOVER, 16),waterprop=rep(County$WATERPROP, 16), popdata=rep(County$POPDATA, 16), tempdata=tempdata_vec,
                  NAO=rep(NAO_March,each=869))

#center and scale covariates covariates
df_scale<-scale(df, center=TRUE, scale=TRUE)
df_scale<-as.data.frame(df_scale)

mh_gp_mu<-example_sim(300000, simset, s, early_arr_max, df_scale)

save(mh_gp_mu, file = 'sigma_nogp_8_all_300K.Rdata')