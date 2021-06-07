#Test run on cluster
load("~/SpatialMaxID/Results12/climate_run_data.Rdata")#from get_climate_data.R
library(stablemix)

source("//storage//home//d//dzw5320//SpatialMaxID//Results8//fun1.R")#change
load("~/work/SpatialMaxID/Results11/krig-data.Rdata")#change



example_sim <- function(nmcmc, simset, s, y, df) {
 
  s <- s
  
  # set.seed(100)
  # zall <- rstabmix(
  #   simset$n,
  #   s,
  #   nbasis = simset$L,
  #   alpha = simset$alpha,
  #   delta = simset$alpha,
  #   theta = simset$theta,
  #   gvar = simset$gvar,
  #   gscl = simset$gscl,
  #   return_all = T,
  #   type = "br"
  # )
  # # lZ <- log(lZ)
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


s <- as.matrix(County[,c(5,4)])
colnames(s)<-c("", "")


early_arr_max<-matrix(NA, nrow=50, ncol=869)


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

tempdata<-County[,c(10:59)]
tempdata<-as.matrix(tempdata)

tempdata_vec<-tempdata[,1]

for(i in c(2:50)){
  
  tempdata_vec<-c(tempdata_vec, tempdata[,i])
  
}





df <- data.frame( lat = rep(s[,1], 50), lon=rep(s[,2], 50), elev=rep(County$ELEVATION, 50),forestcov=rep(County$FORESTCOVER, 50),waterprop=rep(County$WATERPROP, 50), popdata=rep(County$POPDATA, 50), tempdata=tempdata_vec ,
                   NAO=rep(NAO_climate,each=869))

df_scale<-scale(df, center=TRUE, scale=TRUE)
df_scale<-as.data.frame(df_scale)

mh_gp_mu<-example_sim(300000, simset, s, early_arr_max, df_scale)

save(mh_gp_mu, file = '~//work//SpatialMaxID//Results13//samp_data.Rdata')

