#--------------------------------------------------------------------------------------------
#
#--------------------------------------------------------------------------------------------
#devtools::install_github("PalEON-Project/stepps-cal")
library(stepps)
library(dplyr)
library(DT)
library(neotoma)
library(fields)
library(sp)
library(rgdal)
library(abind)
library(rstan)
#--------------------------------------------------------------------------------------------
#
#--------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
setwd('~/workflow_stepps_prediction/prediction/')
help.fun.loc <- 'utils/'
data.loc <- 'data/'
plot.loc <- 'plots/'
code.loc <- 'R/'

########################################################################################################################
#read the huge posterior file of the vegetation model
########################################################################################################################
fit <- read_stan_csv('~/workflow_stepps_calibration/vegetation/output/township_120knots_final.csv')
post <- rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
param.names <-colnames(post[,1,]) #find parameter names of posterior
param.greek <- sapply(param.names,function(x) strsplit(x,'[[]')[[1]][1]) #only take greek letter
param.e.bayes <- param.greek%in%c('eta','rho') #only take eta and rho
est.e.bayes <- colMeans(post[,1,param.e.bayes])

eta <- est.e.bayes[grep('eta',names(est.e.bayes))]
rho <- est.e.bayes[grep('rho',names(est.e.bayes))]
########################################################################################################################
# load knot coords 
########################################################################################################################
load('~/workflow_stepps_calibration/vegetation/data_township/township_data_13_taxa_6796_cells_120_knots.rdata')
######################################################################################################################
# make low resolution coordinates
######################################################################################################################
coord.east.unique <- sort(unique(veg_coords[,1]))
coord.north.unique <- sort(unique(veg_coords[,2]))

coord.east.agg <- seq(coord.east.unique[2],coord.east.unique[length(coord.east.unique)],24000)
coord.north.agg <- seq(coord.north.unique[2],coord.north.unique[length(coord.north.unique)],24000)

coord.agg <- expand.grid(coord.east.agg,coord.north.agg)

veg_coords <- data.frame(veg_coords)
colnames(coord.agg) <- c('east','north')
coord.agg <- data.frame(coord.agg)

dist.coords <- rdist(coord.agg,veg_coords)
min.dist.coords <- apply(dist.coords,1,min)
coord.agg.final <- coord.agg[min.dist.coords==0,] 
########################################################################################################################
veg_box <- bbox_tran(coord.agg.final, '~ east + north',
                       '+init=epsg:3175', 
                       '+init=epsg:3175')
  
reconst_grid <- build_grid(veg_box, resolution = 24000, proj = '+init=epsg:3175')


##########################################################################################################################
#load pollen data
##########################################################################################################################
y <- readRDS('~/workflow_stepps_prediction/pollen_data/pollen_1700_1800.RDS')
y <- round(y)
pollen_coords <- readRDS('~/workflow_stepps_prediction/pollen_data/coordinates_1700_1800.RDS') 
pollen_coords <- matrix(ncol=2,as.numeric(pollen_coords))

###################################################################################################################################
#
###################################################################################################################################
target_taxa <- sort(c('Maple','Birch','Hickory','Chestnut','Other conifer','Other hardwood',
                      'Beech','Ash','Tamarack','Spruce','Pine','Oak','Hemlock'))

calib_trans <- cbind(pollen_coords,y)
calib_trans <- as.data.frame(calib_trans)
colnames(calib_trans) <- c('lon','lat',target_taxa)

veg_trans <- cbind(coord.agg.final,matrix(nrow=754,ncol=13,runif(754*13)))
veg_trans <- as.data.frame(veg_trans)
colnames(veg_trans)[3:ncol(veg_trans)] <- target_taxa


veg_table <- to_stepps_shape(veg_trans,   '~ east + north',      '+init=epsg:3175')
pol_table <- to_stepps_shape(calib_trans, '~ lon + lat', '+init=epsg:4326')


#change this there is a better location for this file!!
source(paste(code.loc,'prep_input_modified.R',sep=''))
neus_agg <- prep_input(veg = veg_table, 
                       pollen = pol_table, 
                       target_taxa = target_taxa,
                       grid   = reconst_grid,hood = 7e+05)

neus_agg$d_pot[,1] <- neus_agg$d_pot[,1]/10^6  
neus_agg$d <- neus_agg$d/10^6
d <- neus_agg$d
idx_cores <- neus_agg$idx_cores
d_pot <- neus_agg$d_pot
N_cores <- neus_agg$N_cores
N <- neus_agg$N_cells
N_cells <- N

T = 1
res <- 3 #correct

  
  
 
  #--------------------------------------------------------------------------------------------
  # have to load parameters from calibration model
  #--------------------------------------------------------------------------------------------
source(paste(help.fun.loc,'pred_helper_funs.r',sep='')) #make sure I use medians of a
source(paste(help.fun.loc,'process_funs.r',sep='')) #make sure I use medians of a
source(paste(code.loc,'build_cal_main.r',sep=''))

runs1 <- list(runs[[4]])
for (run in runs1) {
  kernel <- run$kernel
  num_a <- run$one_a
  one_psi <- run$one_psi
  handle <- strsplit(run$suff_fit,'_A')[[1]][1]
  
  #look at that again...
  for (i in c(0:2,6,7)) { # would have to change this if not run_pl
    if(kernel=='pl'){
      if(num_a==FALSE){
        fname = paste('~/stepps_13_taxa/output/output_cal_pl_Ka_Kgamma_EPs_',i,'.csv',sep='') #sprintf('%s/%s.csv', path_out, suff_fit)
      }
      if((num_a==TRUE)&(i<5)){
        fname = paste('~/stepps_13_taxa//output/output_cal_pl_',i,'.csv',sep='') #sprintf('%s/%s.csv', path_out, suff_fit)
      }
      # if(run$a_const==TRUE){
      #   fname = paste('~/stepps_13_taxa//output/output_cal_pl_Ka_Kgamma_EPs_',i,'.csv',sep='') #sprintf('%s/%s.csv', path_out, suff_fit)
      # }  
    }
    if(kernel=='gaussian'){
      if(one_psi==FALSE){
        fname = paste('~/stepps_13_taxa/output/output_cal_g_Kpsi_Kgamma_EPs_',i,'.csv',sep='') #sprintf('%s/%s.csv', path_out, suff_fit)
      }
      if(one_psi==TRUE){
        fname = paste('~/stepps_13_taxa/output/output_cal_g_',i,'.csv',sep='') #sprintf('%s/%s.csv', path_out, suff_fit)
      }
      }
    fit <- read_stan_csv(fname)
    if (i==0)  {post = rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)}
    else{post_new <-  rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
    post <- abind(post,post_new,along=1)
    }
  }
}    


param.names <-colnames(post[,1,]) #find parameter names of posterior
param.greek <- sapply(param.names,function(x) strsplit(x,'[[]')[[1]][1]) #only take greek letter
param.e.bayes <- param.greek%in%c('gamma','phi') #only take eta and rho
est.e.bayes <- colMeans(post[,1,param.e.bayes])
gamma.all <- est.e.bayes[grep('gamma',names(est.e.bayes))]
gamma <- gamma.all
phi.all <- est.e.bayes[grep('phi',names(est.e.bayes))]
phi <- phi.all

#this is certainly not correct yet
w <- build_weight_matrix(post = post,d = t(neus_agg$d),idx_cores = neus_agg$idx_cores,
                           N = N,N_cores =  N_cores,run = run) #check N
  
  
  
sum_w_pot <- build_sumw_pot(post = post, K=K,N_pot = neus_agg$N_pot,
                              d_pot =  neus_agg$d_pot,run = run)
  
weights_pot <- build_w_pot(post = post, K=K,N_pot = neus_agg$N_pot,d_pot =  neus_agg$d_pot,run = run)
#cbind(weights_pot[,2],neus_agg$d_pot[,2])
gamma.added <- (weights_pot[1,] + weights_pot[2,]/2)/sum_w_pot
gamma.new <- gamma + gamma.added
gamma <- gamma.new
#--------------------------------------------------------------------------------------------
lag <- 0  
#################################################################################################################
# knot distance matrix
#################################################################################################################
d_inter <- rdist(coord.agg.final,knot_coords)/10^6
d_knots <- rdist(knot_coords,knot_coords)/10^6  
  
  
  
  
#-------------------------------------------------------------------------------------------
#d <- rdist(coord.agg.final,coord.agg.final) # check this again
stan_rdump(
  c('K', 'N','T', 'N_knots', 'N_cores','lag', 
    'y', 'res',
    'sum_w_pot',
    'rho','eta','gamma','phi',
    'idx_cores',
    'd','d_knots','d_inter','w'
  ), 
  file=paste('~/workflow_stepps_prediction/prediction/data/prediction_test1.dump',sep=""))


save(K, N,T, N_knots, N_cores,lag, 
     y, res,
     sum_w_pot,
     rho,eta,gamma,phi,
     idx_cores,
     d,d_knots,d_inter,w,
     coord.agg.final, 
     file=paste('~/workflow_stepps_prediction/prediction/data/prediction_test1.rdata',sep=""))
