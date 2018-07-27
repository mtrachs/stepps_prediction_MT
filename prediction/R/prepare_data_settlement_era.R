#--------------------------------------------------------------------------------------------
#
#--------------------------------------------------------------------------------------------
#devtools::install_github("PalEON-Project/stepps-cal")
library(stepps)
library(dplyr)
library(DT)
library(neotoma)
library(sp)
library(fields)
library(rgdal)
library(abind)

#-----------------------------------------------------------------------------------------------
setwd('~/workflow_stepps_calibration/prediction/')
help.fun.loc <- 'utils/'
data.loc <- 'data/'
plot.loc <- 'plots/'
#-----------------------------------------------------------------------------------------------


#load pollen data
stepps_input <- load('~/workflow_stepps_prediction/prediction/data/elicitation_neus_certainty_median.RData')
y_pollen <- y
#pollen_coords

#load data produced to calibrate vegetation model
#load('~/workflow_stepps_calibration/vegetation/data/veg_data_15_taxa_754_cells_260_knots.rdata')#this is not the file we should use
veg_agg <- y
#coord.agg.final
#careful y is vegetation data

#reconst_grid

veg_box <- bbox_tran(coord.agg.final, '~ east + north',
                     '+init=epsg:3175', 
                     '+init=epsg:3175')

reconst_grid <- build_grid(veg_box, resolution = 24000, proj = '+init=epsg:3175')




#---------------------------------------------------------------------------------------------
#
#---------------------------------------------------------------------------------------------
calib_trans <- cbind(pollen_coords,y_pollen)
veg_trans <- cbind(coord.agg.final,veg_agg)

veg_table <- to_stepps_shape(veg_trans,   '~ east + north',      '+init=epsg:3175')
pol_table <- to_stepps_shape(calib_trans, '~ lon + lat', '+init=epsg:4326')

target_taxa <- sort(c('Maple','Birch','Hickory','Chestnut','Other conifer','Other hardwood',
                      'Beech','Ash','Tamarack','Spruce','Pine','Poplar','Oak','Hemlock','Elm'))

source('~/workflow_stepps_calibration/calibration/calibration_helper_funs/prep_input_modified.R')
neus_agg <- prep_input(veg = veg_table, 
                        pollen = pol_table, 
                        target_taxa = target_taxa,
                        grid   = reconst_grid,hood = 7e+05)

neus_agg$d_pot[,1] <- neus_agg$d_pot[,1]/10^6  
neus_agg$d <- neus_agg$d/10^6
T = 1
y <- y_pollen
res <- 24

#--------------------------------------------------------------------------------------------
# have to load parameters from calibration model
#--------------------------------------------------------------------------------------------
library(rstan)
source(paste(help.fun.loc,'pred_helper_funs.r',sep=''))
source('~/workflow_stepps_calibration/calibration/utils/process_funs.r')
# gamma
# a
# b
# phi 
# for settlement era, I do not see if it needs a or b, probably not
# in a first sep I will load phi and gamma from power law kernel Ka_Kgamma
source('R/build_cal_main.r')
run <- run_pl 
kernel <- run$kernel
num_a <- run$one_a

for (i in 0:4) { # would have to change this if not run_pl
  if((kernel=='pl')&(num_a==FALSE)){
  fname = paste('~/stepps_median/output_cal_pl_Ka_Kgamma_EPs_',i,'.csv',sep='') #sprintf('%s/%s.csv', path_out, suff_fit)
  }
  if((kernel=='pl')&(num_a==TRUE)){
  fname = paste('~/workflow_stepps_calibration/results/data/stepps_median/output_cal_pl_',i,'.csv',sep='') #sprintf('%s/%s.csv', path_out, suff_fit)
  }
  if((kernel=='gaussian')&(num_a==TRUE)){
  fname = paste('~/stepps_median/output_cal_g_Kpsi_Kgamma_EPs_',i,'.csv',sep='') #sprintf('%s/%s.csv', path_out, suff_fit)
  }
  fit <- read_stan_csv(fname)
  if (i==0)  {post = rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)}
  else{post_new <-  rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
  post <- abind(post,post_new,along=1)
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
                         N = N,N_cores =  N_cores,run = run)

w <- w/sum(w)

sum_w_pot <- build_sumw_pot(post = post, K=K,N_pot = neus_agg$N_pot,
                            d_pot =  neus_agg$d_pot,run = run)

w1 <- array(dim=c(15,62,754))
for(i in 1:15) {w1[i,,]<-w}
w <- w1
#-------------------------------------------------------------------------------------------
#have to rescale gamma
# simplest is to build the kernel and add the weight of the 8 closest locations...
# 
weights_pot <- build_w_pot(post = post, K=K,N_pot = N_pot,d_pot =  d_pot,run = run)
cbind(weights_pot[,2],d_pot[,2])
gamma.added <- (weights_pot[1,] + weights_pot[2,]/2)/sum_w_pot
gamma.new <- gamma + gamma.added
gamma <- gamma.new
#--------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------
# eta (spatial correlation, correaltion strength)
# rho (spatial correlation decorrealtion length)
fit <- read_stan_csv('~/workflow_stepps_calibration/vegetation/data//stepps_veg_baseline.csv')
post <- rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
param.names <-colnames(post[,1,]) #find parameter names of posterior
param.greek <- sapply(param.names,function(x) strsplit(x,'[[]')[[1]][1]) #only take greek letter
param.e.bayes <- param.greek%in%c('eta','rho') #only take eta and rho
est.e.bayes <- colMeans(post[,1,param.e.bayes])

eta.all <- est.e.bayes[grep('eta',names(est.e.bayes))]
eta <- eta.all#[-length(eta.all)]
rho.all <- est.e.bayes[grep('rho',names(est.e.bayes))]
rho <- rho.all#[-length(rho.all)]



#-------------------------------------------------------------------------------------------
#d <- rdist(coord.agg.final,coord.agg.final) # check this again
d <- neus_agg$d
idx_cores <- neus_agg$idx_cores
#---------------------------------------------------------------------------------------------
#store data in .dump file
stan_rdump(
  c('K', 'N','T', 'N_knots', 'N_cores', 
    'y', 'res',
    'sum_w_pot',
    'rho','eta','gamma','phi',
    'idx_cores',
    'd','d_knots','d_inter','w'
    ), 
file=paste('~/workflow_stepps_calibration/prediction/data/prediction_',K,'_taxa_',N,'_cells_',N_knots,'_knots_pl.dump',sep=""))

