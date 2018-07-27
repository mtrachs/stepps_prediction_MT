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
K
idx_cores
N_cores
d
y
N <- N_cores
taxa <- colnames(y)
#------------------------------------------------------------------------------------------------
#produce knots
#------------------------------------------------------------------------------------------------
clust_n <- 230
knot_coords = kmeans(veg_coords, clust_n)$centers 
knot_coords = unname(knot_coords)
N_knots     = dim(knot_coords)[1]

#-------------------------------------------------------------------------------------------------------------------------------
# new d_inter
#-------------------------------------------------------------------------------------------------------------------------------
d_inter <- rdist(veg_coords,knot_coords)/10^6
d_knots <- rdist(knot_coords,knot_coords)/10^6
#--------------------------------------------------------------------------------------------
# have to load parameters from calibration model
#--------------------------------------------------------------------------------------------
library(rstan)
source(paste(help.fun.loc,'pred_helper_funs.r',sep='')) #make sure I use medians of a
source(paste(help.fun.loc,'process_funs.r',sep='')) #make sure I use medians of a
# gamma
# a
# b
# phi 
# for settlement era, I do not see if it needs a or b, probably not
# in a first sep I will load phi and gamma from power law kernel Ka_Kgamma
source('R/build_cal_main.r') # this is strange....
for (run in runs) {
kernel <- run$kernel
num_a <- run$one_a
one_psi <- run$one_psi
handle <- strsplit(run$suff_fit,'_A')[[1]][1]

#look at that again...
for (i in 0:9) { # would have to change this if not run_pl
  if(kernel=='pl'){
    if(num_a==FALSE){
  fname = paste('~/stepps_median/output_cal_pl_Ka_Kgamma_EPs_',i,'.csv',sep='') #sprintf('%s/%s.csv', path_out, suff_fit)
  }
  if((num_a==TRUE)&(i<5)){
  fname = paste('~/workflow_stepps_calibration/results/data/stepps_median/output_cal_pl_',i,'.csv',sep='') #sprintf('%s/%s.csv', path_out, suff_fit)
  }
  }
  if(kernel=='gaussian'){
    if(one_psi==FALSE){
  fname = paste('~/stepps_median/output_cal_g_Kpsi_Kgamma_EPs_',i,'.csv',sep='') #sprintf('%s/%s.csv', path_out, suff_fit)
  }
  if(one_psi==TRUE){
    fname = paste('~/stepps_median/output_neus_',i,'.csv',sep='') #sprintf('%s/%s.csv', path_out, suff_fit)
  }}
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
w <- build_weight_matrix(post = post,d = t(d),idx_cores = idx_cores,
                         N = N_cells,N_cores =  N_cores,run = run)

w <- w/sum(w)

sum_w_pot <- build_sumw_pot(post = post, K=K,N_pot = N_pot,
                            d_pot =  d_pot,run = run)

# don't remeber what this was doing
# w1 <- array(dim=c(K,N_cores,N_cells))
# for(i in 1:15) {w1[i,,]<-w}
# w <- w1
#-------------------------------------------------------------------------------------------
# #have to rescale gamma
# # simplest is to build the kernel and add the weight of the 8 closest locations...
# # 
# weights_pot <- build_w_pot(post = post, K=K,N_pot = N_pot,d_pot =  d_pot,run = run)
# cbind(weights_pot[,2],d_pot[,2])
# gamma.added <- (weights_pot[1,] + weights_pot[2,]/2)/sum_w_pot
# gamma.new <- gamma + gamma.added
# gamma <- gamma.new
#--------------------------------------------------------------------------------------------
res <- 1

#-------------------------------------------------------------------------------------------
# eta (spatial correlation, correaltion strength)
# rho (spatial correlation decorrealtion length)
fit <- read_stan_csv('~/workflow_stepps_calibration/vegetation/output/pls_correct.csv')
post <- rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
param.names <-colnames(post[,1,]) #find parameter names of posterior
param.greek <- sapply(param.names,function(x) strsplit(x,'[[]')[[1]][1]) #only take greek letter
param.e.bayes <- param.greek%in%c('eta','rho') #only take eta and rho
est.e.bayes <- colMeans(post[,1,param.e.bayes])

eta.all <- est.e.bayes[grep('eta',names(est.e.bayes))]
eta <- eta.all#[-length(eta.all)]
rho.all <- est.e.bayes[grep('rho',names(est.e.bayes))]
rho <- rho.all#[-length(rho.all)]

taxa <- taxa[taxa!='Chestnut']
names(eta) <- taxa[-length(taxa)]
names(rho) <- taxa[-length(taxa)]

#--------------------------------------------------------------------------------------------------------------------
#make up a value for Chestnut
taxa1 <- sort(c(taxa,'Chestnut'))
eta1 <- rep(NA,length(taxa))
names(eta1) <- taxa1[-length(taxa1)]
eta1[names(eta1)%in%names(eta)] <- eta

rho1 <- rep(NA,length(taxa))
names(rho1) <- taxa1[-length(taxa1)]
rho1[names(rho1)%in%names(rho)] <- rho

eta1[names(eta1)=='Chestnut'] <- mean(eta)
rho1[names(rho1)=='Chestnut'] <- mean(rho)

#assign for use in .dump file
rho <- rho1
eta <- eta1
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
file=paste('~/workflow_stepps_calibration/vegetation/data/prediction_',K,'_taxa_',N_cells,'_cells_',N_knots,'_knots_',handle,'.dump',sep=""))


save(K, N,T, N_knots, N_cores, 
       y, res,
       sum_w_pot,
       rho,eta,gamma,phi,
       idx_cores,
       d,d_knots,d_inter,w
, 
file=paste('~/workflow_stepps_calibration/vegetation/data/prediction_',K,'_taxa_',N_cells,'_cells_',N_knots,'_knots_',handle,'.rdata',sep=""))
}
