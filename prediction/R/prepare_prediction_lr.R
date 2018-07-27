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

#-----------------------------------------------------------------------------------------------
setwd('~/workflow_stepps_prediction/prediction/')
help.fun.loc <- 'utils/'
data.loc <- 'data/'
plot.loc <- 'plots/'

########################################################################################################################
#read the huge posterior file of the vegetation model
########################################################################################################################
load('~/workflow_stepps_calibration/vegetation/data_township/township_data_13_taxa_6796_cells_120_knots.rdata')
#fit <- read_stan_csv('~/workflow_stepps_calibration/vegetation/output/veg_township_output_v2171_run1.csv')
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
#########################################################################################################################
load('~/workflow_stepps_calibration/calibration/data/elicitation_neus_certainty_median_80_sites_only_abies_new_species.RData')

  #pollen_coords
  # K
  # idx_cores
  # N_cores
  # d
  # y
  #N <- N_cells
  
  #------------------------------------------------------------------------------------------------
  #produce knots
  #------------------------------------------------------------------------------------------------
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
  runs1 <- list(runs[[3]],runs[[4]])
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
        }}
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


#######################################################################################################################
#load pollen data and coordinates of pollen to make a new idx_files
#######################################################################################################################
y <- readRDS('~/workflow_stepps_prediction/pollen_data/pollen_1700_1800.RDS')
y <- round(y)
pollen_coords <- readRDS('~/workflow_stepps_prediction/pollen_data/coordinates_1700_1800.RDS') 
pollen_coords <- matrix(ncol=2,as.numeric(pollen_coords))

sputm <- SpatialPoints(pollen_coords, proj4string=CRS("+init=epsg:4326"))
spgeo <- spTransform(sputm, CRS("+init=epsg:3175"))
pollen_coord_us <- spgeo@coords

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


# dist.agg.raw <- rdist(veg_coords,coord.agg.final)
# 
# #find point to which distance is mnimal (data is aggregated to data grid point)
# min_dist <-  apply(dist.agg.raw,1,function(x)  which.min(x))
# sort(table(min_dist))# see if the minimum distance is about an equal number of veg_cells contributes to one knot

#######################################################################################################################
N_cells <- nrow(coord.agg.final)
N <- N_cells

dist.cores <- paldist2(pollen_coord_us,coord.agg.final,dist.method = 'euclidean')

idx_cores<- apply(dist.cores,1,function(x) which.min(x))
########################################################################################################################
#make d_pot
########################################################################################################################
gr <- seq(24,700,24)
gr.tot <- c(-rev(gr),0,gr)
grid.tot <- expand.grid(gr.tot,gr.tot)
total.distance <- paldist2(matrix(ncol=2,c(0,0)),grid.tot,dist.method = 'euclidean')
numbers <- table(total.distance)
real.distances <- as.numeric(names(numbers))
d_pot <- cbind(real.distances,numbers)
d_pot <- d_pot[-1,]
d_pot <- d_pot[d_pot[,1]<700,]
d_pot[,1] <-d_pot[,1]/10^3
########################################################################################################################




#these objects need changing (probably best to use functions that makes STEPPS data, no they might re-arrange things)

N_cores <- nrow(pollen_coords)
d <- t(dist.cores)
d <- d/10^6


w <- build_weight_matrix(post = post,d = t(d),idx_cores = idx_cores,
                         N = N_cells,N_cores =  N_cores,run = run)

sum_w_pot <- build_sumw_pot(post = post, K=K,N_pot = N_pot,
                            d_pot =  d_pot,run = run)

if(length(gamma)==1) {
  gamma <- rep(gamma,K)
 w1 <- array(dim=c(K,N_cores,N_cells))
 for(i in 1:K) {w1[i,,]<-w}
 w <- w1
}

#################################################################################################################
# knot distance matrix
#################################################################################################################
d_inter <- rdist(coord.agg.final,knot_coords)/10^6
d_knots <- rdist(knot_coords,knot_coords)/10^6
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
res <- 3 # 
ages    = seq(150,1950,100)
T       = length(ages) 
lag     = unname(as.matrix(dist(matrix(ages), upper=TRUE)))

#-------------------------------------------------------------------------------------------
num_sites <- nrow(pollen_coord_us)

#--------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------
#store data in .dump file
stan_rdump(
  c('K', 'N','T', 'N_knots', 'N_cores','lag', 
    'y', 'res',
    'sum_w_pot',
    'rho','eta','gamma','phi',
    'idx_cores',
    'd','d_knots','d_inter','w'
    ), 
file=paste('~/workflow_stepps_prediction/prediction/data/prediction_',K,'_taxa_',N_cells,'_cells_',N_knots,'_knots_',handle,'_',num_sites,'_sites_full_pred.dump',sep=""))


save(K, N,T, N_knots, N_cores,lag, 
       y, res,
       sum_w_pot,
       rho,eta,gamma,phi,
       idx_cores,
       d,d_knots,d_inter,w
, 
file=paste('~/workflow_stepps_prediction/prediction/data/prediction_',K,'_taxa_',N_cells,'_cells_',N_knots,'_knots_',handle,'_',num_sites,'_sites_full_pred.rdata',sep=""))



#################################################################################################################################
#
#################################################################################################################################

#y <- y[seq(1,nrow(y),19),]
T <- 1

stan_rdump(
  c('K', 'N','T', 'N_knots', 'N_cores','lag', 
    'y', 'res',
    'sum_w_pot',
    'rho','eta','gamma','phi',
    'idx_cores',
    'd','d_knots','d_inter','w'
  ), 
  file=paste('~/workflow_stepps_prediction/prediction/data/prediction_test.dump',sep=""))


save(K, N,T, N_knots, N_cores,lag, 
     y, res,
     sum_w_pot,
     rho,eta,gamma,phi,
     idx_cores,
     d,d_knots,d_inter,w,
     coord.agg.final, 
     file=paste('~/workflow_stepps_prediction/prediction/data/prediction_test.rdata',sep=""))
