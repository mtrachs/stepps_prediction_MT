#--------------------------------------------------------------------------------------------------------------------------------------------
#In this script vegetation is prepared to estimate the spatial model
#--------------------------------------------------------------------------------------------------------------------------------------------
library(fields)
library(rstan)
library(stepps)

#-------------------------------------------------------------------------------------------------------------------
setwd('~/workflow_stepps_prediction/vegetation/')
help.fun.loc <- 'helper_funs/'
data.loc <- 'data/'
plot.loc <- 'plots/'
#-------------------------------------------------------------------------------------------------------------------------------------------
#load data that was used for calibration (is the same vegetation data)
#-------------------------------------------------------------------------------------------------------------------------------------------
pls.counts <- read.csv('data/wiki_outputs/plss_composition_alb_v0.9-10.csv')
veg_coords <- pls.counts[,c('x','y')]

#merge species into other hardwood and other conifer
#colnames(pls.counts)

pls_table <- readr::read_csv('~/workflow_stepps_calibration/calibration/data/veg_trans_edited.csv')
pls_trans <- translate_taxa(pls.counts, pls_table ,id_cols = colnames(pls.counts)[1:3])

y <- pls_trans[,-c(1:3)]




#--------------------------------------------------------------------------------------------------------------------------------------------------
# extract coordinates of knots through k-means 
# use k-means becasue it estiamtes a center of the coordinates
#-----------------------------------------------------------------------------------------------------------------------------------------------
clust_n <- 260
knot_coords = kmeans(veg_coords, clust_n)$centers 
knot_coords = unname(knot_coords)


#--------------------------------------------------------------------------------------------------------------------
#plot coordinates and compute a distance matrix
#--------------------------------------------------------------------------------------------------------------------
plot(veg_coords,pch =15,cex = 0.25)
points(knot_coords,pch = 15,col='red')

distances <- stats::dist(knot_coords)
distances1 <- as.matrix(distances)

#find distance between nearest neighbours
min.dist <- apply(distances1,2,function(x) min(x[x>0]))
s.m.d <- summary(min.dist/10^3)

#--------------------------------------------------------------------------------------------------------------------
hist(distances/10^3)
hist(min.dist/10^3)
#--------------------------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------------------------
#look at data from NEUS
#--------------------------------------------------------------------------------------------------------------------
library(fields)


#-------------------------------------------------------------------------------------------------------------------------------------------
#load data that was used for calibration (is the same vegetation data)
#-------------------------------------------------------------------------------------------------------------------------------------------
load('~/stepps_data/elicitation_neus_certainty_median.RData')


#--------------------------------------------------------------------------------------------------------------------------------------------------
# extract coordinates of knots through k-means 
# use k-means becasue ti estiamtes a center of the coordinates
#-----------------------------------------------------------------------------------------------------------------------------------------------
clust_n <- 230
knot_coords_neus = kmeans(veg_coords, clust_n)$centers 
knot_coords_neus = unname(knot_coords_neus)


#--------------------------------------------------------------------------------------------------------------------
#plot coordinates and compute a distance matrix
#--------------------------------------------------------------------------------------------------------------------
plot(veg_coords,pch =15,cex = 0.25)
points(knot_coords_neus,pch = 15,col='red')

distances.neus <- stats::dist(knot_coords_neus)
distances.neus1 <- as.matrix(distances.neus)

#find distance between nearest neighbours
min.dist.neus <- apply(distances.neus1,2,function(x) min(x[x>0]))
s.m.d.neus <- summary(min.dist.neus/10^3)

hist(distances.neus/10^3)
hist(min.dist.neus/10^3)

#compare distances umw and neus
cbind(s.m.d,s.m.d.neus)

#230 knots seems to give about the same distribution...
#--------------------------------------------------------------------------------------------------------------------
#distances are slightly smaller in the NEUS use fewer knots?



