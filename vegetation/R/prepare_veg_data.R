#--------------------------------------------------------------------------------------------------------------------------------------------
#In this script vegetation is prepared to estimate the spatial model
#--------------------------------------------------------------------------------------------------------------------------------------------
library(fields)
library(rstan)

#-------------------------------------------------------------------------------------------------------------------
setwd('~/workflow_stepps_prediction/vegetation//')
help.fun.loc <- 'helper_funs/'
data.loc <- 'data/'
plot.loc <- 'plots/'
#-------------------------------------------------------------------------------------------------------------------------------------------
#load data that was used for calibration (is the same vegetation data)
#-------------------------------------------------------------------------------------------------------------------------------------------
load('~/stepps_data/elicitation_neus_certainty_median.RData')
y <- ceiling(100*r) #change this to 100 and run again

#--------------------------------------------------------------------------------------------------------------------------------------------------
# extract coordinates of knots through k-means 
# use k-means becasue ti estiamtes a center of the coordinates
#-----------------------------------------------------------------------------------------------------------------------------------------------
clust_n <- 260
knot_coords = kmeans(veg_coords, clust_n)$centers 
knot_coords = unname(knot_coords)
N_knots     = dim(knot_coords)[1]

#-------------------------------------------------------------------------------------------------------------------------------
# new d_inter
#-------------------------------------------------------------------------------------------------------------------------------
d_inter <- rdist(veg_coords,knot_coords)/10^6
d_knots <- rdist(knot_coords,knot_coords)/10^6
##############################################################################################################################
N = dim(veg_coords)[1]
x = matrix(1, nrow=N, ncol=1)
N_p = N

temp = qr(x)
Q = qr.Q(temp)
R = qr.R(temp)

P = Q %*% t(Q)
# M = diag(N) - P

if (all(P-P[1,1]<1.0e-12)){
  P = P[1,1]
  N_p = 1
}
##########################################################################################################################
## chunk: save data to file
##########################################################################################################################
taxa <- colnames(r)

stan_rdump(c('K', 'N', 'N_knots', 
             'y', 
             'd_knots', 'd_inter',
             'P', 'N_p'), 
           file=paste('data/veg_data_', K, '_taxa_', N, '_cells_', N_knots, '_knots.dump',sep=""))

# K = number of species, N = number of vegetation grid points,N_knots = number of knots
# y = vegetation data (rounded to real number), d_knots = distance matrix among knots, d_inter = distance matrix between knots and vegetation grid cells 
# N_p = N, P = N x N matrix all values the same 

save(K, N, N_knots, 
     y, 
     d, d_knots, d_inter, 
     P, N_p,
     veg_coords, knot_coords, taxa,
     file=paste('data/veg_data_', K, '_taxa_', N, '_cells_', N_knots, '_knots.rdata',sep=""))


#-------------------------------------------------------------------------------------------------------------------------------------------
#make pie maps
#------------------------------------------------------------------------------------------------------------------------------------------

library(maptools)


setwd('~/STEPPS/stepps-calibration-master/')
source('r/utils/dataPlotFuns.r')
source('r/utils/simDataFuns.r')

# load pls and pollen data
# pls      = read.table(file='data/pls_UMW.csv', sep=",", row.names=NULL, header=TRUE)
# pollen   = read.table('data/pollen_se_sum.csv', header=TRUE)
# load('r/dump/cal_data_12taxa_mid_comp_UMW_v0.2.rdata')

# load('r/dump/cal_data_12taxa_mid_mi_sp.rdata')
#load('r/dump/cal_data_12taxa_mid_comp_ALL_v0.3.rdata')

#v = '0_3'
# v = 'mi_sp'

########################################################################################################
## pls pie map
########################################################################################################
# reorder data by pollen abundance
new.order = order(colSums(y), decreasing=TRUE)
taxa <- colnames(y)
taxa = taxa[new.order]
#y = veg_knots[,new.order]


r1 = y[,new.order]
r1 <- r1/rowSums(r1)
#colnames(r) = taxa

# col_list = c("#1F78B4", "#33A02C", "#E31A1C", "#6A3D9A", "#B15928", "#FF7F00", 
#              "#FFFF99", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F", "#CAB2D6")
# 
# stepps_cols = data.frame(taxa=taxa, cols=col_list)
# save(stepps_cols, file='r/stepps_cols.rdata')
# 
# taxa[which(taxa == 'OTHER.CONIFER')] = 'FIR'
# taxa[which(taxa == 'OTHER.HARDWOOD')] = 'OTHER HARDWOOD'
########################################################################################################
## pls pie map
########################################################################################################

ntaxa = length(taxa)

centers   = coord.agg.final
colnames(centers) = c('x', 'y')

xlo = min(centers[,1])
xhi = max(centers[,1])
ylo = min(centers[,2])
yhi = max(centers[,2])


# pls_props  = t(apply(pls_coarse, 1, function(x) if (sum(x) != 0){x/sum(x)} else {x}))

shift=30000

# postscript('r/data/figs/pie_plot_pls_UMW_v0.2.eps', width=8, height=6)
pdf(paste('~/stepps_process_knots/pie_map_aggregated_vegetation.pdf',sep=''), width=12, height=10)
par(mfrow=c(1,1),mar =c(2,2,2,2))
pieMap(proportions = r1, 
       centers  = centers,
       restrict = FALSE,
       inputRestricted = FALSE,
       xlim   = c(xlo+shift, xhi-shift),
       ylim   = c(ylo+shift, yhi-shift),
       radius = 10500,
       scale  = 1,
       xlab   = 'x',
       ylab   = 'y',
       add_legend = TRUE, 
       main_title='')
dev.off()

