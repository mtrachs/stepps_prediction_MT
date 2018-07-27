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


#merge species into other hardwood and other conifer
#colnames(pls.counts)

pls_table <- readr::read_csv('~/workflow_stepps_calibration/calibration/data/veg_trans_edited_pls.csv')
pls_trans <- translate_taxa(pls.counts, pls_table ,id_cols = colnames(pls.counts)[1:3])


veg_coords <-pls_trans[,c('x','y')] 
#----------------------------------------------------------------------------------------------------------------------
#code to eliminate scarce parts of the data matrix
xs <- sort(unique(pls.counts$x))
ys <- sort(unique(pls.counts$y))

daten <- data.frame(x = c(45000,245000),y=c(1054000,814000))

lm.cord <- lm(y~x,data=daten)
x.data <- data.frame(x = xs[((xs>=min(daten$x))&(xs<=max(daten$x)))])
y.pred <- predict(lm.cord,x.data)

coord.reg <- data.frame(x = x.data,y = y.pred)


remove <- sapply(1:nrow(coord.reg),function(x) {
  which(((veg_coords$x <= coord.reg$x[x]) & 
           (veg_coords$y <= coord.reg$y[x]))| (veg_coords$x <= min(coord.reg$x)))
})

remove1 <- unique(unlist(remove))

pls_trans <- pls_trans[-remove1,]
y <- pls_trans[,-c(1:3)]
veg_coords <- pls_trans[,1:2]
plot(veg_coords,pch = 16)
#--------------------------------------------------------------------------------------------------------------------------------------------------
# extract coordinates of knots through k-means 
# use k-means becasue ti estiamtes a center of the coordinates
#-----------------------------------------------------------------------------------------------------------------------------------------------
clust_n <- 260
knot_coords = kmeans(veg_coords, clust_n)$centers 
knot_coords = unname(knot_coords)
N_knots     = dim(knot_coords)[1]
K <- ncol(y)
d <- rdist(veg_coords,veg_coords)/10^6

#-------------------------------------------------------------------------------------------------------------------
plot(veg_coords,pch = 16,cex = 0.5)
points(knot_coords,col=2,pch =16)
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
taxa <- colnames(y)

stan_rdump(c('K', 'N', 'N_knots', 
             'y', 
             'd_knots', 'd_inter',
             'P', 'N_p'), 
           file=paste('data/pls_data_', K, '_taxa_', N, '_cells_', N_knots, '_knots.dump',sep=""))

# K = number of species, N = number of vegetation grid points,N_knots = number of knots
# y = vegetation data (rounded to real number), d_knots = distance matrix among knots, d_inter = distance matrix between knots and vegetation grid cells 
# N_p = N, P = N x N matrix all values the same 

save(K, N, N_knots, 
     y, 
     d, d_knots, d_inter, 
     P, N_p,
     veg_coords, knot_coords, taxa,
     file=paste('data/pls_data_', K, '_taxa_', N, '_cells_', N_knots, '_knots.rdata',sep=""))


#-------------------------------------------------------------------------------------------------------------------------------------------
#make pie maps
#------------------------------------------------------------------------------------------------------------------------------------------
# 
# library(maptools)
# 
# 
# setwd('~/STEPPS/stepps-calibration-master/')
# source('r/utils/dataPlotFuns.r')
# source('r/utils/simDataFuns.r')

# load pls and pollen data
# pls      = read.table(file='data/pls_UMW.csv', sep=",", row.names=NULL, header=TRUE)
# pollen   = read.table('data/pollen_se_sum.csv', header=TRUE)
# load('r/dump/cal_data_12taxa_mid_comp_UMW_v0.2.rdata')

# load('r/dump/cal_data_12taxa_mid_mi_sp.rdata')
#load('r/dump/cal_data_12taxa_mid_comp_ALL_v0.3.rdata')

#v = '0_3'
# v = 'mi_sp'

# ########################################################################################################
# ## pls pie map
# ########################################################################################################
# # reorder data by pollen abundance
# new.order = order(colSums(y), decreasing=TRUE)
# taxa <- colnames(y)
# taxa = taxa[new.order]
# #y = veg_knots[,new.order]
# 
# 
# r1 = y[,new.order]
# r1 <- r1/rowSums(r1)
# #colnames(r) = taxa
# 
# # col_list = c("#1F78B4", "#33A02C", "#E31A1C", "#6A3D9A", "#B15928", "#FF7F00", 
# #              "#FFFF99", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F", "#CAB2D6")
# # 
# # stepps_cols = data.frame(taxa=taxa, cols=col_list)
# # save(stepps_cols, file='r/stepps_cols.rdata')
# # 
# # taxa[which(taxa == 'OTHER.CONIFER')] = 'FIR'
# # taxa[which(taxa == 'OTHER.HARDWOOD')] = 'OTHER HARDWOOD'
# ########################################################################################################
# ## pls pie map
# ########################################################################################################
# 
# ntaxa = length(taxa)
# 
# centers   = coord.agg.final
# colnames(centers) = c('x', 'y')
# 
# xlo = min(centers[,1])
# xhi = max(centers[,1])
# ylo = min(centers[,2])
# yhi = max(centers[,2])
# 
# 
# # pls_props  = t(apply(pls_coarse, 1, function(x) if (sum(x) != 0){x/sum(x)} else {x}))
# 
# shift=30000
# 
# # postscript('r/data/figs/pie_plot_pls_UMW_v0.2.eps', width=8, height=6)
# pdf(paste('~/stepps_process_knots/pie_map_aggregated_vegetation.pdf',sep=''), width=12, height=10)
# par(mfrow=c(1,1),mar =c(2,2,2,2))
# pieMap(proportions = r1, 
#        centers  = centers,
#        restrict = FALSE,
#        inputRestricted = FALSE,
#        xlim   = c(xlo+shift, xhi-shift),
#        ylim   = c(ylo+shift, yhi-shift),
#        radius = 10500,
#        scale  = 1,
#        xlab   = 'x',
#        ylab   = 'y',
#        add_legend = TRUE, 
#        main_title='')
# dev.off()

