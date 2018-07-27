library(rstan)
library(RColorBrewer)
library(fields)
#-----------------------------------------------------------------------------------------------
setwd('~/workflow_stepps_calibration/prediction/')
help.fun.loc <- 'utils/'
data.loc <- 'data/'
plot.loc <- 'plots/'
#-----------------------------------------------------------------------------------------------
taxa <- sort(c('Maple','Birch','Hickory','Chestnut','Other conifer','Other hardwood',
                      'Beech','Ash','Tamarack','Spruce','Pine','Poplar','Oak','Hemlock','Elm'))

#fit <- read_stan_csv(paste(data.loc,stepps_pred_settlement1.csv',sep='')
#fit <- read_stan_csv(paste(data.loc,stepps_pred_settlement_test_wrong_distance_matrix.csv',sep='')#is indeed pretty bad...
#fit <- read_stan_csv(paste(data.loc,stepps_pred_settlement_pl_short_run.csv',sep='')
#fit <- read_stan_csv(paste(data.loc,stepps_pred_settlement_pl_old_weight_matrix.csv',sep='')
#fit <- read_stan_csv(paste(data.loc,stepps_pred_settlement_pl.csv',sep='',sep=''))
fit <- read_stan_csv(paste(data.loc,'stepps_pred_settlement_ka_kgamma.csv',sep=''))
#fit <- read_stan_csv('~/cmdstan-2.6.2/stepps_pred_settlement_g_Kpsi_Kgamma.csv',sep='')
post <- rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
#post <- post[751:1000,1,]
dim(post)
var_names<- colnames(post[,1,])
par_names <- sapply(var_names, function(x) unlist(strsplit(x,'[[]'))[1])
gauss.process <- grep('g',par_names)
gauss.process <- gauss.process[-1]

#---------------------------------------------------------------------------------------------------------------------
post.exp <- exp(post[,1,gauss.process])
#find indexes that start the fun
ash.index <- seq(1,ncol(post.exp),14)

summary.proportion <- 
  sapply(ash.index,function(x){
    exp.site <- cbind(post.exp[,x:(x+13)],rep(1,nrow(post.exp))) # last one is adding exp(0) = 1 baseline taxon
    prop.site <- exp.site/rowSums(exp.site)
    prop.mean.site <- colMeans(prop.site)
    uncertainty <- apply(prop.site,2,function(x) quantile(x,probs=c(0.025,0.5,0.975)))
    data.frame(mean.proportion = prop.mean.site,median = uncertainty['50%',],lb = uncertainty['2.5%',],ub = uncertainty['97.5%',])
  })
summary.proportion <- t(summary.proportion)

mean.proportion <- matrix(ncol=15,unlist(summary.proportion[,'mean.proportion']),byrow=TRUE)
colnames(mean.proportion) <- taxa


mean.proportion <- matrix(ncol=15,unlist(summary.proportion[,'mean.proportion']),byrow=TRUE)
colnames(mean.proportion) <- taxa


median.proportion <- matrix(ncol=15,unlist(summary.proportion[,'median']),byrow=TRUE)
colnames(median.proportion) <- taxa

lb.proportion <- matrix(ncol=15,unlist(summary.proportion[,'lb']),byrow=TRUE)
colnames(lb.proportion) <- taxa

ub.proportion <- matrix(ncol=15,unlist(summary.proportion[,'ub']),byrow=TRUE)
colnames(ub.proportion) <- taxa

uncertainty <- ub.proportion - lb.proportion




#mean.gauss.process <- colMeans(exp(post[,1,gauss.process]))
# mean.gauss.process <- colMeans(post[,1,gauss.process])
# mean.gauss.process <- matrix(mean.gauss.process,ncol = 14,byrow=TRUE)
# #mean.gauss.process <- cbind(mean.gauss.process,rep(1,nrow(mean.gauss.process)))
# mean.gauss.process <- cbind(mean.gauss.process,rep(0,nrow(mean.gauss.process)))
# 
# #test.prop <- t(apply(mean.gauss.process,1,function(x) (x)/sum((x))))          
# test.prop <- t(apply(mean.gauss.process,1,function(x) exp(x)/sum(exp(x))))          
# colnames(test.prop) <- taxa[1:15]

#transfrom pollen coordinates to us coordinates
source('R/prepare_data_settlement_era.R')
sputm <- SpatialPoints(pol_table@coords, proj4string=CRS("+init=epsg:4326"))
spgeo <- spTransform(sputm, CRS("+init=epsg:3175"))
pollen_coord_us <- spgeo@coords
load('~/workflow_stepps_calibration/vegetation/data/veg_data_15_taxa_754_cells_260_knots.rdata')

# head(test.prop)
# rowSums(test.prop)
#pdf(paste(plot.loc,'stepps_settlement/settlement_vegetation_pl.pdf',sep=''),sep=''),width = 10,height = 5)
pdf(paste(plot.loc,'settlement_vegetation_pl_ka_gamma.pdf',sep=''),width = 10,height = 5)
#pdf(paste(plot.loc,'settlement_vegetation_g_psi_gamma.pdf',sep=''),width = 10,height = 5)
#pdf(paste(plot.loc,'settlement_vegetation_g_psi_gamma_1_250.pdf',sep=''),width = 10,height = 5)
#jpeg(paste(plot.loc,'settlement_vegetation_pl.jpeg',width = 10,height = 5,units='in',res = 300)
par(mfrow=c(1,2),oma=c(1,1,1,2),cex = 0.8,cex.axis = 0.8)
sapply(taxa, function(x){
  

breaks <- c(0,0.01,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,1)
categories <- cut(mean.proportion,breaks)
cat.sort <- sort(levels(categories))
categories <- matrix(categories,ncol=ncol(mean.proportion))
colours <- rev(brewer.pal(10,'RdYlBu'))
colours.plot <- matrix(ncol = ncol(categories),nrow=nrow(categories))

for(i in cat.sort) {
  colours.plot[categories==i] <- colours[cat.sort==i]
}
colnames(colours.plot) <- taxa

east <- sort(unique(coord.agg.final$east))
north <- sort(unique(coord.agg.final$north))

breaks1 <- c(0,0.01,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7)
breaks2 <- c(0,0.01,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,NA)

  image.plot(east, north,z = matrix(ncol= length(north),rep(-1,length(east)*length(north))),zlim=c(0,0.7),col = colours,
           breaks= breaks1,lab.breaks = breaks2,main = paste(x,'STEPPS'),cex.axis = 0.8)
  points(coord.agg.final,col=colours.plot[,x],pch = 15)
  points(pollen_coord_us,col=1,pch = 16,cex = 0.75)

#----------------------------------------------------------------------------------------------------------------------
veg_agg_matrix <- as.matrix(veg_agg)

categories <- cut(veg_agg_matrix,breaks)
cat.sort <- sort(levels(categories))
categories <- matrix(categories,ncol=ncol(mean.proportion))


colours.plot <- matrix(ncol = ncol(categories),nrow=nrow(categories))

for(i in cat.sort) {
  colours.plot[categories==i] <- colours[cat.sort==i]
}
colnames(colours.plot) <- taxa

#sapply(taxa, function(x){
  image.plot(east, north,z = matrix(ncol= length(north),rep(-1,length(east)*length(north))),zlim=c(0,0.7),col = colours,
             breaks= breaks1,lab.breaks = breaks2,main = paste(x, 'Paciorek et al. (2016)'),cex.axis = 0.8)
  points(coord.agg.final,col=colours.plot[,x],pch = 15)
  points(pollen_coord_us,col=1,pch = 16,cex = 0.75)
  })
dev.off()






#pdf(paste(plot.loc,'settlement_vegetation_uncertainty_pl.pdf',sep=''),width = 10,height = 5)
pdf(paste(plot.loc,'settlement_vegetation_uncertainty_pl_ka_kgamma.pdf',sep=''),width = 10,height = 5)
#pdf(paste(plot.loc,'settlement_vegetation_uncertainty_g_psi_gamma.pdf',sep=''),width = 10,height = 5)
#pdf(paste(plot.loc,'settlement_vegetation_uncertainty_g_psi_gamma_1_250.pdf',sep=''),width = 10,height = 5)
#jpeg(paste(plot.loc,'settlement_vegetation_pl.jpeg',width = 10,height = 5,units='in',res = 300)
par(mfrow=c(1,2),oma=c(1,1,1,2),cex = 0.8,cex.axis = 0.8)
sapply(taxa, function(x){
  
  
  breaks <- c(0,0.01,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,1)
  categories <- cut(mean.proportion,breaks)
  cat.sort <- sort(levels(categories))
  categories <- matrix(categories,ncol=ncol(mean.proportion))
  colours <- rev(brewer.pal(10,'RdYlBu'))
  colours.plot <- matrix(ncol = ncol(categories),nrow=nrow(categories))
  
  for(i in cat.sort) {
    colours.plot[categories==i] <- colours[cat.sort==i]
  }
  colnames(colours.plot) <- taxa
  
  east <- sort(unique(coord.agg.final$east))
  north <- sort(unique(coord.agg.final$north))
  
  breaks1 <- c(0,0.01,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7)
  breaks2 <- c(0,0.01,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,NA)
  
  image.plot(east, north,z = matrix(ncol= length(north),rep(-1,length(east)*length(north))),zlim=c(0,0.7),col = colours,
             breaks= breaks1,lab.breaks = breaks2,main = paste(x,'STEPPS'),cex.axis = 0.8)
  points(coord.agg.final,col=colours.plot[,x],pch = 15)
  points(pollen_coord_us,col=1,pch = 16,cex = 0.75)
  
  
  breaks <- c(0,0.01,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,1)
  categories <- cut(uncertainty,breaks)
  cat.sort <- sort(levels(categories))
  categories <- matrix(categories,ncol=ncol(mean.proportion))
  colours <- rev(brewer.pal(10,'RdYlBu'))
  colours.plot <- matrix(ncol = ncol(categories),nrow=nrow(categories))
  
  for(i in cat.sort) {
    colours.plot[categories==i] <- colours[cat.sort==i]
  }
  colnames(colours.plot) <- taxa
  
  east <- sort(unique(coord.agg.final$east))
  north <- sort(unique(coord.agg.final$north))
  
  breaks1 <- c(0,0.01,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7)
  breaks2 <- c(0,0.01,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,NA)
  
  image.plot(east, north,z = matrix(ncol= length(north),rep(-1,length(east)*length(north))),zlim=c(0,0.7),col = colours,
             breaks= breaks1,lab.breaks = breaks2,main = paste(x,'97.5pct - 2.5pct'),cex.axis = 0.8)
  points(coord.agg.final,col=colours.plot[,x],pch = 15)
  points(pollen_coord_us,col=1,pch = 16,cex = 0.75)
})
dev.off()




#-------------------------------------------------------------------------------------
#ri = gi/(1+sum[k=1:N-1]gk)
#rN = 1/(1+sum[k=1:N-1]gk)

#auf letzten Wert des Gaussian process auf 0 setzen ist korrekt 

#have to evaluate these predictions:
#   - CRPS (evaluates an ensemble of predictions, does it also evaluate different species at the same time?)
#   - some distance metric (see what Andria used)
#   - make sure metrics are comparable