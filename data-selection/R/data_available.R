#check how to get radiocarbon data:
library(neotoma)
library(maps)
library(maptools)
#####################################################################################################################
# load data used for age-depth modelling
#####################################################################################################################
folder.loc <- '~/workflow_stepps_prediction/'
setwd(paste(folder.loc,'age-depth',sep=''))
plot.loc <- 'plots/'
data.loc <- 'data/'
bchron.plot.loc <- 'Bchron/plots/sample_depth/'
bchron.output.loc <- 'Bchron/output/'

data.NE <- readRDS(paste(data.loc,'data_NE.RDS',sep=''))
geocontrol.NE <- sapply(data.NE,function(x) x$geo)

sapply(1:length(geocontrol.NE),function(x){
  dim(geocontrol.NE[[x]]$chron.control)[1] >1
})








######################################################################################################################
# setwd
######################################################################################################################
file.loc <- '~/workflow_stepps_prediction/'
setwd(paste(file.loc,'data-selection/',sep=''))
plot.loc <-'plots/' 

#---------------------------------------------------------------------------------------------------------------------
#First, find sites in the northeastern US
#---------------------------------------------------------------------------------------------------------------------
#first the brackeeting ages are defined, we only consider cores that have data in this time period
older.age <- 6000
younger.age <- -70



get_table('GeoPoliticalUnits') 

#load fossil data for NEUS (given by coordinates in loc )
new_england_fossil_age <- get_dataset(datasettype = 'pollen',ageold = older.age,ageyoung = younger.age,
                                  loc=c(-81,39.75,-67,47.55)) #returns 276 reconrds,with age constraints 231 samples

coordinates_neus_age <- sapply(1:length(new_england_fossil_age),function(x){
  lon <- new_england_fossil_age[[x]]$site.data$long
  lat <- new_england_fossil_age[[x]]$site.data$lat
  data.frame(lon=lon,lat=lat)
})  

coordinates_neus_age <- t(coordinates_neus_age)

new_england_fossil <- get_dataset(datasettype = 'pollen',
                                      loc=c(-81,39.75,-67,47.55)) #returns 276 reconrds,with age constraints 231 samples

coordinates_neus <- sapply(1:length(new_england_fossil),function(x){
  lon <- new_england_fossil[[x]]$site.data$long
  lat <- new_england_fossil[[x]]$site.data$lat
  data.frame(lon=lon,lat=lat)
})  



coordinates_neus <- t(coordinates_neus)

pdf(paste(plot.loc,'site_locations.pdf',sep=''),width = 11,height= 10)
  map('state',xlim=c(-82,-67),ylim=c(38,47.55))
  points(coordinates_neus,pch = 16,col=4)
  points(coordinates_neus_age,pch = 16)
  box()
  axis(1)
  axis(2)
  mtext(side=1,lin=2.2,'Longitude')
  mtext(side=2,lin=2.2,'Latitude')
  legend('bottomright',col=c(1,4),pch = rep(16,2),legend = c('Age','No age'))
dev.off()
