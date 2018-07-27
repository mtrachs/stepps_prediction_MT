library(stepps)#
library(dplyr)#
library(DT)#
library(neotoma)#
library(sp)#
library(fields)#
library(readr)
library(rioja)
#-------------------------------------------------------------------------------------------------------------------
setwd('~/workflow_stepps_calibration/calibration/')
help.fun.loc <- 'calibration_helper_funs/'
data.loc <- 'data/'
plot.loc <- 'plots/'
#-------------------------------------------------------------------------------------------------------------------

#load vegetation data using other code

if('veg_mean.RDS'%in%list.files(data.loc)){
  veg_mean <- readRDS(paste(data.loc,'veg_mean.RDS',sep=''))
  coords.neus <- matrix(ncol=2,c(rep(c(-80.5,-67),each =2),rep(c(39.5,47.5),2)))
  coords.neus <- as.data.frame(coords.neus)
  colnames(coords.neus) <- c('Lon','Lat')

  pol_box <- bbox_tran(coords.neus, '~ Lon + Lat',
                       '+init=epsg:4326',
                       '+init=epsg:4326')

  veg_box <- bbox_tran(coords.neus, '~ Lon + Lat',
                       '+init=epsg:4326',
                       '+init=epsg:3175')

  reconst_grid <- build_grid(veg_box, resolution = 8000, proj = '+init=epsg:3175')


}else {
  source(paste(help.fun.loc,'load_vegetation_data_median.R',sep=''))
  veg_mean <- total.dat
  colnames(veg_mean)[1:2] <- c('meters.east','meters.north')

  #veg_mean <- readr::read_csv('data/composition_v0.3.csv')
  coords.neus <- matrix(ncol=2,c(rep(c(-80.5,-67),each =2),rep(c(39.5,47.5),2)))
  coords.neus <- as.data.frame(coords.neus)
  colnames(coords.neus) <- c('Lon','Lat')

  pol_box <- bbox_tran(coords.neus, '~ Lon + Lat',
                       '+init=epsg:4326',
                       '+init=epsg:4326')

  veg_box <- bbox_tran(coords.neus, '~ Lon + Lat',
                       '+init=epsg:4326',
                       '+init=epsg:3175')

  #extract vegetation data that is in the bounding box
  veg_mean <- veg_mean[(veg_mean$meters.east>veg_box[1]) & (veg_mean$meters.north>veg_box[2]),]



  veg_mean <- veg_mean[!is.na(veg_mean$Ash),]
  veg_mean <- replace(veg_mean,veg_mean==0,2e-7) #does this have an effect
  reconst_grid <- build_grid(veg_box, resolution = 8000, proj = '+init=epsg:3175')

  setwd('~/workflow_stepps_calibration/calibration/')
  saveRDS(veg_mean,paste(data.loc,'veg_mean.RDS',sep=''))
}





source(paste(help.fun.loc, 'get_meta_data_cal.R',sep=''))#neotoma::get_dataset(loc = pol_box, datasettype = 'pollen')#perhaps use Andrias Code for that
datasets <- meta 

#--------------------------------------------------
#load site.ids of sites that have data between present and 4000 cal BP

site_ids_us <- read.table('~/workflow_stepps_calibration/expert_elicitation/data/site_ids.txt',header=TRUE)
ind.plot <- read.table('~/workflow_stepps_calibration/expert_elicitation/data/plot_index.txt')

# find site ids that were plotted and ultimately used
site_ids_us <- site_ids_us[[1]][ind.plot[[1]]]


#----------------------------------------------------------------------------------------------------------

if(!'downloads.rds' %in% list.files('data/')) {
  downloads <- neotoma::get_download(datasets)
  saveRDS(downloads, paste(data.loc,'downloads.rds',sep=''))
} else {
  downloads <- readRDS(paste(data.loc,'downloads.rds',sep=''))
}


#find data that has been plotted
downloads.clean <-
  lapply(1:length(downloads), function(x) {
    daten <- NA
    if(downloads[[x]]$dataset$dataset.meta$dataset.id %in% site_ids_us) daten <- downloads[[x]]
    daten
  })

sites.pull.use <- which(sapply(1:length(downloads.clean), function(x) unique(!is.na(downloads.clean[[x]]))))

downloads.clean <- lapply(sites.pull.use, function(x) downloads.clean[[x]])

#try to sort site ids so that sample.id is ascending (is alrady done in fact not necessary)
site.ids.downloads.clean <- sapply(1:length(downloads.clean), function(x) 
  downloads.clean[[x]]$dataset$dataset.meta$dataset.id)


#site.ids.downloads.clean[order(site.ids.downloads.clean)]
downloads.clean1 <- lapply((1:length(downloads.clean))[order(site.ids.downloads.clean)],function(x) downloads.clean[[x]])
#-----------------------------------------------------------------------------------------------------------------------
#load evaluation of elecitation exercise
#source(paste(help.fun.loc,'evaluate_elicitation_certainty.R',sep=''))

#(1:length(sites.use))[sites.use]
#sample.use


#--------------------------------------------------------------------------------------------------------------------
# find sites that are outside the domain of vegetation
#immediately remove the three sites that were not used in elicitation
if('meta_data.RDS'%in%list.files('~/workflow_stepps_calibration/expert_elicitation/data/')){
  meta.data.neus <- readRDS('~/workflow_stepps_calibration/expert_elicitation/data/meta_data.RDS')
}else {
  source('~/workflow_stepps_calibration/expert_elicitation/elicitation_helper_funs/get_meta_data.R')
  saveRDS(meta.data.neus,'~/workflow_stepps_calibration/expert_elicitation/data/meta_data.RDS')
}  

ind.plot <- ind.plot[[1]]
site.real <- meta.data.neus$site[order(meta.data.neus$datasetID)][ind.plot]
long.real <- meta.data.neus$long[order(meta.data.neus$datasetID)][ind.plot]
lat.real <- meta.data.neus$lat[order(meta.data.neus$datasetID)][ind.plot]


#transform to us coordinates
sputm <- SpatialPoints(cbind(long.real,lat.real), proj4string=CRS("+init=epsg:4326"))
spgeo <- spTransform(sputm, CRS("+init=epsg:3175"))

lakes.coord <- spgeo@coords
#calcualte distances
dist.lakes <- paldist2(lakes.coord,cbind(veg_mean$meters.east,veg_mean$meters.north), dist.method="euclidean")
min.dist.lakes <- apply(dist.lakes,1,min)
min.dist.lakes.index <- which(min.dist.lakes<8000)
#----------------------------------------------------------------------------------------------------------------------
#total removal criterion
total.exclude <- ((min.dist.lakes > 8000))#|((long.real> (-70.7))& (lat.real<43))|((long.real> (-70))& (lat.real<44.5))|

#######################################################################################################################
total.use <- which(total.exclude==FALSE)

downloads.clean2 <- lapply(total.use,function(x) downloads.clean1[[x]])

#######################################################################################################################
# load taxon translation table 
#######################################################################################################################
source(paste(help.fun.loc,'taxa_translation_only_abies.R',sep=''))
pol_table <- readr::read_csv(paste(data.loc,'taxon_names_translated_other_hardwood.csv',sep=''))


#perhaps need to comment this
#######################################################################################################################
setwd('~/workflow_stepps_prediction/')
pollen.output.loc <- 'pollen_data/'

site.names <- matrix(length(downloads.clean2))
 
for(x in 1:length(downloads.clean2)){
  counts <- downloads.clean2[[x]]$counts
  colnames(counts) <- gsub(' ','.',colnames(counts))
  colnames(counts) <- gsub('/','.',colnames(counts))
  id <- 1:nrow(counts)
  counts_id <- as.data.frame(cbind(id,counts))
  colnames(counts_id)[1] <-'ID'
  counts_trans <- translate_taxa(counts_id, 
                              pol_table,
                              id_cols = colnames(counts_id)[1])
  
  sample.depth <- downloads.clean2[[x]]$sample.meta$depth
  dataset.id <- downloads.clean2[[x]]$dataset$dataset.meta$dataset.id
  site.name <- downloads.clean2[[x]]$dataset$site.data$site.name
  if(length(sample.depth)!=nrow(counts_trans)){
   counts_save_int <- matrix(nrow=length(sample.depth),ncol=(ncol(counts_trans)),data=0)
   counts_matrix <- matrix(nrow=nrow(counts_trans),unlist(counts_trans))
   colnames(counts_matrix) <- colnames(counts_trans)
   counts_save_int[rowSums(counts)>0,] <- counts_matrix
   colnames(counts_save_int) <- colnames(counts_trans)
   counts_save <- data.frame(depths = sample.depth,counts_save_int[,-1]) 
  }else{
  counts_save <- data.frame(depths = sample.depth,counts_trans[,-1]) 
  }
  site.names[x] <- paste(site.name,'_dataset_id_',dataset.id,sep='')
  write.table(counts_save, paste(pollen.output.loc,'pollen_counts_',site.name,'_dataset_id_',dataset.id,".csv",sep=''),
              sep=',', col.names = TRUE, row.names = FALSE)
  
}

write.table(site.names,paste(pollen.output.loc,"site_names.txt",sep=''))

#######################################################################################################################
# produce some meta data
#######################################################################################################################
meta.data.sites  <-
  sapply(1:length(downloads.clean2),function(z){
    daten <- downloads.clean2[[z]]$dataset
    data.frame(dataset.id = daten$dataset.meta$dataset.id,
             site.name = daten$site.data$site.name,
             lat = daten$site.data$lat,
             lon = daten$site.data$long)
  })

saveRDS(meta.data.sites,paste(pollen.output.loc,"meta_data_sites.RDS",sep=''))
