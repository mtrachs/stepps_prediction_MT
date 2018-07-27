###########################################################################################################################
# load ages and pollen data
###########################################################################################################################
setwd('~/workflow_stepps_prediction/')
chron.loc <- 'age-depth/Bchron/output/'
pollen.loc <- 'pollen_data/'  

chron_prefix <- 'pollen_sample_depth_'
pollen_prefix <-'pollen_counts_'

site.names <- read.table(paste(pollen.loc,'site_names.txt',sep='')) 

file.chronology <- list.files(chron.loc)
file.chronology <- file.chronology[grep(chron_prefix,file.chronology)]

lakes.with.chron <- unlist(strsplit(file.chronology,chron_prefix))[seq(2,2*length(file.chronology),2)]

file.pollen <- list.files(pollen.loc)
file.pollen <- file.pollen[grep('pollen_counts',file.pollen)]
pollen.within.domain <- unlist(strsplit(file.pollen,pollen_prefix))[seq(2,2*length(file.pollen),2)]

lakes.available <- intersect(pollen.within.domain,lakes.with.chron)
lakes.available <- lakes.available[-35]

pollen_chronology <- 
lapply(lakes.available,function(x){
#for(x in lakes.available){ 
  lake_name <- x#site.names[x,1] 
  pollen.data <- read.csv(paste(pollen.loc,pollen_prefix,lake_name,sep=''))
  chron.data <- read.csv(paste(chron.loc,chron_prefix,lake_name,sep=''))
  depths.equal <- all.equal(chron.data$depths-1,pollen.data$depth) # important: these sites have the same depths!!
  if(depths.equal != TRUE) {stop('Depths do not match')}
  
  data.frame(age= round(rowMeans(chron.data[,-1],na.rm=TRUE)),pollen.data[,-1])
}
)

names(pollen_chronology) <- lakes.available

saveRDS(pollen_chronology,paste(pollen.loc,'compiled_pollen.RDS',sep=''))


#########################################################################################################################
# load meta data
#########################################################################################################################
meta.data.sites <- readRDS(paste(pollen.loc,'meta_data_sites.RDS',sep='')) 

meta.dataset.id <- unlist(meta.data.sites['dataset.id',])

lakes.available2 <- unlist(strsplit(lakes.available,'.csv'))
lakes.available2 <- unlist(strsplit(lakes.available2,'dataset_id_'))[seq(2,2*length(lakes.available),2)] 
lakes.available2 <- as.numeric(lakes.available2)



sum(meta.dataset.id%in%lakes.available2)

meta.data.sites <- meta.data.sites[,meta.dataset.id%in%lakes.available2]
meta.data.sites <- meta.data.sites

saveRDS(meta.data.sites,paste(pollen.loc,'meta_data_sites.RDS',sep='')) 
saveRDS(lakes.available2,paste(pollen.loc,'order_dataset_ids.RDS'))
