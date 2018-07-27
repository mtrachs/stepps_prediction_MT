###################################################################################################################3
# Aggregate pollen to a certain age
###################################################################################################################
setwd('~/workflow_stepps_prediction/')
chron.loc <- 'age-depth/Bchron/output/'
pollen.loc <- 'pollen_data/'  


pollen_chronology <- readRDS(paste(pollen.loc,'compiled_pollen.RDS',sep=''))

taxa <- c('Ash','Beech','Birch','Chestnut','Hemlock','Hickory','Maple','Oak','Other.conifer','Other.hardwood','Pine',
          'Spruce','Tamarack')
K <- length(taxa)

fake.counts <- data.frame(t(rep(0,K+1)))
names(fake.counts) <- c('age',taxa) 

# we have to reorder the data
ordering_dataset_ids <- readRDS(paste(pollen.loc,'order_dataset_ids.RDS'))


aggregated.pollen <- 
lapply(order(ordering_dataset_ids),function(z){
#for(z in 1:length(pollen_chronology)){  
  ages.agg <- seq(150,1950,100)
  ages <- pollen_chronology[[z]]$age
  ages.2k.index <- ((ages>50)&(ages<2050))
  ages.2k <- ages[ages.2k.index]  
  if(length(ages.2k)>0){
    
  ages.assign <- sapply(ages.2k, function(x) ages.agg[which.min(abs(ages.agg - x))])
  
  pollen.aggregated <- aggregate(pollen_chronology[[z]][ages.2k.index,-1],list(age=ages.assign),FUN=sum)
  pollen.aggregated <- matrix(ncol=14,unlist(analogue::join(fake.counts,pollen.aggregated)$pollen))
  #pollen.aggregated <- as.data.frame(pollen.aggregated)
  colnames(pollen.aggregated) <- colnames(fake.counts)
  
  pollen.final <- matrix(ncol = 14, nrow=length(ages.agg),data=0)
  pollen.final[,1] <- ages.agg
  pollen.final[ages.agg%in%pollen.aggregated[,'age'],] <- pollen.aggregated 
  colnames(pollen.final) <- colnames(pollen.aggregated)
  pollen.final
  }
})  

names(aggregated.pollen) <- names(pollen_chronology)

saveRDS(aggregated.pollen,paste(pollen.loc,'aggregated_pollen.RDS',sep=''))


########################################################################################################################
# merge all datasets
########################################################################################################################
pollen.def <- aggregated.pollen[[1]] 

for (i in 2:length(aggregated.pollen)){
  pollen.def <- rbind(pollen.def,aggregated.pollen[[i]])  
} 

pollen.def <- pollen.def[,colnames(pollen.def)!='age']

saveRDS(pollen.def,paste(pollen.loc,'pollen_def.RDS',sep=''))


########################################################################################################################
#coordinates of sites used!!
########################################################################################################################
meta.data.sites <- readRDS(paste(pollen.loc,"meta_data_sites.RDS",sep=''))

dataset.ids <- unlist(meta.data.sites['dataset.id',])
index.use <- sapply(1:length(dataset.ids),function(x) length(grep(as.character(dataset.ids[x]),names(aggregated.pollen)))>0)

site.meta <- meta.data.sites[,index.use]

coordinates.sites <- t(site.meta[c('lon','lat'),])

#there are a few sites without data remove theses sites. 
#find datasets without data :)
index.no.data <- sapply(aggregated.pollen,function(x) !is.null(x))
coordinates.sites <- coordinates.sites[index.no.data,]

saveRDS(coordinates.sites,paste(pollen.loc,'coordinates_def.RDS',sep=''))
