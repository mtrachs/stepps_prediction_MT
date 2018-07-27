library(gridExtra)
library(maptools)

us.shp <- readShapeLines('data/map_data/us_alb.shp',
                         proj4string=CRS('+init=epsg:3175'))
us.shp@data$id <- rownames(us.shp@data)
us.fort <- fortify(us.shp, region='id') 


create_figure_path <- function(subDir){
  mainDir <- getwd()
  
  figure_path = file.path(mainDir, subDir)
  
  if (file.exists(subDir)){
    print('Folder already exists: contents may get overwritten!')
  } else {
    dir.create(figure_path)
  }
  
  print(paste('All figures will be saved in: ', figure_path, sep=''))
  
#    return(figure_path)
}

# trace plots
# fit is a stanfit object
# can pass true values as a list if known
trace_plot_pars <- function(post, N_knots, T, N_pars, taxa, suff=suff, save_plots=TRUE, fpath=subDir){
  
#   post = rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
#   post = post[250:500,,]
  n    = dim(post)[3]
  
  #don't plot knots here
  idx_pars = c(seq(1,N_pars),n)
  
  labels = colnames(post[,1,])[idx_pars]
  
  labels[2:(length(labels)-1)] =  paste(labels[2:(length(labels)-1)], taxa[1:(K-1)], sep=' ')

#   avg = summary(fit)$summary[,"mean"]
  
  par(mfrow=c(4,2))
  if (save_plots){
    pdf(paste(fpath, "/trace_pars.pdf", sep=""), width=8, height=12)
    par(mfrow=c(5,2))
  }
  
  
  for (i in idx_pars){
    plot(post[,1,i], type="l", ylab=labels[i], xlab="iter")
#     plot(post[,i], type="l", ylab=labels[i], xlab="iter")
    abline(h=mean(post[,1,i]), col="blue")
    abline(h=quantile(post[,1,i],probs=0.025), col='blue', lty=2)
    abline(h=quantile(post[,1,i],probs=0.975), col='blue', lty=2)
    #     if (dim(post)[2] >= 2){
    #       lines(post[,2,i], col="blue")
    #     }
  }
  
  if (save_plots){
    dev.off()
  }

}

# trace plots
# fit is a stanfit object
# can pass true values as a list if known
trace_plot_mut <- function(post, N_knots, T, N_pars, taxa, mean_type='other', suff=suff, save_plots=TRUE, fpath=subDir){
  
  n    = dim(post)[3]
  
  #idx_pars = seq(3+N_pars+1, 3+N_pars+W*T)
  #idx_pars = seq(N_pars+1, N_pars+W*T)
  idx_pars = seq(N_pars+1, N_pars+W*(T-1))
  
  if (mean_type=='MRF'){
    idx_pars = seq(3, 3+W*T)
  }
#   
#   labels = colnames(post[,1,])[idx_pars]
  
#   labels[4:(length(labels)-1)] =  paste(labels[4:(length(labels)-1)], taxa[1:(K-1)], sep=' ')
  
#   avg = summary(fit)$summary[,"mean"]
  
  
  if (save_plots){
    pdf(paste(fpath, "/trace_mut.pdf", sep=""), width=8, height=12)
    par(mfrow=c(5,2))
  }
  
  for (k in 1:W){
    par(mfrow=c(5,2))
    print(k)
    idx_taxon = seq(k, W*(T-1), by=W)
    #idx_taxon = seq(k, W*T, by=W)
    idx = idx_pars[idx_taxon]
    labels = colnames(post[,1,])[idx]
    draws = post[,1,idx]
    
    for (i in 1:length(idx)){
      plot(draws[,i], type="l", ylab=labels[i], xlab="iter")
      #     plot(post[,i], type="l", ylab=labels[i], xlab="iter")
      abline(h=mean(draws[,i]), col="blue")
      abline(h=quantile(draws[,i],probs=0.025), col='blue', lty=2)
      abline(h=quantile(draws[,i],probs=0.975), col='blue', lty=2)
      #     if (dim(post)[2] >= 2){
      #       lines(post[,2,i], col="blue")
      #     }
    }
  }

  if (save_plots){
    dev.off()
  }
  
}


knot_idx <- function(w, n, t,K){
  K + 1 + (n-1)*T*W + (t-1)*W + w-1
}

# trace plots of knots
# each knot on a different panels, multiple times per pane
trace_plot_knots <- function(fit, N_knots, T, K, N_pars, suff=suff, save_plots=TRUE, fpath=subDir){
  
  post = rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
  labels = colnames(post[,1,])
  niter = dim(post)[1]
  
#   t = 1
#   n = 2
#   w = 1
  
  t = seq(1,T)
  
  
#   labels[6 + (n-1)*T*W + (t-1)*W + w-1]
#   knot = post[,1,6 + (n-1)*T*W + (t-1)*W + w-1]
  
  par(mfrow=c(3,2))
  
  if (save_plots){
    pdf(paste(fpath, "/trace_knots_", suff, ".pdf", sep=""), width=8, height=4)
    par(mfrow=c(1,1))
  }
  
 
  for (w in 1:W){
    for (n in 1:N_knots){
      
      post_knots = post[,1,knot_idx(w,n,t,N_pars)]
  
      ylo = min(post_knots)
      yhi = max(post_knots)
  
      plot(c(0,0), type="n", ylab=paste('Taxon', w, '; Knot', n ,sep=' '), xlab="iter", ylim=c(ylo, yhi), xlim=c(0,niter))
      #plot(c(0,0), type="n", ylab=paste('Knot', n ,sep=' '), xlab="iter", ylim=c(0.0379, 0.038), xlim=c(0,niter))
      cols=c("red","blue", "black")
      for (i in 1:T){
        lines(post_knots[,i], col=cols[i])
      }
    }
  }
  
  if (save_plots){
    dev.off()
  }

}

# trace plots
# fit is a stanfit object
# can pass true values as a list if known
trace_plot_process <- function(r, suff=suff, save_plots, fpath=subDir){
  
  K = dim(r)[2] # n taxa
  
  cols = rep('black', K)#rainbow(K)
  
  par(mfrow=c(4,2))
  if (save_plots){
    if (nchar(suff)>0) suff1 = paste0('_', suff)
    pdf(paste(fpath, "/trace", suff1, ".pdf", sep=""), width=8, height=4)
    par(mfrow=c(1,1))
  }
  
  for (i in 1:10){
    for (k in 1:K){
      plot(r[i,k,], type="l", col=cols[k], ylab=paste(suff, i, 'taxon', k, sep=' '), 
           xlab="iter", ylim=c(min(r[i,k,]),max(r[i,k,])))
    }
  }

  
  if (save_plots){
    dev.off()
  }
}

get_limits <- function(centers){
  xlo = min(centers[,1])
  xhi = max(centers[,1])

  ylo = min(centers[,2])
  yhi = max(centers[,2])

  return(list(xlims=c(xlo,xhi),ylims=c(ylo, yhi)))
}  

plot_pred_maps <- function(r_mean, centers, taxa, t, N, K, T, thresh, limits, type, suff, save_plots, fpath=subDir){
  
  rescale=1000000
  
  if (is.null(taxa)){taxa=seq(1,K)}
  if (type=='prop'){bar_title='Proportions'}else {bar_title='Values'}
  
  prop_dat = data.frame(props=numeric(0), x=integer(0), y=integer(0),time=integer(0), taxon=character())
  for (k in 1:K){
    prop_dat = rbind(prop_dat, data.frame(props = r_mean[,k], 
                                          x     = rep(centers[,1], each=T)*rescale, 
                                          y     = rep(centers[,2], each=T)*rescale, 
                                          time  = rep(t,times=N), 
                                          taxon = rep(taxa[k], N*T)))
  }
  
  if (!is.na(thresh)){
    prop_dat$props[which(prop_dat$props > thresh)] = thresh
  }
    
  p <- ggplot() + geom_tile(data=prop_dat, aes(x=x, y=y, fill=props)) + 
    scale_fill_gradientn(colours= tim.colors(), name=bar_title) + coord_fixed() +
    scale_x_continuous(limits$xlims) + scale_y_continuous(limits$ylims)
  p <- add_map_albers(p, map_data=us.fort, limits)
  p <- p + facet_grid(time~taxon)
  p <- theme_clean(p) + theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))
  print(p)
  Sys.sleep(2)
  if (save_plots){
    ggsave(file=paste(fpath, '/veg_maps_', suff, '.pdf', sep=''), scale=1, width=12,height=2)
#     dev.off()
  }
  return(p)
}


# plot_pred_maps_select(r_mean, centers_veg, taxa=taxa, ages, N, K, T, thresh=0.5, limits, type='prop', suff=suff,  save_plots=save_plots)

plot_pred_maps_select <- function(r_mean, centers, taxa, ages, N, K, T, thresh, limits, type, suff, save_plots, fpath=subDir){
  
  if (is.null(taxa)){taxa=seq(1,K)}
  if (!is.null(taxa)){
    taxa[which(taxa == 'OTHER.HARDWOOD')] = 'OTHER\nHW'
    taxa[which(taxa == 'OTHER.CONIFER')] = 'OTHER\nCON'
    taxa[which(taxa == 'TAMARACK')] = 'LARCH'
    
#     taxa = as.vector(sapply(taxa, simpleCap))
  }
  
  if (type=='prop'){bar_title='Proportions'}else {bar_title='Values'}
  
  # find a better way to do this later
  if (length(ages) > 2){
    idx.keep  = c(1,2,length(ages)/2,T)
    ages.keep = ages[idx.keep]
    T.keep    = length(ages.keep)
  } else if (length(ages) == 1){
    idx.keep = 1
    ages.keep = ages
    T.keep = length(ages.keep)
  }
  
  idx_r_keep = vector(length=0)
  for (i in 1:T.keep){
    idx_orig   = seq(idx.keep[i], N*T, by=T)
    idx_r_keep = c(idx_r_keep, idx_orig)
  }
  
  idx_r_keep = sort(idx_r_keep)
  r_mean_keep = r_mean[idx_r_keep,]
  
  prop_dat = data.frame(props=numeric(0), x=integer(0), y=integer(0),time=integer(0), taxon=character())
  for (k in 1:K){
    prop_dat = rbind(prop_dat, data.frame(props = r_mean_keep[,k], 
                                          x     = rep(centers[,1]*1000000, each=T.keep), 
                                          y     = rep(centers[,2]*1000000, each=T.keep), 
                                          time  = rep(ages.keep, times=N), 
                                          taxon = rep(taxa[k], N*T.keep)))
  }
  
  if (!is.na(thresh)){
    prop_dat$props[which(prop_dat$props > thresh)] = thresh
  }

  p <- ggplot() + geom_tile(data=prop_dat, aes(x=x, y=y, fill=props)) + 
    scale_fill_gradientn(colours=tim.colors(), name=bar_title) + coord_fixed() #+
    #scale_x_continuous(limits$xlims*1000) + scale_y_continuous(limits$ylims*1000)
  p <- add_map_albers(p, map_data=us.fort, limits)
  p <- p + facet_grid(time~taxon)
  p <- theme_clean(p) + theme(strip.text.y = element_text(size = rel(1.3)), strip.text.x = element_text(size = rel(1.0)))
#   print(p)
  Sys.sleep(2)
  if (save_plots){
    fname = paste0(fpath, '/veg_maps_', suff, '.pdf')
    ggsave(file=fname, scale=1, width=12, height=12)
    sys_str = paste("pdfcrop", fname, fname, sep=' ')
    system(sys_str)
    #     dev.off()
  }

  return(p)
}

# plot_pred_maps_binned_select(r_mean, centers_veg, breaks, taxa, ages, N, K, T, limits, suff=suff_figs, save_plots, fpath=subDir)

plot_pred_maps_binned_select <- function(r_mean, centers, breaks, taxa, ages, N, K, T, limits, suff, save_plots, fpath=subDir){
  
  if (is.null(taxa)){taxa=seq(1,K)}
  if (!is.null(taxa)){
    taxa[which(taxa == 'OTHER.HARDWOOD')] = 'OTHER\nHW'
    taxa[which(taxa == 'OTHER.CONIFER')] = 'OTHER\nCON'
    taxa[which(taxa == 'TAMARACK')] = 'LARCH'
    
    #     taxa = as.vector(sapply(taxa, simpleCap))
  }
  
  if (length(ages) > 2){
    idx.keep  = c(1,2,length(ages)/2,T)
    ages.keep = ages[idx.keep]
    T.keep    = length(ages.keep)
  } else if (length(ages) == 1){
    idx.keep  = 1
    ages.keep = ages
    T.keep    = length(ages.keep)
  }
  
  idx_r_keep = vector(length=0)
  for (i in 1:T.keep){
    idx_orig   = seq(idx.keep[i], N*T, by=T)
    idx_r_keep = c(idx_r_keep, idx_orig)
  }
  
  idx_r_keep = sort(idx_r_keep)
  r_mean_keep = r_mean[idx_r_keep,]

  r_mean_binned = matrix(0, nrow=nrow(r_mean_keep), ncol=ncol(r_mean_keep))
  colnames(r_mean_binned) <- colnames(r_mean_keep)
  
  for (i in 1:ncol(r_mean)){
    r_mean_binned[,i] = cut(r_mean_keep[,i], breaks, include.lowest=TRUE, labels=FALSE)
  }
  
  breaklabels = apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1, 
                      function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })
  
  prop_dat = data.frame(props=numeric(0), x=integer(0), y=integer(0),time=integer(0), taxon=character())
  for (k in 1:K){
    prop_dat = rbind(prop_dat, data.frame(props = r_mean_binned[,k], 
                                          x     = rep(centers[,1]*1000000, each=T.keep), 
                                          y     = rep(centers[,2]*1000000, each=T.keep), 
                                          time  = rep(ages.keep, times=N), 
                                          taxon = rep(taxa[k], N*T.keep)))
  }
  
  p <- ggplot() + geom_tile(data=prop_dat, aes(x=x, y=y, fill=factor(props))) + 
    scale_fill_manual(values = tim.colors(length(breaks)), labels=breaklabels, name='Proportions') + 
    coord_fixed() #+ scale_x_continuous(limits$xlims*1000) + scale_y_continuous(limits$ylims*1000)
  p <- add_map_albers(p, map_data=us.fort, limits)
  p <- p + facet_grid(time~taxon)
#   p <- theme_clean(p) + theme(strip.text.x = element_blank(), strip.text.y = element_text(size = rel(1.5)))
  p <- theme_clean(p) + theme(strip.text.y = element_text(size = rel(1.3)), strip.text.x = element_text(size = rel(1.0)))
  #theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))
#   print(p)
  Sys.sleep(2)
  if (save_plots){
    fname = paste0(fpath, '/veg_maps_props_binned_', suff, '.pdf')	
    ggsave(file=fname, scale=1, width=12, height=12)
    sys_str = paste("pdfcrop", fname, fname, sep=' ')
    system(sys_str)
  }
  return(p)
}

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), tolower(substring(s, 2)),
        sep="", collapse=" ")
}

plot_data_maps <- function(y, centers, taxa, t, N, K, T, thresh, limits, suff, save_plots, fpath=subDir){
  
  if (is.null(taxa)){taxa=seq(1,K)}
  
  props_data = t(apply(y, 1, function(x) x/sum(x)))
  #colnames(props_data) = taxa
  
  prop_dat = data.frame(props=numeric(0), x=integer(0), y=integer(0), taxon=character())
  for (k in 1:K){
    prop_dat = rbind(prop_dat, data.frame(props = props_data[,k], 
                                          x     = centers[,1], 
                                          y     = centers[,2], 
                                          taxon = rep(taxa[k], N)))
  }
  
  if (!is.na(thresh)){
    prop_dat$props[which(prop_dat$props > thresh)] = thresh
  }
    
  p <- ggplot() + geom_tile(data=prop_dat, aes(x=x*1000, y=y*1000, fill=props)) + 
    scale_fill_gradientn(colours=tim.colors(), guide='none') + coord_fixed() + 
    scale_x_continuous(limits$xlims*1000) + scale_y_continuous(limits$ylims*1000)
  p <- add_map_albers(p, map_data=us.fort, limits)
  p <- p + facet_grid(~taxon)
  p <- theme_clean(p) + theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))
  print(p)
  Sys.sleep(2)
  if (save_plots){
    fname = paste0(fpath, '/veg_maps_data_', suff, '.pdf')
    ggsave(file=fname, scale=1, width=12)
    sys_str = paste("pdfcrop", fname, fname, sep=' ')
    system(sys_str)
#     dev.off()
  }
   return(p)
}

# plot_both_maps <- function(r_mean, y, centers, taxa, t, N, K, T, thresh, suff, save_plots, fpath=subDir){
#   
#   if (is.null(taxa)){taxa=seq(1,K)}
#   
#   props_data = t(apply(y, 1, function(x) x/sum(x)))
#   #colnames(props_data) = taxa
#   
#   prop_dat = data.frame(props=numeric(0), x=integer(0), y=integer(0), taxon=character())
#   for (k in 1:K){
#     prop_dat = rbind(prop_dat, data.frame(props = props_data[,k], 
#                                           x     = centers[,1], 
#                                           y     = centers[,2], 
#                                           taxon = rep(taxa[k], N)))
#   }
#   
#   if (!is.na(thresh)){
#     prop_dat$props[which(prop_dat$props > thresh)] = thresh
#   }
#   
#   p1 <- ggplot() + geom_raster(data=prop_dat, aes(x=x, y=y, fill=props)) + 
#     scale_fill_gradientn(colours=tim.colors()) + coord_fixed()
#   p1 <- p1 + facet_grid(~taxon)
#   p1 <- theme_clean(p1) + theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))
# #   print(p)
#   
#   
#   prop_preds = data.frame(props=numeric(0), x=integer(0), y=integer(0),time=integer(0), taxon=character())
#   for (k in 1:K){
#     prop_preds = rbind(prop_preds, data.frame(props = r_mean[,k], 
#                                           x     = rep(centers[,1], each=T), 
#                                           y     = rep(centers[,2], each=T), 
#                                           time  = rep(t,times=N), 
#                                           taxon = rep(taxa[k], N*T)))
#   }
#   
#   if (!is.na(thresh)){
#     prop_preds$props[which(prop_preds$props > thresh)] = thresh
#   }
#   
#   p2 <- ggplot() + geom_raster(data=prop_preds, aes(x=x, y=y, fill=props)) + 
#     scale_fill_gradientn(colours=tim.colors()) + coord_fixed()
#   p2 <- p2 + facet_grid(time~taxon)
#   p2 <- theme_clean(p2) + theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))
# #   print(p2)
#   
#   grid.arrange(p1, p2, ncol=1)
#   
#   
# #   Sys.sleep(1)
# #   if (save_plots){
# #     ggsave(file=paste('figures/pred_model/veg_maps_', suff, '.pdf', sep=''), scale=1)
# #   }
#   #   print(p)
# }


plot_pred_maps_binned <- function(r_mean, centers, breaks, taxa, ages, N, K, T, limits, suff, save_plots, fpath=subDir){
  
  rescale=1000000
  
  if (is.null(taxa)){taxa=seq(1,K)}
  
  r_mean_binned = matrix(0, nrow=nrow(r_mean), ncol=ncol(r_mean))
  colnames(r_mean_binned) <- colnames(r_mean)
  
  for (i in 1:ncol(r_mean)){
    r_mean_binned[,i] = cut(r_mean[,i], breaks, include.lowest=TRUE, labels=FALSE)
  }
  
  breaklabels = apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1, 
                      function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })
  
  prop_dat = data.frame(props=numeric(0), x=integer(0), y=integer(0),time=integer(0), taxon=character())
  for (k in 1:K){
    prop_dat = rbind(prop_dat, data.frame(props = r_mean_binned[,k], 
                                          x     = rep(centers[,1], each=T)*rescale, 
                                          y     = rep(centers[,2], each=T)*rescale, 
                                          time  = rep(ages,times=N), 
                                          taxon = rep(taxa[k], N*T)))
  }
  
  p <- ggplot() + geom_tile(data=prop_dat, aes(x=x, y=y, fill=factor(props))) + 
    scale_fill_manual(values = tim.colors(length(breaks)), labels=breaklabels, name='Proportions') + 
    coord_fixed() + scale_x_continuous(limits$xlims) + scale_y_continuous(limits$ylims)
  p <- add_map_albers(p, map_data=us.fort, limits)
  p <- p + facet_grid(time~taxon)
  p <- theme_clean(p) + theme(strip.text.x = element_blank(), strip.text.y = element_text(size = rel(1.5)))
    #theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))
  print(p)
  Sys.sleep(2)
  if (save_plots){
    ggsave(file=paste(fpath, '/veg_maps_binned_', suff, '.pdf', sep=''), scale=1)
#     dev.off()
  }
  return(p)
}

theme_clean <- function(plot_obj){
  plot_obj <- plot_obj + theme(axis.ticks = element_blank(), 
                               axis.text.y = element_blank(), 
                               axis.text.x = element_blank(),
                               axis.title.x = element_blank(),
                               axis.title.y = element_blank(),
                               plot.background = element_rect(fill = "transparent",colour = NA))
  
  return(plot_obj)
}


plot_core_locations <- function(y, centers_polU, centers_pls, ages, limits, fpath){
  
  centers = centers_polU*1000
  
  K = ncol(y)
  N_cores = nrow(centers_polU)
  
  buffer = 0#100000
  
  # xlo = min(centers[,1]) - buffer
  # xhi = max(centers[,1]) + buffer
  # ylo = min(centers[,2]) - buffer
  # yhi = max(centers[,2]) + buffer
  
  xlo = min(centers_pls[,1])*1000 - buffer
  xhi = max(centers_pls[,1])*1000 + buffer
  ylo = min(centers_pls[,2])*1000 - buffer
  yhi = max(centers_pls[,2])*1000 + buffer
  
  core_dat = data.frame(x=integer(0), y=integer(0), age=character(), core_sites=integer(0))
  for (i in 1:length(ages)){
    print(i)
    print((N_cores*(i-1) + 1))
    
    y_sub    = y[(N_cores*(i-1) + 1):(N_cores*i),]
    idx_data = which(rowSums(y_sub) != 0) 
    
    print(idx_data)
    
    core_dat = rbind(core_dat, data.frame(x     = centers_polU[idx_data,1]*1000, 
                                          y     = centers_polU[idx_data,2]*1000, 
                                          age   = rep(ages[i], length(idx_data))
    ))
  }
  
  
  p <- ggplot() + geom_point(data=core_dat, aes(x=x, y=y), colour='#FF6600', shape=19) + 
    coord_fixed() #+ scale_x_continuous(limits$xlims) + scale_y_continuous(limits$ylims)
  p <- add_map_albers(p, map_data=us.fort, limits)#, xlims=c(xlo,xhi), ylims=c(ylo, yhi))
  p <- p + facet_grid(age~.)
  p <- p + theme(strip.text.x = element_blank(),
                 strip.text.y = element_blank())
  p <- p + theme(strip.background = element_blank())
  p <- theme_clean(p) #+ theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))
  
  Sys.sleep(2)
  print(p)
  
  if (save_plots){
    ggsave(file=paste(fpath, '/core_locs_', suff, '.pdf', sep=''), scale=1)
    #     dev.off()
  }
  
  return(p)
}


plot_core_locations_select <- function(y, centers_pol, centers_pls, ages, idx.keep, limits, suff, fpath){
 
  rescale = 1000000
  centers = centers_pol*rescale

  K = ncol(y)
  N_cores = nrow(centers_pol)
  
  buffer = 0#100000

  # xlo = min(centers[,1]) - buffer
  # xhi = max(centers[,1]) + buffer
  # ylo = min(centers[,2]) - buffer
  # yhi = max(centers[,2]) + buffer

  xlo = min(centers_pls[,1])*rescale - buffer
  xhi = max(centers_pls[,1])*rescale + buffer
  ylo = min(centers_pls[,2])*rescale - buffer
  yhi = max(centers_pls[,2])*rescale + buffer
  
#   idx.keep  = c(1,2,length(ages)/2,T)
  ages.keep = ages[idx.keep]
  T.keep    = length(ages.keep)
  
  idx_y = vector(length=0)
  for (i in 1:T.keep){
    idx_orig   = seq(idx.keep[i], N_cores*T, by=T)
    idx_y = c(idx_y, idx_orig)
  }
  
  idx_y = sort(idx_y)
  y_keep = y[idx_y,]
  
  core_dat = data.frame(x=integer(0), y=integer(0), age=character(), core_sites=integer(0))
  for (i in 1:length(ages.keep)){
    print(i)
    print((N_cores*(i-1) + 1))
    
    y_sub    = y_keep[(N_cores*(i-1) + 1):(N_cores*i),]
    idx_data = which(rowSums(y_sub) != 0) 
    
    print(idx_data)
    
    core_dat = rbind(core_dat, data.frame(x     = centers_pol[idx_data,1]*rescale, 
                                          y     = centers_pol[idx_data,2]*rescale, 
                                          age   = rep(ages.keep[i], length(idx_data))
                                          ))
  }
  
  
  p <- ggplot() + geom_point(data=core_dat, aes(x=x, y=y), colour='#FF6600', shape=19) + 
    coord_fixed() #+ scale_x_continuous(limits$xlims) + scale_y_continuous(limits$ylims)
  p <- add_map_albers(p, map_data=us.fort, limits)#, xlims=c(xlo,xhi), ylims=c(ylo, yhi))
  p <- p + facet_grid(age~.)
#   p <- p + theme(strip.text.x = element_blank(),
#                 strip.text.y = element_blank())
#   p <- p + theme(strip.background = element_blank())
  p <- theme_clean(p) + theme(strip.text.y = element_text(size = rel(1.5)), 
                          strip.text.x = element_text(size = rel(1.5)))
  
  Sys.sleep(2)
  print(p)
  
  if (save_plots){
    fname = paste0(fpath, '/core_locs_', suff, '.pdf')	
    ggsave(file=fname, scale=1)
    sys_str = paste("pdfcrop", fname, fname, sep=' ')
    system(sys_str)
  }
  
  return(p)
}




add_map_albers <- function(plot_obj, map_data=us.fort, limits){
  p <- plot_obj + geom_path(data=map_data, aes(x=long, y=lat, group=group),  colour='grey55') + 
    #     scale_x_continuous(limits = c(min(umw.coord$x, na.rm=TRUE), max(umw.coord$x, na.rm=TRUE))) +
    #     scale_y_continuous(limits = c(min(umw.coord$y, na.rm=TRUE), max(umw.coord$y, na.rm=TRUE)))#, colour = "black", size = 1, fill = "white", aes(x=long, y=lat, group = group))
    # #   
#     scale_x_continuous(limits = c(min(dat[,1], na.rm=TRUE), max(dat[,1], na.rm=TRUE))) +
#     scale_y_continuous(limits = c(min(dat[,2], na.rm=TRUE), max(dat[,2], na.rm=TRUE)))
   scale_x_continuous(limits = limits$xlims*1000000) +
   scale_y_continuous(limits = limits$ylims*1000000) #+ coord_map("albers")
  return(p)
  
}




poster_fig <- function(y, y_veg, r_mean, centers_veg, centers_pls, centers_polU, taxa, t, N, K, T, thresh, limits, type, suff, save_plots, fpath=subDir){
  
  if (is.null(taxa)){taxa=seq(1,K)}
  
  props_data = t(apply(y_veg, 1, function(x) x/sum(x)))
  
  K = ncol(y)
  N_cores = nrow(centers_polU)
  
  prop_dat = data.frame(props=numeric(0), x=integer(0), y=integer(0), time=character(0), taxon=character())
#   core_dat = data.frame(x=integer(0), y=integer(0), age=character(), core_sites=integer(0))
  for (i in 1:length(ages)){
    print(i)
    print((N_cores*(i-1) + 1))
    
    y_sub    = y[(N_cores*(i-1) + 1):(N_cores*i),]
    idx_data = which(rowSums(y_sub) != 0) 
    
    print(idx_data)
    
    prop_dat = rbind(prop_dat, data.frame(props = rep(0, length(idx_data)),
                                          x     = centers_polU[idx_data,1]*1000, 
                                          y     = centers_polU[idx_data,2]*1000, 
                                          time  = rep(ages[i], length(idx_data)),
                                          taxon = rep('cores', length(idx_data))))

  }
  
  core_dat = data.frame(x=integer(0), y=integer(0), time=character(0), taxon=character())
  #   core_dat = data.frame(x=integer(0), y=integer(0), age=character(), core_sites=integer(0))
  for (i in 1:length(ages)){
    print(i)
    print((N_cores*(i-1) + 1))
    
    y_sub    = y[(N_cores*(i-1) + 1):(N_cores*i),]
    idx_data = which(rowSums(y_sub) != 0) 
    
    print(idx_data)
    
    core_dat = rbind(core_dat, data.frame(
                                          x     = centers_polU[idx_data,1]*1000, 
                                          y     = centers_polU[idx_data,2]*1000, 
                                          time  = rep(ages[i], length(idx_data)),
                                          taxon = rep('cores', length(idx_data))))
    
  }
  #colnames(props_data) = taxa
  
#   prop_dat = data.frame(props=numeric(0), x=integer(0), y=integer(0), time=character(0), taxon=character())
  for (k in 1:K){
    prop_dat = rbind(prop_dat, data.frame(props = props_data[,k], 
                                          x     = centers_pls[,1]*1000, 
                                          y     = centers_pls[,2]*1000,
                                          time  = rep('pls', times = nrow(centers_pls)),
                                          taxon = rep(taxa[k], nrow(centers_pls))))
  }
  

  
#   prop_dat = data.frame(props=numeric(0), x=integer(0), y=integer(0),time=integer(0), taxon=character())
  for (k in 1:K){
    prop_dat = rbind(prop_dat, data.frame(props = r_mean[,k], 
                                          x     = rep(centers_veg[,1]*1000, each=T), 
                                          y     = rep(centers_veg[,2]*1000, each=T), 
                                          time  = rep(as.character(ages),times=N), 
                                          taxon = rep(taxa[k], N*T)))
  }
  
  if (!is.na(thresh)){
    prop_dat$props[which(prop_dat$props > thresh)] = thresh
  }
  
  p <- ggplot() + geom_tile(data=prop_dat, aes(x=x, y=y, fill=props)) + 
    scale_fill_gradientn(colours=tim.colors()) + coord_fixed() +
     scale_x_continuous(limits$xlims*1000) + scale_y_continuous(limits$ylims*1000)
  p <- add_map_albers(p, us.shp, limits)
  p <- p + facet_grid(time~taxon)
  p <- theme_clean(p) + theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))
  print(p)
  
  p <- p + geom_point(data = core_dat, aes(x=x, y=y))
  
  print(p)
  #ggsave(file='figures/pred/pred_plot_test.pdf', scale=1)
  Sys.sleep(2)
  if (save_plots){
    ggsave(file=paste(fpath, '/veg_maps_', suff, '.pdf', sep=''), scale=1)
    #     dev.off()
  }
  return(p)
}



# plot(centers_veg[,1]*1000, centers_veg[,2]*1000, col='blue', pch=19)
# points(centers_pls[,1]*1000, centers_pls[,2]*1000)
# 
# us.shp <- readShapeLines('r/data/map_data/us_alb.shp',
#                          proj4string=CRS('+init=epsg:3175'))
# plot(us.shp, add=T, lwd=2)


plot_data_maps <- function(y, centers, taxa, t, N, K, T, thresh, limits, suff, save_plots, fpath=subDir){
  
  rescale=1000000
  
  if (is.null(taxa)){taxa=seq(1,K)}
  
  props_data = t(apply(y, 1, function(x) x/sum(x)))
  #colnames(props_data) = taxa
  
  prop_dat = data.frame(props=numeric(0), x=integer(0), y=integer(0), taxon=character())
  for (k in 1:K){
    prop_dat = rbind(prop_dat, data.frame(props = props_data[,k], 
                                          x     = centers[,1]*rescale, 
                                          y     = centers[,2]*rescale, 
                                          taxon = rep(taxa[k], N)))
  }
  
  if (!is.na(thresh)){
    prop_dat$props[which(prop_dat$props > thresh)] = thresh
  }
  
  prop_dat$type = rep('PLS', nrow(prop_dat))
  
  p <- ggplot() + geom_tile(data=prop_dat, aes(x=x, y=y, fill=props)) + 
    scale_fill_gradientn(colours=tim.colors(), guide='none') + coord_fixed() + 
    scale_x_continuous(limits$xlims) + scale_y_continuous(limits$ylims)
  p <- add_map_albers(p, us.shp, limits)
  p <- p + facet_grid(type~taxon)
  p <- theme_clean(p) #+ theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))
  
  p <- p + theme(strip.text.x = element_blank(),
                 strip.text.y = element_blank())
  p <- p + theme(strip.background = element_blank())
  
  print(p)
  Sys.sleep(2)
  if (save_plots){
    ggsave(file=paste(fpath, '/veg_maps_', suff, '.pdf', sep=''), scale=1)
    ggsave(file=paste(fpath, '/veg_maps_', suff, '.eps', sep=''), scale=1)
    #     dev.off()
  }
  return(p)
}

plot_data_maps_binned <- function(y, centers, taxa, t, N, K, T, breaks, limits, suff, save_plots, fpath=subDir){
  
  rescale=1000000
  
  if (is.null(taxa)){taxa=seq(1,K)}
  
  props_data = t(apply(y, 1, function(x) if (sum(x) > 0) {x/sum(x)} else {x}))
  #colnames(props_data) = taxa
  
  props_data_binned = matrix(0, nrow=nrow(props_data), ncol=ncol(props_data))
  colnames(props_data_binned) <- colnames(props_data)
  
  for (i in 1:ncol(props_data)){
    props_data_binned[,i] = cut(props_data[,i], breaks, include.lowest=TRUE, labels=FALSE)
  }
  
  breaklabels = apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1, 
                      function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })
    
  prop_dat = data.frame(props=numeric(0), x=integer(0), y=integer(0), taxon=character())
  for (k in 1:K){
    prop_dat = rbind(prop_dat, data.frame(props = props_data_binned[,k], 
                                          x     = centers[,1]*rescale, 
                                          y     = centers[,2]*rescale, 
                                          taxon = rep(taxa[k], N)))
  }
  
  prop_dat$type = rep('PLS', nrow(prop_dat))
  
  p <- ggplot() + geom_tile(data=prop_dat, aes(x=x, y=y, fill=factor(props))) + 
    scale_fill_manual(values = tim.colors(length(breaks)), labels=breaklabels, name='Proportions') + 
    coord_fixed() + 
    scale_x_continuous(limits$xlims) + scale_y_continuous(limits$ylims)
  p <- add_map_albers(p, us.shp, limits)
  p <- p + facet_grid(type~taxon)
  p <- theme_clean(p) #+ theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))
  
  p <- p + theme(strip.text.x = element_blank(),
                 strip.text.y = element_blank())
  p <- p + theme(strip.background = element_blank())
  
  print(p)
  Sys.sleep(2)
  if (save_plots){
    fname = paste0(fpath, '/veg_maps_data_binned_', suff, '.pdf')	
    ggsave(file=fname, scale=1, width=12)
    sys_str = paste("pdfcrop", fname, fname, sep=' ')
    system(sys_str)
  }
  return(p)
}