#########################################################################################################################
#
#########################################################################################################################
library(rstan)
library(RColorBrewer)
setwd('~/workflow_stepps_prediction/prediction/')
help.fun.loc <- 'utils/'
data.loc <- 'data/'
plot.loc <- 'plots/'
#########################################################################################################################
source(paste(help.fun.loc,'pred_helper_funs.r',sep=''))
source(paste(help.fun.loc,'pred_plot_funs.r',sep=''))


fit <- read_stan_csv(paste(data.loc,'test.csv',sep=''))
post <- rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
var_names<- colnames(post[,1,])
par_names <- sapply(var_names, function(x) unlist(strsplit(x,'[[]'))[1])
post_dat <- list(post = post,par_names = par_names)

r <- build_r_nb(post_dat=post_dat,N = 754,T=1,K=13)
r_mean <- apply(r$r,c(1,2),mean)

load(paste(data.loc,'prediction_test.rdata',sep=''))
limits <- list(xlims=range(coord.agg.final$east),ylims=range(coord.agg.final$north))
colours <- rev(brewer.pal(10,'RdYlBu'))


plot_pred_maps(r_mean=r_mean,
               centers = coord.agg.final,
               taxa = taxa,
               t = 1800,
               N = N,
               K = K,
               T = T,
               thresh = c(0,0.01,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,1),
               limits = limits,
               type = 'prop',
               suff = 'test',
               save_plots =TRUE ,
               fpath = plot.loc)
