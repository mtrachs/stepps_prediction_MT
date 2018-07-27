library(rstan)
library(RColorBrewer)

taxa <- sort(c('Maple','Birch','Hickory','Chestnut','Other conifer','Other hardwood',
               'Beech','Ash','Tamarack','Spruce','Pine','Poplar','Oak','Hemlock','Elm'))

#fit <- read_stan_csv('~/cmdstan-2.6.2/stepps_pred_settlement1.csv')
#fit <- read_stan_csv('~/cmdstan-2.6.2/stepps_pred_settlement_test_wrong_distance_matrix.csv')#is indeed pretty bad...
#fit <- read_stan_csv('~/cmdstan-2.6.2/stepps_pred_settlement_pl.csv')
fit <- read_stan_csv('~/cmdstan-2.6.2/stepps_pred_settlement_pl_old_weight_matrix.csv')
post <- rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
dim(post)
var_names <- colnames(post[,1,])
par_names <- sapply(var_names, function(x) unlist(strsplit(x,'[[]'))[1])
gauss.process <- grep('g',par_names)
gauss.process <- gauss.process[-1]

#find uncertainty


