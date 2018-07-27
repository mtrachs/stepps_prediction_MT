# compute effective sample size form a stanfit object
ess <- function(fit){
  ess = summary(fit)$summary[,"n_eff"]
  return(ess)
}

# compute acceptance rate from a stanfit object
ar <- function(fit){
  post = extract(fit, permuted=FALSE, inc_warmup=FALSE)
  ar = apply(post, "chains", FUN = function(x) nrow(unique(as.data.frame(x)))) / nrow(post) # acceptance rates
  return(ar)
}

# Little function to calculate posterior variances from simulation
colVars <- function (a){
  diff <- a - matrix (colMeans(a), nrow(a), ncol(a), byrow=TRUE)
  vars <- colMeans (diff^2)*nrow(a)/(nrow(a)-1)
  return (vars)
}

# The calculation of Waic!  Returns lppd, p_waic_1, p_waic_2, and waic, which we define
# as 2*(lppd - p_waic_2), as recommmended in BDA
waic <- function (stanfit){
  log_lik <- extract (stanfit, "log_lik")$log_lik
  lppd <- sum (log (colMeans(exp(log_lik))))
  p_waic_1 <- 2*sum (log(colMeans(exp(log_lik))) - colMeans(log_lik))
  p_waic_2 <- sum (colVars(log_lik)) # more stable
  waic_2 <- -2*lppd + 2*p_waic_2
  return (list (waic=waic_2, p_waic=p_waic_2, lppd=lppd, p_waic_1=p_waic_1))
}

# log likelihood vector; averaged over iterations
log_lik <- function(fit){
  log_lik <- extract(fit, "log_lik")$log_lik
  log_lik_iter <- rowSums(log_lik)
  log_lik_mean <- log(colMeans(exp(log_lik)))
  return(list(log_lik=data.frame(log_lik=log_lik_mean), log_lik_iter=data.frame(log_lik_iter),sum_log_lik=sum(log_lik_mean)))
}


aic <- function(fit, npars){
  log_lik <- extract(fit, "log_lik")$log_lik
  log_lik_mean <- log(colMeans(exp(log_lik)))
  aic = 2 * npars - 2 * sum(log_lik_mean)
  return(aic)
}

# dic <- function(fit, npars){
#   log_lik <- extract(stanfit, "log_lik")$log_lik
#   d_bar <- log(colMeans(exp(log_lik)))
#   d_theta_bar <- 
#   aic = 2 * npars - 2 * sum(log_lik_mean)
#   return(aic)
# }


get_phi_stats <- function(phi, taxa){
  phi_bar <- colMeans(phi)
  phi_L   <- apply(phi, 2, function(x){quantile(x, probs = (0.025))})
  phi_U   <- apply(phi, 2, function(x){quantile(x, probs = (0.975))})
  
  phi_stats <- data.frame(mean = phi_bar, L = phi_L, U = phi_U, taxa = taxa)
  
  phi_stats
}

# build the total potential neighborhood weighting
build_w_pot <- function(post, K, N_pot, d_pot, run){
  
  kernel = run$kernel
  
  sum_w = matrix(nrow=nrow(d_pot),rep(NA, K*nrow(d_pot)))
  
  col_names = colnames(post[,1,])
  par_names  = unlist(lapply(col_names, function(x) strsplit(x, "\\[")[[1]][1]))
  
  if (kernel=='gaussian'){
    one_psi = run$one_psi
    if (one_psi){
      psi   = rep(mean(post[,1,which(par_names == 'psi')]), K)
    } else {
      psi   = colMeans(post[,1,which(par_names == 'psi')])
    }
    
    for (k in 1:K)
      sum_w[,k] = (d_pot[,2] * exp(-d_pot[,1]^2/psi[k]^2))
    
  } else if (kernel=='pl'){
    one_a = run$one_a
    if (one_a){
      a = rep(median(post[,1,which(par_names == 'a')]), K)
    } else {
      a = apply(post[,1,which(par_names == 'a')],2,median)
    }
    
    one_b = run$one_b
    if (one_b){
      b = rep(mean(post[,1,which(par_names == 'b')]), K)
    } else {
      b = colMeans(post[,1,which(par_names == 'b')])
    }
    for (k in 1:K)
      sum_w[,k] = (d_pot[,2] * (b[k]-1) * (b[k]-2) / (2 * pi * a[k]  * a[k]) * (1 + d_pot[,1] / a[k])^(-b[k]) )
  }
  return(sum_w)
}

build_sumw_pot <- function(post, K, N_pot, d_pot, run){
  
  kernel = run$kernel
  
  sum_w = rep(NA, K)
  
  col_names = colnames(post[,1,])
  par_names  = unlist(lapply(col_names, function(x) strsplit(x, "\\[")[[1]][1]))
  
  if (kernel=='gaussian'){
    one_psi = run$one_psi
    if (one_psi){
      psi   = rep(median(post[,1,which(par_names == 'psi')]), K)
    } else {
      psi   = apply(post[,1,which(par_names == 'psi')],2,median)
    }
    
    for (k in 1:K)
      sum_w[k] = sum(d_pot[,2] * exp(-d_pot[,1]^2/psi[k]^2))
    
  } else if (kernel=='pl'){
    one_a = run$one_a
    if (one_a){
      a = rep(median(post[,1,which(par_names == 'a')]), K)
    } else {
      a = apply(post[,1,which(par_names == 'a')],2,median)
    }
    
    one_b = run$one_b
    if (one_b){
      b = rep(mean(post[,1,which(par_names == 'b')]), K)
    } else {
      b = colMeans(post[,1,which(par_names == 'b')])
    }
    for (k in 1:K)
      sum_w[k] = sum( d_pot[,2] * (b[k]-1) * (b[k]-2) / (2 * pi * a[k]  * a[k]) * (1 + d_pot[,1] / a[k])^(-b[k]) )
  }
  return(sum_w)
}







sum_hood_props <- function(post, C, N_pot, d_pot, kernel){
  
  d_hood = d_pot[which((d_pot[,1]*1e6 <= 8000) & (d_pot[,1]*1e6 > 1e-10)),]
  
  col_names = colnames(post[,1,])
  par_names  = unlist(lapply(col_names, function(x) strsplit(x, "\\[")[[1]][1]))
 
  gamma   = mean(post[,1,which(par_names == 'gamma')])
  
  if (kernel=='gaussian'){
    if (one_psi){
      psi   = mean(post[,1,which(par_names == 'psi')])
    } else {
      psi   = colMeans(post[,1,which(par_names == 'psi')])
    }
    
    w_hood = d_hood[2] * exp(-d_hood[1]^2/psi^2) 
  } else if (kernel=='pl'){
    a   = mean(post[,1,which(par_names == 'a')])
    b   = mean(post[,1,which(par_names == 'b')])
    w_hood = d_hood[2] * (b-1) * (b-2) / (2 * pi * a  * a) * (1 + d_hood[1] / a)^(-b)
  }
  
  prop_hood = (1 - gamma) * w_hood / C
  
  return(prop_hood)
  
}


power_law <- function(d, a, b) {
  x = (b-1) * (b-2) / (2 * pi * a  * a) * (1 + d / a) ^ (-b)
  return(x)
} 

#predicted pollen based on weighted neighborhoods using estimated pars
pollen_preds <- function(post, N_cores, d, idx_cores, r, sum_w, run){
  
  kernel    = run$kernel
  one_gamma = run$one_gamma
  
  iters   = dim(post)[1]
  K       = ncol(r)
  N_cells = nrow(d)
  
  col_names = colnames(post[,1,])
  par_names  = unlist(lapply(col_names, function(x) strsplit(x, "\\[")[[1]][1]))
   
  phi   = colMeans(post[,1,which(par_names == 'phi')])

  r_new = matrix(NA, nrow=N_cores, ncol=K)
  preds = matrix(NA, nrow=N_cores, ncol=K)
  
  if (one_gamma){
    gamma = rep(mean(post[,1,which(par_names == 'gamma')]), K)
  } else {
    gamma = colMeans(post[,1,which(par_names == 'gamma')])
  }
  
  if (kernel=='gaussian'){
    one_psi = run$one_psi
    if (one_psi){
      psi   = rep(mean(post[,1,which(par_names == 'psi')]), K)
    } else {
      psi   = colMeans(post[,1,which(par_names == 'psi')])
    }
  } else if (kernel=='pl'){
    print("Kernel type : (inverse) power law")
    one_a = run$one_a
    if (one_a){
      a = rep(mean(post[,1,which(par_names == 'a')]), K)
    } else {
      a = colMeans(post[,1,which(par_names == 'a')])
    }
    
    one_b = run$one_b
    if (one_b){
      b = rep(mean(post[,1,which(par_names == 'b')]), K)
    } else {
      b = colMeans(post[,1,which(par_names == 'b')])
    }
#     w = (b-1) * (b-2) / (2 * pi * a  * a) * (1 + d / a) ^ (-b)
  }
  
  for (i in 1:N_cores){
    print(i)
    
    out_sum = rep(0, K)
      for (k in 1:K){
        print(paste0('k = ', k))
        
        if (kernel == 'gaussian'){
          inv_psi2 = 1.0 / psi[k] / psi[k]
        } else if (kernel == 'pl') {
          pl_p1 = (b[k]-1) * (b[k]-2) / (2 * pi * a[k]  * a[k]) 
        }
        
        
        for (j in 1:N_cells){ # changed N_hood to N_locs
          if (j != idx_cores[i]){
            if (kernel == 'gaussian'){
              w = exp( - ( d[j,i]*d[j,i] ) * inv_psi2 )
              out_sum[k] <- out_sum[k] + w * r[j,k]
            } else if (kernel == 'pl'){
              w = pl_p1 * (1 + d[j,i] / a[k]) ^ (-b[k])
              out_sum[k] <- out_sum[k] + w * r[j,k]
            }
          }  
        }
      }
      
    for (k in 1:K){
      r_new[i,k]  = gamma[k]*r[idx_cores[i],k] + (1-gamma[k])*out_sum[k]/sum_w[k]
      preds[i,k]  = phi[k]*r_new[i,k]        
    }
  }
  
  alpha = rowSums(preds)   
  
  #convert to proportions
  preds = t(apply(preds, 1, function(x) x/sum(x)))
  
  return(list(preds=preds, alpha=alpha, sum_w=sum_w))
}

#predicted pollen based on weighted neighborhoods using estimated pars
pollen_preds_distance <- function(post, N_cores, d, idx_cores, r, C, radius){
  
  rescale = 1e6
  iters   = dim(post)[1]
  K       = ncol(r)
  N_cells = nrow(d)
  
  col_names = colnames(post[,1,])
  par_names  = unlist(lapply(col_names, function(x) strsplit(x, "\\[")[[1]][1]))
  
  phi   = colMeans(post[,1,which(par_names == 'phi')])
  gamma = mean(post[,1,which(par_names == 'gamma')])
  if (one_psi){
    psi   = mean(post[,1,which(par_names == 'psi')])
  } else {
    psi   = colMeans(post[,1,which(par_names == 'psi')])
  }
  
  r_local = matrix(NA, nrow=N_cores, ncol=K)
  r_nl = matrix(NA, nrow=N_cores, ncol=K)
  r_int = matrix(NA, nrow=N_cores, ncol=K)
  r_new = matrix(NA, nrow=N_cores, ncol=K)
  preds = matrix(NA, nrow=N_cores, ncol=K)
  preds_tot = matrix(NA, nrow=length(radius), ncol=N_cores)
  preds_int = matrix(NA, nrow=length(radius), ncol=N_cores)
  preds_dist = matrix(NA, nrow=length(radius), ncol=N_cores)
  preds_dist_taxon = array(NA, c(K, N_cores, length(radius)))
  #   w     = matrix(0, nrow=N_cores, ncol=N_cells)
  
  w = exp(-(d*d)/(psi*psi))
  
  for (rad in 1:length(radius)){
    print(paste0('rad = ', rad))
    for (i in 1:N_cores){
      idx_int = which(d[,i]<radius[rad]/rescale)
      
      #       for (j in 1:N_hood){
      #         if (idx_hood[i,j] > 0){
      #           w[i,j] <- exp(-(d[idx_hood[i,j], i])^2/(psi)^2)
      #         } 
      #       }
      
      out_sum = rep(0, K)
      out_sum_int = rep(0, K)
      sum_w <- 0
      
      for (j in 1:N_cells){ # changed N_hood to N_locs
        if (j != idx_cores[i]){
          #         if (idx_hood[i,j] > 0){
          #           out_sum <- out_sum + w[i,j]*r[idx_hood[i,j],]
          out_sum <- out_sum + w[j,i]
          
          if (j %in% idx_int){
            out_sum_int <- out_sum_int + w[j,i]
          }
        }  
      }
      
      sum_w   <- sum(out_sum)
      
      sum_w_int   <- sum(out_sum_int)
      
      #       r_local[i,] = r[idx_cores[i],]
      r_nl[i,]  = gamma + (1-gamma) / C * out_sum
      r_int[i,] = gamma + (1-gamma) / C * out_sum_int
      #       r_new[i,]  = gamma*r[idx_cores[i],] + (1-gamma)*out_sum/sum_w
      #       preds[i,] = phi*r_new[i,]  
      
    }
    
    # do i want to multiply by phi??
    #     preds_loc = rowSums(r_local)
    preds_tot[rad,]  = r_nl[,1]#rowSums(r_nl)
    preds_int[rad,] = r_int[,1]#rowSums(r_int)
    preds_dist[rad,] = preds_int[rad,]/preds_tot[rad,] 
    preds_dist_taxon[,,rad] = r_int/r_nl 
    
  }
  
  #   convert to proportions
  #   preds = t(apply(preds, 1, function(x) x/sum(x)))
  
  return(list(preds_tot=preds_tot, preds_int=preds_int, preds_dist=preds_dist))
}

#predicted pollen based on weighted neighborhoods using estimated pars
# dispersal_decay(post, d_pot, sum_w, radius/rescale, run, taxa)
dispersal_decay <- function(post, dmat, sum_w, radius, run, taxa){
  
  kernel = run$kernel
  one_gamma = run$one_gamma
  if (kernel == 'gaussian'){
    one_psi = run$one_psi
  } else if (kernel == 'pl'){
    one_a = run$one_a
    one_b = run$one_b
  }
  
  rescale = 1e6
  iters   = dim(post)[1]
  K       = ncol(r)
  N_cells = nrow(d)
  
  w = array(NA, c(K, nrow(dmat)))
  
  col_names = colnames(post[,1,])
  par_names  = unlist(lapply(col_names, function(x) strsplit(x, "\\[")[[1]][1]))
  
  phi   = colMeans(post[,1,which(par_names == 'phi')])
  if (one_gamma){
    gamma = rep(mean(post[,1,which(par_names == 'gamma')]), K)
  } else {
    gamma = colMeans(post[,1,which(par_names == 'gamma')])
  }
  
  if (kernel=='gaussian'){
    if (one_psi){
      psi   = rep(mean(post[,1,which(par_names == 'psi')]), K)
    } else {
      psi   = colMeans(post[,1,which(par_names == 'psi')])
    }
    
    for (k in 1:K){
      w[k,] = dmat[,2] * exp(-(dmat[,1] * dmat[,1])/(psi[k] * psi[k]))
    }
  } else if (kernel=='pl'){
    if (one_a){
      a = rep(mean(post[,1,which(par_names == 'a')]), K)
    } else {
      a = colMeans(post[,1,which(par_names == 'a')])
    }
    if (one_b){
      b = rep(mean(post[,1,which(par_names == 'b')]), K)
    } else {
      b = colMeans(post[,1,which(par_names == 'b')])
    }
    for (k in 1:K){
        w[k,] = dmat[,2] * (b[k]-1) * (b[k]-2) / (2 * pi * a[k]  * a[k]) * (1 + dmat[,1] / a[k]) ^ (-b[k])
    }
  }
  
  r_int  = array(NA, c(length(radius), K) ) #vector(length=length(radius), mode='numeric')
#   C_test = sum(w)
  colnames(r_int) = taxa
  
  for (rad in 1:length(radius)){
    for (k in 1:K){
      idx_int    = which(dmat[,1]<radius[rad])#/rescale)
      sum_w_int  = sum(w[k, idx_int])
      r_int[rad, k] = gamma[k] + (1-gamma[k]) / sum_w[k] * sum_w_int     
    }
  }
  
  return(r_int)
}

#build_sumw_pot_ci(post, K, N_pot, d_pot, run)
# build the total potential neighborhood weighting
build_sumw_pot_ci <- function(post, K, N_pot, d_pot, run){
  
  kernel = run$kernel
  
  niter = dim(post)[1]
  sum_w = array(NA, c(K, niter))
  
  col_names = colnames(post[,1,])
  par_names = unlist(lapply(col_names, function(x) strsplit(x, "\\[")[[1]][1]))
  
  if (kernel=='gaussian'){
    one_psi = run$one_psi
    if (one_psi){
      psi   = replicate(K, post[,1,which(par_names == 'psi')])
    } else {
      psi   = post[,1,which(par_names == 'psi')]
    }
    
    for (i in 1:niter)
      for (k in 1:K)
        sum_w[k,i] = sum(d_pot[,2] * exp(-d_pot[,1]^2/psi[i,k]^2))
    
  } else if (kernel=='pl'){
    one_a = run$one_a
    if (one_a){
      a = replicate(K, post[,1,which(par_names == 'a')])
    } else {
      a = post[,1,which(par_names == 'a')]
    }
    
    one_b = run$one_b
    if (one_b){
      b = replicate(K, post[,1,which(par_names == 'b')])
    } else {
      b = post[,1,which(par_names == 'b')]
    }
    for (i in 1:niter)
      for (k in 1:K)
        sum_w[k,i] = sum( d_pot[,2] * (b[i,k]-1) * (b[i,k]-2) / (2 * pi * a[i,k]  * a[i,k]) * 
                            (1 + d_pot[,1] / a[i,k])^(-b[i,k]) )
  }
  return(sum_w)
}

#predicted pollen based on weighted neighborhoods using estimated pars
# dispersal_decay(post, d_pot, sum_w, radius/rescale, run, taxa)
dispersal_decay_ci <- function(post, dmat, sum_w, radius, run, taxa){
  
  kernel = run$kernel
  one_gamma = run$one_gamma
  if (kernel == 'gaussian'){
    one_psi = run$one_psi
  } else if (kernel == 'pl'){
    one_a = run$one_a
    one_b = run$one_b
  }
  
  rescale = 1e6
  niter   = dim(post)[1]
  K       = ncol(r)
  N_cells = nrow(d)
  
  w = array(NA, c(K, nrow(dmat), niter))
  
  col_names = colnames(post[,1,])
  par_names  = unlist(lapply(col_names, function(x) strsplit(x, "\\[")[[1]][1]))
  
  phi   = post[,1,which(par_names == 'phi')]
  if (one_gamma){
    gamma = replicate(K, post[,1,which(par_names == 'gamma')])
  } else {
    gamma = post[,1,which(par_names == 'gamma')]
  }
  
  if (kernel=='gaussian'){
    if (one_psi){
      psi   = replicate(K, post[,1,which(par_names == 'psi')])
    } else {
      psi   = post[,1,which(par_names == 'psi')]
    }
    for (i in 1:niter){
      for (k in 1:K){
        w[k,,i] = dmat[,2] * exp(-(dmat[,1] * dmat[,1])/(psi[i,k] * psi[i,k]))
      }
    }
  } else if (kernel=='pl'){
    if (one_a){
      a = replicate(K, post[,1,which(par_names == 'a')])
    } else {
      a = post[,1,which(par_names == 'a')]
    }
    if (one_b){
      b = replicate(K, post[,1,which(par_names == 'b')])
    } else {
      b = post[,1,which(par_names == 'b')]
    }
    for (i in 1:niter){
      for (k in 1:K){
        w[k,,i] = dmat[,2] * (b[i,k]-1) * (b[i,k]-2) / (2 * pi * a[i,k]  * a[i,k]) * (1 + dmat[,1] / a[i,k]) ^ (-b[i,k])
      }
    }
  }
  
  r_int  = array(NA, c(length(radius), K, niter) ) #vector(length=length(radius), mode='numeric')
  #   C_test = sum(w)
  colnames(r_int) = taxa
  
  for (i in 1:niter){
    for (rad in 1:length(radius)){
      for (k in 1:K){
        idx_int    = which(dmat[,1]<radius[rad])#/rescale)
        sum_w_int  = sum(w[k, idx_int, i])
        r_int[rad,k,i] = gamma[i,k] + (1-gamma[i,k]) / sum_w[k,i] * sum_w_int     
      }
    }
  }
  return(r_int)
}


#scale veg props in focal cell based on phi
phi_scale_veg <- function(post, N_cores, r, idx_cores){
  
  iters = dim(post)[1]
  K     = dim(r)[2]
  
  phi   = summary(fit)$summary[,'mean'][1:K]
  
  r_new   = matrix(NA, nrow=N_cores, ncol=K)
  phi_veg = matrix(NA, nrow=N_cores, ncol=K)
  
  for (i in 1:N_cores){
    for (k in 1:K){
      
      r_new[i,k] = phi[k]*r[idx_cores[i],k]
      
    }
  }
  
  preds = r_new/apply(r_new, 1, sum) 
  
  return(preds)
  
}

get_quants <- function(fit, npars){
  quants <- cbind(summary(fit)$summary[,'mean'][1:npars],
                  summary(fit)$summary[,'2.5%'][1:(npars)],
                  summary(fit)$summary[,'50%'][1:(npars)],
                  summary(fit)$summary[,'97.5%'][1:(npars)],
                  summary(fit)$summary[,'n_eff'][1:(npars)],
                  summary(fit)$summary[,'Rhat'][1:(npars)])
  colnames(quants) = c('mean', '2.5%', '50%', '97.5%', 'n_eff', 'Rhat')
  
  return(quants)
}

compute_props <- function(y, taxa){
  pollen_props <- y/rowSums(y) 
  colnames(pollen_props) <- taxa
  
  return(pollen_props)
}


gaussian <- function(d, psi) {
  x = exp( -( d / psi ) ^ 2 )
  return(x)
} 