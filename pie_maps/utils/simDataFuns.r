#
# Simulate data to test the PLS model
# Implemented covariance functions include the NEGEXP and MATERN
# Choose to implement on Wisconsin (or subset) or on a grid
#

library(mvtnorm)
library(stats)
library(fields)
# library(msm)

# makes a square grid with cells that are 8000 m wide
# keep wisconsin x and y lower limits to keep distance magnitudes 
make_grid <- function(xlimits=c(252000,708000), ylimits=c(676000,1172000), x_ncells=10, y_ncells=10){
  
  xlo = xlimits[1]
  #xhi = xlimits[2]
  ylo = ylimits[1]
  #yhi = ylimits[2]
  
  dx = 8 #(xhi-xlo)/side 
  dy = 8 #(yhi-ylo)/side  
  
  #knots = matrix(nrow=0,ncol=2)
  cells = matrix(nrow=0,ncol=2)
  for (i in 1:x_ncells){
    for (j in 1:y_ncells){
      x = xlo + i*dx
      y = ylo + j*dy
      
      cells= rbind(cells, c(x,y))
#       d   = rdist(matrix(c(x,y), nrow=1), as.matrix(coords))
#       #print(min(d))
#       idx = which.min(d)
#       if (d[idx] <= sqrt(4000^2 + 4000^2)){
#         dnew = rdist(matrix(coords[idx,], nrow=1), as.matrix(coords))
#         if (length(which(dnew<=sqrt(8000^2 + 8000^2))) == 9 ){
#           knots = rbind(knots, coords[idx,])           
#         }
#       }            
    }
  }
  return(cells)
}

#make a subgrid of cells
make_subgrid <- function(cells, n_knots=36){
  
  xlo = min(cells[,1])
  xhi = max(cells[,1])
  ylo = min(cells[,2])
  yhi = max(cells[,2])
  
  Dx = 8000
  Dy = 8000
  
  s = max(floor(sqrt(((yhi-ylo)/Dy+1)*((xhi-xlo)/Dx+1)/n_knots)),2)
  
  dy = s * Dy
  dx = s * Dx
  
  knots = matrix(nrow=0,ncol=2)
  
  Nx = ceiling((xhi - xlo) / dx)
  Ny = ceiling((yhi - ylo) / dy)
  for (i in 1:Nx) {
    x = xlo + i*dx
    
    for (j in 1:Ny) {
      y = ylo + j*dy
      
      d   = rdist(matrix(c(x,y), nrow=1), as.matrix(cells))
      #print(min(d))
      idx = which.min(d)
      if (d[idx] <= sqrt(4000^2 + 4000^2)){
        dnew = rdist(matrix(cells[idx,], nrow=1), as.matrix(cells))
        if (length(which(dnew<=sqrt(8000^2 + 8000^2))) == 9 ){
          knots = rbind(knots, cells[idx,])           
        }
      }
    }
  }
  return(knots)
}

#make a subgrid of cells
make_regular_subgrid <- function(cells, n_knots=36, edges=TRUE){
  
  xlo = min(cells[,1])
  xhi = max(cells[,1])
  ylo = min(cells[,2])
  yhi = max(cells[,2])
  
  if (edges){
    dx = (xhi-xlo)/ceiling(sqrt(n_knots))
    dy = (yhi-ylo)/ceiling(sqrt(n_knots))
  }
    
  knots = matrix(nrow=0,ncol=2)
  colnames(knots) = c("x", "y")
  
  Nx = ceiling((xhi - xlo) / dx)
  Ny = ceiling((yhi - ylo) / dy)
  for (i in 0:Nx) {
    x = xlo + i*dx
    
    for (j in 0:Ny) {
      y = ylo + j*dy
      
      knots = rbind(knots, c(x,y))           
    }
  }
  return(knots)
}

# Simulate data to test the PLS model
# Implemented covariance functions include the NEGEXP and MATERN 
# K = number of species
# pass pars as a list, where each list element contains K-1 elements; 
#   for NEGEXP expect pars to contain: eta, rho, sigma_sq, mu
#   for MATERN expect rho, sigma_sq, mu
sim_data <- function(cells=cells, knots=knots, COV='NEGEXP', pars=pars, K=3, write_dat=TRUE){
  
  W = K-1               #working number of species
  N = nrow(cells)
  N_knots = nrow(knots)
  
  d       = rdist(as.matrix(cells), as.matrix(cells))
  d_knots = rdist(as.matrix(knots), as.matrix(knots))
  d_inter = rdist(as.matrix(cells), as.matrix(knots))
  
  # pull out pars from list
  if (COV=='NEGEXP'){
    eta=pars$eta
  }
  
  rho=pars$rho
  sigma_sq=pars$sigma_sq
  mu=pars$mu
    
  # construct covariance matrix
  Sigma = array(dim=c(W, N, N))
  for (k in 1:W) {
    for (i in 1:N) {
      for (j in i:N) {
        if (COV=='NEGEXP'){
          Sigma[k,i,j] <- eta[k] * exp(-1/rho[k] * d[i,j])
          if (i==j) Sigma[k,i,j] <- Sigma[k,i,j] + sigma_sq[k]
        }
        if (COV=='MATERN'){
          Sigma[k,i,j] <- sigma_sq[k] * (1+sqrt(3)/rho[k]*d[i,j])*exp(-sqrt(3)/rho[k] * d[i,j])
        }
      }
    }
  
    for (i in 1:(N-1)) {
      for (j in (i+1):N) {
        Sigma[k,j,i] <- Sigma[k,i,j];  
      }
    }
  } 

# for (i in 1:K){
#   print(all(eigen(Sigma[i,,])$values >0 ) )
# }
  
  # gaussian process
  g = matrix(nrow=W, ncol=N)
  for (k in 1:W){
    g[k, ] = rmvnorm(1, mu[k]*rep(1, N), Sigma[k,,], method="chol")
  }
  
  # transform process vals to proportions
  r = matrix(nrow=K, ncol=N)
  for (k in 1:W) {
    for (i in 1:N) {
      r[k,i] <- exp(g[k,i]) / (1 + colSums(exp(g))[i])
    }
  }
  
  for (i in 1:N) {
    r[K,i] <- 1 / (1 + colSums(exp(g))[i])
  }
  
  #cell_count=runif(N, min=74, max=282)
  cell_count = rep(100, N) 
  y = matrix(nrow=N, ncol=K)
  for (i in 1:N){
    y[i,] = rmultinom(1, size=cell_count[i], r[,i])
  }
  
  if (write_dat==TRUE){
    fit_data = list(K = K, N=N, N_knots=N_knots, y=y, d=d, d_knots=d_knots, d_inter=d_inter, r=r)
    unname(fit_data)
    dump("fit_data", file=paste("r/dump/fit_data_", tolower(COV),".dump", sep=""))
    
#     eta = eta[1:2]
#     rho = rho[1:2]
#     sigma_sq=sigma_sq[1:2]
#     mu=mu[1:2]
#     g=g[1:2,]
    
    inits = list()
    if (COV=='NEGEXP'){inits$eta=eta}
    inits = c(inits, list(rho=rho, sigma_sq=sigma_sq, mu=mu, g=g))
    
    unname(inits)
    #dump(c("eta", "rho_sq", "sigma_sq", "mu", "g"), file='r/dump/true_inits.dump')
    dump('inits', file=paste("r/dump/inits_", tolower(COV), ".dump", sep=""))
    
#     rout=r[1:2,]
    true = list()
    if (COV=='NEGEXP'){true$eta=eta}
    true = c(true, list(rho=rho, sigma_sq=sigma_sq, mu=mu, g=g, r=r))
    dump('true', file=paste("r/dump/true_", tolower(COV), ".dump", sep=""))
  }
  
  return(list(K=K, N=N, y=y, d=d, d_knots=d_knots, d_inter=d_inter, r=r))
  
}



# # Create data to help find parameters that generate different relative proportions
# # No spatial pattern here though!
# make_data <- function(cells=cells, knots=knots, K=3, suff, write_dat=TRUE){
#   
#   N = nrow(cells)
#   N_knots = nrow(knots)
#   
#   d       = rdist(as.matrix(cells), as.matrix(cells))
#   d_knots = rdist(as.matrix(knots), as.matrix(knots))
#   d_inter = rdist(as.matrix(cells), as.matrix(knots))
#   
#   y = t(rmultinom(N, 100, c(1/3, 1/3, 1/3)))
#   
#   r = t(prop.table(y,1))
#   
#   if (write_dat==TRUE){
#     rig_data = list(K = K, N=N, N_knots=N_knots, y=y, d=d, d_knots=d_knots, d_inter=d_inter, r=r)
#     unname(rig_data)
#     dump("rig_data", file=paste("r/dump/make_data", suff, ".dump",sep=""))
#   }
#   
#   return(list(K=K, N=N, y=y, d=d, d_knots=d_knots, d_inter=d_inter, r=r))
#   
# }


# make some spatial patterns with gaussian bumps!
# need to pass list: npars = list(sigma, nbumps)
# 

gaussian_data <- function(cells=cells, npars, K=3){

  N = nrow(cells)

  xlo = min(cells[,1])
  xhi = max(cells[,1])
  ylo = min(cells[,2])
  yhi = max(cells[,2])
  
  sigma  = npars$sigma
  nbumps = npars$nbumps
  r      = matrix(nrow=K, ncol=N) 
  
  for (k in 1:K){
    
    props = rep(0,N)
    for (j in 1:nbumps[k]){
     x0 = runif(1, min=xlo, max=xhi)
     y0 = runif(1, min=ylo, max=yhi)
     center = matrix(c(x0, y0), nrow=1)
     print(center)
     
     myd   = rdist(center, cells)
     props = props + exp(-(myd)**2/(sigma[k]*0.7)**2)
     #props = props + exp(-(myd)**2/sigma[k]**2)
    }
    
    # a little jitter to cells with low vals
    nlows = length(which(props< median(props)))
    props[props < median(props)] = runif(nlows, min=0, max=1.3*median(props))

    r[k,]=props
    
  }
  
  # scale to get cells proportions
  r = prop.table(r, 2)
  
  return(r)
}    
  
# take simulates proportions and write them to a data file stan can use
props_to_data <- function(r, cells=cells, knots=knots, K=3, suff, write_dat=TRUE){
  
  N = nrow(cells)
  N_knots = nrow(knots)
  
  d       = rdist(as.matrix(cells), as.matrix(cells))
  d_knots = rdist(as.matrix(knots), as.matrix(knots))
  d_inter = rdist(as.matrix(cells), as.matrix(knots))
  
  y = matrix(nrow=N, ncol=K)
  
  for (i in 1:N){  
    y[i,] = rmultinom(1, 100, r[,i])
  }
  
  if (write_dat==TRUE){
    test_data = list(K = K, N=N, N_knots=N_knots, y=y, d=d, d_knots=d_knots, d_inter=d_inter, r=r)
    unname(test_data)
    dump("test_data", file=paste("r/dump/test_data_", suff, ".dump",sep=""))
    dump(c('K', 'N', 'N_knots', 'y', 'd', 'd_knots', 'd_inter', 'r'), file=paste("r/dump/test_dataT_", suff, ".dump",sep=""))
  }
  
  return(list(K=K, N=N, N_knots=N_knots, y=y, d=d, d_knots=d_knots, d_inter=d_inter, r=r))
  
}

knots_in_domain4 <-function(knots, cells, cell_width = 8000){
  
  knots_int = matrix(nrow=0, ncol=2)
  
  for (i in 1:nrow(knots)){
    x = knots[i,1]
    y = knots[i,2]
    
    xright = x + cell_width
    xleft  = x - cell_width
    yup    = y + cell_width
    ydown  = y - cell_width
    
    right = matrix(cbind(xright, y), nrow=1)
    left = matrix(cbind(xleft, y), nrow=1)
    up = matrix(cbind(x, yup), nrow=1)
    down = matrix(cbind(x, ydown), nrow=1)
    
    d = rdist(as.matrix(cells), right)
    nright = length(which(d < cell_width))
    
    d = rdist(as.matrix(cells), left)
    nleft = length(which(d < cell_width))
    
    d = rdist(as.matrix(cells), up)
    nup = length(which(d < cell_width))
    
    d = rdist(as.matrix(cells), down)
    ndown = length(which(d < cell_width))
    
    if ((nup > 0) & (ndown > 0) & (nright > 0) & (nleft > 0)) {
      knots_int = rbind(knots_int, knots[i,])
    }
    
  }
  
  knots_int
}


# take simulates proportions and write them to a data file stan can use
dump_knots_data <- function(dat, cells=cells, knots=knots, suff){
  
  N = unname(nrow(cells))
  N_knots = unname(nrow(knots))
  K = unname(dat$K)
  y = unname(dat$y)
  
  d       = rdist(as.matrix(cells), as.matrix(cells))
  d_knots = rdist(as.matrix(knots), as.matrix(knots))
  d_inter = rdist(as.matrix(cells), as.matrix(knots))
  
  test_data = list(K = dat$K, N=N, N_knots=N_knots, y=dat$y, d=d, d_knots=d_knots, d_inter=d_inter)
  unname(test_data)
  #dump("test_data", file=paste("r/dump/fit_data_", suff, ".dump",sep=""))
  dump(c('K', 'N', 'N_knots', 'y', 'd', 'd_knots', 'd_inter'), 
       file=paste("r/dump/fit_data_", N_knots, 'knots_', N, 'cells_', K, 'taxa', suff, ".dump",sep=""))
}


# make a list of initial conditions
dump_inits <- function(dat, cells=cells, knots=knots, suff){
  
  N = unname(nrow(cells))
  N_knots = unname(nrow(knots))
  K = unname(dat$K)
  y = unname(dat$y)
  
  eta = rep(1, K-1)
  rho = rep(100000, K-1)
  sigma_sq = rep(0.1, K-1) 
  mu = rep(0, K-1)
  
  alpha = matrix(0, nrow=(K-1), ncol = N_knots)
  
  inits = list(eta = eta, rho = rho, sigma_sq = sigma_sq, mu = mu, alpha = alpha)
  unname(inits)
  #dump("test_data", file=paste("r/dump/fit_data_", suff, ".dump",sep=""))
  dump(c('eta', 'rho', 'sigma_sq', 'mu', 'alpha'), 
       file=paste("r/dump/inits_", suff, ".dump",sep=""))
}

# knots_in_domain <- function(knots){
# 
#   #change projections
#   knots_df = as.data.frame(knots)
#   coordinates(knots_df) <- ~ x + y   # or is it x and y?
#   proj4string(knots_df) <- CRS('+init=epsg:3175')
#   knots_p <- spTransform(knots_df, CRS('+proj=longlat +ellps=WGS84'))
#   coordsLL = data.frame(coordinates(knots_p))
# 
#   #get wisconsin sites
#   idx = map.where(database = "state", coordsLL$x, coordsLL$y)
#   wi.idx = which(idx == "wisconsin")
#   in_state = knots[wi.idx,]
# 
#   return(in_state)
# }
# 
# knots_int <- knots_in_domain(knots)
# 
# 
# 
# knots_in_domain2 <-function(knots, cells){
#   
#   knots_int = matrix(nrow=0, ncol=2)
#   
#   for (i in 1:nrow(knots)){
#     d = rdist(as.matrix(cells), matrix(knots[i,], nrow=1))
#     
#     nbs = length(which(d <= sqrt(2*8000^2)))
#     
#     print(nbs)
#     
#     if (nbs > 4){
#       knots_int = rbind(knots_int, knots[i,])
#     }
#     
#   }
#   
#   return(knots_int)
# }
# 
# 
# 
# 
# 
# knots_in_domain3 <-function(knots, cells){
#   
#   knots_int = matrix(nrow=0, ncol=2)
#   
#   for (i in 1:nrow(knots)){
#     x = knots[i,1]
#     y = knots[i,2]
#     
#     nbs = which((abs(cells$x - x) <= 8000) & (abs(cells$y - y) <= 8000))
#     print(nbs)
#     
#     if (length(nbs) > 1){
#       knots_int = rbind(knots_int, knots[i,])
#     }
#     
#   }
#   
#   knots_int
# }

#   y = t(rmultinom(N, 100, c(1/3, 1/3, 1/3)))
# 
#   r = t(prop.table(y,1))
# 
#   if (write_dat==TRUE){
#     rig_data = list(K = K, N=N, N_knots=N_knots, y=y, d=d, d_knots=d_knots, d_inter=d_inter, r=r)
#     unname(rig_data)
#     dump("rig_data", file=paste("r/dump/make_data", suff, ".dump",sep=""))
#   }
# 
#   return(list(K=K, N=N, y=y, d=d, d_knots=d_knots, d_inter=d_inter, r=r))
