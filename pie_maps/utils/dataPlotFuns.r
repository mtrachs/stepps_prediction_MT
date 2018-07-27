# plot settlement era time series data given a site id
se_time_series <- function(id, pollen_se){
  
  rows   = which(pollen_se[,'id'] == id)
  ages   = pollen_se[rows,'ages']
  counts = pollen_se[rows, 9:ncol(pollen_se)]
  props  = t(apply(counts, 1, function(x) x/sum(x)))
  
  sitename = pollen_se[rows[1],'sitename']
  
  if (length(ages) < 2){
    print('Less than 2 ages in bracket!')
    return()
  }
  
  xlo = min(ages) 
  xhi = max(ages)
  
  ylo = min(counts)
  yhi = max(counts)
  
  ylo = min(props)
  yhi = max(props)
  
  #pdf(paste('r/data/figs/se_time_series_id', id, '.pdf', sep=''), width=8, height=6)
  plot(1, type="n", xlim=c(-xhi, -xlo + 0.3*(xhi-xlo)), ylim=c(ylo, yhi), xlab='Age (ybp)', ylab='Proportion', 
       main=paste('Site:', sitename, '; ID:', id, sep=''), xaxt='n')
  at=axis(1, labels=FALSE)
  axis(1, at=at, labels=sapply(abs(at), toString))
  
  #fix this??? ntaxa
  col_list = c("#1F78B4", "#33A02C", "#FF7F00", "#E31A1C", "#6A3D9A", "#B15928", 
               "#FFFF99", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F", "#CAB2D6",
               "magenta", "cyan", "goldenrod", "black")
  #ntaxa = ncol(counts)
  #   cols = rainbow(ntaxa)
  n=16
  for (j in 1:(n-1)){
    #lines(ages, counts[,j], col=cols[j])
    lines(-ages, props[,j], col=col_list[j], lwd=3)#cols[j]) 
  }
  lines(-ages, rowSums(props[,n:ncol(props)]), col=col_list[n], lwd=3)
  
  legend('bottomright', c(colnames(counts)[1:(n-1)], 'rest'), pch=rep(19,n), pt.cex=1.2, col=col_list[1:n], 
         bg='white', title='Taxa groups')
  #dev.off()
}

# this plots observations at locations (x) with color varying with magnitude (y)
pmap2 <-
  function (y, x, xlab = "", ylab = "", col = tim.colors(32), zlim = NULL, 
            shrink = 0.9, cex = 0.7, pch = 19, ...) 
  {
    if (is.null(zlim)) {
      zlim = range(y, na.rm = T) + c(-1e-15, 1e-15)
    }
    out = as.numeric(cut(y, breaks = seq(zlim[1], zlim[2], length = length(col) + 
                                           1)))
    if (shrink < 1) {
      tmp = par()$plt
      tmp[2] = shrink * tmp[2]
      par(plt = tmp)
    }
    plot(x[, 1], x[, 2], col = col[out], xlab = xlab, ylab = ylab, 
         pch = pch, cex = cex, ...)
    if (shrink < 1) {
      image.plot(y, x = x, add = T, legend.only = T, zlim = zlim, 
                 col = col)
    }
  }

# plots pie charts for taxa composition on a map; note this is specialized to the STEPPS1 New England domain in the input args defaults, but should work for a different domain
# proportions should be a n x p matrix with columns being different taxa
# centers should be a n x 2 matrix of locations
pieMap=function(proportions, 
                centers, 
                restrict = FALSE,
                inputRestricted = FALSE,
                xlim = c(-52000,940000),
                ylim = c(676000, 1484000),
                radius = NULL,
                scale = 1,
                xlab = 'x',
                ylab = 'y', 
                add_legend, 
                main_title,
                col_list){#c("#1F78B4", "#33A02C", "#E31A1C", "#6A3D9A", "#B15928", "#FF7F00", 
                             #"#FFFF99", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F", "#CAB2D6",'cyan','grey','white')){
  # plots multiple pie composition charts as a map
  if(!restrict){
    used=1:nrow(centers)
  } else{ # assumes subset of grid
    if(inputRestricted){
      used=1:144
    } else{
      used=49:(256-64)
    }
  }
  centers=as.matrix(centers[used,])
  proportions=as.matrix(proportions[used,])
  if(is.null(xlim)){
    rg=range(centers[,1])
    df=(scale-1)*diff(range(centers[,1]))
    xlim=c(rg[1]-df,rg[2]+df)
  }
  if(is.null(ylim)){
    rg=range(centers[,2])
    df=(scale-1)*diff(range(centers[,2]))
    ylim=c(rg[1]-df,rg[2]+df)
  }
  plot(centers,type='n',xlim=xlim,ylim=ylim,xaxt='n', yaxt='n',ann=FALSE, frame.plot=F, asp=1, main="Pie plots")#, xlab=xlab,ylab=ylab)
  us.shp <- readShapeLines('~/workflow_stepps_calibration/calibration/data/map_data/us_alb.shp',
                           proj4string=CRS('+init=epsg:3175'))
  plot(us.shp, add=T, lwd=2)
  n=length(centers[,1])

#   col_list = c("#1F78B4", "#33A02C", "#E31A1C", "#6A3D9A", "#B15928", "#FF7F00", 
#                "#FFFF99", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F", "#CAB2D6")
  cols     = col_list[1:ncol(proportions)]
#   if (ncol(proportions)<13){
#     cols     = col_list[1:ncol(proportions)]
#   } else {
#     cols = rainbow(ncol(proportions))
#   }
  #if (ncol(proportions) > 12){print('More than 12 categories, aggregate to distinguish!')}
  
  if(is.null(radius)){
    radius=.025*diff(range(centers[,1]))
  }
  minVal=min(proportions)
  # eps=1e-10*max(proportions)
  # if(minVal==0){
  #   warning("Some proportions are zero; proceeding with jittering.")
  #   proportions=proportions+eps
  # }
  if(minVal<0){
    stop("Some proportions are less than zero.")
  }
  if(length(radius)==1){ radius=rep(radius,n)}
  for(i in 1:n){
    if(sum(proportions[i,])>0){
      minVal=min(proportions[i,])
      if(minVal==0){
        warning("Some proportions are zero; proceeding with jittering.")
        eps=1e-10*max(proportions[i,])
        proportions[i,]=proportions[i,]+eps
      }
      pieAdd(as.vector(proportions[i,]),as.vector(centers[i,]),radius=radius[i], col=cols)
    } else{
      #points(centers[i,],pch='X')
    }
  }
  #   par(xpd=NA) 
  #   tmp <- cnvrt.coords(.9,.7, 'tdev')$usr 
  if (add_legend){
    legend.col=c(0,1,1,1,0,1,1,1) 
    #     900000,1500000
    legend('topleft', colnames(proportions), pch=rep(22,n), pt.cex=1.6, cex=1.2, pt.bg=cols, col=rep('black', n),
           bg='white', ncol=2)
    title(main=main_title, cex.main=2)
  }
}

# auxiliary fxn needed by pieMap
pieAdd= function (x, center, labels = names(x), edges = 200, radius = 0.8, density = NULL, 
                  angle = 45, col = NULL, border = NULL, lty = NULL,
                  col_list)# = c("#1F78B4", "#33A02C", "#E31A1C", "#6A3D9A", "#B15928", "#FF7F00", 
                           #    "#FFFF99", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F", "#CAB2D6")) # modified from the pie() function in R
{
  if (!is.numeric(x) || any(is.na(x) | x <= 0)) 
    stop("pie: `x' values must be positive.")
  if (is.null(labels)) 
    labels <- rep("",length(x))
  x <- c(0, cumsum(x)/sum(x))
  dx <- diff(x)
  
  pin <- par("pin")
  nx <- length(dx)
  if (is.null(col)){ 
    col <- if (is.null(density)) 
#     col_list = c("#1F78B4", "#33A02C", "#E31A1C", "#6A3D9A", "#B15928", "#FF7F00", 
#                  "#FFFF99", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F", "#CAB2D6")
    cols     = col_list[1:nx]
#     if (ncol(proportions)<13){
#       cols     = col_list[1:nx]
#     } else {
#       cols = rainbow(ncol(proportions))
#     }
    
    #topo.colors(nx)
    #       cols = brewer.pal(12, 'Paired')
    #       cols = c("white", "black","lightblue", "darkblue", "red","yellow",
    #        "purple","orange","lightgreen","darkgreen", "burlywood", 
    #         "hotpink", "indianred4")
    #       cols = c("#1F78B4", "#33A02C", "#FF7F00", "#E31A1C", "#6A3D9A", "#B15928", 
    #                "#FFFF99", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F", "#CAB2D6")
    #       cols[1:nx]
    #       if (nx > 12){print('More than 12 categories, aggregate to distinguish!')}
  }else par("fg")
  col <- rep(col, length.out = nx)
  border <- rep(border, length.out = nx)
  lty <- rep(lty, length.out = nx)
  angle <- rep(angle, length.out = nx)
  density <- rep(density, length.out = nx)
  for (i in 1:nx) {
    n <- max(2, floor(edges * dx[i]))
    t2p <- 2 * pi * seq(x[i], x[i + 1], length = n)
    xc <- c(cos(t2p), 0) * radius + center[1]
    yc <- c(sin(t2p), 0) * radius + center[2]
    polygon(xc, yc, density = density[i], angle = angle[i], 
            border = border[i], col = col[i], lty = lty[i])
    t2p <- 2 * pi * mean(x[i + 0:1])
    xc <- cos(t2p) * radius + center[1]
    yc <- sin(t2p) * radius + center[2]
    if (!is.na(lab <- labels[i]) && lab != "") {
      lines(c(1, 1.05) * xc, c(1, 1.05) * yc)
      text(1.1 * xc, 1.1 * yc, lab, xpd = TRUE, adj = ifelse(xc < 
                                                               0, 1, 0))
    }
  }
  invisible(NULL)
}


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


#make a subgrid of cells
regular_subgrid <- function(cells, dx, dy){
  
  xlo = min(cells[,1])
  xhi = max(cells[,1])
  ylo = min(cells[,2])
  yhi = max(cells[,2])
  
  knots = matrix(nrow=0,ncol=2)
  colnames(knots) = c("x", "y")
  
  Nx = ceiling((xhi - xlo) / dx)
  Ny = ceiling((yhi - ylo) / dy)
  for (i in 0:Nx) {
    x = xlo + i*dx
    print(x)
    
    for (j in 0:Ny) {
      y = ylo + j*dy
      print(y)
      
      knots = rbind(knots, c(x,y))           
    }
  }
  return(knots)
}