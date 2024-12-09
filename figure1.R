
########################################################################
#12/7/2024
#Goal: Draw figure 1 of Paye & Edge, 2024.


##############################################################3
#second color is surface, third is points
#pal <- c('#66c2a5','#fc8d62','#8da0cb')
pal <- c('#1f78b4', '#a6cee3', '#000000')


#combines surface and line plots into 4-panel figure
ncptraversals <- function(d, a, b, t){
  par(mfrow = c(2, 2), mgp = c(1, 5, 6), mar = c(0, 4, 1, 2), cex.axis = 1, las = 1)
  ncpsurfaceplotfigure_risk(b, t)
  mtext("A", adj = 0)
  ncpsurfaceplotfigure_protective(b, t)
  mtext("B", adj = 0)
  
  # par(mfrow = c(1, 2), mgp = c(3, 1, 0), mar = c(5.1, 4.1, 2.1, 2.1))
  allele <- 'risk'
  
  par(mgp = c(2, 1, 0), mar = c(4, 4, 2, 2))
  riskvsprotectivefigure(d, a, b)
}

#ncp surface plots and traversals
ncpsurfaceplotfigure_risk <- function(b, t){
  #surface plot figure
  
  # d: disease frequency values at fixed intervals
  # a: risk allele frequency values at fixed intervals
  # z: function that defines the surface
  d = seq(0, 0.99, by = 0.03)
  a = seq(0, 0.99, by = 0.03)
  #t = 1
  #b = 0.5
  z <- outer(d,a, function(d,a) (((b-d)^2)*a*t*(d*t-1))/(d*(1+b*(t-1)-d*t)*(-1+d+a+a*b*(t-1)-d*a*t)))
  
  z_new <- pmin(z, 1)
  
  p <- persp(d,a,z_new, theta=20, phi=15, col=pal[2],border = 'white', lwd = 0.2, expand = 0.5, ticktype = "detailed",
             xlab=" Disease frequency", ylab=" Risk allele frequency", zlab="", cex.axis = 0.55, cex.lab = 0.7, xlim = c(0.025, .97), zlim = c(0,1))
  
  points_x <- c(0.1, 0.5)  # Example x-coordinates of points
  points_y <- c(0.01, 0.05)  # Example y-coordinates of points
  points_z <- c(0.091, 0.052)  # Example z-coordinates of points
  
  points_trans <- trans3d(points_x, points_y, points_z, pmat = p)  # Transform points to plot coordinates
  
  points(points_trans$x, points_trans$y, col = pal[3], pch = 16)  # Plot points on the 3D plot
  
  line_x <- c(0.1, 0.5)  # Example x-coordinates of line endpoints
  line_y <- c(0.01, 0.05)  # Example y-coordinates of line endpoints
  line_z <- c(0.091, 0.052)  # Example z-coordinates of line endpoints
  
  line_trans <- trans3d(line_x, line_y, line_z, pmat = p)  # Transform line endpoints to plot coordinates
  
  lines(line_trans$x, line_trans$y, col = pal[1], lwd = 2)  # Draw the line on the plot
  
  
  legend("topleft", legend = "a)", inset= c(-0.2, -.01), bty = "n")
  mtext(expression(lambda), side = 2, line = 1, at = -0.05, las = 2, cex = 1)
  
} 

ncpsurfaceplotfigure_protective <- function(b, t){
  d = seq(0, 0.99, by = 0.03)
  a = seq(0, 0.99, by = 0.03)
  #t = 1
  #b = 0.5
  z <- outer(d,a, function(d,a) ((a*t*(-1+b+d)^2*(-1+t*d))/((b*(-1+t)+t*(-1+d))*(1+a*b*(-1+t)+a*t*(-1+d)-d)*d)))
  z_new <- pmin(z, 1)
  # Plots the surface
  p <- persp(d,a,z_new, theta=-20, phi=15, col=pal[2], border = 'white', lwd = 0.2, expand = 0.5, ticktype = "detailed",
             xlab="Disease frequency", ylab="Protective allele frequency", zlab= "", cex.axis = 0.55, cex.lab = 0.7, xlim = c(0.025, .97), zlim = c(0,1))
  
  points_x <- c(0.1, 0.5)  # Example x-coordinates of points
  points_y <- c(0.01, 0.05)  # Example y-coordinates of points
  points_z <- c(0.001, 0.052)  # Example z-coordinates of points
  
  points_trans <- trans3d(points_x, points_y, points_z, pmat = p)  # Transform points to plot coordinates
  
  points(points_trans$x, points_trans$y, col = pal[3], pch = 16)  # Plot points on the 3D plot
  
  line_x <- c(0.1, 0.5)  # Example x-coordinates of line endpoints
  line_y <- c(0.01, 0.05)  # Example y-coordinates of line endpoints
  line_z <- c(0.001, 0.053)  # Example z-coordinates of line endpoints
  
  line_trans <- trans3d(line_x, line_y, line_z, pmat = p)  # Transform line endpoints to plot coordinates
  
  lines(line_trans$x, line_trans$y, col = pal[1], lwd = 2)  # Draw the line on the plot
  
  legend("topleft", legend = "b)", inset= c(-0.2, -.01), bty = "n")
  mtext(expression(lambda), side = 2, line = .8, at = .05, las = 2, cex = 1)
}

#ncp line plot in haploid case
haploidncpcalculation <- function(d, a, samp.size, df, b, allele, x){
  inflationfactor <- seq(0, x, length.out = 1000)
  
  ncpdata <- vector()
  
  for (t in inflationfactor){
    
    if (allele == 'risk'){
      ncp <- ((b-d)^2*a*t*(d*t-1))/(d*(1+b*(t-1)-d*t)*(-1+d+a+a*b*(t-1)-d*a*t))
      
    }else if (allele == 'protective'){
      ncp <- (a*t*(-1+b+d)^2*(-1+t*d))/((b*(-1+t)+t*(-1+d))*(1+a*b*(-1+t)+a*t*(-1+d)-d)*d)
      
    }
    
    ncpdata[match(t, inflationfactor)] <- ncp 
  }
  plot(inflationfactor*(1/x), ncpdata, xlab = "Case Fraction", ylab = "", type = "l", las = 1, bty = "n", pch =19, xlim = c(0, 1))
  mtext(expression(lambda), side = 2, line = 3, las = 2)
}

#haploid case risk vs protective
riskvsprotectivefigure <- function(d, a, penetrance){
  allele <- 'risk'
  # par(mgp = c(3, 1, 0), mar = c(3, 3, 3, 3) , axes = TRUE)
  haploidncpcalculation(d, a, 10000, 1, penetrance, allele, 1/d)
  #legend("topleft", legend = "c)", inset= c(-0.38, -.1), bty = "n")
  mtext("C", adj = 0)
  allele <- 'protective'
  haploidncpcalculation(d, a, 10000, 1, penetrance, allele, 1/d)
  mtext("D", adj = 0)
  #legend("topleft", legend = "d)", inset= c(-0.38, -.1), bty = "n")
}


pdf("fig1.pdf", width = 7, height = 4.5)

ncptraversals(0.3, 0.1, 1, 1)

dev.off()


