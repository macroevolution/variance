

# Correlation is not appropriate for testing validity of BAMM
#   (or any other hierarchical model)

library(mvtnorm)
 
 
plotSetup <- function(){
	plot.new()
	plot.window(xlim = c(0,10), ylim = c(0,10), asp = 0)
	abline(0,1, lwd = 2, lty = "dotted", col = "gray70")
	axis(1, at = seq(-2, 10, by = 2) )
	axis(2, at = seq(-2, 10, by = 2), las = 1)

}

rr <- 0.98
psize <- 1.5
npoints <- 20

quartz.options(height=14, width=12)
plot.new()
par(oma = c(1,1,3,1), mar=c(6,6,4,1), mfrow=c(3,3))

## ---- Case 1: High sampling error at one level

set.seed(1) 
cmat <- matrix(c(1,rr, 0,rr,1, 0, 0, 0, 1), nrow=3, ncol=3) 
 
xx <- rmvnorm(npoints, sigma = cmat)
xx <- (xx - min(xx))    # Rescale
xx <- 10* xx / max(xx)   # Rescale
  
 
plotSetup()
points(xx[,1], xx[,3], pch=21, bg="red", cex=psize)

mtext(side = 1, text = expression(paste("Estimated ", theta, ", discard-data")), line=3, cex=1.15 ) 
mtext(side = 2, text = expression(paste("Estimated ", theta, ", use-all-data")), line=3, cex=1.15 ) 
 
 
plotSetup()
points(xx[,1], xx[,2], pch=21, bg="red", cex=psize)

mtext(side = 1, text = expression(paste("True ", theta)), line=3, cex=1.15 ) 
mtext(side = 2, text = expression(paste("Estimated ", theta, ", use-all-data")), line=3, cex=1.15 ) 
 
mtext(side = 3, text = "Estimates uncorrelated due to sampling variation", line = 1.5, cex=1.3, font=4) 
 
 
plotSetup()
points(xx[,2], xx[,3], pch=21, bg="red", cex = psize)
mtext(side = 1, text = expression(paste("True ", theta)), line=3, cex=1.15 ) 
mtext(side = 2, text = expression(paste("Estimated ", theta, ", discard-data")), line=3, cex=1.15 ) 
 
 
cmat <- matrix(c(1,0, 0, 0,1, 0, 0, 0, 1), nrow=3, ncol=3) 
xx <- rmvnorm(20, sigma = cmat)
xx <- (xx - min(xx))    # Rescale
xx <- 1.5 * xx / max(xx) + 5  # Rescale
  

plotSetup()
points(xx[,1], xx[,3], pch=21, bg="red", cex = psize)

mtext(side = 1, text = expression(paste("Estimated ", theta, ", discard-data")), line=3, cex=1.15 ) 
mtext(side = 2, text = expression(paste("Estimated ", theta, ", use-all-data")), line=3, cex=1.15 ) 
 
 
plotSetup()
points(xx[,1], xx[,2], pch=21, bg="red", cex = psize)

mtext(side = 1, text = expression(paste("True ", theta)), line=3, cex=1.15 ) 
mtext(side = 2, text = expression(paste("Estimated ", theta, ", use-all-data")), line=3, cex=1.15 ) 
 
mtext(side = 3, text = "Estimates uncorrelated, but highly accurate", line = 1.5, cex=1.3, font=4) 
 
plotSetup()
points(xx[,2], xx[,3], pch=21, bg="red", cex = psize)
mtext(side = 1, text = expression(paste("True ", theta)), line=3, cex=1.15 ) 
mtext(side = 2, text = expression(paste("Estimated ", theta, ", discard-data")), line=3, cex=1.15 ) 
 
 

 
cmat <- matrix(c(1,rr, 0, rr,1, 0, 0, 0, 1), nrow=3, ncol=3) 
xx <- rmvnorm(20, sigma = cmat)
xx <- (xx - min(xx))    # Rescale
xx <- 10 * xx / max(xx)  # Rescale
  

plotSetup()
points(xx[,1], xx[,2], pch=21, bg="red", cex = psize)

mtext(side = 1, text = expression(paste("Estimated ", theta, ", discard-data")), line=3, cex=1.15 ) 
mtext(side = 2, text = expression(paste("Estimated ", theta, ", use-all-data")), line=3, cex=1.15 ) 
 
 
plotSetup()
points(xx[,2], xx[,3], pch=21, bg="red", cex = psize)
mtext(side = 1, text = expression(paste("True ", theta)), line=3, cex=1.15 ) 
mtext(side = 2, text = expression(paste("Estimated ", theta, ", use-all-data")), line=3, cex=1.15 ) 
 
mtext(side = 3, text = "Estimates highly correlated but inaccurate", line = 1.5, cex=1.3, font=4) 
  
plotSetup()
points(xx[,1], xx[,3], pch=21, bg="red", cex = psize)
mtext(side = 1, text = expression(paste("True ", theta)), line=3, cex=1.15 ) 
mtext(side = 2, text = expression(paste("Estimated ", theta, ", discard-data")), line=3, cex=1.15 ) 
 
 




