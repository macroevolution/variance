
# Simple script to illustrate properties of the variance in diversification rate 
#  estimates, as a function of tree size (and overall).
#   
#  The variance in estimated rates associated with even a simple constant-rate 
#    process, is quite high. It is worth going through these exercises to see
#    just how high the estimation error can be. 

# Here, we will simulate trees under a pure-birth process


# Necessary functions for simulating trees & fitting constant-rate model
source("src/diversification_functions1.R")

# Analysis of sampling variance for estimators of net diversification rate

library(diversitree)

# Read MRW2018 tortoise rates; we will use these to parameterize our simulations
#   but the variance is not especially sensitive to speciation and extinction 
#   parameters. Much more sensitive to clade size.

xx <- read.csv("data/MRW-tableS1.csv", stringsAsFactors=F, header=F)
colnames(xx) <- c("genus", "fulltree_speciation", "fulltree_extinction", "fulltree_netdiv", "subtree_speciation", "subtree_extinction", "subtree_netdiv")


mean_lambda <- mean(xx$fulltree_speciation)
mean_mu     <- mean(xx$fulltree_extinction)

 
# We will fix diversification rates for this exercise...
lambda <- mean_lambda
mu     <- mean_mu


# vector of clade sizes
nvec   <- 3:63

REPS <- 500

rmat <- matrix(NA, nrow= length(nvec), ncol = REPS)

for (ii in 1:length(nvec)){
	cat(ii, "\n")
	for (kk in 1:REPS){
		
		
		tmp <- simulateTree(c(lambda, mu), max.taxa = nvec[ii])
	
		rr <- fitCRBD(tmp, lmax=50)
		rmat[ii, kk] <- as.numeric(rr[1] - rr[2])
		
		
	}
 
	
}

# compute variance of the estimator across different clade sizes:

var_ndr <- apply(rmat, 1, var)


# Plot variance of estimator as a function of clade size

plot(log10(var_ndr) ~ nvec)

# The variance is extreme, and much greater, for clades with < 5 
#   species


# compute ratios with respect to the full dataset:
ratio <- var_ndr / var_ndr[length(var_ndr)]


quartz.options(height=6, width=10)
par(mar=c(6,6,1,1), oma=c(1,1,1,1), mfrow=c(1,3))

plot.new()
plot.window(xlim=c(0,64), ylim=c(0,105))
points(nvec, log10(var_ndr), pch=21, cex=0.7, bg="red")


hist(rmat[2,], breaks=50) 


#------------------------------------------------------
# Strictly for fun, we will repeat this analysis 
#  but comparing the variance of MS estimators, versus 
#    the constant-rate birth-death process.
#   To keep the numbers of parameters the same, we will simulate
 #   under a pure-birth process, and constrain our estimates to 
 #    be performed under pure-birth process as well.
#  

pbmat <- matrix(NA, nrow= length(nvecX), ncol = REPS)
msmat <- matrix(NA, nrow= length(nvecX), ncol = REPS)

nvecX <- seq(5, 500, by = 5)

for (ii in 1:length(nvecX)){
 
	cat(ii, "\n")
	for (kk in 1:REPS){
		
		
		tmp <- simulateTree(c(0.1, 0), max.taxa = nvecX[ii])
	 	
	 	age <- max(branching.times(tmp))
	 	
	 	ss <- sum(tmp$edge.length)  
	 
		pbmat[ii, kk] <- as.numeric((nvecX[ii]-2 )/ ss)
		msmat[ii, kk] <- (log(nvecX[ii]) - log(2) ) / age
		
	}
 
	
}

# We can also compute the analytical variances for lambda from the Moran (1951) estimator
#    (discussed by S. Nee, 2001, "Inferring speciation rates from phylogenies", Evolution)
#
# 

var_lambda_sim        <- apply(pbmat, 1, var)
var_lambda_analytical <- 0.1^2 / (nvecX - 2)

fit <- lm(var_lambda_sim ~ var_lambda_analytical)

summary(fit)

# The first few observations are pretty bad compared to the analytical variance.
#  I expect this has to do with how the final waiting time is simulated in 
#   diversitree, or alternatively, the assumption about this wait time made by the
#     the Moran estimator (i'm unsure of this... )

fit2 <- lm(var_lambda_sim[-(1:3)] ~ var_lambda_analytical[-(1:3)])
summary(fit2)
# much better after dropping first 3 (n = 3:5)


var_lambda_ms <- apply(msmat, 1, var)

plot(var_lambda_ms, var_lambda_analytical)

# How about on a log scale?

plot(log(var_lambda_ms), log(var_lambda_analytical), asp=1)
abline(0, 1)

# you can easily see that the variance is vastly higher for the MS estimators
#   and it gets proportionately worse, as you increase tree size. 
#  This makes intuitive sense, because there is simply much more data
#   in the branch lengths. Why discard them?


# ANother view, focusing on 500-taxon clades. Histograms
#   of rate estimates for the 2 approaches.

bks <- seq(-0.004, 0.3, by=0.002)
plot.new()
par(mfrow=c(2,1))
z <- hist(pbmat[100,], breaks=bks, xlim=c(0.05, 0.15), ylim=c(0,100))
hist(msmat[100,], breaks=bks, xlim=c(0.05, 0.15), ylim=c(0,100))


 










