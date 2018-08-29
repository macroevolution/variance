
library(diversitree)
library(phytools)
library(BAMMtools)

# Generate pure-birth tree

set.seed(1)   # Set seed to integer 1, 
				# so that this analysis is repeatable!

v <- tree.bd(c(0.1,0), max.taxa=50)

plot(v)
 

# Picking set of genus-nodes
genera <- c(62, 61, 92, 91, 90, 71, 77, 78, 68, 74, 73, 56)

spset <- list()
for (ii in 1:length(genera)){
	
	spset[[ii]] <- extract.clade(v, node = genera[ii])$tip.label
 	
}

#--------------------------------------------------------------
# Here we read in the BAMM results, both for the complete tree
#  and from the analyses performed separately on each subclade.

# The code below assumes that you are in a directory from which you 
#   can access the data/bamm-pb subdirectory etc


# Read in full BAMM analysis:
edall <- getEventData(v, "data/bamm-pb/pb-all/event_data.txt", burnin=0.1, nsamples=1000)

# Matrices for the marginal posteriors for each genus
#   spmat_full & exmat_full are for the use-all-data posteriors
#   spmat_sub  & exmat_sub  are for the "genera"
spmat_full <- matrix(NA, nrow = length(genera), ncol=1000)
rownames(spmat_full) <- paste("g", genera, sep="")
spmat_sub  <- spmat_full
exmat_full <- spmat_full
exmat_sub  <- spmat_full

 

for (ii in 1:length(genera)){
  
	tree <- extract.clade(v, node = genera[ii])
	treefilename <- paste("pb-", genera[ii], ".tre", sep="")
	newdir <- paste("data/bamm-pb/pb-", genera[ii], sep="")
	
	tree <- read.tree(paste(newdir, "/", treefilename, sep=""))
	
	edtmp <- getEventData(tree, paste(newdir, "/event_data.txt", sep=""), burnin=0.1, nsamples=1000)
	
	crates_sub <- getCladeRates(edtmp)
 
	crates_full <- getCladeRates(edall, node = genera[ii])
	
	spmat_full[ii, ] <- crates_full$lambda
	exmat_full[ii, ] <- crates_full$mu
	spmat_sub[ii, ]  <- crates_sub$lambda
	exmat_sub[ii, ]  <- crates_sub$mu
	
}
 
fx <- function(x) return(quantile(x, c(0.025, 0.5, 0.975)))	

# Get marginal posterior distributions for each "genus", 
#   for net diversification rates (speciation - extinction)

ndrmat_full <- spmat_full - exmat_full
ndrmat_sub  <- spmat_sub - exmat_sub
cimat_sub <- matrix(NA, nrow = length(genera), ncol=3)
cimat_full <- matrix(NA, nrow = length(genera), ncol=3)

# Compute variances for each posterior
var_full <- apply(ndrmat_full, 1, var)
var_sub  <- apply(ndrmat_sub, 1, var)

# Confidence intervals
for (ii in 1:length(genera)){
	cimat_full[ii, ] <- fx(ndrmat_full[ii, ])
	cimat_sub[ii, ]  <- fx(ndrmat_sub[ii, ])
}

# Means for use-all-data (full) and discard-data (sub)
means_full <- apply(ndrmat_full, 1, mean)
means_sub  <- apply(ndrmat_sub, 1, mean)

 

#--------------------------------------------------------
# Figure 4: 
#   Showing marginal rate distributions for 4 focal genera
#     

lmat <- matrix(c(rep(1,4), 2:9), nrow=4, ncol=3, byrow=F)

focal <- c(2,6,8,11)
fnodes <- genera[focal]

dset <- getDescendants(v, node = fnodes[1])
for (ii in 2:length(fnodes)){
	dset <- c(dset, getDescendants(v, fnodes[ii]))
}

ew <- rep(1.5, nrow(v$edge))
ew[v$edge[,2] %in% dset] <- 4
ec <- rep("gray50", nrow(v$edge))
ec[v$edge[,2] %in% dset] <- "darkgreen"


histFx <- function(x, breaks, col = "gray90", bcol="black"){
 
	tmp <- hist(x, breaks=breaks, plot=F)
	for (ii in 1:(length(breaks)-1)){
		xv <- tmp$breaks[c(ii,ii,ii+1,ii+1)]
		yv <- c(0, tmp$density[c(ii,ii)], 0)
		polygon(xv, yv, col=col, border=bcol)
	}
	
}

quartz.options(height=12, width=13)
plot.new()
layout(mat =  lmat)

par(oma=c(1,1,1,1), mar=c(6,6,1,4))
 
plot(v, edge.color = ec, edge.width = ew, show.tip.label=F)

bks <- seq(-2, 2.6, by=0.02)

tmp <- hist(ndrmat_sub[ii,], plot=F, breaks=bks)
par(mar=c(6,6,1,1))
for (ii in focal){
	
	plot.new()
	plot.window(xlim=c(-0.2, 0.4), ylim=c(0, 25))

	histFx(ndrmat_full[ii, ], breaks=bks, col= "gray90")

	lines(x=c(0.1, 0.1), y = c(0, 25), lwd= 3, col = "blue")	
	lines(x=means_full[c(ii,ii)], y = c(0, 25), lwd=3, col = "red", lty="dotted")
	
	mtext(side=1, text = "Net diversification rate", line=3, cex=1)
	mtext(side =2, text="Density", line=3, cex=1)
	axis(1, at = seq(-0.4, 0.8, by=0.2))
	axis(2, at=seq(-5, 25, by=5), labels=T, las=1)
}


for (ii in focal){
	
	plot.new()
	plot.window(xlim=c(-0.2, 0.4), ylim=c(0, 25))
	histFx(ndrmat_sub[ii, ], breaks=bks, col = "gray90")
	lines(x=c(0.1, 0.1), y = c(0, 25), lwd= 3, col = "blue" )	
	lines(x=means_sub[c(ii,ii)], y = c(0, 25), lwd=3, col = "red", lty="dotted")
	axis(1, at = seq(-0.4, 0.8, by=0.2))
	axis(2, at=seq(-5, 25, by=5), labels=T)	
	mtext(side=1, text = "Net diversification rate", line=3, cex=1)
	mtext(side =2, text="Density", line=3, cex=1)
	
}






