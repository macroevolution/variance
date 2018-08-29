

library(BAMMtools)
library(diversitree)


v <- read.tree("data/bamm-tortoise/t-tortoise/tortoise.tre")
 


genus <- character(length(v$tip.label))
for (ii in 1:length(v$tip.label)){
	genus[ii] <- unlist(strsplit(v$tip.label[ii], "_"))[[1]]
}


tx <- table(genus)

good <- names(tx)[tx >= 2]
counts <- tx[good]

#--------------------
# Get BAMM results for each genus from discard-data and use-all-data methods
 
# Full BAMM analysis:
edall <- getEventData(v, "data/bamm-tortoise/t-tortoise/event_data.txt", burnin=0.1, nsamples=1000)

spmat_full <- matrix(NA, nrow = length(good), ncol=1000)
rownames(spmat_full) <- good
spmat_sub  <- spmat_full
exmat_full <- spmat_full
exmat_sub  <- spmat_full

 

for (ii in 1:length(good)){
	
	gen <- good[ii]
	gset <- v$tip.label[grep(gen, v$tip.label)]
	tree <- drop.tip(v, setdiff(v$tip.label, gset))
	treefilename <- paste(gen, ".tre", sep="")
	newdir <- paste("data/bamm-tortoise/t-", gen, sep="")
	
	tree <- read.tree(paste(newdir, "/", treefilename, sep=""))
	
	edtmp <- getEventData(tree, paste(newdir, "/event_data.txt", sep=""), burnin=0.1, nsamples=1000)
	
	crates_sub <- getCladeRates(edtmp)
	
	mnode <- getMRCA(v, gset)
	
	crates_full <- getCladeRates(edall, node = mnode)
	
	spmat_full[ii, ] <- crates_full$lambda
	exmat_full[ii, ] <- crates_full$mu
	spmat_sub[ii, ]  <- crates_sub$lambda
	exmat_sub[ii, ]  <- crates_sub$mu
	
}
 

fx <- function(x) return(quantile(x, c(0.025, 0.5, 0.975)))	

ndrmat_full <- spmat_full - exmat_full
ndrmat_sub  <- spmat_sub - exmat_sub
cimat_sub <- matrix(NA, nrow = length(good), ncol=3)
cimat_full <- matrix(NA, nrow = length(good), ncol=3)
varvec_sub <- numeric(length(good))
varvec_full <- numeric(length(good))

for (ii in 1:length(good)){
	cimat_full[ii, ] <- fx(ndrmat_full[ii, ])
	cimat_sub[ii, ]  <- fx(ndrmat_sub[ii, ])
	varvec_sub[ii] <- var(ndrmat_sub[ii, ])
	varvec_full[ii] <- var(ndrmat_full[ii, ])
}

mean_rates_full <- apply(ndrmat_full, 1, mean)
mean_rates_sub  <- apply(ndrmat_sub, 1, mean)

var_ratio <- varvec_sub / varvec_full
mean(var_ratio)
median(var_ratio)


#------------------------------------------------
# ---- helper function for density plots, figure 5

plotDens <- function(x1, x2, xlim=c(-0.25, 0.25), ylim=c(0,48), col1, col2){
	plot.new()
	plot.window(xlim=xlim, ylim = ylim)
	
	dx <- density(x1)
	x  <- c(dx$x, rev(dx$x))
	y  <- c(dx$y, rep(0, length(dx$y)))

	polygon(x, y, col=col1)
 
	dx <- density(x2)
	x  <- c(dx$x, rev(dx$x))
	y  <- c(dx$y, rep(0, length(dx$y)))
 
	polygon(x, y, col=col2)
	axis(1, at= round(seq(-0.3, 0.3, by=0.1), 1) )
	axis(2, at=seq(-10, 40, by=10), las=1)
 
	
}

col1 <- BAMMtools:::transparentColor("red", 0.5)
col2 <- BAMMtools:::transparentColor("blue", 0.5)

quartz.options(height = 10, width=10)
par(oma=c(5,5,1,1))
par(mar=c(4,3,1,1))
plot.new()
par(mfrow=c(4,4))
 

for (ii in 1:13){
	plotDens(ndrmat_sub[ii,], ndrmat_full[ii,], col1=col1, col2=col2)
	label <- paste(good[ii], ", n = ", counts[ii], sep="")
	mtext(label, side=3, cex=1, font = 3, line = -1)
	if (ii %in% c(1,5, 9, 13))
		mtext("Density", side=2, line=3)
}
	
plot.new()
plot.window(xlim=c(0,10), ylim=c(0,10))
points(x=2, y = 7, pch=21, bg = col1, cex=5)
points(x=2, y = 3, pch=21, bg = col2, cex=5)
text("Discard-data", x = 3, y = 7, pos=4,cex = 1.7)
text("Use-all-data", x = 3, y = 3, pos=4, cex = 1.7)

#----------------------------------
# Figure 6: Comparison of credible intervals by clade
 

quartz.options(height = 6, width=10)
plot.new()
par(mar=c(6,6,1,1))
par(oma=c(1,1,1,1))

plot.window(xlim=c(0,14), ylim = c(-0.3, 0.3))

offset <- 0.15
for (ii in 1:length(good)){
	
	lines(x=c(ii-offset, ii-offset), y=c(cimat_full[ii,c(1,3)]), lwd=2, col="blue")
	lines(x=c(ii+offset, ii+offset), y=c(cimat_sub[ii,c(1,3)]), lwd=3, col="red")
 
}

points(x=1:13 - offset, y = cimat_full[,2], pch=21, cex=1.3, bg = "blue")
points(x=1:13 + offset, y = cimat_sub[,2], pch=21, cex=1.3, bg = "red")

axis(2, at=round(seq(-0.6, 0.6, by=0.1),2), las=1)
axis(1, at=seq(-1, 14, by=1), labels=F)
text(x=1:13, par("usr")[3], labels = good, srt = 45, adj = c(1.2,1.2), xpd = TRUE, cex=.9)

mtext(side=2, text = "Net diversification rate", line=3.5, cex=1.3)

lines(x=c(0, 1), y=c(-0.2, -0.2), lwd=1.3, col= "blue")
lines(x=c(0, 1), y=c(-0.27, -0.27), lwd=3, col= "red")
text(x=1, y=-0.2, "Use-all-data", pos=4, cex=1)
text(x=1, y=-0.27, "Discard-data", pos=4, cex=1)


#-------------------------------------
#----------------------------
# Power analysis 
 
 
# Here we assume that the "true" means for each clade
#    are those estimated from the use-all-data analysis: 
true_means <- apply(ndrmat_full, 1, mean)

# And we will empirically parameterize the variance for the discard-data 
#   level using the observed sampling variation
true_var   <- apply(ndrmat_sub, 1, var)



cor_mat <- matrix(NA, nrow=1000, ncol=2)

for (ii in 1:1000){
	
	# 1. Draw a "discard-data" vector with mean equal to the 
	#     use-all-data means, and with std deviation equal to the observed
	#     standard deviation.
	sim_means <- rnorm(length(true_means), mean = true_means, sd=sqrt(true_var))
	
	# 2. Compute correlation between the rate vectors
	stat <- cor.test(sim_means, true_means)

	cor_mat[ii, 1] <- stat$estimate
	cor_mat[ii, 2] <- stat$p.value

}


sum(cor_mat[,1]  < 0)    / 1000 # How many negative correlations?
sum(cor_mat[,1]  < 0.49) / 1000 # How many less than observed for tortoise data?
sum(cor_mat[,2]  < 0.05) / 1000 # How many "significant"?


 

#----------------------------
# Plotting distribution of correlations and p-values 
#   from the power analysis:

#- helper function for histogram plots
histFx <- function(x, breaks='Sturges', col='black', border='black') {
	z <- hist(x, breaks=breaks, plot=F) 
	for (i in 1:(length(z$breaks)-1))	
		polygon(x=c(rep(z$breaks[i],2), rep(z$breaks[i+1],2)), y=c(0, z$counts[i], z$counts[i], 0), col=col, border=border) 

}


quartz.options(height = 6 , width=10)
plot.new()
par(mfcol=c(1,2))
bks <- seq(-1.05, 1.05, by=0.05)

plot.new()
plot.window(xlim=c(-1,1), ylim=c(0, 100))
histFx(cor_mat[,1], breaks=bks, col="gray80")
lines(x=c(0.49, 0.49), y=c(0, 100), lwd=4, col="red", lty="dashed")
text(x=0.49, y=75, pos=4, labels = "Observed\n (tortoise)", cex=1.2)
axis(1, at=round(seq(-1.2, 1.2,by=0.4),1))
axis(2, at=seq(-25, 125, by=25), labels=F)
mtext("Correlation: use-all-data vs discard-data", cex=1.2, line=3, side=1)
mtext("Frequency", cex=1.2, line=1.3, side=2)


plot.new()
plot.window(xlim=c(0,1), ylim=c(0, 300))
histFx(cor_mat[,2], breaks=seq(0,1,by=0.05), col="gray80")

axis(1, at=round(seq(-1.2, 1.2,by=0.4),1))
axis(2, at=seq(-50, 350, by=50), labels=F)
mtext("p-value", cex=1.2, line=3, side=1)
mtext("Frequency", cex=1.2, line=1.3, side=2)

 
#- 





