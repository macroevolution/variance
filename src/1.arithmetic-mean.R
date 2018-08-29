

clades <- 13

set.seed(1)
 
# Choose means for each of 13 populations 
means  <- runif(clades, min = 0.1, max=0.5)
stdev  <- 1   
subgroup_size <- 4   
 
rmat <- matrix(NA, nrow=clades, ncol=5) 
 
for (ii in 1:clades){
	
	# Sample 63 individuals, as in tortoise dataset
	use_all_data <- rnorm(63, mean = means[ii], sd=stdev)
 
	# Generate dataset for "discard-data" approach
	#  by randomly sampling  
	#  from the use_all_data vector
	discard_data <- sample(use_all_data, size =  subgroup_size)
	
	# Corresponding means at these 2 levels:
	mean_use_all <- mean(use_all_data)
	mean_discard <- mean(discard_data)
	
	
	rmat[ii, 1:3]   <- c(mean_use_all, mean_discard, means[ii])
	
	# Confidence intervals:
	sigma_all <- sd(use_all_data)
	sigma_discard <- sd(discard_data)
	rmat[ii, 5] <- 1.96 * sigma_discard/ sqrt(4)
	rmat[ii, 4] <- 1.96 * sigma_all / sqrt(63)
	
}

#Correlation between use-all-data and discard-data estimates:
cor.test(rmat[,1], rmat[,2])

 
#------------------------
# plot:
quartz.options(height = 6, width=6)
plot.new()
par(mar=c(6,6,1,1))
plot.window(xlim = c(-0.6, 1 ), ylim=c(-0.6, 1), asp=1)

abline(0,1, lwd=1.5, col="black", lty="dotted")
points(rmat[,2], rmat[,1], pch=21, bg= "red", cex = 1.5)
axis(1, at=seq(-1,1,by=0.5))
axis(2, at=seq(-1,1,by=0.5), las=1)
mtext("Arithmetic mean, discard-data", side=1, line = 3, cex=1.5)
mtext("Arithmetic mean, use-all-data", side=2, line = 3, cex=1.5)
text(x=-0.6, y = 0.75, labels = "r = 0.09", cex=1.3, font=3, pos=4)

# ---------
# Here is a replotting showing 95% CIs on each point

plot.new()
par(mar=c(6,6,1,1))
plot.window(xlim = c(-1, 1.6), ylim=c(-1, 1.6), asp=1)

#abline(0,1, lwd=1.5, col="black", lty="dotted")

upper_all <- rmat[,1] + rmat[,4]
lower_all <- rmat[,1] - rmat[,4]
upper_sub <- rmat[,2] + rmat[,5]
lower_sub <- rmat[,2] - rmat[,5]

segments(x0=rmat[,2], y0=lower_all, y1=upper_all)
segments(x0=lower_sub, x1=upper_sub, y0=rmat[,1])
 

points(rmat[,2], rmat[,1], pch=21, bg= "red", cex = 1.5)
axis(1, at=seq(-1.5,2,by=0.5))
axis(2, at=seq(-1.5,2,by=0.5), las=1)
mtext("Arithmetic mean, discard-data", side=1, line = 3, cex=1.5)
mtext("Arithmetic mean, use-all-data", side=2, line = 3, cex=1.5)
text(x=-1.5, y = 2, labels = "r = 0.09", cex=1.3, font=3, pos=4)



#-------------------------
# Repeat for large sample

clades <- 13

REPS <- 1000
cor_mat <- matrix(NA, nrow=REPS, ncol=2)

for (kk in 1:REPS){
	
	means  <- runif(clades, min = 0.1, max=0.5)
	
	rmat <- matrix(NA, nrow=clades, ncol=3) 
 
	for (ii in 1:clades){
	
		# Sample 63 individuals, as in tortoise dataset
		use_all_data <- rnorm(63, mean = means[ii], sd=stdev)
	
		# Generate dataset for "discard-data" approach
		#  by randomly sampling 3 individuals 
		#  from the use_all_data vector
		discard_data <- sample(use_all_data, size =  subgroup_size) 
	
		mean_use_all <- mean(use_all_data)
		mean_discard <- mean(discard_data)
	
		rmat[ii, ]   <- c(mean_use_all, mean_discard, means[ii])
	
	}
 	stat <- cor.test(rmat[,1], rmat[,2])
	cor_mat[kk, ] <- c(stat$estimate, stat$p.value)
}
 
 
sum(cor_mat[,1] < 0) / nrow(cor_mat) # what fraction less than 0?

sum(cor_mat[,1] < 0.49) / nrow(cor_mat) # fraction less than tortoise observed?

sum(cor_mat[,2] < 0.05) /nrow(cor_mat)  # fraction significant?



#----------------------------------------------
# The variance of the estimators is not equal!

# Given a true population mean, what is the variance
#   of the estimator, given a sample of 63 individuals?

REPS <- 10000 
mu_bigsample <- rep(NA, REPS)
mu_smallsample <- rep(NA, REPS)

for (ii in 1:REPS){
	mu_bigsample[ii]   <- mean(rnorm(63, sd=stdev))
	mu_smallsample[ii] <- mean(rnorm(subgroup_size, sd=stdev))
}

var(mu_bigsample)
var(mu_smallsample)

# How does this accord with the analytical expectation?

# For the big sample (n = 63):

(1 / 63) * (stdev^2)


# For the discard-data subsample (n = 3):

(1 / 3) * (stdev^2)


#----------------------------------------------
# Variance x sample-size curve

nvec <- 2:63
#analytical_var <- (1 / nvec^2) * (nvec * stdev^2) = (1 / nvec) * stdev^2


nvec <- 2:63
analytical_var <- (1 / nvec) 

ratio <- analytical_var / analytical_var[62]


quartz.options(height=6, width=6)
plot.new()
par(oma=c(1,1,1,1), mar=c(6,6,1,1))
plot.window(xlim = c(0,64), ylim=c(0, 0.6))

points(nvec, analytical_var, pch=21, bg="red", cex=1)

axis(1, at=seq(-10, 60, by=10))
axis(2, at=seq(-0.1, 0.6, by=0.1), las=1)

mtext("Sample size to compute estimator \n (arithmetic mean)", side=1, line=4, cex=1.3)
mtext("Normalized variance (analytical)", side=2, line=3.5, cex=1.3)







