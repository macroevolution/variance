
 
library(diversitree)


# Read the tree for full tortoise dataset,
#     Pull out the rates for Chelonoidis:
v <- read.tree("data/bamm-tortoise/t-tortoise/tortoise.tre")

chelo <- v$tip.label[grep("Chelonoidis" , v$tip.label)]

chelo_tree <- extract.clade(v, node = getMRCA(v, tip = chelo))


# Make constant-rate birth-death likelihood function
#   using diversitree

lfx <- make.bd(chelo_tree)

# This function can easily be optimized with and without survival 
#   conditioning, since by default, the diversitree make.bd construction
#    allows you to pass conditioning scheme as an argument:

# We will then use the diversitree::find.mle function to maximize the 
#    likelihood and to estimate the corresponding rate parameters at the max

rates_conditioned <- find.mle(lfx, c(0.5, 0.3), condition.surv = T)$par

rates_unconditioned <- find.mle(lfx, c(0.5, 0.3), condition.surv = F)$par

# Net diversification rates, conditioned:

net_con <- as.numeric(rates_conditioned[1] - rates_conditioned[2])
net_uncon <- as.numeric(rates_unconditioned[1] - rates_unconditioned[2])


net_con # versus 0.058 that MRW2018 obtain from BAMM for use-all-data

net_uncon   # versus 0.125 that MRW2018 obtain from BAMM for use-all-data


 


