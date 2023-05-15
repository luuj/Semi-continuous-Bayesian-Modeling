source("gendata.R")
source("jags_nodeath.R")

##### Params for cluster #####
params <- commandArgs(trailingOnly=TRUE)
simNum <- as.numeric(params[[1]])

##### Create dataset #####
set.seed(simNum)
n_k <- 50 # Number of individuals per cluster
n_t <- 12 # Max number of time points per individual
k <- 30 # Number of clusters
dat <- gen_data(n_k=n_k, n_t=n_t, k=k) # Generate data

##### Main #####
# Run algorithm
results <- GS(dat=dat)
save(results, file=paste0("Output/jags", simNum, ".RData"))



