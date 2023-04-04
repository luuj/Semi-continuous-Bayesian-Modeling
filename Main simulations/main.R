source("likelihood.R")
source("gendata.R")
source("jags_re.R")
# source("jags.R")

##### Create dataset #####
set.seed(3)
n_k <- 60 # Number of individuals per cluster
n_t <- 12 # Max number of time points per individual
k <- 30 # Number of clusters
dat <- gen_data(n_k=n_k, n_t=n_t, k=k) # Generate data
write.csv(dat, file = "dat.csv")

##### Main #####
dat <- read.csv("dat.csv")
truth_coef <- c(alphay=-1.5, betay=1.3, gammay=0.8,
                alpham=1, betam=3, shape=2,
                alphad=-2.5, betad=-0.5, gammad=-1.3)

# Run algorithm
results <- GS(dat=dat)

# Results
results
par <- c("alphay", "betay", "gammay","alpham", "betam", "shape", "alphad", "betad", "gammad",
        "V1", "V1_sig", "V2_sig", "Vmu_sig")
traceplot(results, mfrow=c(3,3), ask=F)










##### Params for cluster #####
# params <- commandArgs(trailingOnly=TRUE)
# simNum <- as.numeric(params[[1]])
# write.csv(results, paste0("out", simNum, ".csv"))
# load(file="Output/jags_re.RData")
# save(results, file="Output/jags.RData")

