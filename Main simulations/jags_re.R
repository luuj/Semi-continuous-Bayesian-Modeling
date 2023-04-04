library("R2jags")

# Binary estimation - Probability of accruing a cost
m1 <- function(){
  # Priors:
  alphay ~ dnorm(0, 0.0001) # mean, precision
  betay ~ dnorm(0, 0.0001)
  gammay ~ dnorm(0, 0.0001)
  alpham ~ dnorm(0, 0.0001)
  betam ~ dnorm(0, 0.0001)
  shape ~ dunif(0, 100)
  alphad ~ dnorm(0, 0.0001)
  betad ~ dnorm(0, 0.0001)
  gammad ~ dnorm(0, 0.0001)
  
  # Priors for cluster random effects
  V1 ~ dmnorm(a0, (1 / (V1_sig * V1_sig)) * A0[,])
  V1_num ~ dnorm(0, 0.0016)
  V1_denom ~ dnorm(0,1)
  V1_sig <- abs(V1_num / V1_denom)

  V2 ~ dmnorm(a0, (1 / (V2_sig * V2_sig)) * A0[,])
  V2_num ~ dnorm(0, 0.0016)
  V2_denom ~ dnorm(0,1)
  V2_sig <- abs(V2_num / V2_denom)
  
  Vmu ~ dmnorm(a1, (1 / (Vmu_sig * Vmu_sig)) * A1[,])
  Vmu_num ~ dnorm(0, 0.0016)
  Vmu_denom ~ dnorm(0,1)
  Vmu_sig <- abs(Vmu_num / Vmu_denom)
  
  # Likelihood data model:
  for (i in 1:N) {
    lp1[i] <- alphay + betay * x[i] + gammay * y_prev[i] + V1[re[i]]
    y1[i] ~ dbern(ilogit(lp1[i]))
    
    lp3[i] <- alphad + betad * x[i] + gammad * y_prev[i] + V2[re[i]]
    y2[i] ~ dbern(ilogit(lp3[i]))
  }
  
  for (i in 1:N2){
    lp2[i] <- alpham + betam * x2[i] + Vmu[re2[i]]
    cost[i] ~ dgamma(shape, shape / exp(lp2[i]))
  }
}

# Returns estimated MCMC results
GS <- function(dat.in){
  dat.dist <- dat.in %>% distinct(cost, .keep_all=T) %>% filter(cost>0)
  Nre <- length(unique(dat$clst))
  Nre2 <- length(unique(dat.dist$clst))
  
  model.data <- list(
    y1 = dat.in$y1, # Outcomes
    y2 = dat.in$y2,
    y_prev = dat.in$y1_prev,
    x = dat.in$trt, # Covariates
    N = nrow(dat.in), # Full dataset
    cost = dat.dist$cost,
    x2 = dat.dist$trt,
    N2 = nrow(dat.dist), # Positive costs dataset 
    re = as.numeric(dat$clst), # Cluster RE
    re2 = as.numeric(dat.dist$clst), # Cluster RE positive cost
    a0 = rep(0,Nre), # Initial mean and variance matrices
    A0 = diag(Nre),
    a1 = rep(0,Nre2),
    A1 = diag(Nre2)
  )
  
  par <- c("alphay", "betay", "gammay","alpham", "betam", "shape", "alphad", "betad", "gammad",
           "V1", "V2", "Vmu", "V1_sig", "V2_sig", "Vmu_sig")
  
  inits <- function () {
    list(V1 = rnorm(Nre, 0, 1),
         V1_num = runif(1, 0, 5),
         V1_denom = runif(1, 0, 1),
         V2 = rnorm(Nre, 0, 1),
         V2_num = runif(1, 0, 5),
         V2_denom = runif(1, 0, 1),
         Vmu = rnorm(Nre2, 0, 1),
         Vmu_num = runif(1, 0, 5),
         Vmu_denom = runif(1, 0, 1))}
  
  mod <- R2jags::jags.parallel(data = model.data, inits = inits,
                               parameters.to.save = par, model.file = m1,
                               n.chains=3, n.iter=6000)
  
  return(mod)
}



