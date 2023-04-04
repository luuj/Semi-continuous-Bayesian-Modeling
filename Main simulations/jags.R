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
  
  # Likelihood data model:
  for (i in 1:N) {
    lp1[i] <- alphay + betay * x[i] + gammay * y_prev[i]
    y1[i] ~ dbern(ilogit(lp1[i]))
    
    lp3[i] <- alphad + betad * x[i] + gammad * y_prev[i]
    y2[i] ~ dbern(ilogit(lp3[i]))
  }
  
  for (i in 1:N2){
    lp2[i] <- alpham + betam * x2[i]
    cost[i] ~ dgamma(shape, shape / exp(lp2[i]))
  }
}

# Returns estimated MCMC results
GS <- function(dat.in){
  dat.dist <- dat.in %>% distinct(cost, .keep_all=T) %>% filter(cost>0)
  
  model.data <- list(
    y1 = dat.in$y1,
    y2 = dat.in$y2,
    y_prev = dat.in$y1_prev,
    x = dat.in$trt,
    N = nrow(dat.in),
    cost = dat.dist$cost,
    x2 = dat.dist$trt,
    N2 = nrow(dat.dist)
  )

  par <- c("alphay", "betay", "gammay","alpham", "betam", "shape", "alphad", "betad", "gammad")
  mod <- R2jags::jags.parallel(data = model.data,
                       parameters.to.save = par, model.file = m1,
                       n.chains=3, n.iter=3000)
  
  return(mod)
}













