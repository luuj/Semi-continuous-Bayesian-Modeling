source("likelihood.R")

# Gibbs sampler
GS <- function(dat, par, jump_sigma, n_iter, burn_in){
  # Total # of iterations
  total_iter <- n_iter + burn_in
  
  # Coefficient estimate storage
  theta_hat <- matrix(0, total_iter, length(par))
  colnames(theta_hat) <- c("beta1", "beta2", "betat", "beta1y", "beta2y", "alpha1", "alpha2", "alphat", "shape", "rate1", "rate2")
  theta_hat[1,] <- par
  
  # Iterate through all coefficients
  for(i in 1:(total_iter-1)){
    for (j in colnames(theta_hat)){
      # Propose new coefficient for current variable
      prev_theta <- prop_theta <- theta_hat[i,]
      prop_theta[j] <- prop_theta[j] + rnorm(1,0,jump_sigma[1])
      
      # Use appropriate likelihood function for MH algorithm
      if (j %in% c("shape", "rate1", "rate2")){
        if (prop_theta[j] > 0)
          theta_hat[i,] <- MH(dat, prop_theta, prev_theta, lc) # Coefficients must be >0 for gamma distribution
      }else{
        theta_hat[i,] <- MH(dat, prop_theta, prev_theta, lb)
      }
    }

    theta_hat[i+1,] <- theta_hat[i,]
  }
  return(theta_hat)
}

# Metropolis-hastings update step
MH <- function(dat, prop, prev, log_lik){
  # Calculate log posterior probability of proposed
  log_prop <- log_lik(prop, dat)
  # Calculate log posterior probability of previous
  log_prev <- log_lik(prev, dat)

  # Calculate ratio
  log_ratio <- log_prop - log_prev

  # Flip coin to choose acceptance
  if (log(runif(1)) < log_ratio){
    return(prop)
  } else{
    return(prev)
  }
}

# Run test
set.seed(3)
dat <- gen_data(n=100)
truth <- c(beta1=1.3, beta2=-0.5, betat=1, beta1y=0.8, beta2y=-1.3, 
           alpha1=-1, alpha2=-0.5, alphat=1, 
           shape=5, rate1=0.544, rate2=10.92)
jump_sigma <- c(1,1,1,1,1,1,1,1,1,1,1)
coef_index <- c(1,2,4,5,6,7,9,10,11)
init_par <- rep(0,11)
n_iter <- 10000

# Results
coef_names <- c("beta1","beta2","beta1y","beta2y","alpha1","alpha2","shape","rate1","rate2")

results <- GS(dat=dat, par=truth, jump_sigma=jump_sigma, n_iter=n_iter, burn_in=0)
results.mean <- colMeans(results)

results2 <- GS(dat=dat, par=truth, jump_sigma=jump_sigma, n_iter=n_iter, burn_in=0)
results2.mean <- colMeans(results2)

results3 <- GS(dat=dat, par=truth, jump_sigma=jump_sigma, n_iter=n_iter, burn_in=0)
results3.mean <- colMeans(results3)

tibble(Variable=coef_names, truth=truth[coef_index], Estimate1=results.mean[coef_index], 
       Estimate2=results2.mean[coef_index], Estimate3=results3.mean[coef_index])


# Look at diagnostics / convergence rate
# Trace plots
par(mfrow=c(3,3))
j <- 1
for(i in coef_index){
  plot(1:n_iter, results[,i], type="l", xlab="Iteration", ylab="Beta_hat", main=coef_names[j])
  lines(1:n_iter, results2[,i], col="red")
  lines(1:n_iter, results3[,i], col="blue")
  j <- j+1
}

# Potential scale reduction - rule of thumb is 1.1
library(coda)
ptr <- mcmc.list(mcmc(results[,coef_index]), mcmc(results2[,coef_index]), mcmc(results3[,coef_index]))
gelman.diag(ptr)




# Convert to C++

# Applied to FAS cluster
# Speed up by optimizing which parts are run depending on coefficient - 249s to 24 seconds for 100 iterations
# Can run overnight - 10000 iterations - multiple chains










