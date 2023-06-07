source("likelihood.R")
library(stringr)

## Gibbs sampler
# dat = Dataset
# par = Main coefficients
# cre = Cluster-level random effects
# cre-sigma = Variance vector for cluster-level random effects
# jump_sigma = Variance for how far to jump for the proposal distribution
# n_iter = Number of iterations to run for the gibbs sampler
# burn_in = Number of initial iterations discarded
GS <- function(dat, par, cre=NULL, cre_sigma=NULL, 
               jump_sigma, n_iter, burn_in){
  # Total # of iterations
  total_iter <- n_iter + burn_in

  # Coefficient estimate storage
  theta_hat <- matrix(0, total_iter, length(par))
  colnames(theta_hat) <- c("beta1", "beta2", "betat", "beta1y", "beta2y", 
                           "alpha1", "alpha2", "alphat", 
                           "alpham", "betam", "shape")

  re_hat <- matrix(0, total_iter, length(cre) + 4)
  colnames(re_hat) <- c(paste0("V1_", 1:(length(cre)/3) ),
                        paste0("V2_", 1:(length(cre)/3) ),
                        paste0("Vmu_", 1:(length(cre)/3) ),
                        c("sigma_1c", "sigma_2c", "sigma_3c", "rho_c"))
  theta_hat[1,] <- par
  re_hat[1,] <- c(cre,cre_sigma)

  # Iterate through all coefficients
  for(i in 1:(total_iter-1)){
    for (j in colnames(theta_hat)){
      # Propose new coefficient for current variable
      prev_theta <- prop_theta <- theta_hat[i,]
      prop_theta[j] <- prop_theta[j] + rnorm(1,0,jump_sigma[which(names(par)==j)])
      
      # Choose what to do depending on variable
      if (j %in% c("shape", "rate1", "rate2")){
        if (prop_theta[j] > 0) # Coefficients must be >0 for gamma distribution
          theta_hat[i,] <- MH(dat, prop_theta, prev_theta, lc, re_hat[i,]) # Continuous likelihood
      }else{
        theta_hat[i,] <- MH(dat, prop_theta, prev_theta, lb, re_hat[i,]) # Binary likelihood
      }
    }
    
    # Update random effects every x iterations
    if ((i %% 20) == 1){
      for (j in colnames(re_hat)){
        prev_re <- prop_re <- re_hat[i,]
        prop_re[j] <- prop_re[j] + rnorm(1,0,0.2)
  
        # Choose what to do depending on variable
        if (str_detect(j,"V1")){
          re_hat[i,] <- MH_re(dat, theta_hat[i,], lV, prop_re, prev_re, prop_re[j], 1)
        }else if (str_detect(j,"V2")){
          re_hat[i,] <- MH_re(dat, theta_hat[i,], lV, prop_re, prev_re, prop_re[j], 2)
        }else if (str_detect(j,"Vmu")){
          re_hat[i,] <- MH_re(dat, theta_hat[i,], lVmu, prop_re, prev_re, prop_re[j], 3)
        }
      }
    }

    theta_hat[i+1,] <- theta_hat[i,]
    re_hat[i+1,] <- re_hat[i,]
  }
  return(list(theta=theta_hat,re=re_hat))
}

# Metropolis-hastings update step for theta parameters
MH <- function(dat, prop, prev, log_lik, re){
  # Calculate log posterior probability of proposed
  log_prop <- log_lik(prop, dat, re)
  # Calculate log posterior probability of previous
  log_prev <- log_lik(prev, dat, re)

  # Calculate ratio
  log_ratio <- log_prop - log_prev

  # Flip coin to choose acceptance
  if (log(runif(1)) < log_ratio){
    return(prop)
  } else{
    return(prev)
  }
}

# Metropolis-hastings update step for random effects
MH_re <- function(dat, theta, log_lik, prop_re, prev_re, curr, type){
  # Calculate log posterior probability of proposed
  log_prop <- log_lik(theta, dat, prop_re, curr, type)
  # Calculate log posterior probability of previous
  log_prev <- log_lik(theta, dat, prev_re, curr, type)
  
  # Calculate ratio
  log_ratio <- log_prop - log_prev
  
  # Flip coin to choose acceptance
  if (log(runif(1)) < log_ratio){
    return(prop_re)
  } else{
    return(prev_re)
  }
}












