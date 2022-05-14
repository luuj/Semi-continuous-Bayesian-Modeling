# BLR implementation
# Inverse-logit
expit <- function(xb){
  exp(xb) / (1+exp(xb))
}

# Logistic log likelihood
log_lik <- function(theta, y, x){
  p <- expit(x%*%theta) # probability of success
  sum(dbinom(y, size=1, prob=p, log=T))
}

# Log normal prior
log_prior <- function(theta, mu, sigma){
  sum(dnorm(theta, mu, sigma, log = T))
}

# Posterior distribution for logistic
log_post <- function(theta, y, x, prior_mu, prior_sigma) {
  return(log_lik(theta, y, x) + log_prior(theta, prior_mu, prior_sigma))
}

# Change prior so that mu does not change in the prior


# Metropolis-Hastings algorithm
MH <- function(y, x, theta_prior, prior_sigma, jump_sigma, num_iter, burn_in){
  # Total # of iterations
  total_iter <- num_iter + burn_in + 1
  x <- cbind(1,x) # append 1's to x for intercept term
  
  # Acceptance probability storage
  accept <- matrix(0, 2, length(theta_prior))
  
  # Coefficient estimate storage
  theta_hat <- matrix(NA, total_iter, length(theta_prior)) 
  theta_hat[1,] <- theta_prior
  
  # Probability vector for random scan
  prob <- c(0.4,0.1,0.1,0.4)
  num_par <- length(prob)
  
  # Run algorithm
  for(i in 2:total_iter){
    # Pick which coefficient to run
    j <- sample(1:num_par, 1, prob=prob)
    
    # Generate a proposal beta using a normal distribution
    grads <- gradients(theta_hat[i-1,]+0.4, y, x, j)
    proposal_theta_j <- rnorm(1, theta_hat[i-1,j] - grads$D1/grads$D2, -(2.4^2)/grads$D2)
    
    #proposal_theta_j <- theta_hat[i-1,j] + rnorm(1,0,jump_sigma[j])
    # Store previous and proposal betas
    prev_theta <- prop_theta <- theta_hat[i-1,]
    prop_theta[j] <- proposal_theta_j
    
    # Calculate log posterior probability of proposed
    log_prop <- log_post(prop_theta, y, x, prev_theta, prior_sigma)
    
    # Calculate log posterior probability of previous
    log_prev <- log_post(prev_theta, y, x, prev_theta, prior_sigma)
    
    # Calculate ratio
    log_ratio <- log_prop - log_prev
    
    # Flip coin to choose acceptance
    if (log(runif(1)) < log_ratio){
      theta_hat[i,] <- prop_theta
      accept[1,j] <- accept[1,j] + 1
    } else{
      theta_hat[i,] <- theta_hat[i-1,]
    }
    
    accept[2,j] <- accept[2,j] + 1
  }
  
  list(coef = theta_hat[(burn_in+2):total_iter,], accept=accept)
}




## Gradients for MH normal optimization
gradients <- function(theta, y, x, j){
  eta <- x %*% as.vector(theta)
  pi  <- exp(eta) / (1 + exp(eta))
  D1 <- as.numeric(t(y - pi) %*% x[,j])
  D2 <- -as.numeric(t(pi * (1-pi)) %*% x[,j]^2)

  list(D1=D1,D2=D2)
}



  