# BLR implementation
# Inverse-logit
expit <- function(xb){
  1 / (1+exp(-xb))
}

logit <- function(xb){
  log(xb/(1-xb))
}

# Logistic log likelihood
log_lik.u <- function(theta, y, x){
  p <- expit(x%*%theta) # probability of success
  sum(dbinom(y, size=1, prob=p, log=T))
}

# Posterior distribution for logistic - flat prior
log_post.u <- function(theta, y, x) {
  return(log_lik.u(theta, y, x))
}

# Metropolis-Hastings algorithm
# y - outcome values
# x - covariates
# beta_init - starting beta values for MH algorithm
# jump_sigma - variance for proposal distribution
# n_iter - number of MH steps
# burn_in - number of iterations to remove for burn-in
MH.u <- function(y, x, beta_init, jump_sigma, n_iter, burn_in){
  # Format data
  y <- as.matrix(y)
  x <- as.matrix(x)
  
  # Total # of iterations
  total_iter <- n_iter + burn_in + 1
  
  # Format x
  x <- cbind(1,x) 
  
  # Acceptance probability storage
  accept <- matrix(0, 2, length(beta_init))
  
  # Coefficient estimate storage
  theta_hat <- matrix(NA, total_iter, length(beta_init)) 
  theta_hat[1,] <- beta_init
  
  # Probability vector for random scan
  prob <- c(0.4,0.1,0.1,0.4)
  num_par <- length(prob)
  
  # Run algorithm
  for(i in 2:total_iter){
    # Pick which coefficient to run
    j <- sample(1:num_par, 1, prob=prob)
    
    # Generate a proposal beta using a normal distribution
    proposal_theta_j <- theta_hat[i-1,j] + rnorm(1,0,jump_sigma[j])
    
    # Store previous and proposal betas
    prev_theta <- prop_theta <- theta_hat[i-1,]
    prop_theta[j] <- proposal_theta_j
    
    # Calculate log posterior probability of proposed
    log_prop <- log_post.u(prop_theta, y, x)
    
    # Calculate log posterior probability of previous
    log_prev <- log_post.u(prev_theta, y, x)
    
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
# gradients <- function(theta, y, x, j){
#   eta <- x %*% as.vector(theta)
#   pi  <- exp(eta) / (1 + exp(eta))
#   D1 <- as.numeric(t(y - pi) %*% x[,j])
#   D2 <- -as.numeric(t(pi * (1-pi)) %*% x[,j]^2)
# 
#   list(D1=D1,D2=D2)
# }

#grads <- gradients(theta_hat[i-1,]+0.4, y, x, j)
#proposal_theta_j <- rnorm(1, theta_hat[i-1,j] - grads$D1/grads$D2, -(2.4^2)/grads$D2)

  