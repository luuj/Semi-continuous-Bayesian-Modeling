# BLR clustered implementation
library(NormalGamma)

# Conditional log likelihood for beta
log_lik.beta <- function(theta, y, x, vk){
  p <- expit(x%*%theta + vk) # probability of success
  sum(dbinom(y, size=1, prob=p, log=T))
}

# Conditional log likelihood for vk
log_lik.vk <- function(theta, y, x, vk, sigma_v){
  beta_post(theta,y,x,vk) * dnorm(vk,0,sigma_v)
}

# Conditional log likelihood for sigmav
log_lik.sv <- function(y, sigma_v, a, b){
  dnormgam(c(0,sigma_v,a,b),y)
}

# Metropolis-Hastings algorithm
# y - outcome values
# x - covariates
# beta.init - starting beta values for MH algorithm
# jump_sigma - variance for proposal distribution
# num_iter - number of MH steps
# burn_in - number of iterations to remove for burn-in
MH.c <- function(y, x, beta_init, gamma_init, jump_sigma, num_iter, burn_in){
  # Total # of iterations
  total_iter <- num_iter + burn_in + 1
  x <- cbind(1,x) # append 1's to x for intercept term
  
  # Coefficient estimate storage - betas, vk, sigma_v
  num_beta <- length(beta_init)
  theta_hat <- matrix(NA, total_iter, num_beta) 
  theta_hat[1,] <- c(beta_init,gamma_init)
  
  # Run algorithm
  for(i in 2:total_iter){
    for(j in 1:num_beta){
        # Generate a proposal theta using a normal distribution
        proposal_theta_j <- theta_hat[i-1,j] + rnorm(1,0,jump_sigma[j])
        
        # Store previous and proposal theta
        prev_theta <- prop_theta <- theta_hat[i-1,]
        prop_theta[j] <- proposal_theta_j
        
        if (j <= num_beta-2){
          # Update beta
          
          
        }else if (j == num_beta-1){
          # Update Vk
        }else if (j == num_beta){
          # Update sigma_v
        }
        
    }
    
   
    

    
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
    
  }
  
  theta_hat[(burn_in+2):total_iter,]
}






  