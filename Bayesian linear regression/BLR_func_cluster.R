# BLR clustered implementation
library(NormalGamma)

# Conditional log likelihood for beta
log_lik.beta <- function(beta, y, x, vk){
  p <- expit(x%*%beta + vk) # probability of success
  sum(dbinom(y, size=1, prob=p, log=T))
}

# Conditional log likelihood for vk
log_lik.vk <- function(beta, y, x, vk, sigma_v){
  log_lik.beta(beta,y,x,vk) * dnorm(vk,0,sigma_v)
}

# Conditional log likelihood for sigmav
log_lik.sv <- function(sigma_v, y, a, b){
  dnormgam(c(0,sigma_v,a,b),y)
}

# Metropolis-Hastings algorithm
# y - outcome values
# x - covariates
# beta.init - starting beta values for MH algorithm
# jump_sigma - variance for proposal distribution
# num_cluster - number of clusters
# num_iter - number of MH steps
# burn_in - number of iterations to remove for burn-in
MH.c <- function(y, x, theta_init, gamma_init, jump_sigma, num_iter, burn_in){
  # Total # of iterations
  total_iter <- num_iter + burn_in + 1
  
  # Record cluster information
  num_cluster <- length(unique(x$cluster))
  num_beta <- num_theta-1-num_cluster
  index_cluster <- x$cluster
  
  # Format x
  x <- cbind(1,x) # append 1's to x for intercept term
  x <- x %>% select(-cluster)
  
  # Format data
  y <- as.matrix(y)
  x <- as.matrix(x)
  
  # Coefficient estimate storage - betas, vk, sigma_v
  num_theta <- length(theta_init)
  theta_hat <- matrix(NA, total_iter, num_theta) 
  theta_hat[1,] <- theta_init

  # Run algorithm
  for(i in 2:total_iter){
    for(j in 1:num_theta){
      # Generate a proposal theta using a normal distribution
      proposal_theta_j <- theta_hat[i-1,j] + rnorm(1,0,jump_sigma[j])
      
      # Store previous and proposal theta
      prev_theta <- prop_theta <- theta_hat[i-1,]
      prop_theta[j] <- proposal_theta_j
      
      if (j <= num_beta){
        ## Update beta terms
        vk <- prev_theta[num_theta-1]
        beta_prop <- prop_theta[1:num_theta-2]
        beta_prev <- prev_theta[1:num_theta-2]
        
        # Calculate log posterior probability of proposed and previous
        log_prop <- log_lik.beta(beta_prop, y, x, vk)
        log_prev <- log_lik.beta(beta_prev, y, x, vk)
      }else if (j > num_beta && j < num_theta){
        ## Update Vk terms
        beta <- prop_theta[1:num_theta-2]
        vk_prop <- prop_theta[num_theta-1]
        vk_prev <- prev_theta[num_theta-1]
        sv <- prev_theta[num_theta]
        
        # Calculate log posterior probability of proposed and previous
        log_prop <- log_lik.vk(beta, y, x, theta_hat[i-1,num_theta-1], theta_hat[i-1,num_theta])
        log_prev <- log_lik.vk(beta, y, x, theta_hat[i-1,num_theta-1], theta_hat[i-1,num_theta])
      }else if (j == num_theta){
        ## Update sigma_v
        log_prop <- log_lik.sv(prop_theta, y, gamma_init[1], gamma_init[2])
        log_prev <- log_lik.sv(prev_theta, y, gamma_init[1], gamma_init[2])
      }
      
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
  }
  
  theta_hat[(burn_in+2):total_iter,]
}






  