#### Metropolis-Hastings algorithm - Clustered ####

# Conditional log likelihood for beta
log_lik.beta <- function(beta, y, x, vk){
  p <- expit(x%*%beta + vk) # probability of success
  sum(dbinom(y, size=1, prob=p, log=T))
}

# Conditional log likelihood for vk
log_lik.vk <- function(beta, y, x, vk, sigma_v){
  log_lik.beta(beta,y,x,vk) + dnorm(vk,0,sigma_v,log = T)
}

# Conditional log likelihood for sigmav
log_lik.sv <- function(sigma_v, vk, num_cluster, a, b){
  dgamma(sigma_v, shape=a + num_cluster/2, rate=b + sum(vk^2)/2, log = T)
}

# Metropolis-Hastings algorithm
# y - outcome values
# x - covariate matrix
# theta_init - starting beta, vk, and sigma_v values for MH algorithm
# gamma_init - starting gamma prior distribution parameters
# jump_sigma - variance for proposal distribution
# n_iter - number of MH steps
# burn_in - number of iterations to remove for burn-in
MH.c <- function(y, x, theta_init, gamma_init, jump_sigma, n_iter, burn_in){
  # Total # of iterations
  num_iter <- n_iter + burn_in + 1

  # Record count information
  num_total <- length(theta_init)
  num_cluster <- length(unique(x[,"cluster"]))
  num_beta <- num_total-1-num_cluster

  # Store index information
  index_vk <- (num_beta+1):(num_total-1)
  index_beta <- 1:num_beta
  index_sigma <- num_total
  index_cluster <- x[,"cluster"]

  # Format data
  x <- x[,colnames(x)!="cluster"]

  # Coefficient estimate storage - betas, vk, sigma_v
  theta_hat <- matrix(NA, num_iter, num_total)
  theta_hat[1,] <- theta_init
  
  # Run algorithm
  for(i in 2:num_iter){
    # Generate a proposal theta using a normal distribution
    proposal_theta_j <- theta_hat[i-1,] + rnorm(num_total,0,jump_sigma)
    for(j in 1:num_total){
      # Store previous and proposal theta
      prev_theta <- prop_theta <- theta_hat[i,]
      prev_theta[j:num_total] <- prop_theta[j:num_total] <- theta_hat[i-1,j:num_total]
      prop_theta[j] <- proposal_theta_j[j]

      if (j %in% index_beta){
        ## Update beta terms
        vk <- prev_theta[index_vk]
        beta_prop <- prop_theta[index_beta]
        beta_prev <- prev_theta[index_beta]

        # Calculate log posterior probability of proposed and previous
        log_prop <- log_lik.beta(beta_prop, y, x, vk)
        log_prev <- log_lik.beta(beta_prev, y, x, vk)
      }else if (j %in% index_vk){
        ## Update Vk terms
        beta <- prop_theta[index_beta]
        vk_prop <- prop_theta[j]
        vk_prev <- prev_theta[j]
        sv <- prev_theta[index_sigma]

        # Subset y and x for cluster k
        curr_cluster <- which (j==index_vk)
        curr_x <- x[index_cluster==curr_cluster,]
        curr_y <- y[index_cluster==curr_cluster]

        # Calculate log posterior probability of proposed and previous
        log_prop <- log_lik.vk(beta, curr_y, curr_x, vk_prop, sv)
        log_prev <- log_lik.vk(beta, curr_y, curr_x, vk_prev, sv)
      }else if (j == index_sigma){
        ## Update sigma_v
        sv_prop <- prop_theta[j]
        sv_prev <- prev_theta[j]
        curr_vk <- prop_theta[index_vk]

        log_prop <- log_lik.sv(sv_prop, curr_vk, num_cluster, gamma_init[1], gamma_init[2])
        log_prev <- log_lik.sv(sv_prev, curr_vk, num_cluster, gamma_init[1], gamma_init[2])
      }

      # Calculate ratio
      log_ratio <- log_prop - log_prev

      # Flip coin to choose acceptance
      if (log(runif(1)) < log_ratio){
        theta_hat[i,j] <- prop_theta[j]
      } else{
        theta_hat[i,j] <- theta_hat[i-1,j]
      }
    }
  }

  theta_hat[(burn_in+2):num_iter,]
}






  