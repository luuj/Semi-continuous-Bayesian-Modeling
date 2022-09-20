library(gtools)
library(dplyr)

# Likelihood calculator
# Wrapper function for optim
loglik_wrapper <- function(par, data){
  beta <- data.frame(beta1=par[1])
  alpha <- data.frame(alpha1=par[2])
  return(-loglik(data,beta,alpha))
}

# Calculate full likelihood
loglik <- function(data, beta, alpha){
  ll <- dat %>% group_by(id) %>% group_modify(~ loglik_i(.,beta, alpha))
  sum(ll$ll)
}

# Calculate likelihood for one individual
loglik_i <- function(data, beta, alpha){
  ll <- 0
  n <- nrow(data)
  base <- data[1,]

  # Calculate values of pi1 and pi2, given previous y1
  pi1 <- pi1(beta$beta1, base$trt, alpha$alpha1)
  
  ll <- sum(data$y1*c(log(pi1)) + (1-data$y1)*c(log(1-pi1)))
  data.frame(ll)
}

# Functions to calculate equations 5-7
pi1 <- function(beta1, X, alpha1){
  inv.logit(alpha1 + X%*%beta1)
}

### Generate data
gen_data <- function(n=1000, n_obs=10){
  # Create base data frame
  id <- rep(1:n, each=n_obs)
  t <- rep(70:(70+n_obs-1),n)
  trt <- rep(rbinom(n,1,0.5),each=n_obs)
  pattern <- rep(c(1, rep(0,n_obs-1)),n)
  dat <- tibble(id,t,trt,y1=NA)
  
  # Known coefficients
  beta <- tibble(beta1=2)
  alpha <- tibble(alpha1=-2)
  
  genIndividual <- function(dat.in){
    base <- dat.in[1,]
    pi1 <- pi1(beta$beta1, base$trt, alpha$alpha1)
    
    for(i in 1:n_obs){
      dat.in[i,]$y1 <- rbinom(1,1,pi1)
    # Stop if event
     if (dat.in[i,]$y1 == 1){
        break
     }
    }
    
    return(dat.in)
  }
  
  dat %>% group_by(id) %>% group_modify(~genIndividual(.))
}

dat <- gen_data() %>% na.omit()

max_like <- optim(par=c(2,-2),fn=loglik_wrapper, data=dat)



# Remove 0



