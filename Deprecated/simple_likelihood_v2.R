library(gtools)
library(dplyr)

# Likelihood calculator
# Wrapper function for optim
loglik_wrapper <- function(par, data){
  beta <- data.frame(beta1=par[1],beta2=par[2],beta2y=par[3])
  alpha <- data.frame(alpha1=par[4],alpha2=par[5])
  return(-loglik(data,beta,alpha))
}

# Calculate full likelihood
loglik <- function(data, beta, alpha){
  coefs <- dat %>% group_by(id) %>% group_modify(~ getPi(.,beta,alpha))

  prob_0011 <- data$y1*data$y2*log(coefs$pi12) # Prev is 00, curr is 11
  prob_0010 <- (data$y1 * (1-data$y2))*log(coefs$pi1-coefs$pi12)# Prev is 00, curr is 10
  prob_0001 <- ((1-data$y1) * data$y2)*log(coefs$pi2_0-coefs$pi12) # Prev is 00, curr is 01
  prob_0000 <- ((1-data$y1) * (1-data$y2))*log(1-coefs$pi1-coefs$pi2_0+coefs$pi12) # Prev is 00, curr is 00
  prob_1011 <- data$y2 * log(coefs$pi2_1) # Prev is 10, curr is 11
  prob_1010 <- (1-data$y2)*log(1-coefs$pi2_1) # Prev is 10, curr is 10

  sum((1-data$y1_prev)*(prob_0011+prob_0010+prob_0001+prob_0000) + (data$y1_prev)*(prob_1011+prob_1010))
}

# Calculate likelihood for one individual
getPi <- function(data, beta, alpha){
  base <- data[1,]
  n_row <- nrow(data)

  # Calculate values of pi1 and pi2, given previous y1
  pi1 <- pi1(beta$beta1, base$trt, alpha$alpha1)
  pi2_0 <- pi2(beta$beta2, base$trt, 0, beta$beta2y, alpha$alpha2)
  pi2_1 <- pi2(beta$beta2, base$trt, 1, beta$beta2y, alpha$alpha2)
  pi12 <- pi12(theta,pi1,pi2_0)
  
  return(data.frame(pi1=rep(pi1,each=n_row),pi2_0=rep(pi2_0,each=n_row),pi2_1=rep(pi2_1,each=n_row),pi12=rep(pi12,each=n_row)))
}

# Functions to calculate equations 5-7
pi1 <- function(beta1, X, alpha1){
  inv.logit(alpha1 + X%*%beta1)
}

pi2 <- function(beta2, X, y1, beta2y, alpha2){
  inv.logit(alpha2 + X%*%beta2 + y1*beta2y)
}

pi12 <- function(theta, pi1, pi2){
  return(pi1*pi2)
}

### Generate data
gen_data <- function(n=1000, n_obs=10){
  # Create base data frame
  id <- rep(1:n, each=n_obs)
  t <- rep(70:(70+n_obs-1),n)
  trt <- rep(rbinom(n,1,0.5),each=n_obs)
  pattern <- rep(c("P00", rep(0,n_obs-1)),n)
  dat <- tibble(id,t,trt,y1=0,y2=0,pattern,y1_prev=0)
  
  # Known coefficients
  beta <- tibble(beta1=1.3, beta2=-0.5, beta2y=-1.3)
  alpha <- tibble(alpha1=-1, alpha2=-0.5)
  
  genIndividual <- function(dat.in){
    base <- dat.in[1,]
    
    # Calculate values of pi1 and pi2, given previous y1
    pi1 <- pi1(beta$beta1, base$trt, alpha$alpha1)
    pi2_0 <- pi2(beta$beta2, base$trt, 0, beta$beta2y, alpha$alpha2)
    pi2_1 <- pi2(beta$beta2, base$trt, 1, beta$beta2y, alpha$alpha2)
    pi12 <- pi12(1,pi1,pi2_0)
    
    prob_0011 <- pi12 # Prev is 00, curr is 11
    prob_0010 <- pi1-pi12 # Prev is 00, curr is 10
    prob_0001 <- pi2_0-pi12 # Prev is 00, curr is 01
    prob_0000 <- 1-pi1-pi2_0+pi12 # Prev is 00, curr is 00
    prob_1011 <- pi2_1 # Prev is 10, curr is 11
    prob_1010 <- 1-pi2_1 # Prev is 10, curr is 01
    
    for(i in 2:n_obs){
      dat.in[i,]$y1_prev <- dat.in[i-1,]$y1
      cur <- dat.in[i,]

      # Choose outcome
      if (cur$y1_prev == 0){
        out <- which(rmultinom(1,1,c(prob_0000, prob_0001, prob_0010, prob_0011))==1)
        switch(out,
               {dat.in[i,c("y1","y2","pattern")] <- list(0,0,"P00")},
               {dat.in[i,c("y1","y2","pattern")] <- list(0,1,"P01")},
               {dat.in[i,c("y1","y2","pattern")] <- list(1,0,"P10")},
               {dat.in[i,c("y1","y2","pattern")] <- list(1,1,"P11")})
      }
      else{
        out <- which(rmultinom(1,1,c(prob_1010,prob_1011))==1)
        switch(out,
               {dat.in[i,c("y1","y2","pattern")] <- list(1,0,"P10")},
               {dat.in[i,c("y1","y2","pattern")] <- list(1,1,"P11")})
      }
    }
    
    return(dat.in %>% filter(pattern != 0))
  }
  
  dat %>% group_by(id) %>% group_modify(~genIndividual(.))
}

dat <- gen_data() %>% filter(t!=70)
max_like <- optim(par=c(2,-0.5,0.5,-2,2),fn=loglik_wrapper, data=dat)
loglik_wrapper(c(1,-0.5,0.5,-2,2),dat)



