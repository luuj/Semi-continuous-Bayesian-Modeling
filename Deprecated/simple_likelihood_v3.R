library(gtools)
library(dplyr)

### Likelihood calculator
# Wrapper function for optim
loglik_wrapper <- function(par, data){
  beta <- data.frame(beta1=par[1], beta2=par[2], betat=par[3], beta2y=par[4])
  alpha <- data.frame(alpha1=par[5], alpha2=par[6], alpha3=par[7])
  return(-loglik(data,beta,alpha))
}

# Calculate full likelihood
loglik <- function(data, beta, alpha){
  # Calculate pi, the probability of experiencing a terminal/non-terminal event
  coefs <- dat %>% group_by(id) %>% group_modify(~ getPi(.,beta, alpha))
  
  # Probability (prob_abcd) of cd occurring, given that ab is the current state
  prob_0011 <- data$y1*data$y2*log(coefs$pi12)
  prob_0010 <- (data$y1 * (1-data$y2))*log(coefs$pi1-coefs$pi12)
  prob_0001 <- ((1-data$y1) * data$y2)*log(coefs$pi2_0-coefs$pi12) 
  prob_0000 <- ((1-data$y1) * (1-data$y2))*log(1-coefs$pi1-coefs$pi2_0+coefs$pi12) 
  prob_1011 <- data$y2 * log(coefs$pi2_1) 
  prob_1010 <- (1-data$y2)*log(1-coefs$pi2_1) 
  
  sum((1-data$y1_prev)*(prob_0011+prob_0010+prob_0001+prob_0000) + (data$y1_prev)*(prob_1011+prob_1010))
}

# Calculate probabilities of experiencing events
getPi <- function(data, beta, alpha){
  base <- data[1,]
  n_row <- nrow(data)

  theta <- theta(beta$betat, base$trt, alpha$alpha3)
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

theta <- function(betat, X, alpha3){
  inv.logit(alpha3 + X%*%betat)
}

pi12 <- function(theta, pi1, pi2){
  # Note pi1 and pi2 are functions of the previous y1
  if (theta == 1){
    return(pi1*pi2)
  }
  
  a = (pi1 + pi2) * (theta-1)
  (1/(2*(theta-1))) * (1 + a - sqrt( (1 + a)^2 - 4*theta*(theta-1)*pi1*pi2))
}




### Generate data
gen_data <- function(n=1000, n_obs=10){
  # Create base data frame
  id <- rep(1:n, each=n_obs)
  t <- rep(70:(70+n_obs-1),n)
  trt <- rep(rbinom(n,1,0.5),each=n_obs)
  dat <- tibble(id,t,trt,y1=NA,y2=NA,y1_prev=0)
  
  # Known coefficients
  beta <- tibble(beta1=1.3, beta2=-0.5, betat=1, beta2y=-1.3)
  alpha <- tibble(alpha1=-1, alpha2=-0.5, alpha3=1)
  
  genIndividual <- function(dat.in){
    base <- dat.in[1,]
    
    # Calculate values of pi1 and pi2, given previous y1
    theta <- theta(beta$betat, base$trt, alpha$alpha3)
    pi1 <- pi1(beta$beta1, base$trt, alpha$alpha1)
    pi2_0 <- pi2(beta$beta2, base$trt, 0, beta$beta2y, alpha$alpha2)
    pi2_1 <- pi2(beta$beta2, base$trt, 1, beta$beta2y, alpha$alpha2)
    pi12 <- pi12(theta,pi1,pi2_0)
    
    prob_0011 <- pi12 # Prev is 00, curr is 11
    prob_0010 <- pi1-pi12 # Prev is 00, curr is 10
    prob_0001 <- pi2_0-pi12 # Prev is 00, curr is 01
    prob_0000 <- 1-pi1-pi2_0+pi12 # Prev is 00, curr is 00
    prob_1011 <- pi2_1 # Prev is 10, curr is 11
    prob_1010 <- 1-pi2_1 # Prev is 10, curr is 01
    
    for(i in 1:n_obs){
      # Store previous non-terminal event
      if(i > 1)
        dat.in[i,]$y1_prev <- dat.in[i-1,]$y1
      
      # Choose outcome
      if (dat.in[i,]$y1_prev == 0){
        out <- which(rmultinom(1,1,c(prob_0000, prob_0001, prob_0010, prob_0011))==1)
        switch(out,
               {dat.in[i,c("y1","y2")] <- list(0,0)},
               {dat.in[i,c("y1","y2")] <- list(0,1)},
               {dat.in[i,c("y1","y2")] <- list(1,0)},
               {dat.in[i,c("y1","y2")] <- list(1,1)})
      }
      else{
        out <- which(rmultinom(1,1,c(prob_1010,prob_1011))==1)
        switch(out,
               {dat.in[i,c("y1","y2")] <- list(1,0)},
               {dat.in[i,c("y1","y2")] <- list(1,1)})
      }
      
      # Stop if terminal event
      if (dat.in[i,]$y2 == 1){
        break
      }
    }
    return(dat.in %>% na.omit())
  }
  
  dat %>% group_by(id) %>% group_modify(~genIndividual(.))
}


# Generate data and calculate estimates
dat <- gen_data()
max_like <- optim(par=c(1.3,-0.5,0.5,-1.3,-1,-0.5,0.5), fn=loglik_wrapper, data=dat)















