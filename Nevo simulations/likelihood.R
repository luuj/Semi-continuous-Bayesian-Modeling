library(gtools)
library(dplyr)

# Likelihood calculator
# Wrapper function for optim
loglik_wrapper <- function(data, beta, alpha){
  beta <- data.frame(beta1=beta[1], beta2=beta[2], betat=beta[3], beta1y=beta[4], beta2y=beta[5])
  alpha <- data.frame(alpha1=alpha[1], alpha2=alpha[2], alpha3=alpha[3])
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
  
  for(i in 1:n){
    ll <- ll + loglik_k(beta, alpha, data[i,])
  }
  
  data.frame(ll)
}

# Calculate likelihood for one observation
loglik_k <- function(beta, alpha, cur){
  theta <- theta(beta$betat, cur$trt, cur$t, alpha$alpha3)
  
  # Calculate values of pi1 and pi2, given previous y1
  pi1 <- pi1(beta$beta1, cur$trt, cur$t, cur$y1_prev, beta$beta1y, alpha$alpha1)
  pi2 <- pi2(beta$beta2, cur$trt, cur$t, cur$y1_prev, beta$beta2y, alpha$alpha2)
  pi12 <- pi12(theta,pi1,pi2)

  # Calculate the 8 likelihood contributions, given current and previous y values
  y_00 <- (1-cur$y1)*(1-cur$y2)*log(1-pi1-pi2+pi12)
  y_01 <- (1-cur$y1)*cur$y2*log(pi2-pi12)
  y_10 <- cur$y1*(1-cur$y2)*log(pi1-pi12)
  y_11 <- cur$y1*cur$y2*log(pi12)
  
  return(y_00+y_01+y_10+y_11)
}

# Functions to calculate equations 5-7
pi1 <- function(beta1, X, t, y1, beta1y, alpha1){
  inv.logit(alpha1 + X%*%beta1)
  #inv.logit(alpha1 + X%*%beta1 + y1*beta1y)
}

pi2 <- function(beta2, X, t, y1, beta2y, alpha2){
  inv.logit(alpha2 + X%*%beta2 + y1*beta2y)
}

theta <- function(betat, X, t, alpha3){
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
  pattern <- rep(c("P00", rep(0,n_obs-1)),n)
  dat <- tibble(id,t,trt,y1=0,y2=0,pattern, y1_prev=0)
  
  # Known coefficients
  beta <- tibble(beta1=2, beta2=-0.5, betat=0.5, beta1y=0.5, beta2y=0.5)
  alpha <- tibble(alpha1=0, alpha2=0, alpha3=0)

  genIndividual <- function(dat.in){
    for(i in 2:n_obs){
      dat.in[i,]$y1_prev <- dat.in[i-1,]$y1
      cur <- dat.in[i,]

      # Stop if terminal event
      if (dat.in[i-1,]$y2 == 1){
        break
      }
      
      # Calculate values of pi1 and pi2, given previous y1
      theta <- theta(beta$betat, cur$trt, cur$t, alpha$alpha3)
      pi1 <- pi1(beta$beta1, cur$trt, cur$t, cur$y1_prev, beta$beta1y, alpha$alpha1)
      pi2 <- pi2(beta$beta2, cur$trt, cur$t, cur$y1_prev, beta$beta2y, alpha$alpha2)
      pi12 <- pi12(theta,pi1,pi2)
      
      # Calculate multinomial probabilities
      prob_00 <- 1-pi1-pi2+pi12
      prob_01 <- pi2-pi12
      prob_10 <- pi1-pi12
      prob_11 <- pi12
      
      # Choose outcome
      out <- which(rmultinom(1,1,c(prob_00,prob_01,prob_10,prob_11))==1)
      
      switch(out,
             {dat.in[i,c("y1","y2","pattern")] <- list(0,0,"P00")},
             {dat.in[i,c("y1","y2","pattern")] <- list(0,1,"P01")},
             {dat.in[i,c("y1","y2","pattern")] <- list(1,0,"P10")},
             {dat.in[i,c("y1","y2","pattern")] <- list(1,1,"P11")})
    }

    return(dat.in %>% filter(pattern != 0))
  }

  dat %>% group_by(id) %>% group_modify(~genIndividual(.))
}

dat <- gen_data()
table(dat$trt,dat$pattern)

max_like <- optim(par=c(0,0,0,0,0),fn=loglik,data=dat)
loglik_wrapper(dat,beta=c(2,-0.5,0.5,0.5,0.5), alpha=c(0,0,0))










# # Count censored individuals
# checkCensor <- function(dat.in){
#   tibble(censor=(sum(dat.in$y1)==0 & sum(dat.in$y2)==0))
# }
# censor <- dat %>% group_by(id) %>% group_modify(~ checkCensor(.))
# censor.id <- censor %>% filter(censor == T) %>% pull(id)
# length(censor.id)
# checkFinal <- function(dat.in){
#   n <- nrow(dat.in)
#   tibble(pattern=dat.in[n,"pattern"])
# }
# final <- dat %>% group_by(id) %>% group_modify(~ checkFinal(.))
# table(final$pattern)
