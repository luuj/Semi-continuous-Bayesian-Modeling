library(gtools)
library(dplyr)
library(tidyr)
library(lme4)
library(MASS)
library(ggplot2)

# Likelihood calculator
# Wrapper for binary likelihood
lb <- function(par, dat, re){
  # Setup parameter data structures
  beta <- data.frame(beta1=par[1], beta2=par[2], betat=par[3], beta1y=par[4], beta2y=par[5])
  alpha <- data.frame(alpha1=par[6], alpha2=par[7], alphat=par[8])
  re_size <- (length(re)-4)/3
  V1 <- data.frame(V1= re[1:(re_size)])
  V2 <- data.frame(V2= re[(re_size+1):(re_size*2)])
  
  # Get unique covariates / random effect values
  environment(getPi) <- environment()
  uniqueDat <- dat %>% distinct(trt,V1,V2)
  uniqueDat <- cbind(uniqueDat, pi1_0=0, pi1_1=0, pi2_0=0, pi2_1=0, pi12_0=0, pi12_1=0)
  for (i in 1:nrow(uniqueDat)){
    uniqueDat[i,4:9] <- getPi(uniqueDat[i,])
  }
  dat <- left_join(dat, uniqueDat, by=c("trt","V1"))

  getLL(dat)
}

# Wrapper for continuous likelihood
lc <- function(par, dat, re){                  
  gamma <- data.frame(shape=par[9], rate1=par[10], rate2=par[11])
  conLik <- dat %>% filter(cost > 0 & trt==0) %>% pull(cost) %>% gammaLik(c(gamma$shape,gamma$rate1), .)
  conLik <- conLik + dat %>% filter(cost > 0 & trt==1) %>% pull(cost) %>% gammaLik(c(gamma$shape,gamma$rate2), .)
  return(conLik)
}

# Wrapper for cluster-level random effects
# Specify type = {1=sigma_1c, 2=sigma_2c}
lV <- function(par, dat, re, curr, type=1){
  ll <- lb(par, dat, re)
  sig <- re[length(re)-4+type]
  ll + log(dnorm(curr,0,sig))
}

lVmu <- function(par, dat, re, curr, type=3){
  ll <- lc(par, dat, re)
  sig <- re[length(re)-4+type]
  ll + log(dnorm(curr,0,sig))
}

# Get probability quantities pi1, pi2, pi12
getPi <- function(dat.in){
  # Calculate values of pi1 and pi2, given previous y1
  theta <- theta(beta$betat, dat.in$trt, alpha$alphat)
  pi1_0 <- pi1(beta$beta1, dat.in$trt, 0, beta$beta1y, alpha$alpha1, dat.in$V1)
  pi1_1 <- pi1(beta$beta1, dat.in$trt, 1, beta$beta1y, alpha$alpha1, dat.in$V1)
  
  pi2_0 <- pi2(beta$beta2, dat.in$trt, 0, beta$beta2y, alpha$alpha2, dat.in$V2)
  pi2_1 <- pi2(beta$beta2, dat.in$trt, 1, beta$beta2y, alpha$alpha2, dat.in$V2)
  
  pi12_0 <- pi12(theta,pi1_0,pi2_0)
  pi12_1 <- pi12(theta,pi1_1,pi2_1)
  
  return(c(pi1_0=pi1_0, pi1_1=pi1_1, pi2_0=pi2_0, pi2_1=pi2_1, pi12_0=pi12_0, pi12_1=pi12_1))
}

# Sum up likelihood for binary component
getLL <- function(dat){
  prob_0000 <- ((1-dat$y1) * (1-dat$y2))*as.numeric(log(1-dat$pi1_0-dat$pi2_0+dat$pi12_0))
  prob_0001 <- ((1-dat$y1) * dat$y2)*as.numeric(log(dat$pi2_0-dat$pi12_0))
  prob_0010 <- (dat$y1 * (1-dat$y2))*as.numeric(log(dat$pi1_0-dat$pi12_0))
  prob_0011 <- dat$y1*dat$y2*as.numeric(log(dat$pi12_0))
  prob_1000 <- ((1-dat$y1) * (1-dat$y2))*as.numeric(log(1-dat$pi1_1-dat$pi2_1+dat$pi12_1))
  prob_1001 <- ((1-dat$y1) * dat$y2)*as.numeric(log(dat$pi2_1-dat$pi12_1))
  prob_1010 <- (dat$y1 * (1-dat$y2))*as.numeric(log(dat$pi1_1-dat$pi12_1))
  prob_1011 <- dat$y1*dat$y2*as.numeric(log(dat$pi12_1))
  
  sum((1-dat$y1_prev)*(prob_0011+prob_0010+prob_0001+prob_0000) + (dat$y1_prev)*(prob_1000+prob_1001+prob_1011+prob_1010))
}

# Helper functions to calculate equations 5-7
pi1 <- function(beta1, X, y1, beta1y, alpha1, V1=0){
  inv.logit(alpha1 + X%*%beta1 + y1*beta1y + V1)
}

pi2 <- function(beta2, X, y1, beta2y, alpha2, V2=0){
  inv.logit(alpha2 + X%*%beta2 + y1*beta2y + V2)
}

theta <- function(betat, X, alphat){
  inv.logit(alphat + X%*%betat)
}

pi12 <- function(theta, pi1, pi2){
  if (theta > 0.99 & theta < 1.01){
    return(pi1*pi2)
  }
  
  a = (pi1 + pi2) * (theta-1)
  (1/(2*(theta-1))) * (1 + a - sqrt( (1 + a)^2 - 4*theta*(theta-1)*pi1*pi2))
}

# Function to calculate the continuous portion of the model
mu <- function(betam, X, alpham, Vmu=0){
  exp(alpham + X%*%betam + Vmu)
}

# Calculates the gamma distribution likelihood for the continuous portion
gammaLik <- function(par, x){
  alpha <- par[1]
  beta <- par[2]
  n <- length(x)
  
  (alpha-1)*sum(log(x))-(1/beta)*sum(x) - n*alpha*log(beta) - n*log(gamma(alpha))
}











