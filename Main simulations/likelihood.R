library(gtools)
library(dplyr)
library(lme4)
library(MASS)

# Likelihood calculator
# Wrapper for binary likelihood
lb <- function(par, dat){
  beta <- data.frame(beta1=par[1], beta2=par[2], betat=par[3], beta1y=par[4], beta2y=par[5])
  alpha <- data.frame(alpha1=par[6], alpha2=par[7], alphat=par[8])
  loglik_bin(dat,beta,alpha)
}

# Wrapper for continuous likelihood
lc <- function(par, dat){
  gamma <- data.frame(shape=par[9], rate1=par[10], rate2=par[11])
  loglik_con(dat,gamma)
}

# Calculate full likelihood for binary portion
loglik_bin <- function(dat, beta, alpha){
  # Obtain probabilities for each covariate
  trt <- unique(dat$trt)
  trt_n <- length(trt)
  ti <- dat$trt+1
  coefs <- data.frame(matrix(nrow = trt_n, ncol=6))
  colnames(coefs) <- c("pi1_0", "pi1_1", "pi2_0", "pi2_1", "pi12_0", "pi12_1")

  for (i in 1:trt_n){
    # Calculate values of pi1 and pi2, given previous y1
    x <- trt[i]
    theta <- theta(beta$betat, x, alpha$alphat)
    coefs$pi1_0[i] <- pi1(beta$beta1, x, 0, beta$beta1y, alpha$alpha1)
    coefs$pi1_1[i] <- pi1(beta$beta1, x, 1, beta$beta1y, alpha$alpha1)
    coefs$pi2_0[i] <- pi2(beta$beta2, x, 0, beta$beta2y, alpha$alpha2)
    coefs$pi2_1[i] <- pi2(beta$beta2, x, 1, beta$beta2y, alpha$alpha2)
    coefs$pi12_0[i] <- pi12(theta, coefs$pi1_0[i], coefs$pi2_0[i])
    coefs$pi12_1[i] <- pi12(theta,coefs$pi1_1[i], coefs$pi2_1[i])
  }

  # Probability (prob_abcd) of cd occurring, given that ab is the current state
  prob_0000 <- ((1-dat$y1) * (1-dat$y2))*log(1-coefs$pi1_0[ti]-coefs$pi2_0[ti]+coefs$pi12_0[ti])
  prob_0001 <- ((1-dat$y1) * dat$y2)*log(coefs$pi2_0[ti]-coefs$pi12_0[ti])
  prob_0010 <- (dat$y1 * (1-dat$y2))*log(coefs$pi1_0[ti]-coefs$pi12_0[ti])
  prob_0011 <- dat$y1*dat$y2*log(coefs$pi12_0[ti])
  prob_1000 <- ((1-dat$y1) * (1-dat$y2))*log(1-coefs$pi1_1[ti]-coefs$pi2_1[ti]+coefs$pi12_1[ti])
  prob_1001 <- ((1-dat$y1) * dat$y2)*log(coefs$pi2_1[ti]-coefs$pi12_1[ti])
  prob_1010 <- (dat$y1 * (1-dat$y2))*log(coefs$pi1_1[ti]-coefs$pi12_1[ti])
  prob_1011 <- dat$y1*dat$y2*log(coefs$pi12_1[ti])
  
  # Likelihood for binary portion
  sum((1-dat$y1_prev)*(prob_0011+prob_0010+prob_0001+prob_0000) + (dat$y1_prev)*(prob_1000+prob_1001+prob_1011+prob_1010))
}

# Calculate full likelihood for continuous portion
loglik_con <- function(dat, gamma){
  conLik <- dat %>% filter(cost > 0 & trt==0) %>% pull(cost) %>% gammaLik(c(gamma$shape,gamma$rate1), .)
  conLik <- conLik + dat %>% filter(cost > 0 & trt==1) %>% pull(cost) %>% gammaLik(c(gamma$shape,gamma$rate2), .)
  return(conLik)
}

# Get probability quantities pi1, pi2, pi12
getPi <- function(dat, beta, alpha){
  base <- dat[1,]
  n_row <- nrow(dat)
  
  # Calculate values of pi1 and pi2, given previous y1
  theta <- theta(beta$betat, base$trt, alpha$alphat)
  pi1_0 <- pi1(beta$beta1, base$trt, 0, beta$beta1y, alpha$alpha1)
  pi1_1 <- pi1(beta$beta1, base$trt, 1, beta$beta1y, alpha$alpha1)
  
  pi2_0 <- pi2(beta$beta2, base$trt, 0, beta$beta2y, alpha$alpha2)
  pi2_1 <- pi2(beta$beta2, base$trt, 1, beta$beta2y, alpha$alpha2)
  
  pi12_0 <- pi12(theta,pi1_0,pi2_0)
  pi12_1 <- pi12(theta,pi1_1,pi2_1)
  
  return(data.frame(pi1_0=rep(pi1_0,each=n_row),pi1_1=rep(pi1_1,each=n_row),
                    pi2_0=rep(pi2_0,each=n_row),pi2_1=rep(pi2_1,each=n_row),
                    pi12_0=rep(pi12_0,each=n_row),pi12_1=rep(pi12_1,each=n_row)))
}

# Helper functions to calculate equations 5-7
pi1 <- function(beta1, X, y1, beta1y, alpha1, c=0){
  inv.logit(alpha1 + X%*%beta1 + y1*beta1y + c)
}

pi2 <- function(beta2, X, y1, beta2y, alpha2){
  inv.logit(alpha2 + X%*%beta2 + y1*beta2y)
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
mu <- function(betam, X, alpham, d=0){
  exp(alpham + X%*%betam + d)
}

# Calculates the gamma distribution likelihood for the continuous portion
gammaLik <- function(par, x){
  alpha <- par[1]
  beta <- par[2]
  n <- length(x)
  
  out <- (alpha-1)*sum(log(x))-(1/beta)*sum(x) - n*alpha*log(beta) - n*log(gamma(alpha))
  return(out)
}

# Converts shape and rate terms for each trt group into alpham and betam
getGammaPar <- function(shape, rate1, rate2){
  alpham <- log(rate1*shape)
  betam <- log(rate2*shape) - alpham
  return(c(alpham=alpham,betam=betam))
}

# Converts alpham and betam into rate1 and rate2
convertGammaPar <- function(shape, alpham, betam){
  rate1 <- exp(alpham) / shape
  rate2 <- exp(betam + alpham) / shape
  return(c(rate1=rate1,rate2=rate2))
}





### Generate data
gen_data <- function(n=5000, n_obs=10){
  # Create base data frame
  id <- rep(1:n, each=n_obs)
  t <- rep(70:(70+n_obs-1),n)
  trt <- rep(rbinom(n,1,0.5),each=n_obs)
  dat <- tibble(id,t,trt,y1=NA,y2=NA,y1_prev=0,cost=0,cumcost=0)
  
  # Known coefficients
  beta <- tibble(beta1=1.3, beta2=-0.5, betat=1, beta1y=0.8, beta2y=-1.3, betam=3)
  alpha <- tibble(alpha1=-1, alpha2=-0.5, alphat=0.5, alpham=1)
  
  genIndividual <- function(dat.in){
    base <- dat.in[1,]
    
    # Generate random effect values
    Sigma <- matrix(c(3,0,0,3), nrow=2)
    b <- mvrnorm(n=1, mu=c(0,0), Sigma=Sigma)
    
    # Calculate values of pi1 and pi2, given previous y1
    theta <- theta(beta$betat, base$trt, alpha$alphat)
    pi1_0 <- pi1(beta$beta1, base$trt, 0, beta$beta1y, alpha$alpha1, b[1])
    pi1_1 <- pi1(beta$beta1, base$trt, 1, beta$beta1y, alpha$alpha1, b[1])
    
    pi2_0 <- pi2(beta$beta2, base$trt, 0, beta$beta2y, alpha$alpha2)
    pi2_1 <- pi2(beta$beta2, base$trt, 1, beta$beta2y, alpha$alpha2)
    
    pi12_0 <- pi12(theta,pi1_0,pi2_0)
    pi12_1 <- pi12(theta,pi1_1,pi2_1)
    
    # Calculate multinomial probabilities
    prob_0000 <- 1-pi1_0-pi2_0+pi12_0
    prob_0001 <- pi2_0-pi12_0
    prob_0010 <- pi1_0-pi12_0
    prob_0011 <- pi12_0
    prob_1000 <- 1-pi1_1-pi2_1+pi12_1
    prob_1001 <- pi2_1-pi12_1
    prob_1010 <- pi1_1-pi12_1
    prob_1011 <- pi12_1
    
    base_cost <- mu(beta$betam, base$trt, alpha$alpham, b[2])
    shape <- 5
    
    for(i in 1:n_obs){
      # Store previous non-terminal event
      if(i > 1)
        dat.in[i,]$y1_prev <- dat.in[i-1,]$y1
      
      # Calculate cost accrued based on gamma distribution
      cost <- rgamma(1,shape=shape,scale=base_cost/shape)
      
      # Choose outcome
      if (dat.in[i,]$y1_prev == 0){
        out <- which(rmultinom(1,1,c(prob_0000, prob_0001, prob_0010, prob_0011))==1)
        switch(out,
               {dat.in[i,c("y1","y2","cost")] <- list(0,0,0)},
               {dat.in[i,c("y1","y2","cost")] <- list(0,1,0)},
               {dat.in[i,c("y1","y2","cost")] <- list(1,0,cost)},
               {dat.in[i,c("y1","y2","cost")] <- list(1,1,cost)})
      }
      else{
        out <- which(rmultinom(1,1,c(prob_1000, prob_1001, prob_1010,prob_1011))==1)
        switch(out,
               {dat.in[i,c("y1","y2","cost")] <- list(0,0,0)},
               {dat.in[i,c("y1","y2","cost")] <- list(0,1,0)},
               {dat.in[i,c("y1","y2","cost")] <- list(1,0,cost)},
               {dat.in[i,c("y1","y2","cost")] <- list(1,1,cost)})
      }
      
      # Stop if terminal event
      if (dat.in[i,]$y2 == 1){
        break
      }
    }
    # Calculate cumulative cost
    dat.in$cumcost <- cumsum(dat.in$cost)
    return(dat.in %>% na.omit())
  }
  
  dat %>% group_by(id) %>% group_modify(~genIndividual(.))
}





