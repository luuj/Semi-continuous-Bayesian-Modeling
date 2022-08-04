library(gtools)
library(dplyr)

# Likelihood calculator
loglik_i <- function(beta, data){
  y1_prev <- y2_prev <- ll <- 0
  n <- nrow(data)
  
  for (i in 1:n){
    ll <- ll + loglik_k(beta, data[i,], y1_prev, y2_prev)
    y1_prev <- data[i,]$y1
    y2_prev <- data[i,]$y2
  }
  
  ll
}

# Calculate likelihood for one observation
loglik_k <- function(beta, curr_obs, y1_prev, y2_prev){
  # Extract current obs. values
  trt <- curr_obs$trt
  t <- curr_obs$t
  y1 <- curr_obs$y1
  y2 <- curr_obs$y2
  theta <- theta(beta$betat, trt, t)
  
  # Calculate values of pi1 and pi2, given previous y1
  pi1 <- pi1(beta$beta1, trt, t, y1_prev, beta$beta1y)
  pi2 <- pi2(beta$beta2, trt, t, y1_prev, beta$beta2y)
  pi12 <- pi12(theta,pi1,pi2)

  # Calculate the 8 likelihood contributions, given current and previous y values
  y_00 <- (1-y1)*(1-y2)*log(1-pi1-pi2+pi12)
  y_01 <- (1-y1)*y2*log(pi2-pi12)
  y_10 <- y1*(1-y2)*log(pi1-pi12)
  y_11 <- y1*y2*log(pi12)
  
  return(y_00+y_01+y_10+y_11)
}

# Functions to calculate equations 5-7
pi1 <- function(beta1, X, t, y1, beta1y){
  # Added y1*beta1y term
  inv.logit(alpha1(t) + X%*%beta1 + y1*beta1y)
}

pi2 <- function(beta2, X, t, y1, beta2y){
  inv.logit(alpha2(t) + X%*%beta2 + y1*beta2y)
}

theta <- function(betat, X, t){
  inv.logit(alphat(t) + X%*%betat)
}

pi12 <- function(theta, pi1, pi2){
  # Note pi1 and pi2 are functions of the previous y1
  if (theta == 1){
    return(pi1*pi2)
  }
  
  a = (pi1 + pi2) * (theta-1)
  (1/(2*(theta-1))) * (1 + a - sqrt( (1 + a)^2 - 4*theta*(theta-1)*pi1*pi2))
}

# Functions to calculate time-dependent alphas
alpha1 <- function(t){
  logit(0.2 + 0.02*(t-65) + 0.005*(t-70)^2)
}

temp <- 70:74
plot(temp, inv.logit(alpha1(temp)), type="l")

alpha2 <- function(t){
  logit(0.03 + 0.003*(t-65) + 0.00016*(t-65)^2)
}

alphat <- function(t){
  ifelse(t<=95, 0.9 + 0.07*(t-65) - 0.0032*(t-65)^2, 0)
}


# Data structure
data <- data.frame(id=rep(0,4), t=70:73, y1=c(0,1,1,0), y2=c(0,0,0,1), trt=rep(1,4))
beta <- data.frame(beta1=log(0.7), beta2=log(0.5), betat=0, beta1y=log(1.2), beta2y=log(1.4))

# Run likelihood
loglik_i(beta, data)




# Generate data
gen_data <- function(n=20, n_obs=5){
  trt <- rep(rbinom(n,1,0.5),each=n_obs)
  t <- rep(70:74,n_obs)
  id <- rep(1:n, each=n_obs)
  beta <- data.frame(beta1=5, beta2=3, betat=3, beta1y=log(1.2), beta2y=log(1.4))
  dat <- data.frame(id,t,y1=0,y2=0,trt)

  for (i in 1:n){
    # Generate sequence of outcomes
    for (j in 2:n_obs){
      ci <- (i-1)*n_obs + j
      y1_prev <- dat[ci-1,]$y1
      y2_prev <- dat[ci-1,]$y2
      curr_obs <- dat[ci,]
      trt <- curr_obs$trt
      t <- curr_obs$t

      # Stop if terminal event
      if (y2_prev == 1){
        break
      }
      
      # Calculate values of pi1 and pi2, given previous y1
      theta <- theta(beta$betat, trt, t)
      pi1 <- pi1(beta$beta1, trt, t, y1_prev, beta$beta1y)
      pi2 <- pi2(beta$beta2, trt, t, y1_prev, beta$beta2y)
      pi12 <- pi12(theta,pi1,pi2)
      
      # Calculate multinomial probabilities
      prob_00 <- 1-pi1-pi2+pi12
      prob_01 <- pi2-pi12
      prob_10 <- pi1-pi12
      prob_11 <- pi12
      
      # Choose outcome
      out <- which(rmultinom(1,1,c(prob_00,prob_01,prob_10,prob_11))==1)
      print(c(prob_00,prob_01,prob_10,prob_11,out))
      
      switch(out,
       {dat[ci,"y1"]=0
        dat[ci,"y2"]=0},
       {dat[ci,"y1"]=0
        dat[ci,"y2"]=1},
       {dat[ci,"y1"]=1
        dat[ci,"y2"]=0},
       {dat[ci,"y1"]=1
        dat[ci,"y2"]=1}
      )
    }
  }
  return(dat)
}

gen_data()

# could look at event rates and its effect on beta
# play with alpha curve for good variability to affect baseline
# can start by setting alpha as set values
# can start with freq. maximizing likelihood before bayesian




