library(gtools)

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

loglik_k <- function(beta, curr_obs, y1_prev, y2_prev){
  trt <- curr_obs$trt
  t <- curr_obs$t
  y1 <- curr_obs$y1
  y2 <- curr_obs$y2
  
  pi1 <- pi1(beta$beta1, trt, t)
  pi2_0 <- pi2(beta$beta2, trt, t, 0, beta$beta2y)
  pi2_1 <- pi2(beta$beta2, trt, t, 1, beta$beta2y)
  theta <- theta(beta$betat, trt, t)
  pi12 <- pi12(theta,pi1,pi2_0)
  
  term1 <- (y1 * y2) * log(pi12)
  term2 <- (y1 * (1-y2)) * log(pi1-pi12)
  term3 <- ((1-y1) * y2) * log(pi2_0-pi12)
  term4 <- ((1-y1) * (1-y2)) * log(1-pi1-pi2_0+pi12)
  term5 <- y2 * log(pi2_1)
  term6 <- (1-y2) * log(1-pi2_1)
  
  (1-y1_prev)*(1-y2_prev)*(term1+term2+term3+term4) + (y1_prev)*(1-y2_prev)*(term5+term6)
}

# Functions to calculate equations 5-7
pi1 <- function(beta1, X, t){
  # Equation 5
  inv.logit(alpha1(t) + X%*%beta1)
}

pi2 <- function(beta2, X, t, y1, beta2y){
  # Equation 6
  inv.logit(alpha2(t) + X%*%beta2 + y1*beta2y)
}

theta <- function(betat, X, t){
  # Equation 7
  inv.logit(alphat(t) + X%*%betat)
}

pi12 <- function(theta, pi1, pi2_0){
  if (theta == 1){
    return(pi1*pi2_0)
  }
  
  a = (pi1 + pi2_0) * (theta-1)
  (1/(2*(theta-1))) * (1 + a - sqrt( (1 + a)^2 - 4*theta*(theta-1)*pi1*pi2_0))
}

# Functions to calculate time-dependent alphas
alpha1 <- function(t){
  logit(0.005 + 0.002*(t-65) + 0.0008*(t-70)^2 - 0.0000128*(t-62.5)^3)
}

alpha2 <- function(t){
  logit(0.03 + 0.003*(t-65) + 0.00016*(t-65)^2)
}

alphat <- function(t){
  ifelse(t<=95, 0.9 + 0.07*(t-65) - 0.0032*(t-65)^2, 0)
}
  

# Data structure
beta <- data.frame(beta1=log(0.7), beta2=log(0.5), betat=0, beta2y=log(1.4))
data <- data.frame(id=c(0,0), t=c(70,72), y1=c(0,0), y2=c(0,1), trt=c(1,1))

loglik_i(beta, data)


