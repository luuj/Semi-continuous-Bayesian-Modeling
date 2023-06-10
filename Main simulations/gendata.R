library(gtools)
library(MASS)
library(dplyr)

### Generate data
# n_k = Number of individuals per cluster
# n_t = Number of observations per individual
# k = Number of clusters t
gen_data <- function(n_k=10, n_t=12, k=20){
  # Create base data frame
  n <- n_k*k
  id <- rep(1:n, each=n_t)
  t <- rep(0:(0+n_t-1),n)
  trt <- rep(rbinom(n,1,0.5),each=n_t)
  age <- rep(runif(n,18,90), each=n_t)
  clst <- factor(rep(1:k, each=(n_k*n_t)))
  hometype <- rep(rbinom(k,1,0.5), each=(n_k*n_t))
  perTrtInit <- tapply(trt, clst, mean)
  perTrt <- as.factor(rep(ifelse(perTrtInit < .5,"Less","Greater"), each=(n_k*n_t)))
  dat <- tibble(id,t,clst,trt,age,hometype,perTrt,y1=NA,y2=NA,y1_prev=0,cost=0,cumcost=0)

  # Known coefficients
  beta <- tibble(beta1=c(1.3,0.02,-0.5), beta2=c(-0.5,0.02,-0.5), betat=c(1,0.02,-0.5), betam=c(1.5,0.02,-0.5))
  betay <- tibble(beta1y=0.8, beta2y=-1.3)
  alpha <- tibble(alpha1=-1.5, alpha2=-2.5, alphat=0.5, alpham=1)
  shape <- 2
  
  # beta <- tibble(beta1=c(0.01,0.002,-0.01), beta2=c(-2,0.03,-1), betat=c(1,0.02,-0.5), betam=c(0.01,0.001,-0.01))
  # betay <- tibble(beta1y=0.8, beta2y=-1.3)
  # alpha <- tibble(alpha1=-0.5, alpha2=-2.5, alphat=0.5, alpham=4)
  # shape <- 2
  
  # Generate cluster-level random effects
  nDiff <- sum(perTrtInit < .5)
  vSigma=matrix(c(1,0,0,0,1,0,0,0,0.01), nrow=3)
  V <- mvrnorm(n=k-nDiff, mu=rep(0,3), Sigma=vSigma)
  colnames(V) <- c("V1", "V2", "Vmu")
  V <- cbind(index=which(perTrtInit >= .5), V)
  
  # Generate second level of cluster-level random effects
  vSigma2=matrix(c(1,0,0,0,1,0,0,0,0.01), nrow=3)
  V2 <- mvrnorm(n=nDiff, mu=rep(0,3), Sigma=vSigma2)
  colnames(V2) <- c("V1", "V2", "Vmu")
  V2 <- cbind(index=which(perTrtInit < .5), V2)
  
  # Combine the two RE matrices
  V <- rbind(V,V2)
  V <- V[order(V[,"index"]),-1]
  dat <- cbind(dat,V[dat$clst,])
  
  genIndividual <- function(dat.in){
    base <- dat.in[1,]
    x <- as.matrix(base[,c("trt","age","hometype")])
    
    # Calculate values of pi1 and pi2, given previous y1
    theta <- theta(beta$betat, x, alpha$alphat)
    pi1_0 <- pi1(beta$beta1, x, 0, betay$beta1y, alpha$alpha1, base$V1)
    pi1_1 <- pi1(beta$beta1, x, 1, betay$beta1y, alpha$alpha1, base$V1)

    pi2_0 <- pi2(beta$beta2, x, 0, betay$beta2y, alpha$alpha2, base$V2)
    pi2_1 <- pi2(beta$beta2, x, 1, betay$beta2y, alpha$alpha2, base$V2)
    
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
    
    base_cost <- mu(beta$betam, x, alpha$alpham, base$Vmu)
    
    for(i in 1:n_t){
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
  dat %>% group_by(id) %>% group_modify(~genIndividual(.)) %>% ungroup()
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







