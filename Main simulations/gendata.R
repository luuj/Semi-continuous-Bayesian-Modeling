### Generate data
# n_k = Number of individuals per cluster
# n_t = Number of observations per individual
# k = Number of clusters t
gen_data <- function(n_k=80, n_t=12, k=150, 
                     vSigma=matrix(c(1,0,0,0,1,0,0,0,0.01), nrow=3),
                     coefs=c(beta1=1.3,beta2=-0.5,betat=1,
                             beta1y=0.8,beta2y=-1.3,betam=3,
                             alpha1=-1.5,alpha2=-2.5,alphat=0.5,alpham=1,shape=2)){
  # Create base data frame
  n <- n_k*k
  id <- rep(1:n, each=n_t)
  t <- rep(0:(0+n_t-1),n)
  trt <- rep(rbinom(n,1,0.5),each=n_t)
  clst <- factor(rep(1:k, each=(n_k*n_t)))
  dat <- tibble(id,t,trt,clst,y1=NA,y2=NA,y1_prev=0,cost=0,cumcost=0)
  
  # Known coefficients
  beta <- tibble(beta1=coefs[1], beta2=coefs[2], betat=coefs[3], beta1y=coefs[4], beta2y=coefs[5], betam=coefs[6])
  alpha <- tibble(alpha1=coefs[7], alpha2=coefs[8], alphat=coefs[9], alpham=coefs[10])
  shape <- coefs[11]
  
  # Generate cluster-level random effects
  V <- mvrnorm(n=k, mu=rep(0,3), Sigma=vSigma)
  colnames(V) <- c("V1", "V2", "Vmu")
  dat <- cbind(dat,V[dat$clst,])
  
  genIndividual <- function(dat.in){
    base <- dat.in[1,]
    
    # Calculate values of pi1 and pi2, given previous y1
    theta <- theta(beta$betat, base$trt, alpha$alphat)
    pi1_0 <- pi1(beta$beta1, base$trt, 0, beta$beta1y, alpha$alpha1, base$V1)
    pi1_1 <- pi1(beta$beta1, base$trt, 1, beta$beta1y, alpha$alpha1, base$V1)
    
    pi2_0 <- pi2(beta$beta2, base$trt, 0, beta$beta2y, alpha$alpha2, base$V2)
    pi2_1 <- pi2(beta$beta2, base$trt, 1, beta$beta2y, alpha$alpha2, base$V2)
    
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
    
    base_cost <- mu(beta$betam, base$trt, alpha$alpham, base$Vmu)
    
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



