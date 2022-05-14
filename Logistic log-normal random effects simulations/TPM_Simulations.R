library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(mvtnorm)
expit <- function(x){plogis(x)}



## Parameters
# n = number of individuals simulated
# time_periods = number of time periods each individual is observed
# beta = coefficients for beta
# gamma = coefficients for gamma
# re_sigma1, re_sigma2, re_rho = parameters for random effects
# ln_sigma = parameter for log-normal distribution
genData <- function(n=10000, time_periods = 3, 
                    beta, gamma,
                    re_sigma1, re_sigma2, re_rho,
                    ln_sigma,
                    seed=123){
  set.seed(seed)
  
  # Binary time invariant covariate - takes value 0 or 1 with equal probability
  X1 <- sample(0:1, n, replace=TRUE)
  
  # Time variant covariate - how many months an individual is observed
  X2 <- rep(time_periods, n)
  
  # Random effect values
  Sigma <- matrix(c(re_sigma1^2,re_rho*re_sigma1*re_sigma2, re_rho*re_sigma1*re_sigma2, re_sigma2^2), nrow = 2)
  rand_vals <- rmvnorm(n, sigma = Sigma)
  c <- rand_vals[,1]
  d <- rand_vals[,2]
  
  # Duplicate values for data frame 
  X1.vals <- rep(X1, X2)
  X2.vals <- rep(seq(1:time_periods), n)
  c.vals <- rep(c, X2)
  d.vals <- rep(d, X2)
  id <- rep(1:n, X2)
  dat <- data.frame(id=id, time=X2.vals, gender=X1.vals, cost=0, re_c=c.vals, re_d=d.vals)
  
  # Logit parameters
  ## Note the time parameter - for positive beta, probability of accruing cost increases
  dat$pi_it <- expit(beta[1] + dat$gender*beta[2] + dat$time*beta[3] + c.vals)
  dat$U_it <- rbinom(nrow(dat), 1, dat$pi_it)
  
  # Log-normal parameters
  dat$mu_it <- gamma[1] + dat$gender*gamma[2] + dat$time*gamma[3] + d.vals 
  
  dat$Z_it <- rnorm(nrow(dat), dat$mu_it, ln_sigma)
  dat$Y_it <- exp(dat$Z_it)
  
  # Set value to 0 depending on eta_ij
  dat$cost[dat$U_it==1] <- round(dat$Y_it[dat$U_it==1],2)
  
  # Calculate cumulative sum
  dat <- dat %>% group_by(id) %>% mutate(csum=cumsum(cost)) %>% ungroup()
  
  return(dat)
}

dat <- genData(beta = c(0,0,-1), gamma=c(0,0,-1.5),
               re_sigma1 = 1, re_sigma2 = 0, re_rho = 0.2,
               ln_sigma = 1)


# Data summaries per time point + all time points
wide_dat <- dat %>% select(id,time,cost) %>% pivot_wider(names_from = time, values_from = cost, names_prefix = "time")
summary(dat$cost)
summary(wide_dat$time1)
summary(wide_dat$time2)
summary(wide_dat$time3)

# Percentage of those with zero cost
table(dat$cost==0)[2] / nrow(dat)

# Average probability
summary(dat$pi_it)

# Filter for only positive cost
pos_dat <- dat %>% filter(cost>0)
summary(pos_dat$cost)
hist(pos_dat$cost)

# Sampled subset of 100 people
sub_dat <- dat %>% filter(id %in% sample(dat$id, 100))
p <- ggplot(data=sub_dat, aes(x=time, y=cost, group=id))
p + geom_line() 

cumu_p <- ggplot(data=sub_dat, aes(x=time, y=csum, group=id))
cumu_p + geom_line() + xlab("Time") + ylab("Cumulative cost") + ggtitle("Cumulative cost trajectories (n=100)")



