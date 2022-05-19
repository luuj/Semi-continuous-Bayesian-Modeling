library(dplyr)
library(lme4)
source("BLR_func.R")
source("BLR_func_cluster.R")

## Generate simple data set
# param = beta coefficients
# n = number of observations
# n_cluster = number of clusters
# sigma_v = variance of random effect distribution

# Keep design matrix outside of data generator
age <- sample(20:70,n,replace=TRUE)
trt <- rbinom(n,1,0.5)
bmi <- rnorm(n,28,4)
x <- data.frame(1,trt,age,bmi)

# J is the number of clusters
# n_j is the number of observations per cluster
# Could add in cluster specific covariate like sex/trt
gen_data <- function(param, x, n, n_cluster, sigma_v){
  # Generate covariates
  cluster <- rep(1:n_cluster, each=n/n_cluster)
  V_j <- 0
  
  # Generate outcome
  if (n_cluster > 1){
    V_j <- rep(rnorm(n_cluster,0,sigma_v), each=n/n_cluster)
    pi <- expit(as.matrix(x)%*%as.matrix(param) + V_j) # Random intercept
  }else{
    pi <- expit(as.matrix(x)%*%as.matrix(param))
  }
  
  y <- rbinom(n,1,pi)
  
  data.frame(y,trt,age,bmi,cluster,V_j)
}





## Unclustered simulation
set.seed(3)
dat <- gen_data(param=c(14.24,0.22,-0.07,-0.44), n=1000, n_cluster=1, sigma_v=3)
dat <- dat %>% select(-cluster, -vk)

# Univariate logistic regression
lr <- glm(y ~ trt + age + bmi, family="binomial", data=dat)
summary(lr)

# Initialize parameters
num_iter <- 70000
burn_in <- 5000
beta.init <- rep(0,ncol(dat))
jump_sigma <- c(0.25,0.2,0.01,0.01)
y <- dat %>% select(y)
x <- dat %>% select(-y)

# Run test
set.seed(3)
run1 <- MH.u(y, x, beta.init, jump_sigma, num_iter, burn_in)
run2 <- MH.u(y, x, beta.init+runif(1,-0.5,0.5), jump_sigma, num_iter, burn_in)
run3 <- MH.u(y, x, beta.init+runif(1,-1,1), jump_sigma, num_iter, burn_in)

# Compare coefficient estimates
coef1 <- run1$coef
coef2 <- run2$coef
coef3 <- run3$coef
data.frame(LR=coefficients(lr), MCMC1=apply(coef1,2,mean), 
           MCMC2=apply(coef2,2,mean), MCMC3=apply(coef3,2,mean))

# Trace plots
par(mfrow=c(2,2))
plot(1:num_iter, coef1[,1], type="l", xlab="Iteration", ylab="Beta_hat", main="Intercept")
lines(1:num_iter, coef2[,1], col="red")
lines(1:num_iter, coef3[,1], col="blue")
plot(1:num_iter, coef1[,2], type="l", xlab="Iteration", ylab="Beta_hat", main="Trt")
lines(1:num_iter, coef2[,2], col="red")
lines(1:num_iter, coef3[,2], col="blue")
plot(1:num_iter, coef1[,3], type="l", xlab="Iteration", ylab="Beta_hat", main="Age")
lines(1:num_iter, coef2[,3], col="red",)
lines(1:num_iter, coef3[,3], col="blue")
plot(1:num_iter, coef1[,4], type="l", xlab="Iteration", ylab="Beta_hat", main="BMI")
lines(1:num_iter, coef2[,4], col="red")
lines(1:num_iter, coef3[,4], col="blue")

# Acceptance probabilities
run1$accept[1,]/run1$accept[2,]

# Potential scale reduction - rule of thumb is 1.1
ptr <- mcmc.list(mcmc(coef1), mcmc(coef2), mcmc(coef3))
gelman.diag(ptr)





### Clustered simulation
set.seed(3)
num_cluster <- 5
dat.cl <- gen_data(param=c(14.24,0.22,-0.07,-0.44), n=1000, n_cluster=num_cluster, sigma_v=0.5)
glmer(y ~ trt + age + bmi + (1 | cluster), family="binomial", data=dat.cl)

# Initialize parameters
dat.cl <- dat.cl %>% select(-vk)

num_iter <- 70000
burn_in <- 5000
theta.init <- rep(0,ncol(dat.cl)+num_cluster)
gamma.init <- c(1,2)
jump_sigma <- c(0.25,0.2,0.01,0.01,rep(0.2, num_cluster),0.2)
y <- dat.cl %>% select(y)
x <- dat.cl %>% select(-y)

# Run test
source("BLR_func_cluster.R")
set.seed(3)
run1 <- MH.c(y, x, theta.init, gamma.init, jump_sigma, num_iter, burn_in)

rgamma(10000,0.5,0.5)

# Try to analytically derive the full conditional is just a gamma
# Complete the square to
