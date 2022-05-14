library(dplyr)
library(coda)
library(lme4)
source("BLR_func.R")
source("BLR_func_cluster.R")



### Unclustered simulation
## Generate simple data set
# param = beta coefficients
# n = number of observations
# n_cluster = number of clusters
# sigma_v = variance of random effect distribution
gen_data <- function(param, n, n_cluster, sigma_v){
  # Generate covariates
  age <- sample(20:70,n,replace=TRUE)
  trt <- rbinom(n,1,0.5)
  bmi <- rnorm(n,28,4)
  x <- data.frame(1,trt,age,bmi)
  cluster <- rep(1:n_cluster, each=n/n_cluster)
  
  # Generate outcome
  if (n_cluster > 1){
    vk <- rnorm(n_cluster,0,sigma_v)
    pi <- expit(as.matrix(x)%*%as.matrix(param) + rep(vk, each=n/n_cluster)) # Random intercept
  }else{
    pi <- expit(as.matrix(x)%*%as.matrix(param))
  }
  
  y <- rbinom(n,1,pi)
  
  data.frame(y,trt,age,bmi,cluster)
}

set.seed(3)
dat <- gen_data(param=c(14.24,0.22,-0.07,-0.44), n=1000, n_cluster=1, sigma_v=3)
dat <- dat %>% select(-cluster)

# Univariate logistic regression
lr <- glm(y ~ trt + age + bmi, family="binomial", data=dat)
summary(lr)

# Initialize parameters
num_iter <- 30000
burn_in <- 5000
num_var <- ncol(dat)
theta.init <- rep(0,num_var)
y <- as.matrix(dat %>% select(y))
x <- as.matrix(dat %>% select(-y))
sigma.prior <- rep(1,num_var)
sigma.jump <- c(0.2,0.2,0.004,0.01)

# Run test
set.seed(3)
run1 <- MH(y, x, theta.init, sigma.prior, sigma.jump, num_iter, burn_in)
run2 <- MH(y, x, theta.init+runif(1,-0.5,0.5), sigma.prior, sigma.jump, num_iter, burn_in)
run3 <- MH(y, x, theta.init+runif(1,-0.2,0.2), sigma.prior, sigma.jump, num_iter, burn_in)

#my.theta <- seq(from=-5, to 20, length=1000)
#my.logPost <- rep(NA, 1000)
#for(i in 1:1000) my.logPost[i] <- log_post(c(my.theta[i], 0, 0, 0),
#                                           y, x, ,sigma.prior)

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

par(mfrow=c(2,2))
hist(c(coef1[,1],coef2[,1],coef3[,1]))
hist(c(coef1[,2],coef2[,2],coef3[,2]))
hist(c(coef1[,3],coef2[,3],coef3[,3]))
hist(c(coef1[,4],coef2[,4],coef3[,4]))





### Clustered simulation
set.seed(3)
dat.cl <- gen_data(param=c(14.24,0.22,-0.07,-0.44), n=1000, n_cluster=10, sigma_v=3)
glmer(y ~ trt + age + bmi + (1 | cluster), family="binomial", data=dat.cl)



