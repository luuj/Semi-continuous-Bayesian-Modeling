### Bayesian Linear Regression - Main Script ####
library(dplyr)
library(lme4)
library(Rcpp)
require(RcppArmadillo)
sourceCpp("BLR_func.cpp")
source("Helper Files/BLR_Data_Gen.R")





## Unclustered simulation
set.seed(4)
dat <- gen_data(param=c(14.24,0.22,-0.07,-0.44), x=design_matrix, n_j=n, J=1, sigma_v=1)
dat <- dat %>% select(-cluster, -V_j)

# Univariate logistic regression
lr <- glm(y ~ trt + age + bmi, family="binomial", data=dat)
summary(lr)

# Initialize parameters
n_iter <- 100000
burn_in <- 10000
beta.init <- rep(0,ncol(dat))
jump_sigma <- c(0.25,0.2,0.01,0.01)
y <- as.matrix(dat %>% select(y))
x <- cbind(1,as.matrix(dat %>% select(-y)))

# Run algorithm
run1 <- MHu(y, x, beta.init, jump_sigma, n_iter, burn_in)
run2 <- MHu(y, x, beta.init+runif(1,-0.5,0.5), jump_sigma, n_iter, burn_in)
run3 <- MHu(y, x, beta.init+runif(1,-1,1), jump_sigma, n_iter, burn_in)

# Compare coefficient estimates
data.frame(LR=coefficients(lr), MCMC1=apply(run1,2,mean), 
           MCMC2=apply(run2,2,mean), MCMC3=apply(run3,2,mean))

# Trace plots
par(mfrow=c(2,2))
plot(1:n_iter, run1[,1], type="l", xlab="Iteration", ylab="Beta_hat", main="Intercept")
lines(1:n_iter, run2[,1], col="red")
lines(1:n_iter, run3[,1], col="blue")
plot(1:n_iter, run1[,2], type="l", xlab="Iteration", ylab="Beta_hat", main="Trt")
lines(1:n_iter, run2[,2], col="red")
lines(1:n_iter, run3[,2], col="blue")
plot(1:n_iter, run1[,3], type="l", xlab="Iteration", ylab="Beta_hat", main="Age")
lines(1:n_iter, run2[,3], col="red",)
lines(1:n_iter, run3[,3], col="blue")
plot(1:n_iter, run1[,4], type="l", xlab="Iteration", ylab="Beta_hat", main="BMI")
lines(1:n_iter, run2[,4], col="red")
lines(1:n_iter, run3[,4], col="blue")

# Potential scale reduction - rule of thumb is 1.1
library(coda)
ptr <- mcmc.list(mcmc(run1), mcmc(run2), mcmc(run3))
gelman.diag(ptr)





### Clustered simulation
set.seed(4)
dat.cl <- gen_data(param=c(14.24,0.22,-0.07,-0.44), x=design_matrix, n_j=n_j, J=J, sigma_v=0.5)
stored.vk <- unique(dat.cl$V_j)
dat.cl <- dat.cl %>% select(-V_j)

# GLMM
glmm <- glmer(y ~ trt + age + bmi + (1 | cluster), family="binomial", data=dat.cl)
summary(glmm)

# Initialize parameters
n_iter <- 70000
burn_in <- 3000
theta.init <- c(rep(0,ncol(dat.cl)+J-1),0.5)
gamma.init <- c(1,1)
jump_sigma <- c(0.25,0.2,0.01,0.01,rep(0.2, J),0.2)
y <- as.matrix(dat.cl %>% select(y))
x <- cbind(1,as.matrix(dat.cl %>% select(-y)))


# Run test
set.seed(3)
run1 <- MH.c(y, x, theta.init, gamma.init, jump_sigma, n_iter, burn_in)

# Coefficient comparison
data.frame(GLMM=c(fixef(glmm),stored.vk,0.5), MCMC1=apply(run1,2,mean))

# Trace plots
par(mfrow=c(2,2))
plot(1:n_iter, run1[,1], type="l", xlab="Iteration", ylab="Beta_hat", main="Intercept")
plot(1:n_iter, run1[,2], type="l", xlab="Iteration", ylab="Beta_hat", main="Trt")
plot(1:n_iter, run1[,3], type="l", xlab="Iteration", ylab="Beta_hat", main="Age")
plot(1:n_iter, run1[,4], type="l", xlab="Iteration", ylab="Beta_hat", main="BMI")

par(mfrow=c(2,2))
plot(1:n_iter, run1[,5], type="l", xlab="Iteration", ylab="Vk_hat", main="Vk1")
plot(1:n_iter, run1[,6], type="l", xlab="Iteration", ylab="Vk_hat", main="Vk2")
plot(1:n_iter, run1[,7], type="l", xlab="Iteration", ylab="Vk_hat", main="Vk3")
plot(1:n_iter, run1[,8], type="l", xlab="Iteration", ylab="Vk_hat", main="Vk4")






source("BLR_func.R")
set.seed(3)
run1 <- MH.c(y, x, theta.init, gamma.init, jump_sigma, n_iter, burn_in)



