library(sn)
library(e1071)
library(flexsurv)

# sample size 200
set.seed(122615)
n<-200
nsims <- 1000
totalN<-n*nsims

simnum <- c(rep(-1,totalN))
for(j in 1:nsims) {
  s <- (n*(j-1)+1)
  f <- n*j
  simnum[s:f] <- c(rep(j,n))
}

x1 <- rbinom(totalN,1,.5)
x2 <- rnorm(totalN, mean=0, sd=.64)
#x3 <- rlnorm(totalN,meanlog=0,sdlog=.5)
x3 <- rpois(totalN,1)
beta0 <- c(rep(6,totalN))
beta1 <- c(rep(0.2,totalN))
beta2 <- c(rep(-0.01,totalN))
beta3 <- c(rep(0.05,totalN))
alpha0 <- c(rep(3,totalN))
alpha1 <- c(rep(-4,totalN))
alpha2 <- c(rep(3.5,totalN))
alpha3 <- c(rep(2.5,totalN))
nu <- exp(beta0+beta1*x1+beta2*x2+beta3*x3) # marginal mean
pi <- exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)/(1+exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)) # P(y_i>0)
sigma<-1.2
eta<-2.5
Q <- 1/sqrt(eta)
mu <- beta0+beta1*x1+beta2*x2+beta3*x3-log(pi)-(sigma*log(Q^2))/Q-log(gamma(1/Q^2+sigma/Q))+log(gamma(1/Q^2))
y <- rgengamma(totalN, mu=mu, sigma=sigma, Q=Q) 
# add zero values
posval<-rbinom(totalN,1,pi)
y[posval==0]<-0
ygtz<-c(rep(1,totalN))
ygtz[y==0]<-0

mean(y)
max(y)
skewness(y)
mydatak2 <- data.frame(simnum, x1, x2, x3,beta0,beta1,beta2, beta3,alpha0, alpha1, alpha2, alpha3,
                       nu,pi,mu,y,ygtz)
write.csv(mydatak2, "C:/Users/Valerie/Documents/BIOS General/Dissertation/Paper 3/Simulation Data/zeros20GG/GG samp size 200 20p zeros 021215.csv")


# sample size 1000
set.seed(122615)
n<-1000
nsims <- 1000
totalN<-n*nsims

simnum <- c(rep(-1,totalN))
for(j in 1:nsims) {
  s <- (n*(j-1)+1)
  f <- n*j
  simnum[s:f] <- c(rep(j,n))
}

x1 <- rbinom(totalN,1,.5)
x2 <- rnorm(totalN, mean=0, sd=.64)
#x3 <- rlnorm(totalN,meanlog=0,sdlog=.5)
x3 <- rpois(totalN,1)
beta0 <- c(rep(6,totalN))
beta1 <- c(rep(0.2,totalN))
beta2 <- c(rep(-0.01,totalN))
beta3 <- c(rep(0.05,totalN))
alpha0 <- c(rep(3,totalN))
alpha1 <- c(rep(-4,totalN))
alpha2 <- c(rep(3.5,totalN))
alpha3 <- c(rep(2.5,totalN))
nu <- exp(beta0+beta1*x1+beta2*x2+beta3*x3) # marginal mean
pi <- exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)/(1+exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)) # P(y_i>0)
sigma<-1.2
eta<-2.5
Q <- 1/sqrt(eta)
mu <- beta0+beta1*x1+beta2*x2+beta3*x3-log(pi)-(sigma*log(Q^2))/Q-log(gamma(1/Q^2+sigma/Q))+log(gamma(1/Q^2))
y <- rgengamma(totalN, mu=mu, sigma=sigma, Q=Q) 
# add zero values
posval<-rbinom(totalN,1,pi)
y[posval==0]<-0
ygtz<-c(rep(1,totalN))
ygtz[y==0]<-0

mean(y)
max(y)
skewness(y)
mydatak2 <- data.frame(simnum, x1, x2, x3,beta0,beta1,beta2, beta3,alpha0, alpha1, alpha2, alpha3,
                       nu,pi,mu,y,ygtz)
write.csv(mydatak2, "C:/Users/Valerie/Documents/BIOS General/Dissertation/Paper 3/Simulation Data/zeros20GG/GG samp size 1000 20p zeros 021215.csv")



# sample size 10000
set.seed(122615)
n<-10000
nsims <- 1000
totalN<-n*nsims

simnum <- c(rep(-1,totalN))
for(j in 1:nsims) {
  s <- (n*(j-1)+1)
  f <- n*j
  simnum[s:f] <- c(rep(j,n))
}

x1 <- rbinom(totalN,1,.5)
x2 <- rnorm(totalN, mean=0, sd=.64)
#x3 <- rlnorm(totalN,meanlog=0,sdlog=.5)
x3 <- rpois(totalN,1)
beta0 <- c(rep(6,totalN))
beta1 <- c(rep(0.2,totalN))
beta2 <- c(rep(-0.01,totalN))
beta3 <- c(rep(0.05,totalN))
alpha0 <- c(rep(3,totalN))
alpha1 <- c(rep(-4,totalN))
alpha2 <- c(rep(3.5,totalN))
alpha3 <- c(rep(2.5,totalN))
nu <- exp(beta0+beta1*x1+beta2*x2+beta3*x3) # marginal mean
pi <- exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)/(1+exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)) # P(y_i>0)
sigma<-1.2
eta<-2.5
Q <- 1/sqrt(eta)
mu <- beta0+beta1*x1+beta2*x2+beta3*x3-log(pi)-(sigma*log(Q^2))/Q-log(gamma(1/Q^2+sigma/Q))+log(gamma(1/Q^2))
y <- rgengamma(totalN, mu=mu, sigma=sigma, Q=Q) 
# add zero values
posval<-rbinom(totalN,1,pi)
y[posval==0]<-0
ygtz<-c(rep(1,totalN))
ygtz[y==0]<-0

mean(y)
max(y)
skewness(y)
mydatak2 <- data.frame(simnum, x1, x2, x3,y,ygtz)
write.csv(mydatak2, "C:/Users/Valerie/Documents/BIOS General/Dissertation/Paper 3/Simulation Data/zeros20GG/GG samp size 10000 20p zeros 021215.csv")

# sample size 50000

# first 500
set.seed(122615)
n<-50000
nsims <- 500
totalN<-n*nsims

simnum <- c(rep(-1,totalN))
for(j in 1:nsims) {
  s <- (n*(j-1)+1)
  f <- n*j
  simnum[s:f] <- c(rep(j,n))
}

x1 <- rbinom(totalN,1,.5)
x2 <- rnorm(totalN, mean=0, sd=.64)
#x3 <- rlnorm(totalN,meanlog=0,sdlog=.5)
x3 <- rpois(totalN,1)
beta0 <- c(rep(6,totalN))
beta1 <- c(rep(0.2,totalN))
beta2 <- c(rep(-0.01,totalN))
beta3 <- c(rep(0.05,totalN))
alpha0 <- c(rep(3,totalN))
alpha1 <- c(rep(-4,totalN))
alpha2 <- c(rep(3.5,totalN))
alpha3 <- c(rep(2.5,totalN))
nu <- exp(beta0+beta1*x1+beta2*x2+beta3*x3) # marginal mean
pi <- exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)/(1+exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)) # P(y_i>0)
sigma<-1.2
eta<-2.5
Q <- 1/sqrt(eta)
mu <- beta0+beta1*x1+beta2*x2+beta3*x3-log(pi)-(sigma*log(Q^2))/Q-log(gamma(1/Q^2+sigma/Q))+log(gamma(1/Q^2))
y <- rgengamma(totalN, mu=mu, sigma=sigma, Q=Q) 
# add zero values
posval<-rbinom(totalN,1,pi)
y[posval==0]<-0
ygtz<-c(rep(1,totalN))
ygtz[y==0]<-0

mean(y)
max(y)
skewness(y)
mydatak2 <- data.frame(simnum, x1, x2, x3,y,ygtz)
write.csv(mydatak2, "C:/Users/Valerie/Documents/BIOS General/Dissertation/Paper 3/Simulation Data/For Paper 3/zeros20GGeta25/GG samp size 50000 20p zeros part 1 062315.csv")

rm(list=ls())
# 2nd 500
set.seed(062315)
n<-50000
nsims <- 500
totalN<-n*nsims

simnum <- c(rep(-1,totalN))
for(j in 1:nsims) {
  s <- (n*(j-1)+1)
  f <- n*j
  simnum[s:f] <- c(rep(j,n))
}

x1 <- rbinom(totalN,1,.5)
x2 <- rnorm(totalN, mean=0, sd=.64)
#x3 <- rlnorm(totalN,meanlog=0,sdlog=.5)
x3 <- rpois(totalN,1)
beta0 <- c(rep(6,totalN))
beta1 <- c(rep(0.2,totalN))
beta2 <- c(rep(-0.01,totalN))
beta3 <- c(rep(0.05,totalN))
alpha0 <- c(rep(3,totalN))
alpha1 <- c(rep(-4,totalN))
alpha2 <- c(rep(3.5,totalN))
alpha3 <- c(rep(2.5,totalN))
nu <- exp(beta0+beta1*x1+beta2*x2+beta3*x3) # marginal mean
pi <- exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)/(1+exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)) # P(y_i>0)
sigma<-1.2
eta<-2.5
Q <- 1/sqrt(eta)
mu <- beta0+beta1*x1+beta2*x2+beta3*x3-log(pi)-(sigma*log(Q^2))/Q-log(gamma(1/Q^2+sigma/Q))+log(gamma(1/Q^2))
y <- rgengamma(totalN, mu=mu, sigma=sigma, Q=Q) 
# add zero values
posval<-rbinom(totalN,1,pi)
y[posval==0]<-0
ygtz<-c(rep(1,totalN))
ygtz[y==0]<-0

mean(y)
max(y)
skewness(y)
mydatak2 <- data.frame(simnum, x1, x2, x3,y,ygtz)
write.csv(mydatak2, "C:/Users/Valerie/Documents/BIOS General/Dissertation/Paper 3/Simulation Data/For Paper 3/zeros20GGeta25/GG samp size 50000 20p zeros part 2 062315.csv")


############## SET BETA1=0 TO GET TYPE 1 ERRORS ###############

# sample size 200
set.seed(122615)
n<-200
nsims <- 1000
totalN<-n*nsims

simnum <- c(rep(-1,totalN))
for(j in 1:nsims) {
  s <- (n*(j-1)+1)
  f <- n*j
  simnum[s:f] <- c(rep(j,n))
}

x1 <- rbinom(totalN,1,.5)
x2 <- rnorm(totalN, mean=0, sd=.64)
#x3 <- rlnorm(totalN,meanlog=0,sdlog=.5)
x3 <- rpois(totalN,1)
beta0 <- c(rep(6,totalN))
beta1 <- c(rep(0,totalN))
beta2 <- c(rep(-0.01,totalN))
beta3 <- c(rep(0.05,totalN))
alpha0 <- c(rep(3,totalN))
alpha1 <- c(rep(-4,totalN))
alpha2 <- c(rep(3.5,totalN))
alpha3 <- c(rep(2.5,totalN))
nu <- exp(beta0+beta1*x1+beta2*x2+beta3*x3) # marginal mean
pi <- exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)/(1+exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)) # P(y_i>0)
sigma<-1.2
eta<-2.5
Q <- 1/sqrt(eta)
mu <- beta0+beta1*x1+beta2*x2+beta3*x3-log(pi)-(sigma*log(Q^2))/Q-log(gamma(1/Q^2+sigma/Q))+log(gamma(1/Q^2))
y <- rgengamma(totalN, mu=mu, sigma=sigma, Q=Q) 
# add zero values
posval<-rbinom(totalN,1,pi)
y[posval==0]<-0
ygtz<-c(rep(1,totalN))
ygtz[y==0]<-0

mean(y)
max(y)
skewness(y)
mydatak2 <- data.frame(simnum, x1, x2, x3,beta0,beta1,beta2, beta3,alpha0, alpha1, alpha2, alpha3,
                       nu,pi,mu,y,ygtz)
write.csv(mydatak2, "C:/Users/Valerie/Documents/BIOS General/Dissertation/Paper 3/Simulation Data/zeros20GGeta25/GG samp size 200 20p zeros type 1 error 021215.csv")


# sample size 1000
set.seed(122615)
n<-1000
nsims <- 1000
totalN<-n*nsims

simnum <- c(rep(-1,totalN))
for(j in 1:nsims) {
  s <- (n*(j-1)+1)
  f <- n*j
  simnum[s:f] <- c(rep(j,n))
}

x1 <- rbinom(totalN,1,.5)
x2 <- rnorm(totalN, mean=0, sd=.64)
#x3 <- rlnorm(totalN,meanlog=0,sdlog=.5)
x3 <- rpois(totalN,1)
beta0 <- c(rep(6,totalN))
beta1 <- c(rep(0,totalN))
beta2 <- c(rep(-0.01,totalN))
beta3 <- c(rep(0.05,totalN))
alpha0 <- c(rep(3,totalN))
alpha1 <- c(rep(-4,totalN))
alpha2 <- c(rep(3.5,totalN))
alpha3 <- c(rep(2.5,totalN))
nu <- exp(beta0+beta1*x1+beta2*x2+beta3*x3) # marginal mean
pi <- exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)/(1+exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)) # P(y_i>0)
sigma<-1.2
eta<-2.5
Q <- 1/sqrt(eta)
mu <- beta0+beta1*x1+beta2*x2+beta3*x3-log(pi)-(sigma*log(Q^2))/Q-log(gamma(1/Q^2+sigma/Q))+log(gamma(1/Q^2))
y <- rgengamma(totalN, mu=mu, sigma=sigma, Q=Q) 
# add zero values
posval<-rbinom(totalN,1,pi)
y[posval==0]<-0
ygtz<-c(rep(1,totalN))
ygtz[y==0]<-0

mean(y)
max(y)
skewness(y)
mydatak2 <- data.frame(simnum, x1, x2, x3,beta0,beta1,beta2, beta3,alpha0, alpha1, alpha2, alpha3,
                       nu,pi,mu,y,ygtz)
write.csv(mydatak2, "C:/Users/Valerie/Documents/BIOS General/Dissertation/Paper 3/Simulation Data/zeros20GGeta25/GG samp size 1000 20p zeros type 1 error 021215.csv")



# sample size 10000
set.seed(122615)
n<-10000
nsims <- 1000
totalN<-n*nsims

simnum <- c(rep(-1,totalN))
for(j in 1:nsims) {
  s <- (n*(j-1)+1)
  f <- n*j
  simnum[s:f] <- c(rep(j,n))
}

x1 <- rbinom(totalN,1,.5)
x2 <- rnorm(totalN, mean=0, sd=.64)
#x3 <- rlnorm(totalN,meanlog=0,sdlog=.5)
x3 <- rpois(totalN,1)
beta0 <- c(rep(6,totalN))
beta1 <- c(rep(0,totalN))
beta2 <- c(rep(-0.01,totalN))
beta3 <- c(rep(0.05,totalN))
alpha0 <- c(rep(3,totalN))
alpha1 <- c(rep(-4,totalN))
alpha2 <- c(rep(3.5,totalN))
alpha3 <- c(rep(2.5,totalN))
nu <- exp(beta0+beta1*x1+beta2*x2+beta3*x3) # marginal mean
pi <- exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)/(1+exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)) # P(y_i>0)
sigma<-1.2
eta<-2.5
Q <- 1/sqrt(eta)
mu <- beta0+beta1*x1+beta2*x2+beta3*x3-log(pi)-(sigma*log(Q^2))/Q-log(gamma(1/Q^2+sigma/Q))+log(gamma(1/Q^2))
y <- rgengamma(totalN, mu=mu, sigma=sigma, Q=Q) 
# add zero values
posval<-rbinom(totalN,1,pi)
y[posval==0]<-0
ygtz<-c(rep(1,totalN))
ygtz[y==0]<-0

mean(y)
max(y)
skewness(y)
mydatak2 <- data.frame(simnum, x1, x2, x3,y,ygtz)
write.csv(mydatak2, "C:/Users/Valerie/Documents/BIOS General/Dissertation/Paper 3/Simulation Data/zeros20GGeta25/GG samp size 10000 20p zeros type 1 error 021215.csv")

# sample size 50000
# 1st 500
set.seed(122615)
n<-50000
nsims <- 500
totalN<-n*nsims

simnum <- c(rep(-1,totalN))
for(j in 1:nsims) {
  s <- (n*(j-1)+1)
  f <- n*j
  simnum[s:f] <- c(rep(j,n))
}

x1 <- rbinom(totalN,1,.5)
x2 <- rnorm(totalN, mean=0, sd=.64)
#x3 <- rlnorm(totalN,meanlog=0,sdlog=.5)
x3 <- rpois(totalN,1)
beta0 <- c(rep(6,totalN))
beta1 <- c(rep(0,totalN))
beta2 <- c(rep(-0.01,totalN))
beta3 <- c(rep(0.05,totalN))
alpha0 <- c(rep(3,totalN))
alpha1 <- c(rep(-4,totalN))
alpha2 <- c(rep(3.5,totalN))
alpha3 <- c(rep(2.5,totalN))
nu <- exp(beta0+beta1*x1+beta2*x2+beta3*x3) # marginal mean
pi <- exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)/(1+exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)) # P(y_i>0)
sigma<-1.2
eta<-2.5
Q <- 1/sqrt(eta)
mu <- beta0+beta1*x1+beta2*x2+beta3*x3-log(pi)-(sigma*log(Q^2))/Q-log(gamma(1/Q^2+sigma/Q))+log(gamma(1/Q^2))
y <- rgengamma(totalN, mu=mu, sigma=sigma, Q=Q) 
# add zero values
posval<-rbinom(totalN,1,pi)
y[posval==0]<-0
ygtz<-c(rep(1,totalN))
ygtz[y==0]<-0

mean(y)
max(y)
skewness(y)
mydatak2 <- data.frame(simnum, x1, x2, x3,y,ygtz)
write.csv(mydatak2, "C:/Users/Valerie/Documents/BIOS General/Dissertation/Paper 3/Simulation Data/For Paper 3/zeros20GGeta25/GG samp size 50000 20p zeros type 1 error part 1 062315.csv")

# 2nd 500
set.seed(062315)
n<-50000
nsims <- 500
totalN<-n*nsims

simnum <- c(rep(-1,totalN))
for(j in 1:nsims) {
  s <- (n*(j-1)+1)
  f <- n*j
  simnum[s:f] <- c(rep(j,n))
}

x1 <- rbinom(totalN,1,.5)
x2 <- rnorm(totalN, mean=0, sd=.64)
#x3 <- rlnorm(totalN,meanlog=0,sdlog=.5)
x3 <- rpois(totalN,1)
beta0 <- c(rep(6,totalN))
beta1 <- c(rep(0,totalN))
beta2 <- c(rep(-0.01,totalN))
beta3 <- c(rep(0.05,totalN))
alpha0 <- c(rep(3,totalN))
alpha1 <- c(rep(-4,totalN))
alpha2 <- c(rep(3.5,totalN))
alpha3 <- c(rep(2.5,totalN))
nu <- exp(beta0+beta1*x1+beta2*x2+beta3*x3) # marginal mean
pi <- exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)/(1+exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)) # P(y_i>0)
sigma<-1.2
eta<-2.5
Q <- 1/sqrt(eta)
mu <- beta0+beta1*x1+beta2*x2+beta3*x3-log(pi)-(sigma*log(Q^2))/Q-log(gamma(1/Q^2+sigma/Q))+log(gamma(1/Q^2))
y <- rgengamma(totalN, mu=mu, sigma=sigma, Q=Q) 
# add zero values
posval<-rbinom(totalN,1,pi)
y[posval==0]<-0
ygtz<-c(rep(1,totalN))
ygtz[y==0]<-0

mean(y)
max(y)
skewness(y)
mydatak2 <- data.frame(simnum, x1, x2, x3,y,ygtz)
write.csv(mydatak2, "C:/Users/Valerie/Documents/BIOS General/Dissertation/Paper 3/Simulation Data/For Paper 3/zeros20GGeta25/GG samp size 50000 20p zeros type 1 error part 2 062315.csv")
