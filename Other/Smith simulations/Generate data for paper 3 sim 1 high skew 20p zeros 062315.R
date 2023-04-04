library(sn)
library(e1071)

set.seed(122615)

# LSN data with kappa=0.5, sample size=200
# This has pi=0.296229, 20.4% zeros
# mean of y is 7752.807
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
kappa <- c(rep(5,totalN))  # using kappa=alpha in standard SN parameterization to avoid confusion with parameters we're using
delta <- kappa/sqrt(1+kappa^2)
omega <- c(rep(1.2,totalN))
omega2 <- omega*omega
nu <- exp(beta0+beta1*x1+beta2*x2+beta3*x3) # marginal mean
pi <- exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)/(1+exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)) # P(y_i>0)
xi <- beta0+beta1*x1+beta2*x2+beta3*x3-log(2)-log(pi)-log(pnorm(omega*delta))-omega2/2  # SN location parameter
ypre <- rsn(totalN, xi, omega, kappa)  # ypre ~ skew normal
y <- exp(ypre)  # y ~ log skew normal

# add zero values
posval<-rbinom(totalN,1,pi)
y[posval==0]<-0
ygtz<-c(rep(1,totalN))
ygtz[y==0]<-0

plot(density(ypre), type="l", col="red")
plot(density(y), type="l", col="red")
skewness(ypre)

#skewtheor <- (4-pi)/2*(delta*sqrt(2/pi))^3/(1-2*delta^2/pi)^(3/2)
#skewtheor[1]

# write to file to import into SAS
mydatak2 <- data.frame(simnum, x1, x2, x3,beta0,beta1,beta2, beta3,alpha0, alpha1, alpha2, alpha3,
	kappa,delta,omega,nu,pi,xi,y,yzero,ygtz)
write.csv(mydatak2, "C:/Users/Valerie/Documents/BIOS General/Dissertation/Paper 3/Simulation Data/zeros20highskew/LSN high skew samp size 200 20p zeros 020915.csv")



# Sample size 1000
set.seed(122615)

# LSN data with kappa=0.5, sample size=200
# This has pi=0.296229, 20.4% zeros
# mean of y is 7752.807
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
kappa <- c(rep(5,totalN))  # using kappa=alpha in standard SN parameterization to avoid confusion with parameters we're using
delta <- kappa/sqrt(1+kappa^2)
omega <- c(rep(1.2,totalN))
omega2 <- omega*omega
nu <- exp(beta0+beta1*x1+beta2*x2+beta3*x3) # marginal mean
pi <- exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)/(1+exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)) # P(y_i>0)
xi <- beta0+beta1*x1+beta2*x2+beta3*x3-log(2)-log(pi)-log(pnorm(omega*delta))-omega2/2  # SN location parameter
ypre <- rsn(totalN, xi, omega, kappa)  # ypre ~ skew normal
y <- exp(ypre)  # y ~ log skew normal

# add zero values
posval<-rbinom(totalN,1,pi)
y[posval==0]<-0
ygtz<-c(rep(1,totalN))
ygtz[y==0]<-0

plot(density(ypre), type="l", col="red")
plot(density(y), type="l", col="red")
skewness(ypre)

#skewtheor <- (4-pi)/2*(delta*sqrt(2/pi))^3/(1-2*delta^2/pi)^(3/2)
#skewtheor[1]

# write to file to import into SAS
mydatak2 <- data.frame(simnum, x1, x2, x3,beta0,beta1,beta2, beta3,alpha0, alpha1, alpha2, alpha3,
                       kappa,delta,omega,nu,pi,xi,y,yzero,ygtz)
write.csv(mydatak2, "C:/Users/Valerie/Documents/BIOS General/Dissertation/Paper 3/Simulation Data/zeros20highskew/LSN high skew samp size 1000 20p zeros 020915.csv")


# Sample size 10000
set.seed(122615)

# LSN data with kappa=0.5, sample size=200
# This has pi=0.296229, 20.4% zeros
# mean of y is 7752.807
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
#beta0 <- c(rep(7,totalN))
#beta1 <- c(rep(0.2,totalN))
#beta2 <- c(rep(-0.05,totalN))
#beta3 <- c(rep(0.05,totalN))
alpha0 <- c(rep(3,totalN))
alpha1 <- c(rep(-4,totalN))
alpha2 <- c(rep(3.5,totalN))
alpha3 <- c(rep(2.5,totalN))
kappa <- c(rep(5,totalN))  # using kappa=alpha in standard SN parameterization to avoid confusion with parameters we're using
delta <- kappa/sqrt(1+kappa^2)
omega <- c(rep(1.2,totalN))
omega2 <- omega*omega
nu <- exp(beta0+beta1*x1+beta2*x2+beta3*x3) # marginal mean
pi <- exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)/(1+exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)) # P(y_i>0)
xi <- beta0+beta1*x1+beta2*x2+beta3*x3-log(2)-log(pi)-log(pnorm(omega*delta))-omega2/2  # SN location parameter
ypre <- rsn(totalN, xi, omega, kappa)  # ypre ~ skew normal
y <- exp(ypre)  # y ~ log skew normal

# add zero values
posval<-rbinom(totalN,1,pi)
y[posval==0]<-0
ygtz<-c(rep(1,totalN))
ygtz[y==0]<-0

plot(density(ypre), type="l", col="red")
plot(density(y), type="l", col="red")
skewness(ypre)

#skewtheor <- (4-pi)/2*(delta*sqrt(2/pi))^3/(1-2*delta^2/pi)^(3/2)
#skewtheor[1]

# write to file to import into SAS
mydatak2 <- data.frame(simnum, x1, x2, x3,y,yzero,ygtz)
write.csv(mydatak2, "C:/Users/Valerie/Documents/BIOS General/Dissertation/Paper 3/Simulation Data/zeros20highskew/LSN high skew samp size 10000 20p zeros 020915.csv")


# Sample size 50000
# first 250
set.seed(122615)
memory.limit(16000)

# LSN data with kappa=0.5, sample size=50000

n<-50000
nsims <- 250
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
#beta0 <- c(rep(7,totalN))
#beta1 <- c(rep(0.2,totalN))
#beta2 <- c(rep(-0.05,totalN))
#beta3 <- c(rep(0.05,totalN))
alpha0 <- c(rep(3,totalN))
alpha1 <- c(rep(-4,totalN))
alpha2 <- c(rep(3.5,totalN))
alpha3 <- c(rep(2.5,totalN))
kappa <- c(rep(5,totalN))  # using kappa=alpha in standard SN parameterization to avoid confusion with parameters we're using
delta <- kappa/sqrt(1+kappa^2)
omega <- c(rep(1.2,totalN))
omega2 <- omega*omega
nu <- exp(beta0+beta1*x1+beta2*x2+beta3*x3) # marginal mean
pi <- exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)/(1+exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)) # P(y_i>0)
xi <- beta0+beta1*x1+beta2*x2+beta3*x3-log(2)-log(pi)-log(pnorm(omega*delta))-omega2/2  # SN location parameter
ypre <- rsn(totalN, xi, omega, kappa)  # ypre ~ skew normal
y <- exp(ypre)  # y ~ log skew normal

# add zero values
posval<-rbinom(totalN,1,pi)
y[posval==0]<-0
ygtz<-c(rep(1,totalN))
ygtz[y==0]<-0

plot(density(ypre), type="l", col="red")
plot(density(y), type="l", col="red")
skewness(ypre)
mean(y)
max(y)

#skewtheor <- (4-pi)/2*(delta*sqrt(2/pi))^3/(1-2*delta^2/pi)^(3/2)
#skewtheor[1]

# write to file to import into SAS
mydatak2 <- data.frame(simnum, x1, x2, x3,y)
write.csv(mydatak2, "C:/Users/Valerie/Documents/BIOS General/Dissertation/Paper 3/Simulation Data/For Paper 3/zeros20highskew/LSN high skew samp size 50000 20p zeros part 1 062315.csv")


# second 250
rm(list=ls())
set.seed(062315)


# LSN data with kappa=0.5, sample size=50000

n<-50000
nsims <- 250
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
#beta0 <- c(rep(7,totalN))
#beta1 <- c(rep(0.2,totalN))
#beta2 <- c(rep(-0.05,totalN))
#beta3 <- c(rep(0.05,totalN))
alpha0 <- c(rep(3,totalN))
alpha1 <- c(rep(-4,totalN))
alpha2 <- c(rep(3.5,totalN))
alpha3 <- c(rep(2.5,totalN))
kappa <- c(rep(5,totalN))  # using kappa=alpha in standard SN parameterization to avoid confusion with parameters we're using
delta <- kappa/sqrt(1+kappa^2)
omega <- c(rep(1.2,totalN))
omega2 <- omega*omega
nu <- exp(beta0+beta1*x1+beta2*x2+beta3*x3) # marginal mean
pi <- exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)/(1+exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)) # P(y_i>0)
xi <- beta0+beta1*x1+beta2*x2+beta3*x3-log(2)-log(pi)-log(pnorm(omega*delta))-omega2/2  # SN location parameter
ypre <- rsn(totalN, xi, omega, kappa)  # ypre ~ skew normal
y <- exp(ypre)  # y ~ log skew normal

# add zero values
posval<-rbinom(totalN,1,pi)
y[posval==0]<-0
ygtz<-c(rep(1,totalN))
ygtz[y==0]<-0

plot(density(ypre), type="l", col="red")
plot(density(y), type="l", col="red")
skewness(ypre)

#skewtheor <- (4-pi)/2*(delta*sqrt(2/pi))^3/(1-2*delta^2/pi)^(3/2)
#skewtheor[1]

# write to file to import into SAS
mydatak2 <- data.frame(simnum, x1, x2, x3,y)
write.csv(mydatak2, "C:/Users/Valerie/Documents/BIOS General/Dissertation/Paper 3/Simulation Data/For Paper 3/zeros20highskew/LSN high skew samp size 50000 20p zeros part 2 062315.csv")


# third 250
rm(list=ls())
set.seed(7654321)


# LSN data with kappa=0.5, sample size=50000

n<-50000
nsims <- 250
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
#beta0 <- c(rep(7,totalN))
#beta1 <- c(rep(0.2,totalN))
#beta2 <- c(rep(-0.05,totalN))
#beta3 <- c(rep(0.05,totalN))
alpha0 <- c(rep(3,totalN))
alpha1 <- c(rep(-4,totalN))
alpha2 <- c(rep(3.5,totalN))
alpha3 <- c(rep(2.5,totalN))
kappa <- c(rep(5,totalN))  # using kappa=alpha in standard SN parameterization to avoid confusion with parameters we're using
delta <- kappa/sqrt(1+kappa^2)
omega <- c(rep(1.2,totalN))
omega2 <- omega*omega
nu <- exp(beta0+beta1*x1+beta2*x2+beta3*x3) # marginal mean
pi <- exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)/(1+exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)) # P(y_i>0)
xi <- beta0+beta1*x1+beta2*x2+beta3*x3-log(2)-log(pi)-log(pnorm(omega*delta))-omega2/2  # SN location parameter
ypre <- rsn(totalN, xi, omega, kappa)  # ypre ~ skew normal
y <- exp(ypre)  # y ~ log skew normal

# add zero values
posval<-rbinom(totalN,1,pi)
y[posval==0]<-0
ygtz<-c(rep(1,totalN))
ygtz[y==0]<-0

plot(density(ypre), type="l", col="red")
plot(density(y), type="l", col="red")
skewness(ypre)

#skewtheor <- (4-pi)/2*(delta*sqrt(2/pi))^3/(1-2*delta^2/pi)^(3/2)
#skewtheor[1]

# write to file to import into SAS
mydatak2 <- data.frame(simnum, x1, x2, x3,y)
write.csv(mydatak2, "C:/Users/Valerie/Documents/BIOS General/Dissertation/Paper 3/Simulation Data/For Paper 3/zeros20highskew/LSN high skew samp size 50000 20p zeros part 3 062315.csv")

# last 250
rm(list=ls())
set.seed(54321895)


# LSN data with kappa=0.5, sample size=50000

n<-50000
nsims <- 250
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
#beta0 <- c(rep(7,totalN))
#beta1 <- c(rep(0.2,totalN))
#beta2 <- c(rep(-0.05,totalN))
#beta3 <- c(rep(0.05,totalN))
alpha0 <- c(rep(3,totalN))
alpha1 <- c(rep(-4,totalN))
alpha2 <- c(rep(3.5,totalN))
alpha3 <- c(rep(2.5,totalN))
kappa <- c(rep(5,totalN))  # using kappa=alpha in standard SN parameterization to avoid confusion with parameters we're using
delta <- kappa/sqrt(1+kappa^2)
omega <- c(rep(1.2,totalN))
omega2 <- omega*omega
nu <- exp(beta0+beta1*x1+beta2*x2+beta3*x3) # marginal mean
pi <- exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)/(1+exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)) # P(y_i>0)
xi <- beta0+beta1*x1+beta2*x2+beta3*x3-log(2)-log(pi)-log(pnorm(omega*delta))-omega2/2  # SN location parameter
ypre <- rsn(totalN, xi, omega, kappa)  # ypre ~ skew normal
y <- exp(ypre)  # y ~ log skew normal

# add zero values
posval<-rbinom(totalN,1,pi)
y[posval==0]<-0
ygtz<-c(rep(1,totalN))
ygtz[y==0]<-0

plot(density(ypre), type="l", col="red")
plot(density(y), type="l", col="red")
skewness(ypre)

#skewtheor <- (4-pi)/2*(delta*sqrt(2/pi))^3/(1-2*delta^2/pi)^(3/2)
#skewtheor[1]

# write to file to import into SAS
mydatak2 <- data.frame(simnum, x1, x2, x3,y)
write.csv(mydatak2, "C:/Users/Valerie/Documents/BIOS General/Dissertation/Paper 3/Simulation Data/For Paper 3/zeros20highskew/LSN high skew samp size 50000 20p zeros part 4 062315.csv")


#### Generate with beta1=0 to get type 1 error ####

set.seed(122615)

# LSN data with kappa=0.5, sample size=200
# This has pi=0.296229, 20.4% zeros
# mean of y is 7752.807
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
kappa <- c(rep(5,totalN))  # using kappa=alpha in standard SN parameterization to avoid confusion with parameters we're using
delta <- kappa/sqrt(1+kappa^2)
omega <- c(rep(1.2,totalN))
omega2 <- omega*omega
nu <- exp(beta0+beta1*x1+beta2*x2+beta3*x3) # marginal mean
pi <- exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)/(1+exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)) # P(y_i>0)
xi <- beta0+beta1*x1+beta2*x2+beta3*x3-log(2)-log(pi)-log(pnorm(omega*delta))-omega2/2  # SN location parameter
ypre <- rsn(totalN, xi, omega, kappa)  # ypre ~ skew normal
y <- exp(ypre)  # y ~ log skew normal

# add zero values
posval<-rbinom(totalN,1,pi)
y[posval==0]<-0
ygtz<-c(rep(1,totalN))
ygtz[y==0]<-0

plot(density(ypre), type="l", col="red")
plot(density(y), type="l", col="red")
skewness(ypre)

#skewtheor <- (4-pi)/2*(delta*sqrt(2/pi))^3/(1-2*delta^2/pi)^(3/2)
#skewtheor[1]

# write to file to import into SAS
mydatak2 <- data.frame(simnum, x1, x2, x3,beta0,beta1,beta2, beta3,alpha0, alpha1, alpha2, alpha3,
                       kappa,delta,omega,nu,pi,xi,y,yzero,ygtz)
write.csv(mydatak2, "C:/Users/Valerie/Documents/BIOS General/Dissertation/Paper 3/Simulation Data/zeros20highskew/LSN high skew samp size 200 20p zeros type 1 error 020915.csv")



# Sample size 1000
set.seed(122615)

# LSN data with kappa=0.5, sample size=200
# This has pi=0.296229, 20.4% zeros
# mean of y is 7752.807
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
kappa <- c(rep(5,totalN))  # using kappa=alpha in standard SN parameterization to avoid confusion with parameters we're using
delta <- kappa/sqrt(1+kappa^2)
omega <- c(rep(1.2,totalN))
omega2 <- omega*omega
nu <- exp(beta0+beta1*x1+beta2*x2+beta3*x3) # marginal mean
pi <- exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)/(1+exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)) # P(y_i>0)
xi <- beta0+beta1*x1+beta2*x2+beta3*x3-log(2)-log(pi)-log(pnorm(omega*delta))-omega2/2  # SN location parameter
ypre <- rsn(totalN, xi, omega, kappa)  # ypre ~ skew normal
y <- exp(ypre)  # y ~ log skew normal

# add zero values
posval<-rbinom(totalN,1,pi)
y[posval==0]<-0
ygtz<-c(rep(1,totalN))
ygtz[y==0]<-0

plot(density(ypre), type="l", col="red")
plot(density(y), type="l", col="red")
skewness(ypre)

#skewtheor <- (4-pi)/2*(delta*sqrt(2/pi))^3/(1-2*delta^2/pi)^(3/2)
#skewtheor[1]

# write to file to import into SAS
mydatak2 <- data.frame(simnum, x1, x2, x3,beta0,beta1,beta2, beta3,alpha0, alpha1, alpha2, alpha3,
                       kappa,delta,omega,nu,pi,xi,y,yzero,ygtz)
write.csv(mydatak2, "C:/Users/Valerie/Documents/BIOS General/Dissertation/Paper 3/Simulation Data/zeros20highskew/LSN high skew samp size 1000 20p zeros type 1 error 020915.csv")


# Sample size 10000
set.seed(122615)

# LSN data with kappa=0.5, sample size=200
# This has pi=0.296229, 20.4% zeros
# mean of y is 7752.807
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
#beta0 <- c(rep(7,totalN))
#beta1 <- c(rep(0.2,totalN))
#beta2 <- c(rep(-0.05,totalN))
#beta3 <- c(rep(0.05,totalN))
alpha0 <- c(rep(3,totalN))
alpha1 <- c(rep(-4,totalN))
alpha2 <- c(rep(3.5,totalN))
alpha3 <- c(rep(2.5,totalN))
kappa <- c(rep(5,totalN))  # using kappa=alpha in standard SN parameterization to avoid confusion with parameters we're using
delta <- kappa/sqrt(1+kappa^2)
omega <- c(rep(1.2,totalN))
omega2 <- omega*omega
nu <- exp(beta0+beta1*x1+beta2*x2+beta3*x3) # marginal mean
pi <- exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)/(1+exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)) # P(y_i>0)
xi <- beta0+beta1*x1+beta2*x2+beta3*x3-log(2)-log(pi)-log(pnorm(omega*delta))-omega2/2  # SN location parameter
ypre <- rsn(totalN, xi, omega, kappa)  # ypre ~ skew normal
y <- exp(ypre)  # y ~ log skew normal

# add zero values
posval<-rbinom(totalN,1,pi)
y[posval==0]<-0
ygtz<-c(rep(1,totalN))
ygtz[y==0]<-0

plot(density(ypre), type="l", col="red")
plot(density(y), type="l", col="red")
skewness(ypre)

#skewtheor <- (4-pi)/2*(delta*sqrt(2/pi))^3/(1-2*delta^2/pi)^(3/2)
#skewtheor[1]

# write to file to import into SAS
mydatak2 <- data.frame(simnum, x1, x2, x3,y,yzero,ygtz)
write.csv(mydatak2, "C:/Users/Valerie/Documents/BIOS General/Dissertation/Paper 3/Simulation Data/zeros20highskew/LSN high skew samp size 10000 20p zeros type 1 error 020915.csv")

# Sample size 50000

# first 500 
rm(list=ls())
set.seed(122615)

# LSN data with kappa=0.5, sample size=50000

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
#beta0 <- c(rep(7,totalN))
#beta1 <- c(rep(0.2,totalN))
#beta2 <- c(rep(-0.05,totalN))
#beta3 <- c(rep(0.05,totalN))
alpha0 <- c(rep(3,totalN))
alpha1 <- c(rep(-4,totalN))
alpha2 <- c(rep(3.5,totalN))
alpha3 <- c(rep(2.5,totalN))
kappa <- c(rep(5,totalN))  # using kappa=alpha in standard SN parameterization to avoid confusion with parameters we're using
delta <- kappa/sqrt(1+kappa^2)
omega <- c(rep(1.2,totalN))
omega2 <- omega*omega
nu <- exp(beta0+beta1*x1+beta2*x2+beta3*x3) # marginal mean
pi <- exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)/(1+exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)) # P(y_i>0)
xi <- beta0+beta1*x1+beta2*x2+beta3*x3-log(2)-log(pi)-log(pnorm(omega*delta))-omega2/2  # SN location parameter
ypre <- rsn(totalN, xi, omega, kappa)  # ypre ~ skew normal
y <- exp(ypre)  # y ~ log skew normal

# add zero values
posval<-rbinom(totalN,1,pi)
y[posval==0]<-0
ygtz<-c(rep(1,totalN))
ygtz[y==0]<-0

plot(density(ypre), type="l", col="red")
plot(density(y), type="l", col="red")
skewness(ypre)

#skewtheor <- (4-pi)/2*(delta*sqrt(2/pi))^3/(1-2*delta^2/pi)^(3/2)
#skewtheor[1]

# write to file to import into SAS
mydatak2 <- data.frame(simnum, x1, x2, x3,y)
write.csv(mydatak2, "C:/Users/Valerie/Documents/BIOS General/Dissertation/Paper 3/Simulation Data/For Paper 3/zeros20highskew/LSN high skew samp size 50000 20p zeros type 1 error part 1 062315.csv")

# Sample size 50000

# second 500 
rm(list=ls())
set.seed(062315)

# LSN data with kappa=0.5, sample size=50000

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
#beta0 <- c(rep(7,totalN))
#beta1 <- c(rep(0.2,totalN))
#beta2 <- c(rep(-0.05,totalN))
#beta3 <- c(rep(0.05,totalN))
alpha0 <- c(rep(3,totalN))
alpha1 <- c(rep(-4,totalN))
alpha2 <- c(rep(3.5,totalN))
alpha3 <- c(rep(2.5,totalN))
kappa <- c(rep(5,totalN))  # using kappa=alpha in standard SN parameterization to avoid confusion with parameters we're using
delta <- kappa/sqrt(1+kappa^2)
omega <- c(rep(1.2,totalN))
omega2 <- omega*omega
nu <- exp(beta0+beta1*x1+beta2*x2+beta3*x3) # marginal mean
pi <- exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)/(1+exp(alpha0+alpha1*x1+alpha2*x2+alpha3*x3)) # P(y_i>0)
xi <- beta0+beta1*x1+beta2*x2+beta3*x3-log(2)-log(pi)-log(pnorm(omega*delta))-omega2/2  # SN location parameter
ypre <- rsn(totalN, xi, omega, kappa)  # ypre ~ skew normal
y <- exp(ypre)  # y ~ log skew normal

# add zero values
posval<-rbinom(totalN,1,pi)
y[posval==0]<-0
ygtz<-c(rep(1,totalN))
ygtz[y==0]<-0

plot(density(ypre), type="l", col="red")
plot(density(y), type="l", col="red")
skewness(ypre)

#skewtheor <- (4-pi)/2*(delta*sqrt(2/pi))^3/(1-2*delta^2/pi)^(3/2)
#skewtheor[1]
memory.limit(16000)
# write to file to import into SAS
mydatak2 <- data.frame(simnum, x1, x2, x3,y)
write.csv(mydatak2, "C:/Users/Valerie/Documents/BIOS General/Dissertation/Paper 3/Simulation Data/For Paper 3/zeros20highskew/LSN high skew samp size 50000 20p zeros type 1 error part 2 062315.csv")
