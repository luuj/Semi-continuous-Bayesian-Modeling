library(sn) # Skew normal library
library(flexsurv) #Generalized gamma library
set.seed(123)

## Define math functions
expit <- function(x){plogis(x)}

## Generates semi-continuous data using the log-skew-normal distribution
# zero: Determines the percentage of zeroes introduced in the outcome. 
#       Can take values {10, 20, 40}
# kappa: Determines the shape of the log-skew-normal distribution
# sigma: Determines the scale of the log-skew-normal distribution
genSCD_LSN <- function(n=200, zero=20, kappa=5, sigma=1.2){
   # Generate x1 ~ Bernoulli(0.5)
   #          x2 ~ Normal(0,1)
   #          x3 ~ Poisson(1)
   x0 <- rep(1, n)
   x1 <- rbinom(n, 1, 0.5)
   x2 <- rnorm(n)
   x3 <- rpois(n, 1)
   x <- matrix(c(x0,x1,x2,x3),
               nrow=n, ncol=4)
   
   # Marginal mean structure
   # exp(6 + 0.2x1 - 0.01x2 + 0.05x3)
   beta <- c(6, 0.2, -0.01, 0.05)
   nu <- exp(x%*%beta)
   
   # Zero structure
   # 10% zeros: logit(pi) = 3 - 2.4x1 + 1.5x2 + 2.0x3
   # 20% zeros: logit(pi) = 3 - 4.0x1 + 3.5x2 + 2.5x3
   # 40% zeros: logit(pi) = 3 - 7.0x1 + 5.0x2 + 2.0x3
   if (zero==10)
      alpha <- c(3, -2.4, 1.5, 2)
   else if (zero==20)
      alpha <- c(3, -4, 3.5, 2.5)
   else if (zero==40)
      alpha <- c(3, -7, 5, 2)
   else{
      print("Invalid zero value. Can only take values {10,20,40}")
      return(-1)
   }
   
   # Re-express the LSN likelihood as a function of beta
   delta <- kappa/sqrt(1+kappa^2)
   pi <- expit(x%*%alpha)
   xi <- x%*%beta - log(2) - log(pi) - log(pnorm(sigma*delta)) - (sigma^2)/2 
   
   # Generate the outcome
   ypre <- rsn(n, xi, sigma, kappa) # Skew normal
   attributes(ypre) <- NULL # Remove unnecessary output
   y <- exp(ypre) # Log skew normal
   
   # Introduce zeros into the outcome
   posval <- rbinom(n, 1, pi)
   y[posval==0] <- 0

   return(y)
}

genSCD_LSN(n=100)








## Generates semi-continuous data using the generalized gamma distribution
# zero: Determines the percentage of zeroes introduced in the outcome. 
#       Can take values {10, 20, 40}
# kappa: Determines the shape of the generalized gamma distribution
# sigma: Determines the scale of the generalized gamma distribution
genSCD_GG <- function(n=200, zero=20, kappa=0.63, sigma=1.2){
   # Generate x1 ~ Bernoulli(0.5)
   #          x2 ~ Normal(0,1)
   #          x3 ~ Poisson(1)
   x0 <- rep(1, n)
   x1 <- rbinom(n, 1, 0.5)
   x2 <- rnorm(n)
   x3 <- rpois(n, 1)
   x <- matrix(c(x0,x1,x2,x3),
               nrow=n, ncol=4)
   
   # Marginal mean structure
   # exp(6 + 0.2x1 - 0.01x2 + 0.05x3)
   beta <- c(6, 0.2, -0.01, 0.05)
   nu <- exp(x%*%beta)
   
   # Zero structure
   # 10% zeros: logit(pi) = 3 - 2.4x1 + 1.5x2 + 2.0x3
   # 20% zeros: logit(pi) = 3 - 4.0x1 + 3.5x2 + 2.5x3
   # 40% zeros: logit(pi) = 3 - 7.0x1 + 5.0x2 + 2.0x3
   if (zero==10)
      alpha <- c(3, -2.4, 1.5, 2)
   else if (zero==20)
      alpha <- c(3, -4, 3.5, 2.5)
   else if (zero==40)
      alpha <- c(3, -7, 5, 2)
   else{
      print("Invalid zero value. Can only take values {10,20,40}")
      return(-1)
   }
   
   # Re-express the GG likelihood as a function of beta
   pi <- expit(x%*%alpha)
   mu <- x%*%beta - log(pi) - (sigma*log(kappa^2))/kappa - 
      log(gamma(1/kappa^2 + sigma/kappa)) +
      log(gamma(1/kappa^2))
   
   # Generate the outcome
   y <- rgengamma(n, mu=mu, sigma=sigma, Q=kappa) 

   # Introduce zeros into the outcome
   posval <- rbinom(n, 1, pi)
   y[posval==0] <- 0
   
   return(y)
}

genSCD_GG()







# Experimental LSN
lsn.data <- genSCD_LSN(n=10000, zero=20, kappa=5, sigma=1.2)
mean(lsn.data)
summary(lsn.data)
sort(lsn.data, decreasing = TRUE)[1:50]
plot <- lsn.data[lsn.data<5000 & lsn.data >0]
hist(plot, xlab = "Cost", main="Cost histogram (values > 0)")

lsn.data.low <- genSCD_LSN(n=10000, zero=20, kappa=0.5, sigma=1.2)
mean(lsn.data.low)
summary(lsn.data.low)
sort(lsn.data.low, decreasing = TRUE)[1:50]

lsn.data.high <- genSCD_LSN(n=10000, zero=20, kappa=10, sigma=1.2)
mean(lsn.data.high)
summary(lsn.data.high)
sort(lsn.data.high, decreasing = TRUE)[1:50]

# Experimental GG
gg.data <- genSCD_GG(n=10000, zero=20, kappa = 0.63, sigma=1.2)
mean(gg.data)
summary(gg.data)
sort(gg.data, decreasing = TRUE)[1:50]

gg.data.high <- genSCD_GG(n=10000, zero=20, kappa = 5, sigma=1.2)
mean(gg.data.high)
summary(gg.data.high)
sort(gg.data.high, decreasing = TRUE)[1:50]
