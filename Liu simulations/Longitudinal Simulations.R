library(dplyr)
library(ggplot2)
library(gridExtra)
expit <- function(x){plogis(x)}

## Data generating function
# n = number of individuals
# a_sigma = sd of random intercept a_i
# b_sigma = sd of random intercept b_i
# e_sigma = sd of error term e_ij
# delta = level of association between a_i and b_i
generateData <- function(n, a_sigma, b_sigma, e_sigma, delta, seed=12){
   set.seed(seed)
   
   # Binary time invariant covariate - takes value 0 or 1 with equal probability
   Z1 <- sample(0:1, n, replace=TRUE)
   
   # Time variant covariate - how many months an individual is observed
   Z2 <- sample(6:12, n, replace=TRUE)
   Z2.vals <- unlist(sapply(Z2, function(x){
      return(0:x)
   }))
   
   # Alpha and beta fixed effects coefficients
   alpha <- c(-1, 1, 0.1)
   beta <- c(-1, -0.5, 0.1)
   
   # Generate a and b - random intercepts
   a_val <- rnorm(n, sd=a_sigma)
   b_val <- rnorm(n, sd=b_sigma)
   a <- rep(a_val, Z2+1)
   b <- rep(b_val, Z2+1)
   
   # Setup data frame
   Z1.vals <- rep(Z1, Z2+1)
   id <- rep(1:n, Z2+1)
   dat <- data.frame(id=id, time=Z2.vals, gender=Z1.vals)
   
   # Generate error term
   e <- rnorm(nrow(dat), sd=e_sigma)
   
   # Generate model values
   dat$eta_ij <- alpha[1] + dat$gender*alpha[2] + dat$time*alpha[3] + a
   dat$m_ij <- beta[1] + dat$gender*beta[2] + dat$time*beta[3] + delta*a + b + e
   dat$cost <- round(exp(dat$m_ij)*100,2)
   
   # Probability of monthly cost being positive
   acceptProb <- runif(nrow(dat)) > expit(dat$eta_ij)

   # Set value to 0 depending on eta_ij
   dat$cost[!acceptProb] <- 0
   
   # Incorporate death into model
   gamma <- 1
   delta2 <- 0.5
   delta3 <- 0.7
   shape <- 1.5
   scale <- 0.05
   
   U <- runif(n)
   deathTime <- (-log(U) / (scale*exp(Z1*gamma + delta2*a_val + delta3*b_val)))^(1/shape)
   deathTime <- ceiling(deathTime)

   # Calculate cumulative sum
   dat <- data.frame(dat %>% group_by(id) %>% mutate(csum=cumsum(cost)))
   
   return(dat)
}

# Original parameters
dat1 <- generateData(n=100, a_sigma=1, b_sigma=1, e_sigma=1.5, delta=0.5)

# Change a_sigma
dat2 <- generateData(n=100, a_sigma=0.1, b_sigma=1, e_sigma=1.5, delta=0.5)
dat3 <- generateData(n=100, a_sigma=0.5, b_sigma=1, e_sigma=1.5, delta=0.5)
dat4 <- generateData(n=100, a_sigma=2, b_sigma=1, e_sigma=1.5, delta=0.5)

# Change b_sigma
dat5 <- generateData(n=100, a_sigma=1, b_sigma=0.1, e_sigma=1.5, delta=0.5)
dat6 <- generateData(n=100, a_sigma=1, b_sigma=0.5, e_sigma=1.5, delta=0.5)
dat7 <- generateData(n=100, a_sigma=1, b_sigma=2, e_sigma=1.5, delta=0.5)

# Change delta
dat8 <- generateData(n=100, a_sigma=1, b_sigma=1, e_sigma=1.5, delta=2)
dat9 <- generateData(n=100, a_sigma=2, b_sigma=2, e_sigma=1.5, delta=1)
dat10 <- generateData(n=100, a_sigma=2, b_sigma=2, e_sigma=1.5, delta=0.1)

# Change a_sigma and b_sigma
dat11 <- generateData(n=100, a_sigma=0.1, b_sigma=0.1, e_sigma=1.5, delta=0.5)
dat12 <- generateData(n=100, a_sigma=0.1, b_sigma=2, e_sigma=1.5, delta=0.5)
dat13 <- generateData(n=100, a_sigma=2, b_sigma=0.1, e_sigma=1.5, delta=0.5)
dat14 <- generateData(n=100, a_sigma=2, b_sigma=2, e_sigma=1.5, delta=0.5)



# Plot some trajectories
y_ax <- 7000
p <- ggplot(data=dat1, aes(x=time, y=csum, group=id))
p1 <- p + geom_line()  + ylim(c(0,y_ax)) + ggtitle("original")

# Increasing a_sigma
p <- ggplot(data=dat2, aes(x=time, y=csum, group=id))
p2 <- p + geom_line()  + ylim(c(0,y_ax)) + ggtitle("a_sig=0.1")

p <- ggplot(data=dat3, aes(x=time, y=csum, group=id))
p3 <- p + geom_line()  + ylim(c(0,y_ax)) + ggtitle("a_sig=0.5")

p <- ggplot(data=dat4, aes(x=time, y=csum, group=id))
p4 <- p + geom_line()  + ylim(c(0,y_ax)) + ggtitle("a_sig=2")

grid.arrange(p1, p2, p3, p4)



# Changing b_sigma
p <- ggplot(data=dat5, aes(x=time, y=csum, group=id))
p5 <- p + geom_line() + ylim(c(0,7500)) + ggtitle("b_sig=0.1")

p <- ggplot(data=dat6, aes(x=time, y=csum, group=id))
p6 <- p + geom_line() + ylim(c(0,7500)) + ggtitle("b_sig=0.5")

p <- ggplot(data=dat7, aes(x=time, y=csum, group=id))
p7 <- p + geom_line() + ylim(c(0,7500)) + ggtitle("b_sig=2")

grid.arrange(p1, p5, p6, p7)



# Changing delta
p <- ggplot(data=dat8, aes(x=time, y=csum, group=id))
p8 <- p + geom_line() + ylim(c(0,7500)) + ggtitle("delta=0.1")

p <- ggplot(data=dat9, aes(x=time, y=csum, group=id))
p9 <- p + geom_line() + ylim(c(0,7500)) + ggtitle("delta=1")

p <- ggplot(data=dat10, aes(x=time, y=csum, group=id))
p10 <- p + geom_line() + ylim(c(0,7500)) + ggtitle("delta=2")

grid.arrange(p1, p8, p9, p10)



# Changing a_sigma and b_sigma
p <- ggplot(data=dat11, aes(x=time, y=csum, group=id))
p11 <- p + geom_line() + ylim(c(0,7500)) + ggtitle("a_sig=0.1, b_sig=0.1")

p <- ggplot(data=dat12, aes(x=time, y=csum, group=id))
p12 <- p + geom_line() + ylim(c(0,7500)) + ggtitle("a_sig=0.1, b_sig=2")

p <- ggplot(data=dat13, aes(x=time, y=csum, group=id))
p13 <- p + geom_line() + ylim(c(0,7500)) + ggtitle("a_sig=2, b_sig=0.1")

p <- ggplot(data=dat14, aes(x=time, y=csum, group=id))
p14 <- p + geom_line() + ylim(c(0,7500)) + ggtitle("a_sig=2, b_sig=2")

grid.arrange(p1, p11, p12, p13, p14)
