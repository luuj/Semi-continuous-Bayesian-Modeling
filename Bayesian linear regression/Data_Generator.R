# Covariate matrix and number of observations
n_j <- 200
J <- 5
n <- n_j*J
design_matrix <- data.frame(1,trt=rbinom(n,1,0.5),
                age=sample(20:70,n,replace=TRUE),
                bmi=rnorm(n,28,4))

## Generate simple data set
# param = beta coefficients
# x = covariate matrix
# n_j is the number of observations per cluster
# J = number of clusters
# sigma_v = variance of random effect distribution
gen_data <- function(param, x, n_j, J, sigma_v){
  # Generate covariates
  cluster <- rep(1:J, each=n_j)
  V_j <- 0
  
  # Generate outcome
  if (J > 1){
    V_j <- rep(rnorm(J,0,sigma_v), each=n_j)
    pi <- expit(as.matrix(x)%*%as.matrix(param) + V_j) # Random intercept
  }else{
    pi <- expit(as.matrix(x)%*%as.matrix(param))
  }
  
  y <- rbinom(n_j*J,1,pi)
  
  data.frame(y,x[,-1],cluster,V_j)
}

