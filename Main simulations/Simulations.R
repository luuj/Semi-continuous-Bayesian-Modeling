source("jags_re.R")
library("dplyr")
desktop <- "/Users/jonathanluu/Desktop/Simulations/"

# Load in all the data files
jagFileCount <- 150
loadRData <- function(fileName,num){
  load(fileName, .GlobalEnv)
  assign(paste0("jags",i), results, .GlobalEnv)
  rm(results, envir=.GlobalEnv)
}

for(i in 1:jagFileCount){
  tryCatch({loadRData(paste0(desktop,"Scenario 4 - Death/jags",i,".Rdata"), i)},error=function(cond){})
}

# Results
rhats <- NULL
means <- NULL
c025 <- NULL
c975 <- NULL

# datname=jags file name, storage=rhats or mean, type=rhats or mean
makeSummary <- function(datname, storage, type){  
  dat <- get(datname)
  if (is.null(storage)){
    return(data.frame(jags1=dat$BUGSoutput$summary[,type]))
  }
  
  out <- cbind(storage, dat$BUGSoutput$summary[,type])
  colnames(out) <- c(names(storage), datname)
  return(out)
}

for (i in 1:jagFileCount){
  filename <- paste0("jags",i)
  rhats <- makeSummary(filename, rhats, "Rhat")
  means <- makeSummary(filename, means, "mean")
  c025 <- makeSummary(filename, c025, "2.5%")
  c975 <- makeSummary(filename, c975, "97.5%")
}

# Traceplots
traceplot(jags1, mfrow=c(3,3), ask=F)

# Average over 150 posterior means
apply(rhats,1,mean)
apply(means,1,mean)

# True values of coefficients
truth <- c(alphay=-1.5, betay=1.3, betay2=0.02, betay3=-0.5, gammay=0.8,
           alpham=1, betam=1.5, betam2=0.02, betam3=-0.5, shape=2,
           alphad=-2.5, betad=-0.5, betad2=0.02, betad3=-0.5, gammad=-1.3,
           V1_sig=.9, V2_sig=0.9, Vmu_sig=0.1)
parnames <- c("alphay","betay","betay2","betay3","gammay",
              "alpham","betam","betam2","betam3","shape",
              "alphad","betad","betad2","betad3","gammad",
              "V1_sig","V2_sig","Vmu_sig")

# Coverage
results <- data.frame(matrix(nrow=jagFileCount, ncol=18))
colnames(results) <- parnames
for (i in 1:length(parnames)){
  for(j in 1:jagFileCount){
    results[j,i] <- between(truth[i], c025[parnames[i],j], c975[parnames[i],j])
  }
}

colMeans(results)

# Maximum PSR over all parameters
sort((apply(rhats,2,max)))
View(rhats)



# Ideas for exploratory analysis
# Summary table - How many nursing homes in dataset, distribution of them by state, average number of individuals in each home
# Cost distributions of a couple nursing homes in different areas (urban vs rural)
# Average length of stay at a nursing home before death
# Cost accrual in last year of life

# Method that ignores mortality - conditional on probability of dying
# Scenario 4 - method ignores mortality
# Scenario 5 - normal method


