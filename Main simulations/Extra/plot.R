library(ggplot2)
library(ggpubr)
source("../gendata.R")

clusterSummary <- function(covar, fn, dat, includeTrt = T){
  pars <- as.list(match.call()[-1]) # Convert covar into variable
  out <- aggregate(eval(pars$covar) ~ id+clst+trt, dat, max)
  colnames(out) <- c("id", "clst", "trt", "max")
  if (includeTrt)
    return(aggregate(max~clst+trt, out, fn))
  else
    return(aggregate(max~clst, out, fn))
}

set.seed(10)
dat <- gen_data(n_k=70, n_t=12, k=250)

dthPlot <- clusterSummary(t, mean, dat)
dthPlot$trt<- factor(dthPlot$trt, labels=c("Group 0", "Group 1"))
d <- dthPlot %>% ggplot(aes(max, color=trt)) + geom_histogram(fill="white", bins = 11) +
  ylab("Number of nursing homes") + ggtitle("Average time to death distribution by cluster") + 
  facet_wrap(~trt) + xlab("Time to death (months)") + theme(legend.position = "none")

costPlot <- clusterSummary(cumcost, mean, dat)
costPlot$trt<- factor(costPlot$trt, labels=c("Group 0", "Group 1"))
c <- costPlot %>% ggplot(aes(max, color=trt)) + geom_histogram(fill="white", bins = 100) +
  ylab("Number of nursing homes") + ggtitle("Average total cost accrued distribution by cluster") +
  facet_wrap(~trt) + xlab("Total cost accrued (dollars)") + theme(legend.position = "none")

dthprop <- clusterSummary(y2, mean, dat)
costprop <- clusterSummary(y1, mean, dat)
scatplot <- left_join(dthprop, costprop, by=c("clst","trt"))
colnames(scatplot) <- c("cluster", "trt", "death", "cost")
scatplot$trt <- factor(scatplot$trt, labels=c("Group 0", "Group 1"))
# Scatter plot for proportion of at least one cost and death for each cluster
s <- scatplot %>% ggplot(aes(death, cost, color=trt)) + geom_point() + ylab("Proportion with any cost accrued") +
  xlab("Proportion of deaths") + ggtitle("Cost vs. death proportions by cluster") + theme(legend.position = "none")

e <- ggarrange(d, c, s,
          labels = c("A", "B", "C"),
          ncol = 1, nrow = 3)

