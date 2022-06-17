#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

// [[Rcpp::export]]
vec expit(vec xb){
  return (1/(1+exp(-xb)));
}

// [[Rcpp::export]]
vec logit(vec xb){
  return log(xb/(1-xb));
}

/*** R
log_lik <- function(theta, y, x){
  p <- expit(x%*%theta) # probability of success
  sum(dbinom(y, size=1, prob=p, log=T))
}
*/

// [[Rcpp::export]]
mat MHu(mat y, mat x, vec beta_init, vec jump_sigma, int n_iter, int burn_in){
  int total_iter = n_iter + burn_in + 1;
  
  //Theta storage
  mat theta_hat = zeros(total_iter, beta_init.n_elem);
  theta_hat.row(0) = beta_init.t();
  
  // Probability vector for random scan
  vec prob = {0.4,0.1,0.1,0.4};
  Rcpp::Function sample = Rcpp::Environment("package:base")["sample"];
  
  vec j = sample(prob.n_elem,1,false,prob);
  
  // Run algorithm
  for(int i=1; i<=total_iter; i++){
  
  }
  
  return theta_hat;
}