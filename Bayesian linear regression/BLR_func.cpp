#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

Rcpp::Function log_lik = Rcpp::Environment::global_env()["log_lik"];

// [[Rcpp::export]]
vec expit(vec xb){
  return (1/(1+exp(-xb)));
}

// [[Rcpp::export]]
vec logit(vec xb){
  return log(xb/(1-xb));
}

// [[Rcpp::export]]
mat MHu(uword total_iter, vec j_all, vec flips, vec jumps, mat theta_hat, vec jump_sigma, mat y, mat x){
  /*Initialize variables*/
  uword j;
  double proposal_theta_j, log_ratio;
  Rcpp::NumericVector log_prop, log_prev;
  
  for (uword i=1; i<total_iter; i++){
    j = j_all(i-1);
    
    /*Generate proposal*/
    proposal_theta_j = theta_hat(i-1, j) + jumps(i-1);
    mat prev_theta = theta_hat.row(i-1);
    mat prop_theta = prev_theta;
    prop_theta(j) = proposal_theta_j;
    
    /*Calculate acceptance*/
    log_prop = log_lik(prop_theta, y, x);
    log_prev = log_lik(prev_theta, y, x);
    log_ratio = log_prop(0) - log_prev(0);
    
    if (flips(i-1) < log_ratio){
      theta_hat.row(i) = prop_theta;
    }else{
      theta_hat.row(i) = theta_hat.row(i-1);
    }
  }
  
  return theta_hat;
}