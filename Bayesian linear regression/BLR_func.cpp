#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

/*Convert vec to NumericVector*/
template <typename T>
Rcpp::NumericVector arma2vec(const T& x) {
  return Rcpp::NumericVector(x.begin(), x.end()); 
}

/*Log likelihood function*/
// [[Rcpp::export]]
double loglik(const mat& theta, const mat& y, const mat& x){
  /*Calculate probabilities*/
  vec lp = x*theta.t();
  Rcpp::NumericVector p = arma2vec(lp);
  vec results = Rcpp::plogis(p,true,false);

  /*Calculate likelihood*/
  int n = x.n_rows;
  double ll = 0;
  for (int i=0; i<n; i++){
    lp = x.row(i) * theta.t();
    ll += (y(i) * lp - log(1+exp(lp))).eval()(0,0);
  }
  
  return(ll);
}

/*Metropolis Hastings algorithm*/
// [[Rcpp::export]]
void MHu(const uword& total_iter, const vec& j_all, const vec& flips, const vec& jumps, 
         mat& theta_hat, const mat& y, const mat& x){
  /*Initialize variables*/
  uword j;
  mat prop_theta;
  double proposal_theta_j, log_ratio, log_prop;
  double log_prev = 0;

  for (uword i=1; i<total_iter; i++){
    j = j_all(i-1);
    
    /*Generate proposal*/
    proposal_theta_j = theta_hat(i-1, j) + jumps(i-1);
    prop_theta = theta_hat.row(i-1);
    prop_theta(j) = proposal_theta_j;
    
    /*Calculate acceptance ratio*/
    if (i==1){
      log_prev = loglik(theta_hat.row(i-1),y,x);
    }
    log_prop = loglik(prop_theta, y, x);
    log_ratio = log_prop - log_prev;
    
    if (flips(i-1) < log_ratio){
      theta_hat.row(i) = prop_theta;
      log_prev = log_prop;
    }else{
      theta_hat.row(i) = theta_hat.row(i-1);
    }
  }
}