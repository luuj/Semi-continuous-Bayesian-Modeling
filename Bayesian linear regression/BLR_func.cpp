/*C++ Metropolis Hastings Code*/
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

/*Convert vec to NumericVector*/
template <typename T>
Rcpp::NumericVector arma2vec(const T& x) {
  return Rcpp::NumericVector(x.begin(), x.end()); 
}





/***** Unclustered MH *****/
/*Log likelihood function*/
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

/* Unclustered MH algorithm
 * y - outcome values
 * x - covariates
 * beta_init - starting beta values for MH algorithm
 * jump_sigma - variance for proposal distribution
 * n_iter - number of MH steps
 * burn_in - number of iterations to remove for burn-in*/
// [[Rcpp::export]]
mat MHu(const mat& y, const mat& x, const vec& beta_init, 
        const vec& jump_sigma, 
        const int& n_iter, const int& burn_in){
  /*MH storage variables*/
  int total_iter = n_iter + burn_in + 1;
  double proposal_theta_j, log_ratio, log_prop, log_prev = 0;
  mat prop_theta, theta_hat = zeros(total_iter,beta_init.n_elem); 
  theta_hat.row(0) = beta_init.t();

  /*Random scan variables*/
  int j;
  Rcpp::NumericVector prob = {0.4,0.1,0.1,0.4};
  Rcpp::NumericVector pool = {0,1,2,3};
  Rcpp::NumericVector jall = Rcpp::sample(pool,total_iter,true,prob);
  
  /*Run algorithm*/
  for (int i=1; i<total_iter; i++){
    j = jall(i-1);
    
    /*Generate proposal*/
    proposal_theta_j = theta_hat(i-1, j) + R::rnorm(0,jump_sigma(j));
    prop_theta = theta_hat.row(i-1);
    prop_theta(j) = proposal_theta_j;
    
    /*Calculate acceptance ratio*/
    if (i==1){
      log_prev = loglik(theta_hat.row(i-1),y,x);
    }
    log_prop = loglik(prop_theta, y, x);
    log_ratio = log_prop - log_prev;
    
    if(log(R::runif(0,1)) < log_ratio){
      theta_hat.row(i) = prop_theta;
      log_prev = log_prop;
    }else{
      theta_hat.row(i) = theta_hat.row(i-1);
    }
  }
  
  return(theta_hat);
}




/***** Clustered MH *****/
/*Log likelihood function for beta*/
double loglik_beta(const mat& theta, const mat& y, const mat& x, const mat& vk){
  /*Calculate probabilities*/
  vec lp = x*theta.t() + vk;
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

/*Log likelihood function for Vk*/
/*Note: Vk is a single number. y and x are reduced to match the corresponding Vk*/
double loglik_vk(const mat& theta, const mat& y, const mat& x, const vec& vk, const vec& sigma_vk){
  double ll = loglik_beta(theta,y,x,vk);
  ll += R::dnorm(vk(0),0,sigma_vk(0),true);
  
  return(ll);
}

/*Log likelihood function for sigma*/
/*num_cluster: total number of clusters*/
/*a and b: initial gamma values*/
double loglik_sv(const vec& vk, const vec& sigma_vk, int num_cluster, double a, double b){
  double vk_sum = accu(vk);
  double ll = R::dgamma(sigma_vk(0), a+num_cluster/2, b+vk_sum/2, true);
  return(ll);
}

/* Clustered MH algorithm
 * y - outcome values
 * x - covariate matrix
 * theta_init - starting beta, vk, and sigma_v values for MH algorithm
 * gamma_init - starting gamma prior distribution parameters
 * jump_sigma - variance for proposal distribution
 * n_iter - number of MH steps
 * burn_in - number of iterations to remove for burn-in*/
// [[Rcpp::export]]
mat MHc(const mat& y, const mat& x, const vec& beta_init, 
        const vec& jump_sigma, 
        const int& n_iter, const int& burn_in){
  /*MH storage variables*/
  int total_iter = n_iter + burn_in + 1;
  
  
  return y;
}



