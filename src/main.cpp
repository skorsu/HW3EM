#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec dexp_arma(arma::vec y, double rate){
  arma::vec result(y.size(), arma::fill::value(rate));
  result %= arma::exp(-result % y);
  return result;
}

// [[Rcpp::export]]
arma::vec dhat_calc(arma::vec y, double p, double lambda, double mu){
  arma::vec result = p * dexp_arma(y, lambda);
  result %= arma::pow(result + ((1 - p) * dexp_arma(y, mu)), -1);
  return result;
}

// [[Rcpp::export]]
Rcpp::List EM_rcpp(arma::vec y, double p0, double lambda0, double mu0, 
                   double eps){
  
  unsigned int iter = 1;
  double p = p0;
  double lambda = lambda0;
  double mu = mu0;
  arma::vec theta = {p, lambda, mu}; 
  
  arma::vec dhat = dhat_calc(y, p, lambda, mu);
  
  double p_new = arma::mean(dhat);
  double lambda_new = arma::accu(dhat)/arma::accu(dhat % y);
  double mu_new = arma::accu(1 - dhat)/arma::accu((1 - dhat) % y);
  
  arma::vec theta_new = {p_new, lambda_new, mu_new};
  
  while(arma::norm(theta - theta_new, 2.0) >= eps){
    iter += 1;
    p = p_new;
    lambda = lambda_new;
    mu = mu_new;
    theta = {p, lambda, mu}; 
    
    p_new = arma::mean(dhat);
    lambda_new = arma::accu(dhat)/arma::accu(dhat % y);
    mu_new = arma::accu(1 - dhat)/arma::accu((1 - dhat) % y);
    
    theta_new = {p_new, lambda_new, mu_new};
    
  }
  
  Rcpp::List result;
  result["iter"] = iter;
  result["p"] = p;
  result["lambda"] = lambda;
  result["mu"] = mu;
  return result;
  
}