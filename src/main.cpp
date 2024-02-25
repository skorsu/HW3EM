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
arma::mat iX(arma::vec y, double p, double lb, double mu){
  
  arma::mat result(3, 3, arma::fill::zeros);
  
  arma::vec e1 = mu * y;
  arma::vec e2 = lb * y;
  arma::vec e5 = arma::exp(-e2);
  arma::vec e6 = arma::exp(-e1);
  double e7 = 1 - p;
  arma::vec e10 = (lb * p * e5) + (mu * e7 * e6);
  arma::vec e11 = 1 - e2;
  arma::vec e12 = 1 - e1;
  arma::vec e14 = (lb * e5) - (mu * e6);
  arma::vec e15 = arma::pow(e10, 2);
  arma::vec e16 = p * e11;
  arma::vec e21 = e12 * e7;
  double e22 = arma::accu(-(e7 * (e16 % e12 % e5 % e6)/e15));
  double e23 = arma::accu(e11 % (1 - p * e14/e10) % e5/e10);
  
  result.row(0).col(0).fill(arma::accu(-(arma::pow(e14, 2)/e15)));
  result.row(1).col(1).fill(arma::accu(-(p * (e11 % (e16 % e5/e10 + y) + y) % e5/e10)));
  result.row(2).col(2).fill(arma::accu(-(((e21 % e6/e10 + y) % e12 + y) * e7 % e6/e10)));
  
  result.row(0).col(1).fill(e23);
  result.row(1).col(0).fill(e23);
  
  result.row(0).col(2).fill(arma::accu(-((e7 * e14/e10 + 1) % e12 % e6/e10)));
  result.row(2).col(0).fill(arma::accu(-((e7 * e14/e10 + 1) % e12 % e6/e10)));
  
  result.row(1).col(2).fill(e22);
  result.row(2).col(1).fill(e22);
  
  return -result;
  
}

// [[Rcpp::export]]
arma::mat iY(arma::vec y, double p, double lb, double mu){
  
  arma::mat result(3, 3, arma::fill::zeros);
  arma::vec db = dhat_calc(y, p, lb, mu);
  
  arma::vec e1 = 1 - db;
  result.row(0).col(0).fill(arma::accu(-(e1/std::pow(1 - p, 2) + db/std::pow(p, 2))));
  result.row(1).col(1).fill(arma::accu(-(db/std::pow(lb, 2))));
  result.row(2).col(2).fill(arma::accu(-(e1/std::pow(mu, 2))));
  
  return -result;
  
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

// [[Rcpp::export]]
arma::mat BootResult(arma::vec y, double p0, double lambda0, double mu0, 
                     double eps, arma::vec EMinit, unsigned int M){
  
  arma::mat result(M, 3, arma::fill::zeros);
  result.row(0) = EMinit.t();
  arma::imat Bindex = arma::randi<arma::imat>(y.size(), M, arma::distr_param(0, y.size() - 1));
  Rcpp::List rcppEM; 
  for(int i = 1; i < M; ++i){
    arma::uvec bindex = arma::conv_to<arma::uvec>::from(Bindex.col(i));
    // Run EM
    rcppEM = EM_rcpp(y.rows(bindex), p0, lambda0, mu0, eps);
    double p = rcppEM["p"];
    double lb = rcppEM["lambda"];
    double mu = rcppEM["mu"];
    result.row(i).col(0) = p;
    result.row(i).col(1) = lb;
    result.row(i).col(2) = mu;
  }
  
  return result;
  
}
