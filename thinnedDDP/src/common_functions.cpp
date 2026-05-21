// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
  
#include "common_functions.h"
#include <RcppArmadilloExtensions/sample.h>
  
// sample fun
int sample_i(arma::vec ids, arma::vec prob) {
  int out ;
  out = Rcpp::RcppArmadillo::sample(ids, 1, false, prob)[0] ;
  return(out) ;
}


// compute stick-breaking weights starting from the vector of beta r.v. Beta(a_k,b_k)
arma::vec stick_breaking(arma::vec beta_var, bool logarithm) 
{
  int len = beta_var.n_elem ;
  arma::vec logm1_beta_var = log(1-beta_var) ;
  arma::vec out(len) ;
  
  out(0) = log(beta_var(0)) ;
  for(int k = 1; k < len; k++)
  {
    out(k) = log(beta_var(k)) + arma::accu( logm1_beta_var(arma::span(0, k-1) ) ) ;
  }
  if(!logarithm) { return exp(out) ; } else { return(out) ; }
}


/*
 * sample multivariate normal distribution
 */
arma::vec rmvnorm(arma::vec mu, arma::mat sigma) {
  int d = mu.n_elem ;
  arma::vec Y = arma::randn(d) ;
  return mu + arma::chol(sigma).t() * Y ;
}



