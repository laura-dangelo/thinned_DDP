#ifndef COMMONFUNCTIONS_H
#define COMMONFUNCTIONS_H

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>


int sample_i(arma::vec ids, arma::vec prob) ;

arma::vec stick_breaking(arma::vec beta_var, bool logarithm) ;

arma::vec rmvnorm(arma::vec mu, arma::mat sigma) ;


#endif