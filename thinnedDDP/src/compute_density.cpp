#include "common_functions.h"

// [[Rcpp::export]]
arma::mat compute_density_1iter(arma::vec grid,
                                arma::mat weights,
                                arma::vec means,
                                arma::vec variances
                                )
{
  // weights is a matrix (trunc x G)
  if(weights.n_rows != means.n_elem) {Rcpp::Rcout << "number of rows weights != means" ;}
  if(weights.n_rows != variances.n_elem) {Rcpp::Rcout << "number of rows weights != variances" ;}
  
  int n_points = grid.n_elem ;
  int trunc = means.n_elem ;
  int G = weights.n_cols ;
  arma::mat out_density( n_points, G, arma::fill::zeros ) ;
  arma::vec tmp_dens(trunc) ;
  
  for(int g = 0; g < G; g++){
    for(int i = 0; i < n_points; i++) {
      for(int j = 0; j < trunc; j++) {
        tmp_dens(j) = weights(j,g) * R::dnorm(grid(i), means(j), std::sqrt(variances(j)), false) ;
      }
      out_density(i,g) = arma::accu(tmp_dens) ;
      tmp_dens.fill(log(0.0)) ;
    }
  }

  return(out_density) ;
}

// [[Rcpp::export]]
Rcpp::List compute_density(arma::vec grid,
                          arma::cube weights,
                          arma::mat means,
                          arma::mat variances
)
{
  // weights is a cube (trunc x G x nrep)
  // means and variances are matrices (trunc x nrep)
  if(weights.n_rows != means.n_rows) {Rcpp::Rcout << "number of rows weights != means" ;}
  if(weights.n_rows != variances.n_rows) {Rcpp::Rcout << "number of rows weights != variances" ;}
  if(variances.n_rows != means.n_rows) {Rcpp::Rcout << "number of rows variances != means" ;}
  if(weights.n_slices != means.n_cols) {Rcpp::Rcout << "number of cols weights != means" ;}
  if(weights.n_slices != variances.n_cols) {Rcpp::Rcout << "number of cols weights != variances" ;}
  if(variances.n_cols != means.n_cols) {Rcpp::Rcout << "number of cols variances != means" ;}
  
  int G = weights.n_cols ;
  int nrep = weights.n_slices ;
  int n_points = grid.n_elem ;
  
  arma::cube tmp_density_nrep(n_points, G, nrep, arma::fill::zeros) ;
  arma::mat out_density(n_points, G, arma::fill::zeros) ;
  
  for(int iter = 0; iter < nrep; iter++) {
    tmp_density_nrep.slice(iter) = compute_density_1iter(grid, weights.slice(iter), 
                           means.col(iter), variances.col(iter) ) ;
  }
  out_density = arma::mean(tmp_density_nrep, 2) ;
  
  return Rcpp::List::create(Rcpp::Named("density_mcmc") = tmp_density_nrep,
                            Rcpp::Named("mean") = out_density,
                            Rcpp::Named("seq") = grid ) ;
}





// [[Rcpp::export]]
arma::mat compute_density_1iter1DP(arma::vec grid,
                                arma::vec weights,
                                arma::vec means,
                                arma::vec variances
)
{
  // weights is a vec (trunc)
  if(weights.n_elem != means.n_elem) {Rcpp::Rcout << "number of rows weights != means" ;}
  if(weights.n_elem != variances.n_elem) {Rcpp::Rcout << "number of rows weights != variances" ;}
  
  int n_points = grid.n_elem ;
  int trunc = means.n_elem ;
  arma::vec out_density( n_points, arma::fill::zeros ) ;
  arma::vec tmp_dens(trunc) ;
  
  for(int i = 0; i < n_points; i++) {
    for(int j = 0; j < trunc; j++) {
      tmp_dens(j) = weights(j) * R::dnorm(grid(i), means(j), std::sqrt(variances(j)), false) ;
    }
    out_density(i) = arma::accu(tmp_dens) ;
    tmp_dens.fill(log(0.0)) ;
  }
  
  return(out_density) ;
}

// [[Rcpp::export]]
Rcpp::List compute_density_1DP(arma::vec grid,
                           arma::mat weights,
                           arma::mat means,
                           arma::mat variances
)
{
  // weights is a matrix (trunc x nrep)
  // means and variances are matrices (trunc x nrep)
  if(weights.n_rows != means.n_rows) {Rcpp::Rcout << "number of rows weights != means" ;}
  if(weights.n_rows != variances.n_rows) {Rcpp::Rcout << "number of rows weights != variances" ;}
  if(variances.n_rows != means.n_rows) {Rcpp::Rcout << "number of rows variances != means" ;}
  if(weights.n_cols != means.n_cols) {Rcpp::Rcout << "number of cols weights != means" ;}
  if(weights.n_cols != variances.n_cols) {Rcpp::Rcout << "number of cols weights != variances" ;}
  if(variances.n_cols != means.n_cols) {Rcpp::Rcout << "number of cols variances != means" ;}
  
  int nrep = weights.n_cols ;
  int n_points = grid.n_elem ;
  
  arma::mat tmp_density_nrep(n_points, nrep, arma::fill::zeros) ;
  arma::vec out_density(n_points, arma::fill::zeros) ;
  
  for(int iter = 0; iter < nrep; iter++) {
    tmp_density_nrep.col(iter) = compute_density_1iter1DP(grid, weights.col(iter), 
                           means.col(iter), variances.col(iter) ) ;
  }
  out_density = arma::mean(tmp_density_nrep, 1) ;
  
  return Rcpp::List::create(Rcpp::Named("density_mcmc") = tmp_density_nrep,
                            Rcpp::Named("mean") = out_density,
                            Rcpp::Named("seq") = grid );
}










