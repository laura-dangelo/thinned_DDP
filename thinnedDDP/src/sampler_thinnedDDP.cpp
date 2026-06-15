#include "common_functions.h"

// [[Rcpp::export]]
Rcpp::List sampler_thinnedDDP_arma(int nrep, // number of replications of the Gibbs sampler
                                   int burnin, // number of replications to discard as burn-in
                                   int thinning_factor, // e.g., a thinning factor of 2 saves half of the chains
                                   const arma::vec & y, // input data
                                   const arma::vec & group, // group assignment for each observation in the vector y
                                   int trunc, // maximum number of clusters (truncation)
                                   double mu0, double tau0, // hyperparameters on the N-iG prior on the mean parameter, mu|sigma2 ~ N(mu0, sigma2 / tau0)
                                   double gamma0, double lambda0, // hyperparameters on the N-iG prior on the variance parameter, 1/sigma2 ~ Gamma(gamma0, lambda0)
                                   double alpha, // DP concentration parameter
                                   double a_beta, double b_beta, // hyperparameters of the beta prior on pi_g (group-specific thinning probabilities)
                                   arma::vec mu_start, 
                                   arma::vec sigma2_start,
                                   arma::vec cl_start,
                                   bool progressbar
                                   ) 
{
  /*
  * Blocked Gibbs sampler for the thinned-DDP for grouped data
  */

  int N = y.n_elem ;
  arma::vec unique_groups = arma::unique(group) ;
  int G = unique_groups.n_elem ;
  
  // allocate output matrices
  arma::mat out_mu(trunc, (nrep-burnin)/thinning_factor, arma::fill::zeros) ; // cluster-specific means
  arma::mat out_sigma2(trunc, (nrep-burnin)/thinning_factor, arma::fill::ones) ; // cluster-specific variances
  arma::mat out_cl(N, (nrep-burnin)/thinning_factor) ; // cluster allocation of each observation
  arma::cube out_pi(trunc, G, (nrep-burnin)/thinning_factor, arma::fill::zeros) ; // group-specific cluster allocation probabilities
  arma::vec v_j = Rcpp::rbeta(trunc, 1.0, alpha) ; // beta r.v. of stick-breaking
  arma::vec thinned_vj = v_j ;
  arma::vec log_m1v_j = log(1.0 - v_j) ;
  arma::cube out_ell(trunc, G, (nrep-burnin)/thinning_factor, arma::fill::ones) ; // thinning variables ell_j
  arma::mat out_prob_thinning(G, (nrep-burnin)/thinning_factor, arma::fill::zeros) ; // group-specific thinning probabilities
  
  double p1 ; 
  double log_p1 ; double log_p0 ;
  
  // initialization
  arma::vec tmp_mu = mu_start ;
  arma::vec tmp_sigma2 = sigma2_start ;
  arma::vec tmp_cl = cl_start ;
  arma::vec tmp_prob_thinning(G) ;
  arma::mat tmp_ell(trunc, G, arma::fill::ones) ;
  
  arma::mat tmp_pi(trunc, G, arma::fill::zeros) ;
  for(int g = 0; g < G; g++) { 
    tmp_pi.col(g) = stick_breaking( v_j, false ) ;
  }
  
  // auxiliary variables
  arma::vec prob_j(trunc) ;
  arma::vec ids = arma::linspace(0, trunc-1, trunc) ;
  
  arma::vec tmp_sum_rows_ell(trunc) ;
  
  int njg ;
  int nmjg ;
  arma::vec nj(trunc) ;
  arma::vec nmj(trunc) ;
  arma::mat nj_mat(G, trunc, arma::fill::zeros) ;
  arma::mat nmj_mat(G, trunc, arma::fill::zeros) ;
  arma::vec max_cl_g(G, arma::fill::zeros) ; 
  
  double a_j = 0.0 ; double b_j = 0.0 ;
  arma::vec tmp_nm_ell(G) ;
  double tmp_sum_ell = 0.0 ; 
  
  int iter_thin = 0 ;
    
  // progress bar 
  bool display_progress = progressbar ;
  Progress p(nrep, display_progress) ;
  
  // START MCMC
  for(int iter = 0; iter < nrep ; iter++)
  {
    if( Progress::check_abort() ) { return -1.0 ; }
    p.increment();
    
    // compute number of observations allocated to each cluster (and to clusters >)
    nj_mat.fill(0) ; nmj_mat.fill(0) ;
    for(int g = 0; g < G; g++) {
      arma::uvec idg = arma::find( group == g ) ;  // indices of group g
      arma::vec cl_g = tmp_cl( idg ) ; // associated cluster allocation
      
      for(int j = 0; j < trunc; j++) {
        arma::uvec ind_j = find(cl_g == j) ;
        arma::uvec ind_mj = find(cl_g > j) ;
        
        nj_mat(g, j) = ind_j.n_elem ;
        nmj_mat(g, j) = ind_mj.n_elem ;
      }
      max_cl_g(g) = arma::max(cl_g) ;
    }
    
    
    // sample thinning
    for(int g = 0; g < G; g++) {
      for(int j = trunc-1; j>-1; j--){
        if(nj_mat(g, j) > 0) { 
          tmp_ell(j, g) = 1.0 ; 
        } else {
          log_p1 =  nmj_mat(g, j) * log(1.0 - v_j(j)) + log( tmp_prob_thinning(g)) ;
          log_p0 = log( 1.0 - tmp_prob_thinning(g) ) ; 
       
          p1 = exp(log_p1) / ( exp(log_p1) + exp(log_p0) ) ;

          tmp_ell(j, g) = R::rbinom( 1, p1 ) ;
        }
      }
    }
    
    // sample thinning prob for each group
    for(int g = 0; g < G; g++) {
      tmp_sum_ell = arma::accu(tmp_ell.col(g)) ;
      tmp_prob_thinning(g) = R::rbeta(a_beta + tmp_sum_ell, b_beta + trunc - tmp_sum_ell) ;
    }
    
    // sample sticks
    a_j = 0 ; b_j = 0 ;
    for(int j = 0; j < trunc; j++) {
      njg = arma::accu( nj_mat.col(j) ) ;
      for(int g = 0; g < G; g++) { tmp_nm_ell(g) = nmj_mat(g,j) * tmp_ell(j, g) ; }
      nmjg = arma::accu( tmp_nm_ell ) ;
      a_j = 1.0 + njg ;
      b_j = alpha + nmjg ;
      v_j(j) = R::rbeta(a_j, b_j) ;
    }
    log_m1v_j = log(1.0 - v_j) ;
    
    
    // reconstruct the mixing weights based on thinning and sticks
    for(int g = 0; g < G; g++) {
      thinned_vj = v_j ;
      for(int j = 0; j < trunc; j++) {
        if(tmp_ell(j, g) == 0) { thinned_vj(j) = 0 ; }
      }
      tmp_pi.col(g) = stick_breaking( thinned_vj , false );
    }

    
    // update cluster allocation
    for(int i = 0; i < N; i++)
    {
      int g = group(i) ;
      prob_j.fill(log(0.0)) ;
      
      for(int j = 0; j < trunc; j++) {
        prob_j(j) = R::dnorm( y(i), tmp_mu(j), std::sqrt(tmp_sigma2(j)), true ) + log(tmp_pi(j, g)) ;
      }
      prob_j =  prob_j - max(prob_j) ;
      prob_j = exp(prob_j) ;
      tmp_cl(i) = sample_i(ids, prob_j) ;
    }
  
    
    // update model parameters
    double new_mu0 ; double new_tau0 ;
    double new_gamma0; double new_lambda0 ;
    
    nj.fill(0) ;
    for(int j = 0; j < trunc ; j++)
    {
      arma::uvec cl_j = find( tmp_cl == j ) ;
      nj(j) = cl_j.n_elem ;
      
      if( nj(j) > 0) {
        
        double sumcl_j = arma::accu( y(cl_j) ); // sum of y_i
        
        // new parameters
        new_tau0 = tau0 + nj(j) ;
        new_mu0 = (tau0 * mu0 + sumcl_j)/new_tau0 ;
        new_gamma0 = gamma0 + nj(j) / 2.0 ;
        new_lambda0 = lambda0 + 0.5 * ( (nj(j)-1.0)*arma::var(y(cl_j)) + 
                                          (tau0*nj(j))/new_tau0 * (sumcl_j/nj(j) - mu0) * (sumcl_j/nj(j) - mu0)  ) ;
        
        tmp_sigma2(j) = 1.0/R::rgamma( new_gamma0, 1.0/new_lambda0) ;
        tmp_mu(j) = R::rnorm( new_mu0, sqrt(tmp_sigma2(j)/new_tau0) ) ;
        
      } else {
        tmp_sigma2(j) = 1.0/R::rgamma(gamma0, 1.0/lambda0) ;
        tmp_mu(j) = R::rnorm( mu0, sqrt(tmp_sigma2(j)/tau0) )  ;
      }
    }
    
    if((iter >= burnin) & (iter % thinning_factor == 0)) {
      out_mu.col(iter_thin) = tmp_mu ;
      out_sigma2.col(iter_thin) = tmp_sigma2 ;
      out_cl.col(iter_thin) = tmp_cl ;
      out_pi.slice(iter_thin) = tmp_pi ;
      out_ell.slice(iter_thin) = tmp_ell ;
      out_prob_thinning.col(iter_thin) = tmp_prob_thinning ;
      iter_thin = iter_thin + 1 ;
    }
    
    //// END 
  }
  
  arma::vec hyperpar(4) ;
  hyperpar(0) = mu0 ; hyperpar(1) = tau0 ; 
  hyperpar(2) = gamma0 ; hyperpar(3) = lambda0 ; 
  
  return Rcpp::List::create(Rcpp::Named("mu") = out_mu.t(),
                            Rcpp::Named("sigma2") = out_sigma2.t(),
                            Rcpp::Named("cl") = out_cl.t(),
                            Rcpp::Named("pi") = out_pi,
                            Rcpp::Named("ell") = out_ell,
                            Rcpp::Named("thinning_prob") = out_prob_thinning, 
                            Rcpp::Named("y") = y,
                            Rcpp::Named("group") = group,
                            Rcpp::Named("nrep") = nrep,
                            Rcpp::Named("burnin") = burnin,
                            Rcpp::Named("trunc") = trunc,
                            Rcpp::Named("hyperpar") = hyperpar,
                            Rcpp::Named("alpha") = alpha,
                            Rcpp::Named("mu_start") = mu_start,
                            Rcpp::Named("sigma2_start") = sigma2_start,
                            Rcpp::Named("cl_start") = cl_start
  );
}
