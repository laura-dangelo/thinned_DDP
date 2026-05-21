sampler_thinnedDDP <- function(nrep, burnin, 
                                     y, group, 
                                     trunc = 50, 
                                     m0 = 0, tau0 = 0.1, 
                                     gamma0 = 3, lambda0 = 2,
                                     alpha = 1, 
                                     a_beta = 1, b_beta = 1, 
                                     mu_start, sigma2_start, cl_start, 
                                     progressbar = TRUE)
{ 
  start = Sys.time()
  out = sampler_thinnedDDP_arma(nrep, burnin, y, 
                                  group, trunc, 
                                  m0, tau0, 
                                  gamma0, lambda0,
                                  alpha, 
                                  a_beta, b_beta, 
                                  mu_start, sigma2_start, cl_start, 
                                  progressbar)
  end = Sys.time()
  out$time = end-start
  return(out)
}
