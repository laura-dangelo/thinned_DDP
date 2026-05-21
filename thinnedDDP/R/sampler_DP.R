sampler_DP <- function(nrep, burnin, 
                             y, 
                             trunc = 50, 
                             m0 = 0, tau0 = 0.1, 
                             gamma0 = 3, lambda0 = 2,
                             alpha = 1, 
                             mu_start, sigma2_start, cl_start, 
                             progressbar = TRUE)
{ 
  start = Sys.time()
  out = sampler_DP_arma(nrep, burnin, y, 
                               trunc, 
                               m0, tau0, 
                               gamma0, lambda0,
                               alpha, 
                               mu_start, sigma2_start, cl_start, 
                               progressbar)
  end = Sys.time()
  out$time = end-start
  return(out)
}
