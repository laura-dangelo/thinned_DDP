
# # # not for use - part of the next function
# input: sequence of beta random variables
# output: SB weights
stick_breaking_construction = function(betas)
{
  trunc = length(betas)
  out = numeric(trunc)
  out[1] = betas[1]
  out[2:trunc] = sapply(2:trunc, function(x) betas[x]*prod((1-betas)[1:(x-1)]) )
  out[trunc] = 1-sum(out[1:(trunc-1)])
  return(out)
}

# input: alpha parameter of the beta r.v.'s
# output: SB weights
stick_breaking = function(alpha, trunc) {
  betas = rbeta(trunc, 1, alpha)
  stick_breaking_construction(betas)
}

# input:
## - sample_size: number of sampled values
## - alpha: DP concentration parameter
## - rP0: distribution to sample the i.i.d. atoms
## - trunc: truncation value
sim_DP = function(sample_size, alpha, 
                  rP0=NULL, atoms=NULL, 
                  trunc, seed = NULL, ...) {
  
  out = list()
  if(is.null(seed)) { seed = .Random.seed[1] }
  out$seed = seed
  
  set.seed(seed) 
  if(!is.null(atoms)) {
    loc = atoms[1:trunc] } else if(!is.null(rP0)) {
      loc = rP0(trunc, ...) } else { errorCondition("no P0 nor atoms") }
  
  weights = stick_breaking(alpha, trunc)
  
  out$weights = weights # pi_1, pi_2, pi_3, ..., pi_trunc
  out$loc = loc # theta_1, theta_2, ..., theta_trunc
  
  out$cl = sample(1:trunc, sample_size, replace = TRUE, prob = weights)
  out$sim = loc[out$cl]
  return(out)
}

# sim_DP(10, 1, rgamma, 100, shape = 1, scale = 1)


# input:
## - sample_sizes: number of sampled values for each group
## - thinning (non-stochastic) matrix (G x trunc)
## - alpha: DP concentration parameter
## - mu_k: vector length trunc with cluster centers (means)
## - var_k: vector length trunc with cluster centers (variances)
sim_thinned_mixtureN = function(sample_sizes, 
                                weights, 
                                mu_k, var_k, seed = NULL)
{
  out = list()
  G = ncol(weights) 
  trunc = nrow(weights)
  
  if(is.null(seed)) { seed = .Random.seed[1] }
  out$seed = seed
  set.seed(seed)
  
  out$cl = unlist(sapply(1:G, function(x) sample( c(1:trunc), sample_sizes[x], replace = TRUE, prob = weights[,x]) ))

  out$y = rnorm(length(out$cl), mean = mu_k[out$cl], sd = sqrt(var_k[out$cl]) )
  out$group = unlist(sapply(1:G, function(x) rep(x, sample_sizes[x])))
  
  out$weights = weights
  out$mu_k = mu_k
  out$var_k = var_k
  return(out)
}
