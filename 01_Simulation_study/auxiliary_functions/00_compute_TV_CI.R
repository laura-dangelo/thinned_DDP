true_density = function(data, seq)
{
  G = length(unique(data$group))
  out = matrix(0, G, length(seq))
  
  for(g in 1:G) {
    for(i in 1:length(seq)){
      out[g,i] = sum(data$weights[,g] * dnorm(seq[i], mean = data$mu_k, sd = sqrt(data$var_k)) )
    }
  }
  return(out)
}


TV_distance = function(seq, dens1, dens2) {
  sum(.5 * abs(dens1-dens2))
}


compute_densityCI_thinned = function(density_est, data) {
  seqq = density_est$seq
  
  density_CI_thinnedDDP = matrix(rep(1,length(seqq)), length(seqq), 1)
  density_CI_thinnedDDP = cbind(density_CI_thinnedDDP, seqq)
  density_CI_thinnedDDP = cbind(density_CI_thinnedDDP, apply(density_est$density_mcmc[,1,], 1, function(x) emp.hpd(x)[[1]] ))
  density_CI_thinnedDDP = cbind(density_CI_thinnedDDP, apply(density_est$density_mcmc[,1,], 1, function(x) median(x) ))
  density_CI_thinnedDDP = cbind(density_CI_thinnedDDP, apply(density_est$density_mcmc[,1,], 1, function(x) emp.hpd(x)[[2]] ))
  
  for(gg in 2:length(unique(data$group)) ) {
    tmp = matrix(rep(gg,length(seqq)), length(seqq), 1)
    tmp = cbind(tmp, seqq)
    tmp = cbind(tmp, apply(density_est$density_mcmc[,gg,], 1, function(x) emp.hpd(x)[[1]] ))
    tmp = cbind(tmp, apply(density_est$density_mcmc[,gg,], 1, function(x) median(x) ))
    tmp = cbind(tmp, apply(density_est$density_mcmc[,gg,], 1, function(x) emp.hpd(x)[[2]] ))
    
    density_CI_thinnedDDP = rbind(density_CI_thinnedDDP, tmp)
  }

  density_CI_thinnedDDP = data.frame(density_CI_thinnedDDP)
  colnames(density_CI_thinnedDDP) = c("Group", "Seq", "Q1", "Median", "Q3")
  return(density_CI_thinnedDDP)  
}


compute_densityCI_pool = function(density_est, data) {
  seqq = density_est$seq
  
  density_CI_pool = matrix(seqq, length(seqq), 1)
  density_CI_pool = cbind(density_CI_pool, apply(density_est$density_mcmc, 1, function(x) emp.hpd(x)[[1]] ))
  density_CI_pool = cbind(density_CI_pool, apply(density_est$density_mcmc, 1, function(x) median(x) ))
  density_CI_pool = cbind(density_CI_pool, apply(density_est$density_mcmc, 1, function(x) emp.hpd(x)[[2]] ))
  
  density_CI_pool = data.frame(density_CI_pool)
  colnames(density_CI_pool) = c("Seq", "Q1", "Median", "Q3")
  return(density_CI_pool)
}

compute_densityCI_nopool = function(density_est, data) {
  seqq = density_est[[1]]$seq
  
  density_CI_nopool = matrix(rep(1,length(seqq)), length(seqq), 1)
  density_CI_nopool = cbind(density_CI_nopool, seqq)
  density_CI_nopool = cbind(density_CI_nopool, apply(density_est[[1]]$density_mcmc, 1, function(x) emp.hpd(x)[[1]] ))
  density_CI_nopool = cbind(density_CI_nopool, apply(density_est[[1]]$density_mcmc, 1, function(x) median(x) ))
  density_CI_nopool = cbind(density_CI_nopool, apply(density_est[[1]]$density_mcmc, 1, function(x) emp.hpd(x)[[2]] ))
  
  for(gg in 2:length(unique(data$group))) {
    tmp = matrix(rep(gg,length(seqq)), length(seqq), 1)
    tmp = cbind(tmp, seqq)
    tmp = cbind(tmp, apply(density_est[[gg]]$density_mcmc, 1, function(x) emp.hpd(x)[[1]] ))
    tmp = cbind(tmp, apply(density_est[[gg]]$density_mcmc, 1, function(x) median(x) ))
    tmp = cbind(tmp, apply(density_est[[gg]]$density_mcmc, 1, function(x) emp.hpd(x)[[2]] ))
    
    density_CI_nopool = rbind(density_CI_nopool, tmp)
  }
  rm(tmp)
  str(density_CI_nopool)
  
  density_CI_nopool = data.frame(density_CI_nopool)
  colnames(density_CI_nopool) = c("Hosp", "Seq", "Q1", "Median", "Q3")
  return(density_CI_nopool)
}

