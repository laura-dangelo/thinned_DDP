library(TeachingDemos)

n_groups = c(2, 10)
ssg = c(10, 30)
n_ss = 4 # number of different configuration of sample size
n_datasets = 50

tot_datasets = n_datasets * n_ss * length(n_groups)



for(repl in 1:n_datasets) {
  for(i in 1:length(n_groups) ){
    for(j in 1:n_ss){ 
      
      name_file_open = paste0("01_Simulation_study/data/data_", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      data = readRDS(file = name_file_open)
      
      nameopen = paste0("01_Simulation_study/results/true_density_", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      truth = readRDS(nameopen)
      
      # #-----# #-----# #-----# #-----#
      # #-----#   Thinned DDP   #-----#
      # #-----# #-----# #-----# #-----#
      # 
      # nameopen = paste0("01_Simulation_study/results/density_est_thinnedDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      # density_est_thinnedDDP = readRDS(nameopen)
      # 
      # seqq = density_est_thinnedDDP$seq
      # density_CI_thinnedDDP = matrix(rep(1,length(seqq)), length(seqq), 1)
      # density_CI_thinnedDDP = cbind(density_CI_thinnedDDP, seqq)
      # density_CI_thinnedDDP = cbind(density_CI_thinnedDDP, apply(density_est_thinnedDDP$density_mcmc[,1,], 1, function(x) emp.hpd(x)[[1]] ))
      # density_CI_thinnedDDP = cbind(density_CI_thinnedDDP, apply(density_est_thinnedDDP$density_mcmc[,1,], 1, function(x) mean(x) ))
      # density_CI_thinnedDDP = cbind(density_CI_thinnedDDP, apply(density_est_thinnedDDP$density_mcmc[,1,], 1, function(x) emp.hpd(x)[[2]] ))
      # density_CI_thinnedDDP = cbind(density_CI_thinnedDDP, truth[1,])
      # 
      # 
      # for(gg in 2:length(unique(data$group))) {
      #   tmp = matrix(rep(gg,length(seqq)), length(seqq), 1)
      #   tmp = cbind(tmp, seqq)
      #   tmp = cbind(tmp, apply(density_est_thinnedDDP$density_mcmc[,gg,], 1, function(x) emp.hpd(x)[[1]] ))
      #   tmp = cbind(tmp, apply(density_est_thinnedDDP$density_mcmc[,gg,], 1, function(x) mean(x) ))
      #   tmp = cbind(tmp, apply(density_est_thinnedDDP$density_mcmc[,gg,], 1, function(x) emp.hpd(x)[[2]] ))
      #   tmp = cbind(tmp, truth[gg,])
      # 
      #   density_CI_thinnedDDP = rbind(density_CI_thinnedDDP, tmp)
      # }
      # rm(tmp)
      # 
      # density_CI_thinnedDDP = data.frame(density_CI_thinnedDDP)
      # colnames(density_CI_thinnedDDP) = c("Group", "Seq", "lower", "mean", "upper", "true")
      # 
      # namesave = paste0("01_Simulation_study/results/density_CI_thinnedDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      # saveRDS(density_CI_thinnedDDP, file = namesave)
      # 
      # rm(density_CI_thinnedDDP)
      # rm(density_est_thinnedDDP)


      #-----#  #-----#  #-----#  #-----#  #-----#
      #-----#        UNIFORM PRIOR        #-----#
      #-----#  #-----#  #-----#  #-----#  #-----#
      
      namesave = paste0("03_Sensitivity_study/results_sensitivity/density_est_thinnedDDP_unif", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      density_est_thinnedDDP_unif = readRDS(nameopen)

      seqq = density_est_thinnedDDP_unif$seq
      density_CI_thinnedDDP_unif = matrix(rep(1,length(seqq)), length(seqq), 1)
      density_CI_thinnedDDP_unif = cbind(density_CI_thinnedDDP_unif, seqq)
      density_CI_thinnedDDP_unif = cbind(density_CI_thinnedDDP_unif, apply(density_est_thinnedDDP_unif$density_mcmc[,1,], 1, function(x) emp.hpd(x)[[1]] ))
      density_CI_thinnedDDP_unif = cbind(density_CI_thinnedDDP_unif, apply(density_est_thinnedDDP_unif$density_mcmc[,1,], 1, function(x) mean(x) ))
      density_CI_thinnedDDP_unif = cbind(density_CI_thinnedDDP_unif, apply(density_est_thinnedDDP_unif$density_mcmc[,1,], 1, function(x) emp.hpd(x)[[2]] ))
      density_CI_thinnedDDP_unif = cbind(density_CI_thinnedDDP_unif, truth[1,])
      
      
      for(gg in 2:length(unique(data$group))) {
        tmp = matrix(rep(gg,length(seqq)), length(seqq), 1)
        tmp = cbind(tmp, seqq)
        tmp = cbind(tmp, apply(density_est_thinnedDDP_unif$density_mcmc[,gg,], 1, function(x) emp.hpd(x)[[1]] ))
        tmp = cbind(tmp, apply(density_est_thinnedDDP_unif$density_mcmc[,gg,], 1, function(x) mean(x) ))
        tmp = cbind(tmp, apply(density_est_thinnedDDP_unif$density_mcmc[,gg,], 1, function(x) emp.hpd(x)[[2]] ))
        tmp = cbind(tmp, truth[gg,])
        
        density_CI_thinnedDDP_unif = rbind(density_CI_thinnedDDP_unif, tmp)
      }
      rm(tmp)
      density_CI_thinnedDDP_unif = data.frame(density_CI_thinnedDDP_unif)
      colnames(density_CI_thinnedDDP_unif) = c("Seq", "lower", "mean", "upper")

      namesave = paste0("03_Sensitivity_study/results/density_CI_thinnedDDP_unif", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      saveRDS(density_CI_thinnedDDP_unif, file = namesave)

      rm(density_CI_thinnedDDP_unif)
      rm(density_est_thinnedDDP_unif)

      
      
      #-----#  #-----#  #-----#  #-----#  #-----#
      #-----#         BETA(10,10)         #-----#
      #-----#  #-----#  #-----#  #-----#  #-----#
      
      namesave = paste0("03_Sensitivity_study/results_sensitivity/density_est_thinnedDDP_1010", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      density_est_thinnedDDP_1010 = readRDS(nameopen)
      
      seqq = density_est_thinnedDDP_1010$seq
      density_CI_thinnedDDP_1010 = matrix(rep(1,length(seqq)), length(seqq), 1)
      density_CI_thinnedDDP_1010 = cbind(density_CI_thinnedDDP_1010, seqq)
      density_CI_thinnedDDP_1010 = cbind(density_CI_thinnedDDP_1010, apply(density_est_thinnedDDP_1010$density_mcmc[,1,], 1, function(x) emp.hpd(x)[[1]] ))
      density_CI_thinnedDDP_1010 = cbind(density_CI_thinnedDDP_1010, apply(density_est_thinnedDDP_1010$density_mcmc[,1,], 1, function(x) mean(x) ))
      density_CI_thinnedDDP_1010 = cbind(density_CI_thinnedDDP_1010, apply(density_est_thinnedDDP_1010$density_mcmc[,1,], 1, function(x) emp.hpd(x)[[2]] ))
      density_CI_thinnedDDP_1010 = cbind(density_CI_thinnedDDP_1010, truth[1,])
      
      
      for(gg in 2:length(unique(data$group))) {
        tmp = matrix(rep(gg,length(seqq)), length(seqq), 1)
        tmp = cbind(tmp, seqq)
        tmp = cbind(tmp, apply(density_est_thinnedDDP_1010$density_mcmc[,gg,], 1, function(x) emp.hpd(x)[[1]] ))
        tmp = cbind(tmp, apply(density_est_thinnedDDP_1010$density_mcmc[,gg,], 1, function(x) mean(x) ))
        tmp = cbind(tmp, apply(density_est_thinnedDDP_1010$density_mcmc[,gg,], 1, function(x) emp.hpd(x)[[2]] ))
        tmp = cbind(tmp, truth[gg,])
        
        density_CI_thinnedDDP_1010 = rbind(density_CI_thinnedDDP_1010, tmp)
      }
      rm(tmp)
      density_CI_thinnedDDP_1010 = data.frame(density_CI_thinnedDDP_1010)
      colnames(density_CI_thinnedDDP_1010) = c("Seq", "lower", "mean", "upper")
      
      namesave = paste0("03_Sensitivity_study/results/density_CI_thinnedDDP_1010", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      saveRDS(density_CI_thinnedDDP_1010, file = namesave)
      
      rm(density_CI_thinnedDDP_1010)
      rm(density_est_thinnedDDP_1010)
      
      
      
      
      
      #-----#  #-----#  #-----#  #-----#  #-----#
      #-----#         BETA(0.1,0.1)         #-----#
      #-----#  #-----#  #-----#  #-----#  #-----#
      
      namesave = paste0("03_Sensitivity_study/results_sensitivity/density_est_thinnedDDP_0101", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      density_est_thinnedDDP_0101 = readRDS(nameopen)
      
      seqq = density_est_thinnedDDP_0101$seq
      density_CI_thinnedDDP_0101 = matrix(rep(1,length(seqq)), length(seqq), 1)
      density_CI_thinnedDDP_0101 = cbind(density_CI_thinnedDDP_0101, seqq)
      density_CI_thinnedDDP_0101 = cbind(density_CI_thinnedDDP_0101, apply(density_est_thinnedDDP_0101$density_mcmc[,1,], 1, function(x) emp.hpd(x)[[1]] ))
      density_CI_thinnedDDP_0101 = cbind(density_CI_thinnedDDP_0101, apply(density_est_thinnedDDP_0101$density_mcmc[,1,], 1, function(x) mean(x) ))
      density_CI_thinnedDDP_0101 = cbind(density_CI_thinnedDDP_0101, apply(density_est_thinnedDDP_0101$density_mcmc[,1,], 1, function(x) emp.hpd(x)[[2]] ))
      density_CI_thinnedDDP_0101 = cbind(density_CI_thinnedDDP_0101, truth[1,])
      
      
      for(gg in 2:length(unique(data$group))) {
        tmp = matrix(rep(gg,length(seqq)), length(seqq), 1)
        tmp = cbind(tmp, seqq)
        tmp = cbind(tmp, apply(density_est_thinnedDDP_0101$density_mcmc[,gg,], 1, function(x) emp.hpd(x)[[1]] ))
        tmp = cbind(tmp, apply(density_est_thinnedDDP_0101$density_mcmc[,gg,], 1, function(x) mean(x) ))
        tmp = cbind(tmp, apply(density_est_thinnedDDP_0101$density_mcmc[,gg,], 1, function(x) emp.hpd(x)[[2]] ))
        tmp = cbind(tmp, truth[gg,])
        
        density_CI_thinnedDDP_0101 = rbind(density_CI_thinnedDDP_0101, tmp)
      }
      rm(tmp)
      density_CI_thinnedDDP_0101 = data.frame(density_CI_thinnedDDP_0101)
      colnames(density_CI_thinnedDDP_0101) = c("Seq", "lower", "mean", "upper")
      
      namesave = paste0("03_Sensitivity_study/results/density_CI_thinnedDDP_0101", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      saveRDS(density_CI_thinnedDDP_0101, file = namesave)
      
      rm(density_CI_thinnedDDP_0101)
      rm(density_est_thinnedDDP_0101)

      
      
    }
  }
}











