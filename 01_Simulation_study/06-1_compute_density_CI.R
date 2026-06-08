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
      
      #-----# #-----# #-----# #-----#
      #-----#   Thinned DDP   #-----#
      #-----# #-----# #-----# #-----#

      nameopen = paste0("01_Simulation_study/results/density_est_thinnedDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      density_est_thinnedDDP = readRDS(nameopen)

      seqq = density_est_thinnedDDP$seq
      density_CI_thinnedDDP = matrix(rep(1,length(seqq)), length(seqq), 1)
      density_CI_thinnedDDP = cbind(density_CI_thinnedDDP, seqq)
      density_CI_thinnedDDP = cbind(density_CI_thinnedDDP, apply(density_est_thinnedDDP$density_mcmc[,1,], 1, function(x) emp.hpd(x)[[1]] ))
      density_CI_thinnedDDP = cbind(density_CI_thinnedDDP, apply(density_est_thinnedDDP$density_mcmc[,1,], 1, function(x) mean(x) ))
      density_CI_thinnedDDP = cbind(density_CI_thinnedDDP, apply(density_est_thinnedDDP$density_mcmc[,1,], 1, function(x) emp.hpd(x)[[2]] ))
      density_CI_thinnedDDP = cbind(density_CI_thinnedDDP, truth[1,])


      for(gg in 2:length(unique(data$group))) {
        tmp = matrix(rep(gg,length(seqq)), length(seqq), 1)
        tmp = cbind(tmp, seqq)
        tmp = cbind(tmp, apply(density_est_thinnedDDP$density_mcmc[,gg,], 1, function(x) emp.hpd(x)[[1]] ))
        tmp = cbind(tmp, apply(density_est_thinnedDDP$density_mcmc[,gg,], 1, function(x) mean(x) ))
        tmp = cbind(tmp, apply(density_est_thinnedDDP$density_mcmc[,gg,], 1, function(x) emp.hpd(x)[[2]] ))
        tmp = cbind(tmp, truth[gg,])

        density_CI_thinnedDDP = rbind(density_CI_thinnedDDP, tmp)
      }
      rm(tmp)

      density_CI_thinnedDDP = data.frame(density_CI_thinnedDDP)
      colnames(density_CI_thinnedDDP) = c("Group", "Seq", "lower", "mean", "upper", "true")

      namesave = paste0("01_Simulation_study/results/density_CI_thinnedDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      saveRDS(density_CI_thinnedDDP, file = namesave)

      rm(density_CI_thinnedDDP)
      rm(density_est_thinnedDDP)


      #-----# #-----# #-----# #-----#
      #-----#      Pool       #-----#
      #-----# #-----# #-----# #-----#

      nameopen = paste0("01_Simulation_study/results/density_est_pool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      density_est_pool = readRDS(nameopen)

      density_CI_pool = matrix(seqq, length(seqq), 1)
      density_CI_pool = cbind(density_CI_pool, apply(density_est_pool$density_mcmc, 1, function(x) emp.hpd(x)[[1]] ))
      density_CI_pool = cbind(density_CI_pool, apply(density_est_pool$density_mcmc, 1, function(x) mean(x) ))
      density_CI_pool = cbind(density_CI_pool, apply(density_est_pool$density_mcmc, 1, function(x) emp.hpd(x)[[2]] ))

      density_CI_pool = data.frame(density_CI_pool)
      colnames(density_CI_pool) = c("Seq", "lower", "mean", "upper")

      namesave = paste0("01_Simulation_study/results/density_CI_pool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      saveRDS(density_CI_pool, file = namesave)

      rm(density_CI_pool)
      rm(density_est_pool)


      #-----# #-----# #-----# #-----#
      #-----#     No pool     #-----#
      #-----# #-----# #-----# #-----#

      nameopen = paste0("01_Simulation_study/results/density_est_nopool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      density_est_nopool = readRDS(nameopen)

      density_CI_nopool = matrix(rep(1,length(seqq)), length(seqq), 1)
      density_CI_nopool = cbind(density_CI_nopool, seqq)
      density_CI_nopool = cbind(density_CI_nopool, apply(density_est_nopool[[1]]$density_mcmc, 1, function(x) emp.hpd(x)[[1]] ))
      density_CI_nopool = cbind(density_CI_nopool, apply(density_est_nopool[[1]]$density_mcmc, 1, function(x) mean(x) ))
      density_CI_nopool = cbind(density_CI_nopool, apply(density_est_nopool[[1]]$density_mcmc, 1, function(x) emp.hpd(x)[[2]] ))

      for(gg in 2:length(unique(data$group))) {
        tmp = matrix(rep(gg,length(seqq)), length(seqq), 1)
        tmp = cbind(tmp, seqq)
        tmp = cbind(tmp, apply(density_est_nopool[[gg]]$density_mcmc, 1, function(x) emp.hpd(x)[[1]] ))
        tmp = cbind(tmp, apply(density_est_nopool[[gg]]$density_mcmc, 1, function(x) mean(x) ))
        tmp = cbind(tmp, apply(density_est_nopool[[gg]]$density_mcmc, 1, function(x) emp.hpd(x)[[2]] ))

        density_CI_nopool = rbind(density_CI_nopool, tmp)
      }
      rm(tmp)

      density_CI_nopool = data.frame(density_CI_nopool)
      colnames(density_CI_nopool) = c("Group", "Seq", "lower", "mean", "upper")

      namesave = paste0("01_Simulation_study/results/density_CI_nopool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      saveRDS(density_CI_nopool, file = namesave)

      rm(density_CI_nopool)
      rm(density_est_nopool)




      #-----# #-----# #-----# #-----#
      #-----#       CAM       #-----#
      #-----# #-----# #-----# #-----#

      nameopen = paste0("01_Simulation_study/results/density_est_CAM", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      density_est_CAM = readRDS(nameopen)

      seqq = density_est_CAM$seq
      density_CI_CAM = matrix(rep(1,length(seqq)), length(seqq), 1)
      density_CI_CAM = cbind(density_CI_CAM, seqq)
      density_CI_CAM = cbind(density_CI_CAM, apply(density_est_CAM$density_mcmc[,1,], 1, function(x) emp.hpd(x)[[1]] ))
      density_CI_CAM = cbind(density_CI_CAM, apply(density_est_CAM$density_mcmc[,1,], 1, function(x) mean(x) ))
      density_CI_CAM = cbind(density_CI_CAM, apply(density_est_CAM$density_mcmc[,1,], 1, function(x) emp.hpd(x)[[2]] ))
      density_CI_CAM = cbind(density_CI_CAM, truth[1,])


      for(gg in 2:length(unique(data$group))) {
        tmp = matrix(rep(gg,length(seqq)), length(seqq), 1)
        tmp = cbind(tmp, seqq)
        tmp = cbind(tmp, apply(density_est_CAM$density_mcmc[,gg,], 1, function(x) emp.hpd(x)[[1]] ))
        tmp = cbind(tmp, apply(density_est_CAM$density_mcmc[,gg,], 1, function(x) mean(x) ))
        tmp = cbind(tmp, apply(density_est_CAM$density_mcmc[,gg,], 1, function(x) emp.hpd(x)[[2]] ))
        tmp = cbind(tmp, truth[gg,])

        density_CI_CAM = rbind(density_CI_CAM, tmp)
      }
      rm(tmp)

      density_CI_CAM = data.frame(density_CI_CAM)
      colnames(density_CI_CAM) = c("Group", "Seq", "lower", "mean", "upper", "true")

      namesave = paste0("01_Simulation_study/results/density_CI_CAM", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      saveRDS(density_CI_CAM, file = namesave)

      rm(density_CI_CAM)
      rm(density_est_CAM)





      #-----# #-----# #-----# #-----#
      #-----#      GM-DDP     #-----#
      #-----# #-----# #-----# #-----#

      nameopen = paste0("01_Simulation_study/results/density_est_gmDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      density_est_gmDDP = readRDS(nameopen)

      seqq = density_est_gmDDP$seq
      density_CI_gmDDP = matrix(rep(1,length(seqq)), length(seqq), 1)
      density_CI_gmDDP = cbind(density_CI_gmDDP, seqq)
      density_CI_gmDDP = cbind(density_CI_gmDDP, apply(density_est_gmDDP$density_mcmc[,1,], 1, function(x) emp.hpd(x)[[1]] ))
      density_CI_gmDDP = cbind(density_CI_gmDDP, apply(density_est_gmDDP$density_mcmc[,1,], 1, function(x) mean(x) ))
      density_CI_gmDDP = cbind(density_CI_gmDDP, apply(density_est_gmDDP$density_mcmc[,1,], 1, function(x) emp.hpd(x)[[2]] ))
      density_CI_gmDDP = cbind(density_CI_gmDDP, truth[1,])


      for(gg in 2:length(unique(data$group))) {
        tmp = matrix(rep(gg,length(seqq)), length(seqq), 1)
        tmp = cbind(tmp, seqq)
        tmp = cbind(tmp, apply(density_est_gmDDP$density_mcmc[,gg,], 1, function(x) emp.hpd(x)[[1]] ))
        tmp = cbind(tmp, apply(density_est_gmDDP$density_mcmc[,gg,], 1, function(x) mean(x) ))
        tmp = cbind(tmp, apply(density_est_gmDDP$density_mcmc[,gg,], 1, function(x) emp.hpd(x)[[2]] ))
        tmp = cbind(tmp, truth[gg,])

        density_CI_gmDDP = rbind(density_CI_gmDDP, tmp)
      }
      rm(tmp)

      density_CI_gmDDP = data.frame(density_CI_gmDDP)
      colnames(density_CI_gmDDP) = c("Group", "Seq", "lower", "mean", "upper", "true")

      namesave = paste0("01_Simulation_study/results/density_CI_gmDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      saveRDS(density_CI_gmDDP, file = namesave)

      rm(density_CI_gmDDP)
      rm(density_est_gmDDP)







      
      #-----# #-----# #-----# #-----#
      #-----#       HDP       #-----#
      #-----# #-----# #-----# #-----# 
      
      nameopen = paste0("01_Simulation_study/results/density_est_HDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      density_est_HDP = readRDS(nameopen)
      
      seqq = density_est_HDP$seq
      density_CI_HDP = matrix(rep(1,length(seqq)), length(seqq), 1)
      density_CI_HDP = cbind(density_CI_HDP, seqq)
      density_CI_HDP = cbind(density_CI_HDP, apply(density_est_HDP$density_mcmc[,1,], 1, function(x) emp.hpd(x)[[1]] ))
      density_CI_HDP = cbind(density_CI_HDP, apply(density_est_HDP$density_mcmc[,1,], 1, function(x) mean(x) ))
      density_CI_HDP = cbind(density_CI_HDP, apply(density_est_HDP$density_mcmc[,1,], 1, function(x) emp.hpd(x)[[2]] ))
      density_CI_HDP = cbind(density_CI_HDP, truth[1,])
      
      
      for(gg in 2:length(unique(data$group))) {
        tmp = matrix(rep(gg,length(seqq)), length(seqq), 1)
        tmp = cbind(tmp, seqq)
        tmp = cbind(tmp, apply(density_est_HDP$density_mcmc[,gg,], 1, function(x) emp.hpd(x)[[1]] ))
        tmp = cbind(tmp, apply(density_est_HDP$density_mcmc[,gg,], 1, function(x) mean(x) ))
        tmp = cbind(tmp, apply(density_est_HDP$density_mcmc[,gg,], 1, function(x) emp.hpd(x)[[2]] ))
        tmp = cbind(tmp, truth[gg,])
        
        density_CI_HDP = rbind(density_CI_HDP, tmp)
      }
      rm(tmp)
      
      density_CI_HDP = data.frame(density_CI_HDP)
      colnames(density_CI_HDP) = c("Group", "Seq", "lower", "mean", "upper", "true")
      
      namesave = paste0("01_Simulation_study/results/density_CI_HDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      saveRDS(density_CI_HDP, file = namesave)
      
      rm(density_CI_HDP)
      rm(density_est_HDP)
      
      
      
    }
  }
}











