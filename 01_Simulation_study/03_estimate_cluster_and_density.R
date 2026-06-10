#----------------# #----------------# #----------------# 
#             ESTIMATE CLUSTER AND DENSITY             #
#----------------# #----------------# #----------------# 

# This script takes as input the MCMC runs and estimates the partition and posterior group-specific densities.
# The results are then saved into the "01_Simulation_study/results" folder.


library(salso)
library(mclust)


# number of groups and sample sizes (to import and save data)
n_groups = c(2, 10)
ssg = c(10, 30)
n_ss = 4 # number of different configuration of sample size
n_datasets = 50
tot_datasets = n_datasets * n_ss * length(n_groups)
trunc = 50


# start loop for importing the data
for(repl in 1:n_datasets) {                            
  for(i in 1:length(n_groups) ){
    for(j in 1:n_ss){
      
      name_file_open = paste0("01_Simulation_study/data/data_", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      data = readRDS(file = name_file_open)
      
      cat(paste0("Replication ", repl, "\n"))
      cat(paste0("i = ", i, "/", length(n_groups), "\n"))
      cat(paste0("j = ", j, "/", n_ss, "\n"))


      #-----#  #-----#  #-----#  #-----#  #-----#
      #-----#          THINNED DDP        #-----#
      #-----#  #-----#  #-----#  #-----#  #-----#

      nameopen = paste0("01_Simulation_study/results/run_thinnedDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      run_gibbs_thinnedDDP = readRDS(file = nameopen)


      #-----# Estimate partition of observations #-----#
      cl_point_est_thinnedDDP = salso::salso((run_gibbs_thinnedDDP$cl+1), nCores = 3 )
      namesave = paste0("01_Simulation_study/results/est_cl_thinnedDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      saveRDS(cl_point_est_thinnedDDP, file = namesave)

      rand_thinnedDDP = 0
      for(gg in 1:length(unique(data$group))){
          rand_thinnedDDP = rand_thinnedDDP + adjustedRandIndex(cl_point_est_thinnedDDP[data$group==gg], data$cl[data$group==gg])
      }
      rand_thinnedDDP = rand_thinnedDDP/(length(unique(data$group)))

      namesave = paste0("01_Simulation_study/results/rand_thinnedDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      saveRDS(rand_thinnedDDP, file = namesave)


      #-----# #-----# Compute density estimate #-----# #-----#
      seqq = seq(range(data$y)[1]-2, range(data$y)[2]+2, length.out = 300)
      density_est_thinnedDDP = thinnedDDP::compute_density(seqq,
                                                      weight = run_gibbs_thinnedDDP$pi,
                                                      means = t(run_gibbs_thinnedDDP$mu),
                                                      variances = t(run_gibbs_thinnedDDP$sigma2)
                                                      )
      namesave = paste0("01_Simulation_study/results/density_est_thinnedDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      saveRDS(density_est_thinnedDDP, file = namesave)


      rm(run_gibbs_thinnedDDP)
      rm(cl_point_est_thinnedDDP)
      rm(density_est_thinnedDDP)
      rm(rand_thinnedDDP)




      #-----#  #-----#  #-----#  #-----#  #-----#
      #-----#       COMPLETE POOLING      #-----#
      #-----#  #-----#  #-----#  #-----#  #-----#

      nameopen = paste0("01_Simulation_study/results/run_pool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      run_gibbs_pool = readRDS(file = nameopen)

      #-----# Estimate partition of observations #-----#
      cl_point_est_pool = salso::salso((run_gibbs_pool$cl+1), nCores = 3 )
      namesave = paste0("01_Simulation_study/results/est_cl_pool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      saveRDS(cl_point_est_pool, file = namesave)

      rand_pool = 0
      for(gg in 1:length(unique(data$group))){
            rand_pool = rand_pool + adjustedRandIndex(cl_point_est_pool[data$group==gg], data$cl[data$group==gg])
      }
      rand_pool = rand_pool/(length(unique(data$group)))
      namesave = paste0("01_Simulation_study/results/rand_pool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      saveRDS(rand_pool, file = namesave)


      #-----# #-----# Compute density estimate #-----# #-----#
      seqq = seq(range(data$y)[1]-2, range(data$y)[2]+2, length.out = 300)
      density_est_pool = thinnedDDP::compute_density_1DP(seqq,
                                                          weight = t(run_gibbs_pool$pi),
                                                          means = t(run_gibbs_pool$mu),
                                                          variances = t(run_gibbs_pool$sigma2) )
      namesave = paste0("01_Simulation_study/results/density_est_pool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      saveRDS(density_est_pool, file = namesave)


      rm(run_gibbs_pool)
      rm(cl_point_est_pool)
      rm(density_est_pool)
      rm(rand_pool)




      #-----#  #-----#  #-----#  #-----#
      #-----#       NO POOL      #-----#
      #-----#  #-----#  #-----#  #-----#

      nameopen = paste0("01_Simulation_study/results/run_nopool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      run_gibbs_nopool = readRDS(file = nameopen)

      #-----# Estimate partition of observations #-----#
      cl_point_est_nopool = list()
      for(gg in 1:length(unique(data$group))){
        cl_point_est_nopool[[gg]] = salso::salso((run_gibbs_nopool[[gg]]$cl+1), nCores = 3 )
      }
      namesave = paste0("01_Simulation_study/results/est_cl_nopool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      saveRDS(cl_point_est_nopool, file = namesave)

      rand_nopool = 0
      for(gg in 1:length(unique(data$group))){
          rand_nopool = rand_nopool + adjustedRandIndex(cl_point_est_nopool[[gg]], data$cl[data$group==gg])
      }
      rand_nopool = rand_nopool/(length(unique(data$group)))
      namesave = paste0("01_Simulation_study/results/rand_nopool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      saveRDS(rand_nopool, file = namesave)

      #-----# #-----# Compute density estimate #-----# #-----#
      density_est_nopool = list()
      seqq = seq(range(data$y)[1]-2, range(data$y)[2]+2, length.out = 300)
      for(gg in 1:length(unique(data$group))){
            density_est_nopool[[gg]] = thinnedDDP::compute_density_1DP(seqq,
                                                                  weight = t(run_gibbs_nopool[[gg]]$pi),
                                                                  means = t(run_gibbs_nopool[[gg]]$mu),
                                                                  variances = t(run_gibbs_nopool[[gg]]$sigma2) )
      }
      namesave = paste0("01_Simulation_study/results/density_est_nopool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      saveRDS(density_est_nopool, file = namesave)


      rm(run_gibbs_nopool)
      rm(cl_point_est_nopool)
      rm(density_est_nopool)
      rm(rand_nopool)





      #-----#  #-----#  #-----#  #-----#
      #-----#         CAM        #-----#
      #-----#  #-----#  #-----#  #-----#

      nameopen = paste0("01_Simulation_study/results/run_CAM", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      run_gibbs_CAM = readRDS(file = nameopen)

      #-----# Estimate partition of observations #-----#
      cl_point_est_CAM = salso::salso((run_gibbs_CAM$sim$obs_cluster), nCores = 3 )
      namesave = paste0("01_Simulation_study/results/est_cl_CAM", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      saveRDS(cl_point_est_CAM, file = namesave)

      rand_CAM = 0
      for(gg in 1:length(unique(data$group))){
          rand_CAM = rand_CAM + adjustedRandIndex(cl_point_est_CAM[data$group==gg], data$cl[data$group==gg])
      }
      rand_CAM = rand_CAM/(length(unique(data$group)))
      namesave = paste0("01_Simulation_study/results/rand_CAM", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      saveRDS(rand_CAM, file = namesave)


      #-----# #-----# Compute density estimate #-----# #-----#
      seqq = seq(range(data$y)[1]-2, range(data$y)[2]+2, length.out = 300)

      weights_CAM = array(data=NA, dim=c(trunc, length(unique(data$group)), nrow(run_gibbs_CAM$sim$mu)) )
      for(gg in 1:length(unique(data$group))){
        for(iter in 1:nrow(run_gibbs_CAM$sim$mu))
            weights_CAM[,gg,iter] = run_gibbs_CAM$sim$omega[,run_gibbs_CAM$sim$distr_cluster[iter,gg],iter]
      }

      density_est_CAM = thinnedDDP::compute_density(seqq,
                                                   weight = weights_CAM,
                                                   means = t(run_gibbs_CAM$sim$mu),
                                                   variances = t(run_gibbs_CAM$sim$sigma2)
                                                   )

      namesave = paste0("01_Simulation_study/results/density_est_CAM", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      saveRDS(density_est_CAM, file = namesave)

      rm(run_gibbs_CAM)
      rm(cl_point_est_CAM)
      rm(density_est_CAM)
      rm(rand_CAM)




      #-----#  #-----#  #-----#  #-----#
      #-----#       GM-DDP       #-----#
      #-----#  #-----#  #-----#  #-----#

      nameopen = paste0("01_Simulation_study/results/run_gmDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      run_gibbs_gmDDP = readRDS(file = nameopen)

      #-----# Estimate partition of observations #-----#
      cl_point_est_gmDDP = salso::salso((run_gibbs_gmDDP$clust), nCores = 3 )
      namesave = paste0("01_Simulation_study/results/est_cl_gmDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      saveRDS(cl_point_est_gmDDP, file = namesave)

      rand_gmDDP = 0
      for(gg in 1:length(unique(data$group))){
          rand_gmDDP = rand_gmDDP + adjustedRandIndex(cl_point_est_gmDDP[data$group==gg], data$cl[data$group==gg])
      }
      rand_gmDDP = rand_gmDDP/(length(unique(data$group)))
      namesave = paste0("01_Simulation_study/results/rand_gmDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      saveRDS(rand_gmDDP, file = namesave)

      #-----# #-----# Compute density estimate #-----# #-----#

      mean_dens = matrix(NA,300,length(unique(data$group)))
      for(gg in 1:length(unique(data$group))) {
          mean_dens[,gg] = rowMeans(run_gibbs_gmDDP$density[,gg,])
      }

      density_est_gmDDP = list(density_mcmc = run_gibbs_gmDDP$density, mean = mean_dens, seq = seqq)
      namesave = paste0("01_Simulation_study/results/density_est_gmDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      saveRDS(density_est_gmDDP, file = namesave)
      rm(mean_dens)

      rm(run_gibbs_gmDDP)
      rm(cl_point_est_gmDDP)
      rm(density_est_gmDDP)
      rm(rand_gmDDP)







      #-----#  #-----#  #-----#  #-----#
      #-----#         HDP        #-----#
      #-----#  #-----#  #-----#  #-----#

      nameopen = paste0("01_Simulation_study/results/run_HDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      run_gibbs_HDP = readRDS(file = nameopen)
      seqq = seq(range(data$y)[1]-2, range(data$y)[2]+2, length.out = 300)

      #-----# Estimate partition of observations #-----#

      tmp_cluster_alloc_HDP = matrix(0, length(run_gibbs_HDP), length(data$y))
      for(ind in 1:length(run_gibbs_HDP)) { tmp_cluster_alloc_HDP[ind,] = unlist(run_gibbs_HDP[[ind]]$Z) }

      cl_point_est_HDP = salso::salso((tmp_cluster_alloc_HDP), nCores = 3 )
      namesave = paste0("01_Simulation_study/results/est_cl_HDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      saveRDS(cl_point_est_HDP, file = namesave)

      rand_HDP = 0
      for(gg in 1:length(unique(data$group))){
        rand_HDP = rand_HDP + adjustedRandIndex(cl_point_est_HDP[data$group==gg], data$cl[data$group==gg])
      }
      rand_HDP = rand_HDP/(length(unique(data$group)))
      namesave = paste0("01_Simulation_study/results/rand_HDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      saveRDS(rand_HDP, file = namesave)

      #-----# #-----# Compute density estimate #-----# #-----#

      HDP_density = array(NA, dim = c(300, length(unique(data$group)), 2500))
      for(ind in 1:length(run_gibbs_HDP)) {
        for(gg in 1:length(unique(data$group))) {
          HDP_density[,gg,ind] = run_gibbs_HDP[[ind]]$density[gg,]
        }
      }

      mean_dens = matrix(NA,300,length(unique(data$group)))
      for(gg in 1:length(unique(data$group))) {
          mean_dens[,gg] = rowMeans(HDP_density[,gg,])
      }

      density_est_HDP = list(density_mcmc = HDP_density, mean = mean_dens, seq = seqq)
      namesave = paste0("01_Simulation_study/results/density_est_HDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      saveRDS(density_est_HDP, file = namesave)
      rm(mean_dens)

      rm(run_gibbs_HDP)
      rm(cl_point_est_HDP)
      rm(density_est_HDP)
      rm(rand_HDP)


    }
  }
}
