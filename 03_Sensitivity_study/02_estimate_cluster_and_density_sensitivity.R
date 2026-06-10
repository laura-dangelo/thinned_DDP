#----------------# #----------------# #----------------# 
#             ESTIMATE CLUSTER AND DENSITY             #
#----------------# #----------------# #----------------# 

# This script takes as input the MCMC runs and estimates the partition and posterior group-specific densities.
# The results are then saved into the "03_Sensitivity_study/results_sensitivity" folder.


library(salso)
library(mclust)


# number of groups and sample sizes (to import and save data)
n_groups = c(2, 10)
ssg = c(10, 30)
n_ss = 4 # number of different configuration of sample size
n_datasets = 25
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
      #-----#        UNIFORM PRIOR        #-----#
      #-----#  #-----#  #-----#  #-----#  #-----#

      nameopen = paste0("03_Sensitivity_study/results_sensitivity/run_thinnedDDP_unif", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      run_gibbs_thinnedDDP_unif = readRDS(file = nameopen)


      #-----# Estimate partition of observations #-----#
      cl_point_est_thinnedDDP_unif = salso::salso((run_gibbs_thinnedDDP_unif$cl+1), nCores = 3 )
      namesave = paste0("03_Sensitivity_study/results_sensitivity/est_cl_thinnedDDP_unif", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      saveRDS(cl_point_est_thinnedDDP_unif, file = namesave)

      rand_thinnedDDP_unif = 0
      for(gg in 1:length(unique(data$group))){
          rand_thinnedDDP_unif = rand_thinnedDDP_unif + adjustedRandIndex(cl_point_est_thinnedDDP_unif[data$group==gg], data$cl[data$group==gg])
      }
      rand_thinnedDDP_unif = rand_thinnedDDP_unif/(length(unique(data$group)))

      namesave = paste0("03_Sensitivity_study/results_sensitivity/rand_thinnedDDP_unif", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      saveRDS(rand_thinnedDDP_unif, file = namesave)


      #-----# #-----# Compute density estimate #-----# #-----#
      seqq = seq(range(data$y)[1]-2, range(data$y)[2]+2, length.out = 300)
      density_est_thinnedDDP_unif = thinnedDDP::compute_density(seqq,
                                                      weight = run_gibbs_thinnedDDP_unif$pi,
                                                      means = t(run_gibbs_thinnedDDP_unif$mu),
                                                      variances = t(run_gibbs_thinnedDDP_unif$sigma2)
                                                      )
      namesave = paste0("03_Sensitivity_study/results_sensitivity/density_est_thinnedDDP_unif", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      saveRDS(density_est_thinnedDDP_unif, file = namesave)


      rm(run_gibbs_thinnedDDP_unif)
      rm(cl_point_est_thinnedDDP_unif)
      rm(density_est_thinnedDDP_unif)
      rm(rand_thinnedDDP_unif)


      
      
      #-----#  #-----#  #-----#  #-----#  #-----#
      #-----#         BETA(10,10)         #-----#
      #-----#  #-----#  #-----#  #-----#  #-----#
      
      nameopen = paste0("03_Sensitivity_study/results_sensitivity/run_thinnedDDP_1010", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      run_gibbs_thinnedDDP_1010 = readRDS(file = nameopen)
      

      #-----# Estimate partition of observations #-----#
      cl_point_est_thinnedDDP_1010 = salso::salso((run_gibbs_thinnedDDP_1010$cl+1), nCores = 3 )
      namesave = paste0("03_Sensitivity_study/results_sensitivity/est_cl_thinnedDDP_1010", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      saveRDS(cl_point_est_thinnedDDP_1010, file = namesave)

      rand_thinnedDDP_1010 = 0
      for(gg in 1:length(unique(data$group))){
        rand_thinnedDDP_1010 = rand_thinnedDDP_1010 + adjustedRandIndex(cl_point_est_thinnedDDP_1010[data$group==gg], data$cl[data$group==gg])
      }
      rand_thinnedDDP_1010 = rand_thinnedDDP_1010/(length(unique(data$group)))

      namesave = paste0("03_Sensitivity_study/results_sensitivity/rand_thinnedDDP_1010", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      saveRDS(rand_thinnedDDP_1010, file = namesave)

      
      #-----# #-----# Compute density estimate #-----# #-----#
      seqq = seq(range(data$y)[1]-2, range(data$y)[2]+2, length.out = 300)
      density_est_thinnedDDP_1010 = thinnedDDP::compute_density(seqq,
                                                           weight = run_gibbs_thinnedDDP_1010$pi,
                                                           means = t(run_gibbs_thinnedDDP_1010$mu),
                                                           variances = t(run_gibbs_thinnedDDP_1010$sigma2)
      )
      namesave = paste0("03_Sensitivity_study/results_sensitivity/density_est_thinnedDDP_1010", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      saveRDS(density_est_thinnedDDP_1010, file = namesave)
      
      
      rm(run_gibbs_thinnedDDP_1010)
      rm(cl_point_est_thinnedDDP_1010)
      rm(density_est_thinnedDDP_1010)
      rm(rand_thinnedDDP_1010)
      
      
      
      
      
      
      #-----#  #-----#  #-----#  #-----#  #-----#
      #-----#       BETA(0.1,0.1)         #-----#
      #-----#  #-----#  #-----#  #-----#  #-----#
      
      nameopen = paste0("03_Sensitivity_study/results_sensitivity/run_thinnedDDP_0101", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      run_gibbs_thinnedDDP_0101 = readRDS(file = nameopen)
      

      #-----# Estimate partition of observations #-----#
      cl_point_est_thinnedDDP_0101 = salso::salso((run_gibbs_thinnedDDP_0101$cl+1), nCores = 3 )
      namesave = paste0("03_Sensitivity_study/results_sensitivity/est_cl_thinnedDDP_0101", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      saveRDS(cl_point_est_thinnedDDP_0101, file = namesave)

      rand_thinnedDDP_0101 = 0
      for(gg in 1:length(unique(data$group))){
        rand_thinnedDDP_0101 = rand_thinnedDDP_0101 + adjustedRandIndex(cl_point_est_thinnedDDP_0101[data$group==gg], data$cl[data$group==gg])
      }
      rand_thinnedDDP_0101 = rand_thinnedDDP_0101/(length(unique(data$group)))

      namesave = paste0("03_Sensitivity_study/results_sensitivity/rand_thinnedDDP_0101", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      saveRDS(rand_thinnedDDP_0101, file = namesave)

      
      #-----# #-----# Compute density estimate #-----# #-----#
      seqq = seq(range(data$y)[1]-2, range(data$y)[2]+2, length.out = 300)
      density_est_thinnedDDP_0101 = thinnedDDP::compute_density(seqq,
                                                           weight = run_gibbs_thinnedDDP_0101$pi,
                                                           means = t(run_gibbs_thinnedDDP_0101$mu),
                                                           variances = t(run_gibbs_thinnedDDP_0101$sigma2)
      )
      namesave = paste0("03_Sensitivity_study/results_sensitivity/density_est_thinnedDDP_0101", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      saveRDS(density_est_thinnedDDP_0101, file = namesave)
      
      
      rm(run_gibbs_thinnedDDP_0101)
      rm(cl_point_est_thinnedDDP_0101)
      rm(density_est_thinnedDDP_0101)
      rm(rand_thinnedDDP_0101)

    }
  }
}
