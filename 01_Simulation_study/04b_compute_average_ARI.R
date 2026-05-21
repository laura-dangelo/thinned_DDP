library(devtools)
library(salso)
library(mclust)
library(mcclust.ext)

# library(SANple)
# library(BNPmix)

n_groups = c(2, 10)
ssg = c(10, 30)
n_ss = 4 # number of different configuration of sample size
n_datasets = 50

tot_datasets = n_datasets * n_ss * length(n_groups)
# document()


repl=2
i=1
j=1

for(repl in 1: 50 ){#n_datasets) {
  for(i in 1:length(n_groups) ){
    for(j in 1:n_ss){
      
      name_file_open = paste0("../04_Simulation_commonatoms/Data/data_", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
      load(file = name_file_open)
      
      

      #-----#  #-----#  #-----#  #-----#  #-----#
      #-----#     THINNED DDP random pi   #-----#
      #-----#  #-----#  #-----#  #-----#  #-----#
          namesave = paste0("/media/laura/TOSHIBA EXT/Postdoc_Bicocca (updated 11-24)/1.9_Thinned_DDP (updated 17-12-24)/04_Bernoulli_thinning/04_Simulation_commonatoms/Results/run_thinnedDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
          load(namesave)
          
         
          #-----# Estimate partition of observations #-----#
          group_cl_point_est_thinnedDDP = numeric(length(data$y))
          G = length(unique(data$group))
          for(gg in 1:G) {
            est_cl_group = salso::salso((run_gibbs_thinnedDDP$cl[,data$group==gg]+1), nCores = 3 )
            group_cl_point_est_thinnedDDP[data$group==gg] = est_cl_group
          }
            
          namesave = paste0("../04_Simulation_commonatoms/Results/group_spec_est_cl_thinnedDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
          save(group_cl_point_est_thinnedDDP, file = namesave)
          
          
          group_spec_rand_thinnedDDP = numeric(G)
          for(gg in 1:length(unique(data$group))){
            group_spec_rand_thinnedDDP[gg] = adjustedRandIndex(group_cl_point_est_thinnedDDP[data$group==gg], data$cl[data$group==gg])
          }
          
          average_rand_thinnedDDP = mean(group_spec_rand_thinnedDDP)
          namesave = paste0("../04_Simulation_commonatoms/Results/average_rand_thinnedDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
          save(average_rand_thinnedDDP, file = namesave)
          
          
          rm(run_gibbs_thinnedDDP)
          rm(group_cl_point_est_thinnedDDP)
          rm(group_spec_rand_thinnedDDP)
          rm(average_rand_thinnedDDP)
      
      
      #-----#  #-----#  #-----#  #-----#
      #-----#       POOLING      #-----#
      #-----#  #-----#  #-----#  #-----#
      
          namesave = paste0("/media/laura/TOSHIBA EXT/Postdoc_Bicocca (updated 11-24)/1.9_Thinned_DDP (updated 17-12-24)/04_Bernoulli_thinning/04_Simulation_commonatoms/Results/run_pool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
          load(namesave)
          
          
          #-----# Estimate partition of observations #-----#
          group_cl_point_est_pool = numeric(length(data$y))
          G = length(unique(data$group))
          for(gg in 1:G) {
            est_cl_group = salso::salso((run_gibbs_pool$cl[,data$group==gg]+1), nCores = 3 )
            group_cl_point_est_pool[data$group==gg] = est_cl_group
          }
          
          namesave = paste0("../04_Simulation_commonatoms/Results/group_spec_est_cl_pool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
          save(group_cl_point_est_pool, file = namesave)
          
          
          group_spec_rand_pool = numeric(G)
          for(gg in 1:length(unique(data$group))){
            group_spec_rand_pool[gg] = adjustedRandIndex(group_cl_point_est_pool[data$group==gg], data$cl[data$group==gg])
          }
          
          average_rand_pool = mean(group_spec_rand_pool)
          namesave = paste0("../04_Simulation_commonatoms/Results/average_rand_pool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
          save(average_rand_pool, file = namesave)
          
          
          rm(run_gibbs_pool)
          rm(group_cl_point_est_pool)
          rm(group_spec_rand_pool)
          rm(average_rand_pool)
      
      # 
      # 
      # #-----#  #-----#  #-----#  #-----#
      # #-----#       NO POOL      #-----#
      # #-----#  #-----#  #-----#  #-----#
      # tryCatch(
      #   {
      #     run_gibbs_nopool = list()
      #     time_nopool = 0
      #     for(gg in 1:length(unique(data$group))){
      #       set.seed(123)
      #       run_gibbs_nopool[[gg]] = Bthin::blocked_gibbs_DP_burnin(nrep = nrep,
      #                                                               burnin = burnin,
      #                                                               y = data$y[data$group==gg],
      #                                                               trunc = trunc,
      #                                                               m0 = mean(data$y), tau0 = tau0,
      #                                                               gamma0 = gam0, lambda0 = lam0,
      #                                                               alpha = 1,
      #                                                               mu_start = mu_start,
      #                                                               sigma2_start = sigma2_start,
      #                                                               cl_start = cl_start[data$group==gg]-1,
      #                                                               progressbar = T)
      #       time_nopool = time_nopool + run_gibbs_nopool[[gg]]$time
      #     }
      #     namesave = paste0("../04_Simulation_commonatoms/Results/run_nopool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
      #     save(run_gibbs_nopool, file = namesave)
      #     namesave = paste0("../04_Simulation_commonatoms/Results/time_nopool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
      #     save(time_nopool, file = namesave)
      #     
      #     #-----# Estimate partition of observations #-----#
      #     cl_point_est_nopool = list()
      #     for(gg in 1:length(unique(data$group))){
      #       clus_mat = salso::psm((run_gibbs_nopool[[gg]]$cl+1), nCores = 3 )
      #       minv = minVI(clus_mat)
      #       cl_point_est_nopool[[gg]] = minv$cl
      #       # length(unique(cl_point_est_nopool[[gg]])) ### estimate
      #     }
      #     namesave = paste0("../04_Simulation_commonatoms/Results/est_cl_nopool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
      #     save(cl_point_est_nopool, file = namesave)
      #     
      #     rand_nopool = 0
      #     for(gg in 1:length(unique(data$group))){
      #       rand_nopool = rand_nopool + adjustedRandIndex(cl_point_est_nopool[[gg]], data$cl[data$group==gg])
      #     }
      #     rand_nopool = rand_nopool/(length(unique(data$group)))
      #     namesave = paste0("../04_Simulation_commonatoms/Results/rand_nopool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
      #     save(rand_nopool, file = namesave)
      #     
      #     #-----# #-----# Compute density estimate #-----# #-----#
      #     density_est_nopool = list()
      #     seqq = seq(range(data$y)[1]-2, range(data$y)[2]+2, length.out = 300)
      #     for(gg in 1:length(unique(data$group))){
      #       density_est_nopool[[gg]] = Bthin::compute_density_1DP(seqq,
      #                                                             weight = t(run_gibbs_nopool[[gg]]$pi),
      #                                                             means = t(run_gibbs_nopool[[gg]]$mu),
      #                                                             variances = t(run_gibbs_nopool[[gg]]$sigma2) )
      #     }
      #     namesave = paste0("../04_Simulation_commonatoms/Results/density_est_nopool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
      #     save(density_est_nopool, file = namesave)
      #     
      #     
      #     rm(run_gibbs_nopool)
      #     rm(cl_point_est_nopool)
      #     rm(density_est_nopool)
      #     rm(rand_nopool)
      #   }, error = function(e) {
      #     nameerror = paste0("../04_Simulation_commonatoms/Results/ERROR_run_nopool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl, ".Rdata")
      #     errorfile = list("i" = i, "j" = j, "repl" = repl)
      #     save(errorfile, file = nameerror)
      #     print(nameerror)
      #   } )
      # 
      # 
      
      
      
      #-----#  #-----#  #-----#  #-----#
      #-----#         CAM        #-----#
      #-----#  #-----#  #-----#  #-----#
          
          namesave = paste0("/media/laura/TOSHIBA EXT/Postdoc_Bicocca (updated 11-24)/1.9_Thinned_DDP (updated 17-12-24)/04_Bernoulli_thinning/04_Simulation_commonatoms/Results/run_CAM", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
          load(namesave)
          
          
          #-----# Estimate partition of observations #-----#
          group_cl_point_est_CAM = numeric(length(data$y))
          G = length(unique(data$group))
          for(gg in 1:G) {
            est_cl_group = salso::salso((run_gibbs_CAM$sim$obs_cluster[,data$group==gg]), nCores = 3 )
            group_cl_point_est_CAM[data$group==gg] = est_cl_group
          }
          
          namesave = paste0("../04_Simulation_commonatoms/Results/group_spec_est_cl_CAM", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
          save(group_cl_point_est_CAM, file = namesave)
          
          
          group_spec_rand_CAM = numeric(G)
          for(gg in 1:length(unique(data$group))){
            group_spec_rand_CAM[gg] = adjustedRandIndex(group_cl_point_est_CAM[data$group==gg], data$cl[data$group==gg])
          }
          
          average_rand_CAM = mean(group_spec_rand_CAM)
          namesave = paste0("../04_Simulation_commonatoms/Results/average_rand_CAM", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
          save(average_rand_CAM, file = namesave)
          
          
          rm(run_gibbs_CAM)
          rm(group_cl_point_est_CAM)
          rm(group_spec_rand_CAM)
          rm(average_rand_CAM)
          
          
          
      
      
      #-----#  #-----#  #-----#  #-----#
      #-----#       GM-DDP       #-----#
      #-----#  #-----#  #-----#  #-----#
          
          
          namesave = paste0("/media/laura/TOSHIBA EXT/Postdoc_Bicocca (updated 11-24)/1.9_Thinned_DDP (updated 17-12-24)/04_Bernoulli_thinning/04_Simulation_commonatoms/Results/run_gmDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
          load(namesave)
          
          
          #-----# Estimate partition of observations #-----#
          group_cl_point_est_gmDDP = numeric(length(data$y))
          G = length(unique(data$group))
          for(gg in 1:G) {
            est_cl_group = salso::salso((run_gibbs_gmDDP$clust[,data$group==gg]), nCores = 3 )
            group_cl_point_est_gmDDP[data$group==gg] = est_cl_group
          }
          
          namesave = paste0("../04_Simulation_commonatoms/Results/group_spec_est_cl_gmDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
          save(group_cl_point_est_gmDDP, file = namesave)
          
          
          group_spec_rand_gmDDP = numeric(G)
          for(gg in 1:length(unique(data$group))){
            group_spec_rand_gmDDP[gg] = adjustedRandIndex(group_cl_point_est_gmDDP[data$group==gg], data$cl[data$group==gg])
          }
          
          average_rand_gmDDP = mean(group_spec_rand_gmDDP)
          namesave = paste0("../04_Simulation_commonatoms/Results/average_rand_gmDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
          save(average_rand_gmDDP, file = namesave)
          
          
          rm(run_gibbs_gmDDP)
          rm(group_cl_point_est_gmDDP)
          rm(group_spec_rand_gmDDP)
          rm(average_rand_gmDDP)
          
          
      
      
      
      
      
      
      
    }
  }
}
