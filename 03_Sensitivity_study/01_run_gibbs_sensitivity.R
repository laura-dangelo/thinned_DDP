#----------------# #----------------# #----------------# 
#                   SENSITIVITY STUDY                  #
#----------------# #----------------# #----------------# 

# This script runs the Gibbs sampler to fit the thinned DDP with different combinations of the 
# hyperparameters on the Beta prior on pi_g.

# The results are then saved into the "03_Sensitivity_study/results_sensitivity" folder.

library(thinnedDDP)

# number of groups and sample sizes (to import and save data)
n_groups = c(2, 10)
ssg = c(10, 30)
n_ss = 4 
n_datasets = 25
tot_datasets = n_datasets * n_ss * length(n_groups)


# start loop for fitting the models
for(repl in 1:n_datasets) {
  for(i in 1:length(n_groups) ){
    for(j in 1:n_ss){
      
      name_file_open = paste0("01_Simulation_study/data/data_", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      data = readRDS(file = name_file_open)
      
      # set number of iterations
      nrep = 10000
      burnin = (floor(nrep/2))
      
      # define hyperparameters
      tau0 = 0.01
      gam0 = 2.5
      lam0 = (gam0-1)
      
      trunc = 50
      
      # define reasonable starting points of the Markov chains
      set.seed(2)
      cl_start = kmeans(data$y, centers = 4)
      cl_start = cl_start$cluster
      
      mu_start = rep(0, trunc)
      mu_start[1:length(unique(cl_start))] = tapply(data$y, cl_start, mean)
      sigma2_start = rep(var(data$y)/3, trunc)
      sigma2_start[1:length(unique(cl_start))] = tapply(data$y, cl_start, var)
      sigma2_start[1:length(unique(cl_start))] = sigma2_start[1:length(unique(cl_start))]
      sigma2_start[is.na(sigma2_start)] = var(data$y)/3
      
      
      cat(paste0("Replication ", repl, "\n"))
      cat(paste0("i = ", i, "/", length(n_groups), "\n"))
      cat(paste0("j = ", j, "/", n_ss, "\n"))
      
      
      #-----#  #-----#  #-----#  #-----#  #-----#
      #-----#        UNIFORM PRIOR        #-----#
      #-----#  #-----#  #-----#  #-----#  #-----#
      cat("Running thinned-DDP \n")
      tryCatch(
        {
          set.seed(12345)
          run_gibbs_thinnedDDP_unif = thinnedDDP::sampler_thinnedDDP(nrep = nrep,
                                                                burnin = burnin,
                                                                y = data$y,
                                                                group = data$group-1,
                                                                trunc = trunc,
                                                                m0 = mean(data$y), tau0 = tau0,
                                                                gamma0 = gam0, lambda0 = lam0,
                                                                alpha = 1,
                                                                a_beta = 1, b_beta = 1,
                                                                mu_start = mu_start,
                                                                sigma2_start = sigma2_start,
                                                                cl_start = cl_start-1,
                                                                progressbar = T)
          
          namesave = paste0("03_Sensitivity_study/results_sensitivity/run_thinnedDDP_unif", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
          saveRDS(run_gibbs_thinnedDDP_unif, file = namesave)
          

          rm(run_gibbs_thinnedDDP_unif)
          
        }, error = function(e) {
          nameerror = paste0("03_Sensitivity_study/results_sensitivity/ERROR_run_thinnedDDP_unif", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl, ".RDS")
          errorfile = list("i" = i, "j" = j, "repl" = repl)
          saveRDS(errorfile, file = nameerror)
          print(nameerror)
        } )
      
      
  
      #-----#  #-----#  #-----#  #-----#  #-----#
      #-----#         BETA(10,10)         #-----#
      #-----#  #-----#  #-----#  #-----#  #-----#
      cat("Running thinned-DDP \n")
      tryCatch(
        {
          set.seed(12345)
          run_gibbs_thinnedDDP_1010 = thinnedDDP::sampler_thinnedDDP(nrep = nrep,
                                                                     burnin = burnin,
                                                                     y = data$y,
                                                                     group = data$group-1,
                                                                     trunc = trunc,
                                                                     m0 = mean(data$y), tau0 = tau0,
                                                                     gamma0 = gam0, lambda0 = lam0,
                                                                     alpha = 1,
                                                                     a_beta = 10, b_beta = 10,
                                                                     mu_start = mu_start,
                                                                     sigma2_start = sigma2_start,
                                                                     cl_start = cl_start-1,
                                                                     progressbar = T)
          
          namesave = paste0("03_Sensitivity_study/results_sensitivity/run_thinnedDDP_1010", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
          saveRDS(run_gibbs_thinnedDDP_1010, file = namesave)
          
          
          rm(run_gibbs_thinnedDDP_1010)
          
        }, error = function(e) {
          nameerror = paste0("03_Sensitivity_study/results_sensitivity/ERROR_run_thinnedDDP_1010", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl, ".RDS")
          errorfile = list("i" = i, "j" = j, "repl" = repl)
          saveRDS(errorfile, file = nameerror)
          print(nameerror)
        } )
      
      
      
      #-----#  #-----#  #-----#  #-----#  #-----#
      #-----#         BETA(0.1, 0.1)      #-----#
      #-----#  #-----#  #-----#  #-----#  #-----#
      cat("Running thinned-DDP \n")
      tryCatch(
        {
          set.seed(12345)
          run_gibbs_thinnedDDP_0101 = thinnedDDP::sampler_thinnedDDP(nrep = nrep,
                                                                     burnin = burnin,
                                                                     y = data$y,
                                                                     group = data$group-1,
                                                                     trunc = trunc,
                                                                     m0 = mean(data$y), tau0 = tau0,
                                                                     gamma0 = gam0, lambda0 = lam0,
                                                                     alpha = 1,
                                                                     a_beta = 0.1, b_beta = 0.1,
                                                                     mu_start = mu_start,
                                                                     sigma2_start = sigma2_start,
                                                                     cl_start = cl_start-1,
                                                                     progressbar = T)
          
          namesave = paste0("03_Sensitivity_study/results_sensitivity/run_thinnedDDP_0101", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
          saveRDS(run_gibbs_thinnedDDP_0101, file = namesave)
          
          
          rm(run_gibbs_thinnedDDP_0101)
          
        }, error = function(e) {
          nameerror = paste0("03_Sensitivity_study/results_sensitivity/ERROR_run_thinnedDDP_0101", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl, ".RDS")
          errorfile = list("i" = i, "j" = j, "repl" = repl)
          saveRDS(errorfile, file = nameerror)
          print(nameerror)
        } )
      
    }
  }
}
