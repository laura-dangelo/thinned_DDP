#----------------# #----------------# #----------------# 
#                   RUN GIBBS SAMPLER                  #
#----------------# #----------------# #----------------# 

# This script runs the Gibbs sampler to fit the various models to the synthetic data. 
# In particular, the fitted models are:
# - thinned DDP
# - a single DP mixture (complete pooling)
# - independent DP mixtures (a separate model for each group - no pooling)
# - CAM (from the SANple package)
# - GM-DDP (from the BNPmix package)

# The results are then saved into the "01_Simulation_study/results" folder.

library(sanba)
library(BNPmix)
library(thinnedDDP)
source("01_Simulation_study/auxiliary_functions/BGS.R")
source("01_Simulation_study/auxiliary_functions/postestimates.R")


# number of groups and sample sizes (to import and save data)
n_groups = c(2, 10)
ssg = c(10, 30)
n_ss = 4 # number of different configuration of sample size
n_datasets = 50
tot_datasets = n_datasets * n_ss * length(n_groups)

seq_thinning = seq(1, 5000, by = 2)

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
      #-----#          THINNED DDP        #-----#
      #-----#  #-----#  #-----#  #-----#  #-----#
      cat("Running thinned-DDP \n")
      tryCatch(
        {
          set.seed(12345)
          run_gibbs_thinnedDDP = thinnedDDP::sampler_thinnedDDP(nrep = nrep,
                                                                burnin = burnin,
                                                                y = data$y,
                                                                group = data$group-1,
                                                                trunc = trunc,
                                                                m0 = mean(data$y), tau0 = tau0,
                                                                gamma0 = gam0, lambda0 = lam0,
                                                                alpha = 1,
                                                                a_beta = 3, b_beta = 3,
                                                                mu_start = mu_start,
                                                                sigma2_start = sigma2_start,
                                                                cl_start = cl_start-1,
                                                                progressbar = T)

          namesave = paste0("01_Simulation_study/results/run_thinnedDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
          saveRDS(run_gibbs_thinnedDDP, file = namesave)

          time_thinnedDDP = run_gibbs_thinnedDDP$time
          namesave = paste0("01_Simulation_study/results/time_thinnedDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
          saveRDS(time_thinnedDDP, file = namesave)

          rm(run_gibbs_thinnedDDP)

        }, error = function(e) {
          nameerror = paste0("01_Simulation_study/results/ERROR_run_thinnedDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl, ".RDS")
          errorfile = list("i" = i, "j" = j, "repl" = repl)
          saveRDS(errorfile, file = nameerror)
          print(nameerror)
        } )


      #-----#  #-----#  #-----#  #-----#
      #-----#       POOLING      #-----#
      #-----#  #-----#  #-----#  #-----#
      cat("Running pooled DP mixture \n")
      tryCatch(
        {
          set.seed(12345)
          run_gibbs_pool = thinnedDDP::sampler_DP(nrep = nrep,
                                                  burnin = burnin,
                                                  y = data$y,
                                                  trunc = trunc,
                                                  m0 = mean(data$y), tau0 = tau0,
                                                  gamma0 = gam0, lambda0 = lam0,
                                                  alpha = 1,
                                                  mu_start = mu_start,
                                                  sigma2_start = sigma2_start,
                                                  cl_start = cl_start-1,
                                                  progressbar = T)
          namesave = paste0("01_Simulation_study/results/run_pool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
          saveRDS(run_gibbs_pool, file = namesave)

          time_pool = run_gibbs_pool$time
          namesave = paste0("01_Simulation_study/results/time_pool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
          saveRDS(time_pool, file = namesave)

          rm(run_gibbs_pool)

        }, error = function(e) {
          nameerror = paste0("01_Simulation_study/results/ERROR_run_pool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl, ".RDS")
          errorfile = list("i" = i, "j" = j, "repl" = repl)
          saveRDS(errorfile, file = nameerror)
          print(nameerror)
        } )




      #-----#  #-----#  #-----#  #-----#
      #-----#       NO POOL      #-----#
      #-----#  #-----#  #-----#  #-----#
      cat("Running no-pooling DP mixtures \n")
      tryCatch(
        {
          run_gibbs_nopool = list()
          time_nopool = 0
          for(gg in 1:length(unique(data$group))){
            set.seed(12345)
            run_gibbs_nopool[[gg]] = thinnedDDP::sampler_DP(nrep = nrep,
                                                            burnin = burnin,
                                                            y = data$y[data$group==gg],
                                                            trunc = trunc,
                                                            m0 = mean(data$y), tau0 = tau0,
                                                            gamma0 = gam0, lambda0 = lam0,
                                                            alpha = 1,
                                                            mu_start = mu_start,
                                                            sigma2_start = sigma2_start,
                                                            cl_start = cl_start[data$group==gg]-1,
                                                            progressbar = T)
            time_nopool = time_nopool + run_gibbs_nopool[[gg]]$time
          }
          namesave = paste0("01_Simulation_study/results/run_nopool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
          saveRDS(run_gibbs_nopool, file = namesave)

          namesave = paste0("01_Simulation_study/results/time_nopool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
          saveRDS(time_nopool, file = namesave)

          rm(run_gibbs_nopool)

        }, error = function(e) {
          nameerror = paste0("01_Simulation_study/results/ERROR_run_nopool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl, ".RDS")
          errorfile = list("i" = i, "j" = j, "repl" = repl)
          saveRDS(errorfile, file = nameerror)
          print(nameerror)
        } )





      #-----#  #-----#  #-----#  #-----#
      #-----#         CAM        #-----#
      #-----#  #-----#  #-----#  #-----#
      cat("Running CAM \n")
      tryCatch(
        {
          set.seed(12345)
          maxK = 20
          if(length(unique(data$group))==10) {maxK = 40}
          run_gibbs_CAM = fit_CAM(y = data$y, group = data$group, est_method = "MCMC",
                                  mcmc_param = list(nrep = nrep, burn = burnin,
                                                    warmstart = FALSE,
                                                    mu_start = mu_start,
                                                    maxK = maxK,
                                                    sigma2_start = sigma2_start,
                                                    M_start = cl_start-1,
                                                    S_start = 1:n_groups[i] ))

          run_gibbs_CAM$sim$mu = run_gibbs_CAM$sim$mu[seq_thinning,]
          run_gibbs_CAM$sim$sigma2 = run_gibbs_CAM$sim$sigma2[seq_thinning,]
          run_gibbs_CAM$sim$obs_cluster = run_gibbs_CAM$sim$obs_cluster[seq_thinning,]
          run_gibbs_CAM$sim$distr_cluster = run_gibbs_CAM$sim$distr_cluster[seq_thinning,]
          run_gibbs_CAM$sim$pi = run_gibbs_CAM$sim$pi[seq_thinning,]
          run_gibbs_CAM$sim$omega = run_gibbs_CAM$sim$omega[,,seq_thinning]
          run_gibbs_CAM$sim$alpha = run_gibbs_CAM$sim$alpha[seq_thinning]
          run_gibbs_CAM$sim$beta = run_gibbs_CAM$sim$beta[seq_thinning]
          run_gibbs_CAM$sim$maxK = run_gibbs_CAM$sim$maxK[seq_thinning]
          run_gibbs_CAM$sim$maxL = run_gibbs_CAM$sim$maxL[seq_thinning]

          namesave = paste0("01_Simulation_study/results/run_CAM", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
          saveRDS(run_gibbs_CAM, file = namesave)

          time_CAM = run_gibbs_CAM$time
          namesave = paste0("01_Simulation_study/results/time_CAM", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
          saveRDS(time_CAM, file = namesave)

          rm(run_gibbs_CAM)

        }, error = function(e) {
          nameerror = paste0("01_Simulation_study/results/ERROR_run_CAM", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl, ".RDS")
          errorfile = list("i" = i, "j" = j, "repl" = repl)
          saveRDS(errorfile, file = nameerror)
          print(nameerror)
        } )




      #-----#  #-----#  #-----#  #-----#
      #-----#       GM-DDP       #-----#
      #-----#  #-----#  #-----#  #-----#
      cat("Running GM-DDP \n")
      tryCatch(
        {
          seqq = seq(range(data$y)[1]-2, range(data$y)[2]+2, length.out = 300)
          set.seed(12345)
          run_gibbs_gmDDP = BNPmix::DDPdensity(y = data$y,
                                               group = data$group,
                                               mcmc = list(niter = nrep,
                                                           nburn = burnin),
                                               prior = list(m0 = mean(data$y),
                                                            k0 = tau0,
                                                            a0 = 2, # gam0,
                                                            b0 = 1, # lam0,
                                                            wei = 0.5),
                                               output = list(grid = seqq)
          )
          run_gibbs_gmDDP$density = run_gibbs_gmDDP$density[,,seq_thinning]
          run_gibbs_gmDDP$clust = run_gibbs_gmDDP$clust[seq_thinning,]
          run_gibbs_gmDDP$group_log = run_gibbs_gmDDP$group_log[seq_thinning,]
          run_gibbs_gmDDP$wvals = run_gibbs_gmDDP$wvals[seq_thinning,]

          namesave = paste0("01_Simulation_study/results/run_gmDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
          saveRDS(run_gibbs_gmDDP, file = namesave)

          time_gmDDP = run_gibbs_gmDDP$tot_time
          namesave = paste0("01_Simulation_study/results/time_gmDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
          saveRDS(time_gmDDP, file = namesave)

          rm(run_gibbs_gmDDP)

        }, error = function(e) {
          nameerror = paste0("01_Simulation_study/results/ERROR_run_gmDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl, ".RDS")
          errorfile = list("i" = i, "j" = j, "repl" = repl)
          saveRDS(errorfile, file = nameerror)
          print(nameerror)
        } )





      #-----#  #-----#  #-----#  #-----#
      #-----#         HDP        #-----#
      #-----#  #-----#  #-----#  #-----#

      yHDP = list()
      for(g in unique(data$group)){
        yHDP[[g]] = data$y[data$group==g]
      }
      cat("Running HDP \n")

      tryCatch(
        {
          seqq = seq(range(data$y)[1]-2, range(data$y)[2]+2, length.out = 300)
          set.seed(12345)
          tmp_HDP = blocked_gibbs(x = yHDP, L.max = trunc, gam = 1, phi.param = c(0, 1, 1), b0 = 0.1,
                                        N = 1, Burn.in = burnin, M = 5000, est.density = TRUE, y.grid = seqq)


          run_gibbs_HDP = list()
          for(idlist in 1:length(seq_thinning)) {
            run_gibbs_HDP[[idlist]] = tmp_HDP[[seq_thinning[idlist]]]
          }

          namesave = paste0("01_Simulation_study/results/run_HDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
          saveRDS(run_gibbs_HDP, file = namesave)

          time_HDP = run_gibbs_HDP$tot_time
          namesave = paste0("01_Simulation_study/results/time_HDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
          saveRDS(time_HDP, file = namesave)

          rm(run_gibbs_HDP)

        }, error = function(e) {
          nameerror = paste0("01_Simulation_study/results/ERROR_run_HDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl, ".RDS")
          errorfile = list("i" = i, "j" = j, "repl" = repl)
          saveRDS(errorfile, file = nameerror)
          print(nameerror)
        } )


    }
  }
}
