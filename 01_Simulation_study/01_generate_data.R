#----------------# #----------------# #----------------# 
#                GENERATE SYNTHETIC DATA               #
#----------------# #----------------# #----------------# 

# This script generates the synthetic data used in the simulation study.
# The data are organized into G groups, with G = 2, 10.
# They are generated from group-specific mixtures of Gaussian kernels,
# where some of the groups may share the same data-generating distribution.

# The simulated data are then saved in the "01_Simulation_study/data" folder.


source("01_Simulation_study/auxiliary_functions/00_sim_mixN_thinned.R")

n_groups = c(2, 10)
ssg = c(10, 30)
n_ss = 4 # number of different configuration of sample size
n_datasets = 50 # number of replicates




sample_sizes_groups = list()
for(i in 1:length(n_groups))
{
  tmp = matrix(0, n_ss, n_groups[i])
  for(j in 1:n_ss){
    tmp[j,] = sort(rep(ssg*j, n_groups[i]/2))
  }
  sample_sizes_groups[[i]] = tmp
}
sample_sizes_groups
lapply(sample_sizes_groups, rowSums)



mu_k = c(-5, 0, 5, 10) # means of the mixture components
var_k = rep(0.6, length(mu_k)) # variance of the mixture components

weights = matrix(0,10,length(mu_k)) # mixture probabilities
for(i in 1:5){
  weights[2*i-1,] = c(0.5, 0.25, 0.25, 0)
  weights[2*i,] = c(0, 0, 0.4, 0.6)
}


## loop that generates the data for the n_datasets replicates
for(repl in 1:n_datasets) {
  
  for(i in 1:length(n_groups) ){
    for(j in 1:n_ss){
      set.seed(repl*2)
      data = sim_thinned_mixtureN(sample_sizes = c(sample_sizes_groups[[i]][j,]),
                                  weights = t(weights[1:n_groups[i],]),
                                  mu_k = mu_k, var_k = var_k,
                                  seed = repl*2 )
      totalss = length(data$y)
      namefile = paste0("01_Simulation_study/data/data_", n_groups[i], "groups_", totalss, "n_", repl,".RDS")
      saveRDS(data, file = namefile)
    }
  }
}


