
#----------------# #----------------# #----------------# 
#                 EFFECTIVE SAMPLE SIZE                #
#----------------# #----------------# #----------------# 

library(devtools)
library(ggplot2)
library(viridis)
library(coda)

n_groups = c(2, 10)
ssg = c(10, 30)
n_ss = 4 
n_datasets = 50

tot_datasets = n_datasets * n_ss * length(n_groups)


deviance_i = function(y_i, g_i, group, cl, mu, sigma2){
  cluster_sizes = table(cl[group==g_i])
  atoms = as.numeric(attr(cluster_sizes, "dimnames")[[1]])
  comp = cluster_sizes/sum(group==g_i) * dnorm(y_i, mu[atoms], sqrt(sigma2[atoms]))
  out = log(sum(comp))
  out
}

deviance_repl = function(y, group, cl, mu, sigma2){
  tmp = cbind(y, group)
  out = apply(tmp, 1, function(x) deviance_i(x[1], x[2], group, cl, mu, sigma2))
  return(-2*sum(out))
}

if(!file.exists("01_Simulation_study/output_RDS/effective_sample_size.RDS")){
  ess_df = data.frame("Parameter" = character(), "G" = numeric(), "n" = numeric(), "ess" = double())
  
  for(repl in 1:n_datasets) {
    start = Sys.time()
    for(i in 1:length(n_groups) ){
      for(j in 1:n_ss){
        
        G = n_groups[i]
        nameopen = paste0("01_Simulation_study/results/run_thinnedDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
        run_thinnedDDP = readRDS(nameopen)
        
        y = run_thinnedDDP$y
        group = run_thinnedDDP$group + 1
        chain_nclust = apply(run_thinnedDDP$cl, 1, function(x) length(unique(x)))
        chain_cl = run_thinnedDDP$cl +1
        chain_means = run_thinnedDDP$mu
        chain_vars = run_thinnedDDP$sigma2
        rm(run_thinnedDDP)
        
        out_dev = sapply(1:length(chain_nclust), function(r) deviance_repl(y, group, chain_cl[r,], chain_means[r,], chain_vars[r,]) )
        
        newdata = c("Number of clusters", n_groups[i], sum(n_groups[i]/2*ssg*j), effectiveSize(chain_nclust))
        ess_df[nrow(ess_df)+1,] = newdata
        
        newdata = c("Deviance", n_groups[i], sum(n_groups[i]/2*ssg*j), effectiveSize(out_dev))
        ess_df[nrow(ess_df)+1,] = newdata
      }
    }
    end = Sys.time()
    print(paste0("Finished replication ", repl, " in ", (end-start), " sec"))
  }
  str(ess_df)
  
  ess_df$n = as.numeric(ess_df$n)
  ess_df$G = as.numeric(ess_df$G)
  ess_df$ess = as.numeric(ess_df$ess)
  
  ess_df = ess_df[order(ess_df$Parameter, ess_df$G, ess_df$n), ]
  
  ess_df$n = as.factor(ess_df$n)
  ess_df$G = as.factor(ess_df$G)
  ess_df$Parameter = as.factor(ess_df$Parameter)
  
  saveRDS(ess_df, file = "01_Simulation_study/output_RDS/effective_sample_size.RDS")
} else {
  ess_df = readRDS(file = "01_Simulation_study/output_RDS/effective_sample_size.RDS")
}

ggplot(ess_df, aes(x = n, y = ess, fill=Parameter ) ) +
  geom_boxplot(alpha = 0.6) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill='transparent', color=NA), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.minor = element_blank(), #remove minor gridlines
    panel.border = element_rect(color = "darkgray ", fill=NA),
    axis.line.x.bottom = element_line(color="gray"),
    legend.position = "bottom",
    legend.background = element_rect(fill='transparent', color = 'white'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent'), #transparent legend panel
    legend.text = element_text(size=10),
    strip.text = element_text(size=10),
    strip.background = element_rect( fill=NA, color="gray" )
  )+
  scale_fill_manual( values = c("darkgoldenrod1", "cyan4",  "deeppink4", "turquoise") ) +
  xlab("Sample size")  +
  ylab("Effective sample size")+
  facet_wrap( ~ G, scales = "free",
              labeller = labeller(G = c("2" = "2 groups",
                                        "10" = "10 groups") )
  )

ggsave("01_Simulation_study/output_images/09_ess.pdf", width = 8, height = 4)




