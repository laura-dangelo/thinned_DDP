
#----------------# #----------------# #----------------# 
#                 EFFECTIVE SAMPLE SIZE                #
#----------------# #----------------# #----------------# 

library(devtools)
library(ggplot2)
library(viridis)
library(coda)

n_groups = c(2, 10)
ssg = c(10, 30)
n_ss = 4 # number of different configuration of sample size
n_datasets = 50

tot_datasets = n_datasets * n_ss * length(n_groups)

ess_df = data.frame("Parameter" = character(), "G" = numeric(), "n" = numeric(), "ess" = double())


for(repl in 1:n_datasets) {
  start = Sys.time()
  for(i in 1:length(n_groups) ){
    for(j in 1:n_ss){
      
      G = n_groups[i]
      nameopen = paste0("01_Simulation_study/results/run_thinnedDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      run_thinnedDDP = readRDS(nameopen)
      
      tmpell = sapply(1:G, function(g) effectiveSize(t(run_thinnedDDP$ell[,g,])))
      newdata = c("Thinning variables", n_groups[i], sum(n_groups[i]/2*ssg*j), mean(tmpell))
      ess_df[nrow(ess_df)+1,] = newdata
  
      
      tmppi = sapply(1:G, function(g) effectiveSize(t(run_thinnedDDP$pi[,g,])))
      newdata = c("Probabilities", n_groups[i], sum(n_groups[i]/2*ssg*j), mean(tmppi))
      ess_df[nrow(ess_df)+1,] = newdata
      
      
      tmpcl = effectiveSize(run_thinnedDDP$cl)
      newdata = c("Cluster allocation", n_groups[i], sum(n_groups[i]/2*ssg*j), mean(tmpcl))
      ess_df[nrow(ess_df)+1,] = newdata
      
      
      tmpthp = effectiveSize(t(run_thinnedDDP$thinning_prob))
      newdata = c("Thinning probabilities", n_groups[i], sum(n_groups[i]/2*ssg*j), mean(tmpthp))
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

ggplot(ess_df, aes(x = n, y = ess, fill=Parameter ) ) +
  geom_boxplot(alpha = 0.6) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill='transparent', color=NA), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    # panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    panel.border = element_rect(color = "darkgray ", fill=NA),
    # axis.line.y.left = element_line(color="gray"),
    axis.line.x.bottom = element_line(color="gray"),
    #
    legend.position = "bottom",
    legend.background = element_rect(fill='transparent', color = 'white'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent'), #transparent legend panel
    legend.text = element_text(size=10),
    strip.text = element_text(size=10),
    strip.background = element_rect( fill=NA, color="gray" )
  )+
  # scale_fill_manual( values = c(rocket(8)[4], mako(8)[5], inferno(8)[7]) ) +
  scale_fill_manual( values = c("darkgoldenrod1", "cyan4",  "deeppink4", "turquoise") ) +
  xlab("Sample size")  +
  ylab("Effective sample size of the MCMC chain")+
  facet_wrap( ~ G, scales = "free",
              labeller = labeller(G = c("2" = "2 groups",
                                        "10" = "10 groups") )
  )

ggsave("01_Simulation_study/output_images/09_ess.pdf", width = 8, height = 4)




