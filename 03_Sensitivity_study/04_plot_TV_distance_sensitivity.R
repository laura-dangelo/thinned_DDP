
#----------------# #----------------# #----------------# 
#             PLOT TOTAL VARIATION DISTANCE            #
#----------------# #----------------# #----------------# 

library(devtools)
library(salso)
library(mclust)
library(mcclust.ext)
library(ggplot2)
library(viridisLite)


n_groups = c(2, 10)
ssg = c(10, 30)
n_ss = 4 # number of different configuration of sample size
n_datasets = 25

tot_datasets = n_datasets * n_ss * length(n_groups)



true_density = function(data, seq)
{
  G = length(unique(data$group))
  out = matrix(0, G, length(seq))

  for(g in 1:G) {
    for(i in 1:length(seq)){
      out[g,i] = sum(data$weights[,g] * dnorm(seq[i], mean = data$mu_k, sd = sqrt(data$var_k)) )
    }
  }
  return(out)
}

TV_distance = function(seq, dens1, dens2) {
  sum(.5 * abs(dens1-dens2))
}


#-----------#  COMPARE TV DISTANCE WITH COMPLETE AND NO-POOLING MODELS  #-----------#

TVdist_df2 = data.frame("Model" = character(), "G" = numeric(), "n" = numeric(), "TV" = double())

for(repl in 1:n_datasets) {
  for(i in 1:length(n_groups) ){
    for(j in 1:n_ss){

      name_file_open = paste0("01_Simulation_study/data/data_", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      data = readRDS(file = name_file_open)



      #-----#  #-----#  #-----#  #-----#  #-----#
      #-----#         BETA(3,3)         #-----#
      #-----#  #-----#  #-----#  #-----#  #-----#
      nameopen = paste0("01_Simulation_study/results/density_est_thinnedDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      density_est_thinnedDDP = readRDS(nameopen)

      seqq = density_est_thinnedDDP$seq
      truth = true_density(data, seqq)

      namesave = paste0("01_Simulation_study/results/true_density_", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      saveRDS(truth, file = namesave)

      TV_thinnedDDP = 0
      for(gg in 1:length(unique(data$group))){
        TV_thinnedDDP = TV_thinnedDDP + TV_distance(seqq, truth[gg,], density_est_thinnedDDP$mean[,gg])
      }
      TV_thinnedDDP = TV_thinnedDDP/(length(unique(data$group)))

      newdata = c("Beta(3, 3)", n_groups[i], sum(n_groups[i]/2*ssg*j), TV_thinnedDDP)
      TVdist_df2[nrow(TVdist_df2)+1,] = newdata

      
      #-----#  #-----#  #-----#  #-----#  #-----#
      #-----#        UNIFORM PRIOR        #-----#
      #-----#  #-----#  #-----#  #-----#  #-----#
      
      namesave = paste0("03_Sensitivity_study/results_sensitivity/density_est_thinnedDDP_unif", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      density_est_unif = readRDS(namesave)

      TV_pool = 0
      for(gg in 1:length(unique(data$group))){
        TV_pool = TV_pool + TV_distance(seqq, truth[gg,], density_est_unif$mean[,gg])
      }
      TV_pool = TV_pool/(length(unique(data$group)))

      newdata = c("Uniform", n_groups[i], sum(n_groups[i]/2*ssg*j), TV_pool)
      TVdist_df2[nrow(TVdist_df2)+1,] = newdata



      #-----#  #-----#  #-----#  #-----#  #-----#
      #-----#         BETA(10,10)         #-----#
      #-----#  #-----#  #-----#  #-----#  #-----#
      
      namesave = paste0("03_Sensitivity_study/results_sensitivity/density_est_thinnedDDP_1010", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      density_est_thinnedDDP_1010 = readRDS(namesave)

      TV_thinnedDDP_1010 = 0
      for(gg in 1:length(unique(data$group))){
        TV_thinnedDDP_1010 = TV_thinnedDDP_1010 + TV_distance(seqq, truth[gg,], density_est_thinnedDDP_1010$mean[,gg])
      }
      TV_thinnedDDP_1010 = TV_thinnedDDP_1010/(length(unique(data$group)))

      newdata = c("Beta(10, 10)", n_groups[i], sum(n_groups[i]/2*ssg*j), TV_thinnedDDP_1010)
      TVdist_df2[nrow(TVdist_df2)+1,] = newdata
      
      
      
      
      #-----#  #-----#  #-----#  #-----#  #-----#
      #-----#         BETA(0.1,0.1)         #-----#
      #-----#  #-----#  #-----#  #-----#  #-----#
      
      namesave = paste0("03_Sensitivity_study/results_sensitivity/density_est_thinnedDDP_0101", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      density_est_thinnedDDP_0101 = readRDS(namesave)
      
      TV_thinnedDDP_0101 = 0
      for(gg in 1:length(unique(data$group))){
        TV_thinnedDDP_0101 = TV_thinnedDDP_0101 + TV_distance(seqq, truth[gg,], density_est_thinnedDDP_0101$mean[,gg])
      }
      TV_thinnedDDP_0101 = TV_thinnedDDP_0101/(length(unique(data$group)))
      
      newdata = c("Beta(0.1, 0.1)", n_groups[i], sum(n_groups[i]/2*ssg*j), TV_thinnedDDP_0101)
      TVdist_df2[nrow(TVdist_df2)+1,] = newdata

    }
  }
}
str(TVdist_df2)

TVdist_df2$n = as.numeric(TVdist_df2$n)
TVdist_df2$G = as.numeric(TVdist_df2$G)
TVdist_df2$TV = as.numeric(TVdist_df2$TV)

TVdist_df2 = TVdist_df2[order(TVdist_df2$Model, TVdist_df2$G, TVdist_df2$n), ]

TVdist_df2$n = as.factor(TVdist_df2$n)
TVdist_df2$G = as.factor(TVdist_df2$G)
TVdist_df2$Model = as.factor(TVdist_df2$Model)


saveRDS(TVdist_df2, file = "03_Sensitivity_study/output_RDS/TVdist_df_sensitivity.RDS")

TVdist_df2 = readRDS("03_Sensitivity_study/output_RDS/TVdist_df_sensitivity.RDS")



ggplot(TVdist_df2, aes(x = n, y = TV, fill=Model ) ) +
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
  scale_fill_manual( values = c("forestgreen", "royalblue3", "deeppink4","orange" ) ) +
  xlab("Sample size")  +
  ylab("TV distance")+
  # ylim(0.7,5.9)+
  facet_wrap( ~ G, scales = "free",
              labeller = labeller(G = c("2" = "2 groups",
                                        "10" = "10 groups") )
  )


ggsave("03_Sensitivity_study/output_images/05_01_TVdistance_sensitivity.pdf", width = 8, height = 4)



