
#----------------# #----------------# #----------------# 
#             PLOT TOTAL VARIATION DISTANCE            #
#----------------# #----------------# #----------------# 

# This script produces Figure 4 and Figure 6 in the article.
# The script loads the dataframes TVdist_df_pool.RDS and TVdist_df_models.RDS from the output_RDS folder.


library(devtools)
library(salso)
library(mclust)
library(mcclust.ext)
library(ggplot2)
library(viridisLite)


n_groups = c(2, 10)
ssg = c(10, 30)
n_ss = 4 
n_datasets = 50
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

if(!file.exists("01_Simulation_study/output_RDS/TVdist_df_pool.RDS")){
  TVdist_df2 = data.frame("Model" = character(), "G" = numeric(), "n" = numeric(), "TV" = double())
  
  for(repl in 1:n_datasets) {
    for(i in 1:length(n_groups) ){
      for(j in 1:n_ss){
        
        name_file_open = paste0("01_Simulation_study/data/data_", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
        data = readRDS(file = name_file_open)
        
        #-----#          THINNED DDP        #-----#
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
        newdata = c("Thinned DDP", n_groups[i], sum(n_groups[i]/2*ssg*j), TV_thinnedDDP)
        TVdist_df2[nrow(TVdist_df2)+1,] = newdata
        
        
        #-----#       COMPLETE POOLING      #-----#
        nameopen = paste0("01_Simulation_study/results/density_est_pool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
        density_est_pool = readRDS(nameopen)
        
        TV_pool = 0
        for(gg in 1:length(unique(data$group))){
          TV_pool = TV_pool + TV_distance(seqq, truth[gg,], density_est_pool$mean)
        }
        TV_pool = TV_pool/(length(unique(data$group)))
        newdata = c("Complete pooling", n_groups[i], sum(n_groups[i]/2*ssg*j), TV_pool)
        TVdist_df2[nrow(TVdist_df2)+1,] = newdata
        
        
        #-----#       NO POOLING      #-----#
        nameopen = paste0("01_Simulation_study/results/density_est_nopool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
        density_est_nopool = readRDS(nameopen)
        
        TV_nopool = 0
        for(gg in 1:length(unique(data$group))){
          TV_nopool = TV_nopool + TV_distance(seqq, truth[gg,], density_est_nopool[[gg]]$mean)
        }
        TV_nopool = TV_nopool/(length(unique(data$group)))
        newdata = c("No pooling", n_groups[i], sum(n_groups[i]/2*ssg*j), TV_nopool)
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
  
  saveRDS(TVdist_df2, file = "01_Simulation_study/output_RDS/TVdist_df_pool.RDS")
} else {
  TVdist_df2 = readRDS("01_Simulation_study/output_RDS/TVdist_df_pool.RDS")
}




ggplot(TVdist_df2, aes(x = n, y = TV, fill=Model ) ) +
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
  scale_fill_manual( values = c("forestgreen", "royalblue3", "deeppink4" ) ) +
  xlab("Sample size")  +
  ylab("TV distance")+
  facet_wrap( ~ G, scales = "free",
              labeller = labeller(G = c("2" = "2 groups",
                                        "10" = "10 groups") ))

ggsave("01_Simulation_study/output_images/05_01_TVdistance_comparison_pooling.pdf", width = 8, height = 4)





#-----------# #-----------# #-----------# #-----------# #-----------# #-----------# #-----------# 
#-----------# #-----------# #-----------# #-----------# #-----------# #-----------# #-----------# 

#-----------#  COMPARE TV DISTANCE FOR DIFFERENT MODELS #-----------#

if(!file.exists("01_Simulation_study/output_RDS/TVdist_df_models.RDS")) {
  TVdist_df = data.frame("Model" = character(), "G" = numeric(), "n" = numeric(), "TV" = double())
  
  for(repl in 1:n_datasets) {
    for(i in 1:length(n_groups) ){
      for(j in 1:n_ss){
        
        name_file_open = paste0("01_Simulation_study/data/data_", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
        data = readRDS(file = name_file_open)
        
        #-----#          THINNED DDP        #-----#
        nameopen = paste0("01_Simulation_study/results/density_est_thinnedDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
        density_est_thinnedDDP = readRDS(nameopen)
        
        seqq = density_est_thinnedDDP$seq
        
        nameopen = paste0("01_Simulation_study/results/true_density_", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
        truth = readRDS(nameopen)
        
        TV_thinnedDDP = 0
        for(gg in 1:length(unique(data$group))){
          TV_thinnedDDP = TV_thinnedDDP + TV_distance(seqq, truth[gg,], density_est_thinnedDDP$mean[,gg])
        }
        TV_thinnedDDP = TV_thinnedDDP/(length(unique(data$group)))
        newdata = c("Thinned DDP", n_groups[i], sum(n_groups[i]/2*ssg*j), TV_thinnedDDP)
        TVdist_df[nrow(TVdist_df)+1,] = newdata
        
        
        #-----#          CAM        #-----#
        nameopen = paste0("01_Simulation_study/results/density_est_CAM", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
        density_est_CAM = readRDS(file = nameopen)
        
        TV_CAM = 0
        for(gg in 1:length(unique(data$group))){
          TV_CAM = TV_CAM + TV_distance(seqq, truth[gg,], density_est_CAM$mean[,gg])
        }
        TV_CAM = TV_CAM/(length(unique(data$group)))
        newdata = c("CAM", n_groups[i], sum(n_groups[i]/2*ssg*j), TV_CAM)
        TVdist_df[nrow(TVdist_df)+1,] = newdata
        
        
        #-----#          GM-DDP        #-----#
        nameopen = paste0("01_Simulation_study/results/density_est_gmDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
        density_est_gmDDP = readRDS(file = nameopen)
        
        TV_gmDDP = 0
        for(gg in 1:length(unique(data$group))){
          TV_gmDDP = TV_gmDDP + TV_distance(seqq, truth[gg,], density_est_gmDDP$mean[,gg])
        }
        TV_gmDDP = TV_gmDDP/(length(unique(data$group)))
        newdata = c("GM-DDP", n_groups[i], sum(n_groups[i]/2*ssg*j), TV_gmDDP)
        TVdist_df[nrow(TVdist_df)+1,] = newdata
        
        
        #-----#          HDP        #-----#
        nameopen = paste0("01_Simulation_study/results/density_est_HDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
        density_est_HDP = readRDS(file = nameopen)
        
        TV_HDP = 0
        for(gg in 1:length(unique(data$group))){
          TV_HDP = TV_HDP + TV_distance(seqq, truth[gg,], density_est_HDP$mean[,gg])
        }
        TV_HDP = TV_HDP/(length(unique(data$group)))
        newdata = c("HDP", n_groups[i], sum(n_groups[i]/2*ssg*j), TV_HDP)
        TVdist_df[nrow(TVdist_df)+1,] = newdata
      }
    }
  }
  str(TVdist_df)
  
  TVdist_df$n = as.numeric(TVdist_df$n)
  TVdist_df$G = as.numeric(TVdist_df$G)
  TVdist_df$TV = as.numeric(TVdist_df$TV)
  
  TVdist_df = TVdist_df[order(TVdist_df$Model, TVdist_df$G, TVdist_df$n), ]
  
  TVdist_df$n = as.factor(TVdist_df$n)
  TVdist_df$G = as.factor(TVdist_df$G)
  TVdist_df$Model = as.factor(TVdist_df$Model)
  
  saveRDS(TVdist_df, file = "01_Simulation_study/output_RDS/TVdist_df_models.RDS")
} else {
  TVdist_df = readRDS("01_Simulation_study/output_RDS/TVdist_df_models.RDS")
}

ggplot(TVdist_df, aes(x = n, y = TV, fill=Model ) ) +
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
  scale_fill_manual( values = c("darkgoldenrod1", "cyan4", "blue", "deeppink4") ) +
  xlab("Sample size")  +
  ylab("TV distance")+
  facet_wrap( ~ G, scales = "free",
              labeller = labeller(G = c("2" = "2 groups",
                                        "10" = "10 groups") ))

ggsave("01_Simulation_study/output_images/05_02_TVdistance_comparison_models.pdf", width = 8, height = 4)









