
#----------------# #----------------# #----------------# 
#               PLOT ADJUSTED RAND INDEX               #
#----------------# #----------------# #----------------# 

# This script produces Figure 3 and Figure 5 in the article.


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




#-----------#  COMPARE ARI WITH COMPLETE AND NO-POOLING MODELS  #-----------#

if(!file.exists("01_Simulation_study/output_RDS/ARI_df_pool.RDS")){
  rand_df2 = data.frame("Model" = character(), "G" = numeric(), "n" = numeric(), "rand" = double())
  
  for(repl in 1:n_datasets) {
    for(i in 1:length(n_groups) ){
      for(j in 1:n_ss){
        nameopen = paste0("01_Simulation_study/results/rand_thinnedDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
        rand_thinnedDDP = readRDS(nameopen)
        newdata = c("Thinned DDP", n_groups[i], sum(n_groups[i]/2*ssg*j), rand_thinnedDDP)
        rand_df2[nrow(rand_df2)+1,] = newdata
        
        nameopen = paste0("01_Simulation_study/results/rand_pool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
        rand_pool = readRDS(nameopen)
        newdata = c("Complete pooling", n_groups[i], sum(n_groups[i]/2*ssg*j), rand_pool)
        rand_df2[nrow(rand_df2)+1,] = newdata
        
        nameopen = paste0("01_Simulation_study/results/rand_nopool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
        rand_nopool = readRDS(nameopen)
        newdata = c("No pooling", n_groups[i], sum(n_groups[i]/2*ssg*j), rand_nopool)
        rand_df2[nrow(rand_df2)+1,] = newdata
      }
    }
  }
  str(rand_df2)
  
  rand_df2$n = as.numeric(rand_df2$n)
  rand_df2$G = as.numeric(rand_df2$G)
  rand_df2$rand = as.numeric(rand_df2$rand)
  
  rand_df2 = rand_df2[order(rand_df2$Model, rand_df2$G, rand_df2$n), ]
  
  rand_df2$n = as.factor(rand_df2$n)
  rand_df2$G = as.factor(rand_df2$G)
  rand_df2$Model = as.factor(rand_df2$Model)
  
  saveRDS(rand_df2, file = "01_Simulation_study/output_RDS/ARI_df_pool.RDS")
  
} else {
  rand_df2 = readRDS("01_Simulation_study/output_RDS/ARI_df_pool.RDS")
}

ggplot(rand_df2, aes(x = n, y = rand, fill=Model ) ) +
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
  ylim(0,1) +
  ylab("Average ARI")+
  facet_wrap( ~ G, scales = "free",
              labeller = labeller(G = c("2" = "2 groups",
                                        "10" = "10 groups") ) )

ggsave("01_Simulation_study/output_images/04_01_ARI_comparison_pooling.pdf", width = 8, height = 4)





#-----------# #-----------# #-----------# #-----------# #-----------# #-----------# #-----------# 
#-----------# #-----------# #-----------# #-----------# #-----------# #-----------# #-----------# 

#-----------#  COMPARE ARI FOR DIFFERENT MODELS #-----------#

if(!file.exists("01_Simulation_study/output_RDS/ARI_df_models.RDS")){
  rand_df = data.frame("Model" = character(), "G" = numeric(), "n" = numeric(), "rand" = double())
  
  for(repl in 1:n_datasets) {
    for(i in 1:length(n_groups) ){
      for(j in 1:n_ss){
        nameopen = paste0("01_Simulation_study/results/rand_thinnedDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
        rand_thinnedDDP = readRDS(nameopen)
        newdata = c("Thinned DDP", n_groups[i], sum(n_groups[i]/2*ssg*j), rand_thinnedDDP)
        rand_df[nrow(rand_df)+1,] = newdata
        
        nameopen = paste0("01_Simulation_study/results/rand_CAM", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
        rand_CAM = readRDS(nameopen)
        newdata = c("CAM", n_groups[i], sum(n_groups[i]/2*ssg*j), rand_CAM)
        rand_df[nrow(rand_df)+1,] = newdata
        
        nameopen = paste0("01_Simulation_study/results/rand_gmDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
        rand_gmDDP = readRDS(nameopen)
        newdata = c("GM-DDP", n_groups[i], sum(n_groups[i]/2*ssg*j), rand_gmDDP)
        rand_df[nrow(rand_df)+1,] = newdata
        
        nameopen = paste0("01_Simulation_study/results/rand_HDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
        rand_HDP = readRDS(nameopen)
        newdata = c("HDP", n_groups[i], sum(n_groups[i]/2*ssg*j), rand_HDP)
        rand_df[nrow(rand_df)+1,] = newdata
      }
    }
  }
  str(rand_df)
  
  rand_df$n = as.numeric(rand_df$n)
  rand_df$G = as.numeric(rand_df$G)
  rand_df$rand = as.numeric(rand_df$rand)
  
  rand_df = rand_df[order(rand_df$Model, rand_df$G, rand_df$n), ]
  
  rand_df$n = as.factor(rand_df$n)
  rand_df$G = as.factor(rand_df$G)
  rand_df$Model = as.factor(rand_df$Model)
  
  saveRDS(rand_df, file = "01_Simulation_study/output_RDS/ARI_df_models.RDS")
  
} else {
  rand_df = readRDS("01_Simulation_study/output_RDS/ARI_df_models.RDS")
}


ggplot(rand_df, aes(x = n, y = rand, fill=Model ) ) +
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
  scale_fill_manual( values = c("darkgoldenrod1", "cyan4","blue", "deeppink4", "forestgreen", "salmon") ) +
  xlab("Sample size")  +
  ylab("Average ARI")+
  ylim(0,1)+
  facet_wrap( ~ G, scales = "free",
              labeller = labeller(G = c("2" = "2 groups",
                                        "10" = "10 groups") ))

ggsave("01_Simulation_study/output_images/04_02_ARI_comparison_models.pdf", width = 8, height = 4)







