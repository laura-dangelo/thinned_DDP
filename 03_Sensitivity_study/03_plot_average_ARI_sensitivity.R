
#----------------# #----------------# #----------------# 
#               PLOT ADJUSTED RAND INDEX               #
#----------------# #----------------# #----------------# 

# This script produces Figure F2 of the Supplementary Material

library(devtools)
library(salso)
library(mclust)
library(mcclust.ext)
library(ggplot2)
library(viridisLite)

n_groups = c(2, 10)
ssg = c(10, 30)
n_ss = 4 
n_datasets = 25
tot_datasets = n_datasets * n_ss * length(n_groups)


if(!file.exists("03_Sensitivity_study/output_RDS/ARI_sensitivity.RDS")){
  rand_df2 = data.frame("Prior" = character(), "G" = numeric(), "n" = numeric(), "rand" = double())
  
  for(repl in 1:n_datasets) {
    for(i in 1:length(n_groups) ){
      for(j in 1:n_ss){
        
        nameopen = paste0("01_Simulation_study/results/rand_thinnedDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
        rand_thinnedDDP = readRDS(nameopen)
        newdata = c("Beta(3, 3)", n_groups[i], sum(n_groups[i]/2*ssg*j), rand_thinnedDDP)
        rand_df2[nrow(rand_df2)+1,] = newdata
        
        nameopen = paste0("03_Sensitivity_study/results_sensitivity/rand_thinnedDDP_unif", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
        rand_thinnedDDP_unif = readRDS(nameopen)
        newdata = c("Uniform", n_groups[i], sum(n_groups[i]/2*ssg*j), rand_thinnedDDP_unif)
        rand_df2[nrow(rand_df2)+1,] = newdata
        
        nameopen = paste0("03_Sensitivity_study/results_sensitivity/rand_thinnedDDP_1010", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
        rand_thinnedDDP_1010 = readRDS(nameopen)
        newdata = c("Beta(10, 10)", n_groups[i], sum(n_groups[i]/2*ssg*j), rand_thinnedDDP_1010)
        rand_df2[nrow(rand_df2)+1,] = newdata
        
        nameopen = paste0("03_Sensitivity_study/results_sensitivity/rand_thinnedDDP_0101", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
        rand_thinnedDDP_0101 = readRDS(nameopen)
        newdata = c("Beta(0.1, 0.1)", n_groups[i], sum(n_groups[i]/2*ssg*j), rand_thinnedDDP_0101)
        rand_df2[nrow(rand_df2)+1,] = newdata
      }
    }
  }
  str(rand_df2)
  
  rand_df2$n = as.numeric(rand_df2$n)
  rand_df2$G = as.numeric(rand_df2$G)
  rand_df2$rand = as.numeric(rand_df2$rand)
  
  rand_df2 = rand_df2[order(rand_df2$Prior, rand_df2$G, rand_df2$n), ]
  rand_df2$n = as.factor(rand_df2$n)
  rand_df2$G = as.factor(rand_df2$G)
  rand_df2$Prior = as.factor(rand_df2$Prior)
  
  saveRDS(rand_df2, file = "03_Sensitivity_study/output_RDS/ARI_sensitivity.RDS")
  
} else {
  rand_df2 = readRDS("03_Sensitivity_study/output_RDS/ARI_sensitivity.RDS")
}

ggplot(rand_df2, aes(x = n, y = rand, fill=Prior ) ) +
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
  scale_fill_manual( values = c("forestgreen", "royalblue3", "deeppink4","orange" ) ) +
  xlab("Sample size")  +
  ylim(0.7,1) +
  ylab("Average ARI")+
  facet_wrap( ~ G, scales = "free",
              labeller = labeller(G = c("2" = "2 groups",
                                        "10" = "10 groups") )
  )

ggsave("03_Sensitivity_study/output_images/04_01_ARI_comparison_sensitivity.pdf", width = 8, height = 4)


