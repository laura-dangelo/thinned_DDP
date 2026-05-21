library(devtools)
library(salso)
library(mclust)
library(mcclust.ext)
library(ggplot2)
library(viridisLite)

n_groups = c(2, 10)
ssg = c(10, 30)
n_ss = 4 # number of different configuration of sample size
n_datasets = 50

tot_datasets = n_datasets * n_ss * length(n_groups)

# rand_df = data.frame("Model" = character(), "G" = numeric(), "n" = numeric(), "rand" = double())
# 
# 
# 
# for(repl in 1:n_datasets) {
#   for(i in 1:length(n_groups) ){
#     for(j in 1:n_ss){
#       
#       nameopen = paste0("../04_Simulation_commonatoms/Results/rand_thinnedDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
#       load(nameopen)
#       newdata = c("thinnedDDP", n_groups[i], sum(n_groups[i]/2*ssg*j), rand_thinnedDDP)
#       rand_df[nrow(rand_df)+1,] = newdata
# 
#       nameopen = paste0("../04_Simulation_commonatoms/Results/rand_CAM", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
#       load(nameopen)
#       newdata = c("CAM", n_groups[i], sum(n_groups[i]/2*ssg*j), rand_CAM)
#       rand_df[nrow(rand_df)+1,] = newdata
#       
#       nameopen = paste0("../04_Simulation_commonatoms/Results/rand_gmDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
#       load(nameopen)
#       newdata = c("GM-DDP", n_groups[i], sum(n_groups[i]/2*ssg*j), rand_gmDDP)
#       rand_df[nrow(rand_df)+1,] = newdata
#       
#     }
#   }
# }
# str(rand_df)
# 
# rand_df$n = as.numeric(rand_df$n)
# rand_df$G = as.numeric(rand_df$G)
# rand_df$rand = as.numeric(rand_df$rand)
# 
# rand_df = rand_df[order(rand_df$Model, rand_df$G, rand_df$n), ]
# 
# rand_df$n = as.factor(rand_df$n)
# rand_df$G = as.factor(rand_df$G)
# rand_df$Model = as.factor(rand_df$Model)
# 
# saveRDS(rand_df, file = "../04_Simulation_commonatoms/output_RDS/ARI_df_models.RDS")

rand_df = readRDS("../04_Simulation_commonatoms/output_RDS/ARI_df_models.RDS")


ggplot(rand_df, aes(x = n, y = rand, fill=Model ) ) +
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
  scale_fill_manual( values = c("darkgoldenrod1", "cyan4", "deeppink4", "forestgreen", "salmon") ) +
  xlab("Sample size")  +
  ylab("ARI")+
  ylim(min(rand_df$rand),1)+
  facet_wrap( ~ G, scales = "free",
              labeller = labeller(G = c("2" = "2 groups",
                                        "10" = "10 goups") )
  )

ggsave("../04_Simulation_commonatoms/03_02_RandIndex.pdf", width = 8, height = 4)





#-----------# #-----------# #-----------# #-----------# #-----------# #-----------# #-----------# 
#-----------# #-----------# #-----------# #-----------# #-----------# #-----------# #-----------# 



# 
# rand_df2 = data.frame("Model" = character(), "G" = numeric(), "n" = numeric(), "rand" = double())
# 
# for(repl in 1:n_datasets) {
#   for(i in 1:length(n_groups) ){
#     for(j in 1:n_ss){
#       
#       nameopen = paste0("../04_Simulation_commonatoms/Results/rand_thinnedDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
#       load(nameopen)
#       newdata = c("thinnedDDP", n_groups[i], sum(n_groups[i]/2*ssg*j), rand_thinnedDDP)
#       rand_df2[nrow(rand_df2)+1,] = newdata
#       
#       
#       nameopen = paste0("../04_Simulation_commonatoms/Results/rand_pool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
#       load(nameopen)
#       newdata = c("pool", n_groups[i], sum(n_groups[i]/2*ssg*j), rand_pool)
#       rand_df2[nrow(rand_df2)+1,] = newdata
# 
# 
#       nameopen = paste0("../04_Simulation_commonatoms/Results/rand_nopool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
#       load(nameopen)
#       newdata = c("no-pool", n_groups[i], sum(n_groups[i]/2*ssg*j), rand_nopool)
#       rand_df2[nrow(rand_df2)+1,] = newdata
#       
#       
#     }
#   }
# }
# str(rand_df2)
# 
# rand_df2$n = as.numeric(rand_df2$n)
# rand_df2$G = as.numeric(rand_df2$G)
# rand_df2$rand = as.numeric(rand_df2$rand)
# 
# rand_df2 = rand_df2[order(rand_df2$Model, rand_df2$G, rand_df2$n), ]
# 
# rand_df2$n = as.factor(rand_df2$n)
# rand_df2$G = as.factor(rand_df2$G)
# rand_df2$Model = as.factor(rand_df2$Model)
# 
# saveRDS(rand_df2, file = "../04_Simulation_commonatoms/output_RDS/ARI_df_pool.RDS")

rand_df2 = readRDS("../04_Simulation_commonatoms/output_RDS/ARI_df_pool.RDS")


ggplot(rand_df2, aes(x = n, y = rand, fill=Model ) ) +
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
  scale_fill_manual( values = c("forestgreen", "royalblue3", "deeppink4" ) ) +
  xlab("Sample size")  +
  ylab("ARI")+
  ylim(min(rand_df2$rand),1)+
  facet_wrap( ~ G, scales = "free",
              labeller = labeller(G = c("2" = "2 groups",
                                        "10" = "10 goups") )
  )

ggsave("../04_Simulation_commonatoms/03_01_RandIndex.pdf", width = 8, height = 4)

