
#----------------# #----------------# #----------------# 
#           COMPUTE AREA CREDIBLE INTERVALS            #
#----------------# #----------------# #----------------# 


library(TeachingDemos)
library(ggplot2)

n_groups = c(2, 10)
ssg = c(10, 30)
n_ss = 4 # number of different configuration of sample size
n_datasets = 50

tot_datasets = n_datasets * n_ss * length(n_groups)


#-----------#  COMPARE AREA WITH COMPLETE AND NO-POOLING MODELS  #-----------#

area_df2 = data.frame("Model" = character(), "G" = numeric(), "n" = numeric(), "area" = numeric())

for(repl in 1:n_datasets) {
  for(i in 1:length(n_groups) ){
    for(j in 1:n_ss){

      #-----# #-----# #-----# #-----#
      #-----#   Thinned DDP   #-----#
      #-----# #-----# #-----# #-----#

      nameopen = paste0("01_Simulation_study/results/density_CI_thinnedDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      density_CI_thinnedDDP = readRDS(nameopen)

      diff_area = density_CI_thinnedDDP$upper - density_CI_thinnedDDP$lower
      area_CI_thinnedDDP = tapply(diff_area, factor(density_CI_thinnedDDP$Group), sum)
      area_CI_thinnedDDP = mean(area_CI_thinnedDDP)

      newdata = c("Thinned DDP", n_groups[i], sum(n_groups[i]/2*ssg*j), area_CI_thinnedDDP)
      area_df2[nrow(area_df2)+1,] = newdata

      #-----# #-----# #-----# #-----#
      #-----#      Pool       #-----#
      #-----# #-----# #-----# #-----#

      nameopen = paste0("01_Simulation_study/results/density_CI_pool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      density_CI_pool = readRDS(nameopen)

      diff_area = density_CI_pool$upper - density_CI_pool$lower
      area_CI_pool = sum(diff_area)

      newdata = c("Complete pooling", n_groups[i], sum(n_groups[i]/2*ssg*j), area_CI_pool)
      area_df2[nrow(area_df2)+1,] = newdata

      #-----# #-----# #-----# #-----#
      #-----#     No pool     #-----#
      #-----# #-----# #-----# #-----#

      nameopen = paste0("01_Simulation_study/results/density_CI_nopool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      density_CI_nopool = readRDS(nameopen)

      diff_area = density_CI_nopool$upper - density_CI_nopool$lower
      area_CI_nopool = tapply(diff_area, factor(density_CI_nopool$Group), sum)
      area_CI_nopool = mean(area_CI_nopool)

      newdata = c("No pooling", n_groups[i], sum(n_groups[i]/2*ssg*j), area_CI_nopool)
      area_df2[nrow(area_df2)+1,] = newdata

    }
  }
}

str(area_df2)

area_df2$n = as.numeric(area_df2$n)
area_df2$G = as.numeric(area_df2$G)
area_df2$area = as.numeric(area_df2$area)

area_df2 = area_df2[order(area_df2$Model, area_df2$G, area_df2$n), ]

area_df2$n = as.factor(area_df2$n)
area_df2$G = as.factor(area_df2$G)
area_df2$Model = as.factor(area_df2$Model)

saveRDS(area_df2, file = "01_Simulation_study/output_RDS/area_df_pool.RDS")

area_df2 = readRDS("01_Simulation_study/output_RDS/area_df_pool.RDS")
area_df2$area = area_df2$area/300


ggplot(area_df2, aes(x = n, y = area, fill=Model ) ) +
  geom_boxplot(alpha = 0.6) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  theme_bw() +
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
  scale_fill_manual( values = c("forestgreen", "royalblue3", "deeppink4") ) +
  xlab("Sample size")  +
  ylab("Area credible bands")+
  facet_wrap( ~ G, scales = "free",
              labeller = labeller(G = c("2" = "2 groups",
                                        "10" = "10 groups") )
  )


ggsave("01_Simulation_study/output_images/07_01_areaCI.pdf", width = 8, height = 4)





#-----------# #-----------# #-----------# #-----------# #-----------# #-----------# #-----------# 
#-----------# #-----------# #-----------# #-----------# #-----------# #-----------# #-----------# 

#-----------#  COMPARE AREA FOR DIFFERENT MODELS #-----------#

area_df = data.frame("Model" = character(), "G" = numeric(), "n" = numeric(), "area" = numeric())


for(repl in 1:n_datasets) {
  for(i in 1:length(n_groups) ){
    for(j in 1:n_ss){

      #-----# #-----# #-----# #-----#
      #-----#   Thinned DDP   #-----#
      #-----# #-----# #-----# #-----#

      nameopen = paste0("01_Simulation_study/results/density_CI_thinnedDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      density_CI_thinnedDDP = readRDS(nameopen)

      diff_area = density_CI_thinnedDDP$upper - density_CI_thinnedDDP$lower
      area_CI_thinnedDDP = tapply(diff_area, factor(density_CI_thinnedDDP$Group), sum)
      area_CI_thinnedDDP = mean(area_CI_thinnedDDP)

      newdata = c("Thinned DDP", n_groups[i], sum(n_groups[i]/2*ssg*j), area_CI_thinnedDDP)
      area_df[nrow(area_df)+1,] = newdata


      #-----# #-----# #-----# #-----#
      #-----#       CAM       #-----#
      #-----# #-----# #-----# #-----#

      nameopen = paste0("01_Simulation_study/results/density_CI_CAM", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      density_CI_CAM = readRDS(nameopen)

      diff_area = density_CI_CAM$upper - density_CI_CAM$lower
      area_CI_CAM = tapply(diff_area, factor(density_CI_CAM$Group), sum)
      area_CI_CAM = mean(area_CI_CAM)

      newdata = c("CAM", n_groups[i], sum(n_groups[i]/2*ssg*j), area_CI_CAM)
      area_df[nrow(area_df)+1,] = newdata


      #-----# #-----# #-----# #-----#
      #-----#      GM-DDP     #-----#
      #-----# #-----# #-----# #-----#

      nameopen = paste0("01_Simulation_study/results/density_CI_gmDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      density_CI_gmDDP = readRDS(nameopen)

      diff_area = density_CI_gmDDP$upper - density_CI_gmDDP$lower
      area_CI_gmDDP = tapply(diff_area, factor(density_CI_gmDDP$Group), sum)
      area_CI_gmDDP = mean(area_CI_gmDDP)

      newdata = c("GM-DDP", n_groups[i], sum(n_groups[i]/2*ssg*j), area_CI_gmDDP)
      area_df[nrow(area_df)+1,] = newdata
      
      
      #-----# #-----# #-----# #-----#
      #-----#        HDP      #-----#
      #-----# #-----# #-----# #-----#
      
      nameopen = paste0("01_Simulation_study/results/density_CI_HDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      density_CI_HDP = readRDS(nameopen)
      
      diff_area = density_CI_HDP$upper - density_CI_HDP$lower
      area_CI_HDP = tapply(diff_area, factor(density_CI_HDP$Group), sum)
      area_CI_HDP = mean(area_CI_HDP)
      
      newdata = c("HDP", n_groups[i], sum(n_groups[i]/2*ssg*j), area_CI_HDP)
      area_df[nrow(area_df)+1,] = newdata

    }
  }
}

str(area_df)

area_df$n = as.numeric(area_df$n)
area_df$G = as.numeric(area_df$G)
area_df$area = as.numeric(area_df$area)

area_df = area_df[order(area_df$Model, area_df$G, area_df$n), ]

area_df$n = as.factor(area_df$n)
area_df$G = as.factor(area_df$G)
area_df$Model = as.factor(area_df$Model)

saveRDS(area_df, file = "01_Simulation_study/output_RDS/area_df_models.RDS")

area_df = readRDS("01_Simulation_study/output_RDS/area_df_models.RDS")

area_df$area = area_df$area/300

ggplot(area_df, aes(x = n, y = area, fill=Model ) ) +
  geom_boxplot(alpha = 0.6) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  theme_bw() +
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
  scale_fill_manual( values = c("darkgoldenrod1", "cyan4", "blue", "deeppink4") ) +
  xlab("Sample size")  +
  ylab("Area credible bands")+
  facet_wrap( ~ G, scales = "free",
              labeller = labeller(G = c("2" = "2 groups",
                                        "10" = "10 groups") )
  )

ggsave("01_Simulation_study/output_images/07_02_areaCI.pdf", width = 8, height = 4)










