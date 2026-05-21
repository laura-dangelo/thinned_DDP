library(devtools)
library(salso)
library(mclust)
library(mcclust.ext)
library(ggplot2)
library(viridis)

n_groups = c(2, 10)
ssg = c(10, 30)
n_ss = 4 # number of different configuration of sample size
n_datasets = 50

tot_datasets = n_datasets * n_ss * length(n_groups)

time_df = data.frame("Model" = character(), "G" = numeric(), "n" = numeric(), "time" = double())



for(repl in 1:n_datasets) {
  for(i in 1:length(n_groups) ){
    for(j in 1:n_ss){
      
      nameopen = paste0("../04_Simulation_commonatoms/Results/time_thinnedDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
      load(nameopen)
      newdata = c("thinnedDDP", n_groups[i], sum(n_groups[i]/2*ssg*j), time_thinnedDDP)
      time_df[nrow(time_df)+1,] = newdata
      
      
      # nameopen = paste0("../04_Simulation_commonatoms/Results/time_pool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
      # load(nameopen)
      # newdata = c("pool", n_groups[i], sum(n_groups[i]/2*ssg*j), time_pool)
      # time_df[nrow(time_df)+1,] = newdata
      # 
      # 
      # nameopen = paste0("../04_Simulation_commonatoms/Results/time_nopool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
      # load(nameopen)
      # newdata = c("no-pool", n_groups[i], sum(n_groups[i]/2*ssg*j), time_nopool)
      # time_df[nrow(time_df)+1,] = newdata
      
      
      nameopen = paste0("../04_Simulation_commonatoms/Results/time_CAM", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
      load(nameopen)
      newdata = c("CAM", n_groups[i], sum(n_groups[i]/2*ssg*j), time_CAM)
      time_df[nrow(time_df)+1,] = newdata
      
      
      nameopen = paste0("../04_Simulation_commonatoms/Results/time_gmDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
      load(nameopen)
      newdata = c("GM-DDP", n_groups[i], sum(n_groups[i]/2*ssg*j), time_gmDDP)
      time_df[nrow(time_df)+1,] = newdata
      
      
    }
  }
}
str(time_df)

time_df$n = as.numeric(time_df$n)
time_df$G = as.numeric(time_df$G)
time_df$time = as.numeric(time_df$time)

time_df = time_df[order(time_df$Model, time_df$G, time_df$n), ]

time_df$n = as.factor(time_df$n)
time_df$G = as.factor(time_df$G)
time_df$Model = as.factor(time_df$Model)

ggplot(time_df, aes(x = n, y = log(time), fill=Model ) ) +
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
  scale_fill_manual( values = c("darkgoldenrod1", "cyan4", "deeppink4") ) +
  xlab("Sample size")  +
  ylab("log(time)")+
  facet_wrap( ~ G, scales = "free",
              labeller = labeller(G = c("2" = "2 groups",
                                        "10" = "10 goups") )
  )

ggsave("../04_Simulation_commonatoms/08_02_time.pdf", width = 8, height = 4)
