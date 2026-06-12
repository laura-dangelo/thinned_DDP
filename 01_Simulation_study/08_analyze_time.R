
#----------------# #----------------# #----------------# 
#                  PLOT COMPUTING TIME                 #
#----------------# #----------------# #----------------# 

library(devtools)
library(salso)
library(mclust)
library(mcclust.ext)
library(ggplot2)
library(viridis)

n_groups = c(2, 10)
ssg = c(10, 30)
n_ss = 4 
n_datasets = 50

tot_datasets = n_datasets * n_ss * length(n_groups)

time_df = data.frame("Model" = character(), "G" = numeric(), "n" = numeric(), "time" = double())

for(repl in 1:n_datasets) {
  for(i in 1:length(n_groups) ){
    for(j in 1:n_ss){
      
      nameopen = paste0("01_Simulation_study/results/time_thinnedDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      time_thinnedDDP = readRDS(nameopen)
      newdata = c("Thinned DDP", n_groups[i], sum(n_groups[i]/2*ssg*j), time_thinnedDDP)
      time_df[nrow(time_df)+1,] = newdata

      nameopen = paste0("01_Simulation_study/results/time_CAM", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      time_CAM = readRDS(nameopen)
      newdata = c("CAM", n_groups[i], sum(n_groups[i]/2*ssg*j), time_CAM)
      time_df[nrow(time_df)+1,] = newdata
      
      nameopen = paste0("01_Simulation_study/results/time_gmDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".RDS")
      time_gmDDP = readRDS(nameopen)
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
  scale_fill_manual( values = c("darkgoldenrod1", "cyan4",  "deeppink4") ) +
  xlab("Sample size")  +
  ylab("log(time)")+
  facet_wrap( ~ G, scales = "free",
              labeller = labeller(G = c("2" = "2 groups",
                                        "10" = "10 groups") ) )

ggsave("01_Simulation_study/output_images/08_02_time.pdf", width = 8, height = 4)
