library(ggplot2)

n_groups = c(2, 10)
ssg = c(10, 30)
n_ss = 4 # number of different configuration of sample size
n_datasets = 50

tot_datasets = n_datasets * n_ss * length(n_groups)


repl = 10
i = 1
j = 1

name_file_open = paste0("../04_Simulation_commonatoms/Data/data_", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
load(file = name_file_open)

nameopen = paste0("../04_Simulation_commonatoms/Results/density_CI_thinnedDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
load(nameopen)

nameopen = paste0("../04_Simulation_commonatoms/Results/density_CI_CAM", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
load(nameopen)

nameopen = paste0("../04_Simulation_commonatoms/Results/density_CI_gmDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
load(nameopen)

data = data.frame("Group" = data$group, "y" = data$y)


cols = c("Thinned DDP"="deeppink4", "CAM"="darkgoldenrod3")

ggplot(data = density_CI_thinnedDDP, aes(x = Seq, y = mean)) + 
  geom_histogram(data = data, aes(x = y, y=..density..),
                 binwidth = 1.5, colour="gray30", fill="gray80",lwd=0.3, alpha = 0.3) +
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha=0.3, col = "maroon", fill = "maroon", lwd= 0.3) +
  geom_ribbon(data = density_CI_CAM, aes(ymin = lower, ymax = upper), 
              alpha=0.4, col = "gold3", fill = "gold", lwd= 0.3) +
  geom_line(aes(col = "Thinned DDP"), lwd = 1) +
  geom_line(data = density_CI_CAM, aes(x = Seq, y=mean, col="CAM"), lwd = 1 ) +
  geom_line(aes(y = true), lty = "dashed") + 
  scale_colour_manual(name="Model", values = cols) + 
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
  # scale_y_continuous(breaks = c(0.0, 0.1, 0.2,0.3),  limits = c(0,0.3)) +
  xlab("y")  + ylab("Density")+
  facet_wrap( ~reorder(Group, sort(as.numeric(Group))), ncol = 2, dir="h")

ggsave("../04_Simulation_commonatoms/06_02_densityCAM_2G_minSS.pdf", width = 8, height = 3.5)





cols2 = c("Thinned DDP"="deeppink4",  "GM-DDP" = "#114f4f")
ggplot(data = density_CI_thinnedDDP, aes(x = Seq, y = mean)) + 
  geom_histogram(data = data, aes(x = y, y=..density..),
                 binwidth = 1.5, colour="gray30", fill="gray80",lwd=0.3, alpha = 0.3) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.3, 
              col = "maroon", fill = "maroon", lwd= 0.3) +
  geom_ribbon(data = density_CI_gmDDP, aes(ymin = lower, ymax = upper), alpha=0.3, 
              col = "#27A8A8", fill = "#27A8A8", lwd= 0.3) +
  geom_line(aes(col = "Thinned DDP"), lwd = 1) +
  geom_line(aes(y = true), lty = "dashed") + 
  geom_line(data = density_CI_gmDDP, aes(x = Seq, y=mean, col="GM-DDP") , lwd=1) +
  scale_colour_manual(name="Model", values = cols2) + 
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
  xlab("y")  + ylab("Density")+
  facet_wrap( ~reorder(Group, sort(as.numeric(Group))), ncol = 2, dir="h")

ggsave("../04_Simulation_commonatoms/06_02_densityGM_2G_minSS.pdf", width = 8, height = 3.5)




#-----------# #-----------# #-----------# #-----------# #-----------# #-----------# #-----------# 


repl = 10
i = 1
j = 4

name_file_open = paste0("../04_Simulation_commonatoms/Data/data_", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
load(file = name_file_open)

nameopen = paste0("../04_Simulation_commonatoms/Results/density_CI_thinnedDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
load(nameopen)

nameopen = paste0("../04_Simulation_commonatoms/Results/density_CI_CAM", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
load(nameopen)

nameopen = paste0("../04_Simulation_commonatoms/Results/density_CI_gmDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
load(nameopen)

data = data.frame("Group" = data$group, "y" = data$y)


cols = c("Thinned DDP"="deeppink4", "CAM"="darkgoldenrod3")

ggplot(data = density_CI_thinnedDDP, aes(x = Seq, y = mean)) + 
  geom_histogram(data = data, aes(x = y, y=..density..),
                 binwidth = 1.5, colour="gray30", fill="gray80",lwd=0.3, alpha = 0.3) +
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha=0.3, col = "maroon", fill = "maroon", lwd= 0.3) +
  geom_ribbon(data = density_CI_CAM, aes(ymin = lower, ymax = upper), 
              alpha=0.4, col = "gold3", fill = "gold", lwd= 0.3) +
  geom_line(aes(col = "Thinned DDP"), lwd = 1) +
  geom_line(data = density_CI_CAM, aes(x = Seq, y=mean, col="CAM"), lwd = 1 ) +
  geom_line(aes(y = true), lty = "dashed") + 
  scale_colour_manual(name="Model", values = cols) + 
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
  # scale_y_continuous(breaks = c(0.0, 0.1, 0.2,0.3),  limits = c(0,0.3)) +
  xlab("y")  + ylab("Density")+
  facet_wrap( ~reorder(Group, sort(as.numeric(Group))), ncol = 2, dir="h")

ggsave("../04_Simulation_commonatoms/06_02_densityCAM_2G_maxSS.pdf", width = 8, height = 3.5)





cols2 = c("Thinned DDP"="deeppink4",  "GM-DDP" = "#114f4f")
ggplot(data = density_CI_thinnedDDP, aes(x = Seq, y = mean)) + 
  geom_histogram(data = data, aes(x = y, y=..density..),
                 binwidth = 1.5, colour="gray30", fill="gray80",lwd=0.3, alpha = 0.3) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.3, 
              col = "maroon", fill = "maroon", lwd= 0.3) +
  geom_ribbon(data = density_CI_gmDDP, aes(ymin = lower, ymax = upper), alpha=0.3, 
              col = "#27A8A8", fill = "#27A8A8", lwd= 0.3) +
  geom_line(aes(col = "Thinned DDP"), lwd = 1) +
  geom_line(aes(y = true), lty = "dashed") + 
  geom_line(data = density_CI_gmDDP, aes(x = Seq, y=mean, col="GM-DDP") , lwd = 1) +
  scale_colour_manual(name="Model", values = cols2) + 
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
  xlab("y")  + ylab("Density")+
  facet_wrap( ~reorder(Group, sort(as.numeric(Group))), ncol = 2, dir="h")

ggsave("../04_Simulation_commonatoms/06_02_densityGM_2G_maxSS.pdf", width = 8, height = 3.5)



#-----------# #-----------# #-----------# #-----------# #-----------# #-----------# #-----------# 



repl = 10
i = 2
j = 1

name_file_open = paste0("../04_Simulation_commonatoms/Data/data_", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
load(file = name_file_open)

nameopen = paste0("../04_Simulation_commonatoms/Results/density_CI_thinnedDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
load(nameopen)

nameopen = paste0("../04_Simulation_commonatoms/Results/density_CI_CAM", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
load(nameopen)

nameopen = paste0("../04_Simulation_commonatoms/Results/density_CI_gmDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
load(nameopen)

data = data.frame("Group" = data$group, "y" = data$y)


cols = c("Thinned DDP"="deeppink4", "CAM"="darkgoldenrod3")

ggplot(data = density_CI_thinnedDDP, aes(x = Seq, y = mean)) + 
  geom_histogram(data = data, aes(x = y, y=..density..),
                 binwidth = 1.5, colour="gray30", fill="gray80",lwd=0.3, alpha = 0.3) +
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha=0.3, col = "maroon", fill = "maroon", lwd= 0.3) +
  geom_ribbon(data = density_CI_CAM, aes(ymin = lower, ymax = upper), 
              alpha=0.4, col = "gold3", fill = "gold", lwd= 0.3) +
  geom_line(aes(col = "Thinned DDP"), lwd = 1) +
  geom_line(data = density_CI_CAM, aes(x = Seq, y=mean, col="CAM"), lwd = 1 ) +
  geom_line(aes(y = true), lty = "dashed") + 
  scale_colour_manual(name="Model", values = cols) + 
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
  # scale_y_continuous(breaks = c(0.0, 0.1, 0.2,0.3),  limits = c(0,0.3)) +
  xlab("y")  + ylab("Density")+
  facet_wrap( ~reorder(Group, sort(as.numeric(Group))), ncol = 2, dir="h")

ggsave("../04_Simulation_commonatoms/06_02_densityCAM_10G_minSS.pdf", width = 8, height = 10)





cols2 = c("Thinned DDP"="deeppink4",  "GM-DDP" = "#114f4f")
ggplot(data = density_CI_thinnedDDP, aes(x = Seq, y = mean)) + 
  geom_histogram(data = data, aes(x = y, y=..density..),
                 binwidth = 1.5, colour="gray30", fill="gray80",lwd=0.3, alpha = 0.3) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.3, 
              col = "maroon", fill = "maroon", lwd= 0.3) +
  geom_ribbon(data = density_CI_gmDDP, aes(ymin = lower, ymax = upper), alpha=0.3, 
              col = "#27A8A8", fill = "#27A8A8", lwd= 0.3) +
  geom_line(aes(col = "Thinned DDP"), lwd = 1) +
  geom_line(aes(y = true), lty = "dashed") + 
  geom_line(data = density_CI_gmDDP, aes(x = Seq, y=mean, col="GM-DDP"), lwd = 1 ) +
  scale_colour_manual(name="Model", values = cols2) + 
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
  xlab("y")  + ylab("Density")+
  facet_wrap( ~reorder(Group, sort(as.numeric(Group))), ncol = 2, dir="h")

ggsave("../04_Simulation_commonatoms/06_02_densityGM_10G_minSS.pdf", width = 8, height = 10)




#-----------# #-----------# #-----------# #-----------# #-----------# #-----------# #-----------# 


repl = 10
i = 2
j = 4

name_file_open = paste0("../04_Simulation_commonatoms/Data/data_", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
load(file = name_file_open)

nameopen = paste0("../04_Simulation_commonatoms/Results/density_CI_thinnedDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
load(nameopen)

nameopen = paste0("../04_Simulation_commonatoms/Results/density_CI_CAM", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
load(nameopen)

nameopen = paste0("../04_Simulation_commonatoms/Results/density_CI_gmDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
load(nameopen)

data = data.frame("Group" = data$group, "y" = data$y)


cols = c("Thinned DDP"="deeppink4", "CAM"="darkgoldenrod3")

ggplot(data = density_CI_thinnedDDP, aes(x = Seq, y = mean)) + 
  geom_histogram(data = data, aes(x = y, y=..density..),
                 binwidth = 1.5, colour="gray30", fill="gray80",lwd=0.3, alpha = 0.3) +
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha=0.3, col = "maroon", fill = "maroon", lwd= 0.3) +
  geom_ribbon(data = density_CI_CAM, aes(ymin = lower, ymax = upper), 
              alpha=0.4, col = "gold3", fill = "gold", lwd= 0.3) +
  geom_line(aes(col = "Thinned DDP"), lwd = 1) +
  geom_line(data = density_CI_CAM, aes(x = Seq, y=mean, col="CAM"), lwd = 1 ) +
  geom_line(aes(y = true), lty = "dashed") + 
  scale_colour_manual(name="Model", values = cols) + 
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
  # scale_y_continuous(breaks = c(0.0, 0.1, 0.2,0.3),  limits = c(0,0.3)) +
  xlab("y")  + ylab("Density")+
  facet_wrap( ~reorder(Group, sort(as.numeric(Group))), ncol = 2, dir="h")

ggsave("../04_Simulation_commonatoms/06_02_densityCAM_10G_maxSS.pdf", width = 8, height = 10)





cols2 = c("Thinned DDP"="deeppink4",  "GM-DDP" = "#114f4f")
ggplot(data = density_CI_thinnedDDP, aes(x = Seq, y = mean)) + 
  geom_histogram(data = data, aes(x = y, y=..density..),
                 binwidth = 1.5, colour="gray30", fill="gray80",lwd=0.3, alpha = 0.3) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.3, 
              col = "maroon", fill = "maroon", lwd= 0.3) +
  geom_ribbon(data = density_CI_gmDDP, aes(ymin = lower, ymax = upper), alpha=0.3, 
              col = "#27A8A8", fill = "#27A8A8", lwd= 0.3) +
  geom_line(aes(col = "Thinned DDP"), lwd = 1) +
  geom_line(aes(y = true), lty = "dashed") + 
  geom_line(data = density_CI_gmDDP, aes(x = Seq, y=mean, col="GM-DDP"), lwd = 1 ) +
  scale_colour_manual(name="Model", values = cols2) + 
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
  xlab("y")  + ylab("Density")+
  facet_wrap( ~reorder(Group, sort(as.numeric(Group))), ncol = 2, dir="h")

ggsave("../04_Simulation_commonatoms/06_02_densityGM_10G_maxSS.pdf", width = 8, height = 10)



















#-----------# #-----------# #-----------# #-----------# #-----------# #-----------# #-----------# 
#-----------# #-----------# #-----------# #-----------# #-----------# #-----------# #-----------# 




repl = 10
i = 1
j = 1

name_file_open = paste0("../04_Simulation_commonatoms/Data/data_", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
load(file = name_file_open)

nameopen = paste0("../04_Simulation_commonatoms/Results/density_CI_thinnedDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
load(nameopen)

nameopen = paste0("../04_Simulation_commonatoms/Results/density_CI_pool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
load(nameopen)

nameopen = paste0("../04_Simulation_commonatoms/Results/density_CI_nopool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
load(nameopen)

data = data.frame("Group" = data$group, "y" = data$y)


cols = c("Thinned DDP"="deeppink4", "no-pool"="forestgreen")

ggplot(data = density_CI_thinnedDDP, aes(x = Seq, y = mean)) + 
  geom_histogram(data = data, aes(x = y, y=..density..),
                 binwidth = 1.5, colour="gray30", fill="gray80",lwd=0.3, alpha = 0.3) +
  geom_ribbon(data = density_CI_nopool, aes(ymin = lower, ymax = upper), 
              alpha=0.4, col = "forestgreen", fill = "forestgreen", lwd= 0.3) +
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha=0.4, col = "maroon", fill = "maroon", lwd= 0.3) +
  geom_line(aes(col = "Thinned DDP"), lwd = 1) +
  geom_line(data = density_CI_nopool, aes(x = Seq, y=mean, col="no-pool"), lwd = 1 ) +
  geom_line(aes(y = true), lty = "dashed") + 
  scale_colour_manual(name="Model", values = cols) + 
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
  # scale_y_continuous(breaks = c(0.0, 0.1, 0.2,0.3),  limits = c(0,0.3)) +
  xlab("y")  + ylab("Density")+
  facet_wrap( ~reorder(Group, sort(as.numeric(Group))), ncol = 2, dir="h")

ggsave("../04_Simulation_commonatoms/06_01_densitynopool_2G_minSS.pdf", width = 8, height = 3.5)




cols2 = c("Thinned DDP"="deeppink4",  "pool" = "royalblue3")
ggplot(data = density_CI_thinnedDDP, aes(x = Seq, y = mean)) + 
  geom_histogram(data = data, aes(x = y, y=..density..),
                 binwidth = 1.5, colour="gray30", fill="gray80",lwd=0.3, alpha = 0.3) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.4, 
              col = "maroon", fill = "maroon", lwd= 0.4) +
  geom_ribbon(data = density_CI_pool, aes(ymin = lower, ymax = upper), alpha=0.3, 
              col = "royalblue3", fill = "royalblue3", lwd= 0.4) +
  geom_line(aes(col = "Thinned DDP")) +
  geom_line(aes(y = true), lty = "dashed") + 
  geom_line(data = density_CI_pool, aes(x = Seq, y=mean, col="pool") ) +
  scale_colour_manual(name="Model", values = cols2) + 
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
  xlab("y")  + ylab("Density")+
  facet_wrap( ~reorder(Group, sort(as.numeric(Group))), ncol = 2, dir="h")

ggsave("../04_Simulation_commonatoms/06_01_densitypool_2G_minSS.pdf", width = 8, height = 3.5)


#-----------# #-----------# #-----------# #-----------# #-----------# #-----------# #-----------# 




repl = 10
i = 1
j = 4

name_file_open = paste0("../04_Simulation_commonatoms/Data/data_", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
load(file = name_file_open)

nameopen = paste0("../04_Simulation_commonatoms/Results/density_CI_thinnedDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
load(nameopen)

nameopen = paste0("../04_Simulation_commonatoms/Results/density_CI_pool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
load(nameopen)

nameopen = paste0("../04_Simulation_commonatoms/Results/density_CI_nopool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
load(nameopen)

data = data.frame("Group" = data$group, "y" = data$y)


cols = c("Thinned DDP"="deeppink4", "no-pool"="forestgreen")

ggplot(data = density_CI_thinnedDDP, aes(x = Seq, y = mean)) + 
  geom_histogram(data = data, aes(x = y, y=..density..),
                 binwidth = 1.5, colour="gray30", fill="gray80",lwd=0.3, alpha = 0.3) +
  geom_ribbon(data = density_CI_nopool, aes(ymin = lower, ymax = upper), 
              alpha=0.4, col = "forestgreen", fill = "forestgreen", lwd= 0.3) +
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha=0.4, col = "maroon", fill = "maroon", lwd= 0.3) +
  geom_line(aes(col = "Thinned DDP"), lwd = 1) +
  geom_line(data = density_CI_nopool, aes(x = Seq, y=mean, col="no-pool"), lwd = 1 ) +
  geom_line(aes(y = true), lty = "dashed") + 
  scale_colour_manual(name="Model", values = cols) + 
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
  # scale_y_continuous(breaks = c(0.0, 0.1, 0.2,0.3),  limits = c(0,0.3)) +
  xlab("y")  + ylab("Density")+
  facet_wrap( ~reorder(Group, sort(as.numeric(Group))), ncol = 2, dir="h")

ggsave("../04_Simulation_commonatoms/06_01_densitynopool_2G_maxSS.pdf", width = 8, height = 3.5)




cols2 = c("Thinned DDP"="deeppink4",  "pool" = "royalblue3")
ggplot(data = density_CI_thinnedDDP, aes(x = Seq, y = mean)) + 
  geom_histogram(data = data, aes(x = y, y=..density..),
                 binwidth = 1.5, colour="gray30", fill="gray80",lwd=0.3, alpha = 0.3) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.4, 
              col = "maroon", fill = "maroon", lwd= 0.4) +
  geom_ribbon(data = density_CI_pool, aes(ymin = lower, ymax = upper), alpha=0.3, 
              col = "royalblue3", fill = "royalblue3", lwd= 0.4) +
  geom_line(aes(col = "Thinned DDP")) +
  geom_line(aes(y = true), lty = "dashed") + 
  geom_line(data = density_CI_pool, aes(x = Seq, y=mean, col="pool") ) +
  scale_colour_manual(name="Model", values = cols2) + 
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
  xlab("y")  + ylab("Density")+
  facet_wrap( ~reorder(Group, sort(as.numeric(Group))), ncol = 2, dir="h")

ggsave("../04_Simulation_commonatoms/06_01_densitypool_2G_maxSS.pdf", width = 8, height = 3.5)




#-----------# #-----------# #-----------# #-----------# #-----------# #-----------# #-----------# 




repl = 10
i = 2
j = 1

name_file_open = paste0("../04_Simulation_commonatoms/Data/data_", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
load(file = name_file_open)

nameopen = paste0("../04_Simulation_commonatoms/Results/density_CI_thinnedDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
load(nameopen)

nameopen = paste0("../04_Simulation_commonatoms/Results/density_CI_pool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
load(nameopen)

nameopen = paste0("../04_Simulation_commonatoms/Results/density_CI_nopool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
load(nameopen)

data = data.frame("Group" = data$group, "y" = data$y)


cols = c("Thinned DDP"="deeppink4", "no-pool"="forestgreen")

ggplot(data = density_CI_thinnedDDP, aes(x = Seq, y = mean)) + 
  geom_histogram(data = data, aes(x = y, y=..density..),
                 binwidth = 1.5, colour="gray30", fill="gray80",lwd=0.3, alpha = 0.3) +
  geom_ribbon(data = density_CI_nopool, aes(ymin = lower, ymax = upper), 
              alpha=0.4, col = "forestgreen", fill = "forestgreen", lwd= 0.3) +
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha=0.4, col = "maroon", fill = "maroon", lwd= 0.3) +
  geom_line(aes(col = "Thinned DDP"), lwd = 1) +
  geom_line(data = density_CI_nopool, aes(x = Seq, y=mean, col="no-pool"), lwd = 1 ) +
  geom_line(aes(y = true), lty = "dashed") + 
  scale_colour_manual(name="Model", values = cols) + 
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
  # scale_y_continuous(breaks = c(0.0, 0.1, 0.2,0.3),  limits = c(0,0.3)) +
  xlab("y")  + ylab("Density")+
  facet_wrap( ~reorder(Group, sort(as.numeric(Group))), ncol = 2, dir="h")

ggsave("../04_Simulation_commonatoms/06_01_densitynopool_10G_minSS.pdf", width = 8, height = 10)




cols2 = c("Thinned DDP"="deeppink4",  "pool" = "royalblue3")
ggplot(data = density_CI_thinnedDDP, aes(x = Seq, y = mean)) + 
  geom_histogram(data = data, aes(x = y, y=..density..),
                 binwidth = 1.5, colour="gray30", fill="gray80",lwd=0.3, alpha = 0.3) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.4, 
              col = "maroon", fill = "maroon", lwd= 0.4) +
  geom_ribbon(data = density_CI_pool, aes(ymin = lower, ymax = upper), alpha=0.3, 
              col = "royalblue3", fill = "royalblue3", lwd= 0.4) +
  geom_line(aes(col = "Thinned DDP")) +
  geom_line(aes(y = true), lty = "dashed") + 
  geom_line(data = density_CI_pool, aes(x = Seq, y=mean, col="pool") ) +
  scale_colour_manual(name="Model", values = cols2) + 
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
  xlab("y")  + ylab("Density")+
  facet_wrap( ~reorder(Group, sort(as.numeric(Group))), ncol = 2, dir="h")

ggsave("../04_Simulation_commonatoms/06_01_densitypool_10G_minSS.pdf", width = 8, height = 10)



#-----------# #-----------# #-----------# #-----------# #-----------# #-----------# #-----------# 




repl = 10
i = 2
j = 4

name_file_open = paste0("../04_Simulation_commonatoms/Data/data_", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
load(file = name_file_open)

nameopen = paste0("../04_Simulation_commonatoms/Results/density_CI_thinnedDDP", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
load(nameopen)

nameopen = paste0("../04_Simulation_commonatoms/Results/density_CI_pool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
load(nameopen)

nameopen = paste0("../04_Simulation_commonatoms/Results/density_CI_nopool", n_groups[i], "groups_", sum(n_groups[i]/2*ssg*j), "n_", repl,".Rdata")
load(nameopen)

data = data.frame("Group" = data$group, "y" = data$y)


cols = c("Thinned DDP"="deeppink4", "no-pool"="forestgreen")

ggplot(data = density_CI_thinnedDDP, aes(x = Seq, y = mean)) + 
  geom_histogram(data = data, aes(x = y, y=..density..),
                 binwidth = 1.5, colour="gray30", fill="gray80",lwd=0.3, alpha = 0.3) +
  geom_ribbon(data = density_CI_nopool, aes(ymin = lower, ymax = upper), 
              alpha=0.4, col = "forestgreen", fill = "forestgreen", lwd= 0.3) +
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha=0.4, col = "maroon", fill = "maroon", lwd= 0.3) +
  geom_line(aes(col = "Thinned DDP"), lwd = 1) +
  geom_line(data = density_CI_nopool, aes(x = Seq, y=mean, col="no-pool"), lwd = 1 ) +
  geom_line(aes(y = true), lty = "dashed") + 
  scale_colour_manual(name="Model", values = cols) + 
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
  # scale_y_continuous(breaks = c(0.0, 0.1, 0.2,0.3),  limits = c(0,0.3)) +
  xlab("y")  + ylab("Density")+
  facet_wrap( ~reorder(Group, sort(as.numeric(Group))), ncol = 2, dir="h")

ggsave("../04_Simulation_commonatoms/06_01_densitynopool_10G_maxSS.pdf", width = 8, height = 10)




cols2 = c("Thinned DDP"="deeppink4",  "pool" = "royalblue3")
ggplot(data = density_CI_thinnedDDP, aes(x = Seq, y = mean)) + 
  geom_histogram(data = data, aes(x = y, y=..density..),
                 binwidth = 1.5, colour="gray30", fill="gray80",lwd=0.3, alpha = 0.3) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.4, 
              col = "maroon", fill = "maroon", lwd= 0.4) +
  geom_ribbon(data = density_CI_pool, aes(ymin = lower, ymax = upper), alpha=0.3, 
              col = "royalblue3", fill = "royalblue3", lwd= 0.4) +
  geom_line(aes(col = "Thinned DDP")) +
  geom_line(aes(y = true), lty = "dashed") + 
  geom_line(data = density_CI_pool, aes(x = Seq, y=mean, col="pool") ) +
  scale_colour_manual(name="Model", values = cols2) + 
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
  xlab("y")  + ylab("Density")+
  facet_wrap( ~reorder(Group, sort(as.numeric(Group))), ncol = 2, dir="h")

ggsave("../04_Simulation_commonatoms/06_01_densitypool_10G_maxSS.pdf", width = 8, height = 10)










