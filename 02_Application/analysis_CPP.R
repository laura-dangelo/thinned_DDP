
library(salso)
library(ggplot2)
library(cowplot)
library(TeachingDemos)


load("02_Application/CPP.Rdata")
str(CPP)

summary(CPP$gest)
var(CPP$gest)

CPP$hosp = as.numeric(CPP$hosp)



trunc = 300
set.seed(2)
cl_start = kmeans(CPP$gest, centers = 2)
cl_start = cl_start$cluster


mu_start = rep(0, trunc)
mu_start[1:length(unique(cl_start))] = tapply(CPP$gest, cl_start, mean)
sigma2_start = rep(var(CPP$gest)/3, trunc)
sigma2_start[1:length(unique(cl_start))] = tapply(CPP$gest, cl_start, var)
sigma2_start[1:length(unique(cl_start))] = sigma2_start[1:length(unique(cl_start))]
sigma2_start[is.na(sigma2_start)] = var(CPP$gest)/3

tau0 = 0.01
gam0 = 2.5
lam0 = (gam0-1)

nrep = 110000
burnin = 10000

if(!file.exists("02_Application/results/run_CPP.RDS")) {
  run_gibbs_CPP = thinnedDDP::sampler_thinnedDDP(nrep = nrep,
                                                 burnin = burnin,
                                                 thinning_factor = 40,
                                                 y = CPP$gest,
                                                 group = CPP$hosp-1,
                                                 trunc = trunc,
                                                 m0 = mean(CPP$gest), tau0 = tau0,
                                                 gamma0 = gam0, lambda0 = lam0,
                                                 alpha = 1,
                                                 a_beta = 3, b_beta = 3,
                                                 mu_start = mu_start,
                                                 sigma2_start = sigma2_start,
                                                 cl_start = cl_start-1,
                                                 progressbar = T)
  saveRDS(run_gibbs_CPP, file = "02_Application/results/run_CPP.RDS")
} else {
  run_gibbs_CPP = readRDS(file = "02_Application/results/run_CPP.RDS")
}



str(run_gibbs_CPP)

time_CPP = run_gibbs_CPP$time


#-----# Check convergence with traceplots #-----# 
plot(run_gibbs_CPP$mu[,1], type="l" , ylim=c(30,45))
for(j in 2:4) { lines(run_gibbs_CPP$mu[,j], col=j) }



#-----# Estimate partition of observations #-----# 
if(!file.exists("02_Application/results/est_cl_CPP.RDS")){
  cl_point_est_CPP = salso::salso((run_gibbs_CPP$cl+1), nCores = 3 )
  namesave = paste0("02_Application/results/est_cl_CPP.RDS")
  saveRDS(cl_point_est_CPP, file = namesave)
} else {
  cl_point_est_CPP = readRDS(file = "02_Application/results/est_cl_CPP.RDS")
}


length(unique(cl_point_est_CPP))
plot(CPP$gest, CPP$hosp, col = cl_point_est_CPP)


tapply(CPP$gest, cl_point_est_CPP, mean)




#-----# #-----# Compute density estimate #-----# #-----#
if(!file.exists("02_Application/results/density_est_CPP.RDS")){
  seqq = seq(range(CPP$gest)[1]-2, range(CPP$gest)[2]+2, length.out = 300)
  density_est_CPP = thinnedDDP::compute_density(seqq,
                                                weight = run_gibbs_CPP$pi,
                                                means = t(run_gibbs_CPP$mu),
                                                variances = t(run_gibbs_CPP$sigma2)
  )
  namesave = paste0("02_Application/results/density_est_CPP.RDS")
  saveRDS(density_est_CPP, file = namesave)
} else {
  seqq = seq(range(CPP$gest)[1]-2, range(CPP$gest)[2]+2, length.out = 300)
  density_est_CPP = readRDS(file = "02_Application/results/density_est_CPP.RDS")
}

namesave = paste0("02_Application/results/density_CI_CPP.RDS")
if(!file.exists(namesave)){
  seqq = density_est_CPP$seq
  CPPdensity_CI_thinnedDDP = matrix(rep(1,length(seqq)), length(seqq), 1)
  CPPdensity_CI_thinnedDDP = cbind(CPPdensity_CI_thinnedDDP, seqq)
  CPPdensity_CI_thinnedDDP = cbind(CPPdensity_CI_thinnedDDP, apply(density_est_CPP$density_mcmc[,1,], 1, function(x) emp.hpd(x)[[1]] ))
  CPPdensity_CI_thinnedDDP = cbind(CPPdensity_CI_thinnedDDP, apply(density_est_CPP$density_mcmc[,1,], 1, function(x) mean(x) ))
  CPPdensity_CI_thinnedDDP = cbind(CPPdensity_CI_thinnedDDP, apply(density_est_CPP$density_mcmc[,1,], 1, function(x) emp.hpd(x)[[2]] ))
  
  for(gg in 2:length(unique(CPP$hosp))) {
    tmp = matrix(rep(gg,length(seqq)), length(seqq), 1)
    tmp = cbind(tmp, seqq)
    tmp = cbind(tmp, apply(density_est_CPP$density_mcmc[,gg,], 1, function(x) emp.hpd(x)[[1]] ))
    tmp = cbind(tmp, apply(density_est_CPP$density_mcmc[,gg,], 1, function(x) mean(x) ))
    tmp = cbind(tmp, apply(density_est_CPP$density_mcmc[,gg,], 1, function(x) emp.hpd(x)[[2]] ))
    
    CPPdensity_CI_thinnedDDP = rbind(CPPdensity_CI_thinnedDDP, tmp)
  }
  rm(tmp)
  
  CPPdensity_CI_thinnedDDP = data.frame(CPPdensity_CI_thinnedDDP)
  colnames(CPPdensity_CI_thinnedDDP) = c("Group", "Seq", "lower", "mean", "upper")
  
  saveRDS(CPPdensity_CI_thinnedDDP, file = namesave)
} else {
  CPPdensity_CI_thinnedDDP = readRDS(namesave)
}


#-----# #-----# Cluster of hospitals #-----# #-----#

namemat = paste0("02_Application/results/mat_hosp.RDS")
if(!file.exists(namemat)){
  mat_hosp = matrix(0, 12, 12)
  for(iter in 1:7500) {
    for(g1 in 1:12) {
      for(g2 in 1:g1) {
        cl1 = sort(unique(run_gibbs_CPP$cl[iter, run_gibbs_CPP$group==(g1-1)]))
        cl2 =  sort(unique(run_gibbs_CPP$cl[iter, run_gibbs_CPP$group==(g2-1)]))
        if( length(c(setdiff(cl1,cl2),setdiff(cl2,cl1)))<3 ) {
          mat_hosp[g1,g2] = mat_hosp[g1,g2] + 1
          mat_hosp[g2,g1] = mat_hosp[g1,g2] }
      }
    }
  }
  mat_hosp = mat_hosp/7500
  saveRDS(mat_hosp, namemat)
} else {
  mat_hosp = readRDS(namemat)
}


library(mcclust.ext)
minVI(mat_hosp)$cl
which(minVI(mat_hosp)$cl == 1)
which(minVI(mat_hosp)$cl == 2)


# matrix TV distance
TV_distance = function(seq, dens1, dens2) {
  sum(.5 * abs(dens1-dens2))
}

mat_hospTV = matrix(0, 12, 12)
for(iter in 1:7500) {
  for(g1 in 1:12) {
    for(g2 in 1:g1) {
      mat_hospTV[g1,g2] =  mat_hospTV[g1,g2] + TV_distance(seqq, density_est_CPP$density_mcmc[,g1,iter],
                                                           density_est_CPP$density_mcmc[,g2,iter])
      mat_hospTV[g2,g1] = mat_hospTV[g1,g2]
    }
  }
}
mat_hospTV = mat_hospTV/7500

row.names(mat_hospTV) = factor(1:12)
colnames(mat_hospTV) = factor(1:12)

library(pheatmap)
library("RColorBrewer")

hea = pheatmap(mat_hospTV, treeheight_row = 0, treeheight_col = 0, angle_col = 0, 
         color = brewer.pal(n = 9, name = "YlGnBu")[9:1] )


ggsave("02_Application/output_images/CPP_heatmap_clusters_hospitals.pdf", plot=hea, width = 4.5, height = 4)



#-----# #-----# Plots #-----# #-----#
data = data.frame("Group" = CPP$hosp, "y" = CPP$gest,
                  "cl" = as.factor(cl_point_est_CPP), "jit" = runif(length(cl_point_est_CPP), -0.3,0.3))
namesave = paste0("02_Application/results/density_CI_CPP.RDS")
CPPdensity_CI_thinnedDDP = readRDS(namesave)


pp = list()
pp2 = list()

for(g in c(1,4,7,10)){
  pp[[g]] = ggplot(data = CPPdensity_CI_thinnedDDP[CPPdensity_CI_thinnedDDP$Group==g,], 
                   aes(x = Seq, y = mean)) + 
    geom_histogram(data = data[data$Group==g,], aes(x = y, y=..density..),
                   binwidth = 1.5, colour="gray30", fill="gray80",lwd=0.3, alpha = 0.3) +
    geom_ribbon(aes(ymin = lower, ymax = upper), 
                alpha=0.2, col = "maroon", fill = "maroon", lwd= 0.3) +
    geom_line(aes(col = "Thinned DDP"), lwd = 1, col = "deeppink4") +
    theme_minimal() +
    theme(
      plot.margin = unit(c(1,0.1,-0.5,0.1), "pt"),
      panel.background = element_rect(fill='transparent', color=NA), #transparent panel bg
      plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
      # panel.grid.major = element_blank(), #remove major gridlines
      panel.grid.minor = element_blank(), #remove minor gridlines
      panel.border = element_rect(color = "darkgray ", fill=NA),
      # axis.line.y.left = element_line(color="gray"),
      axis.line.x.bottom = element_line(color="gray"),
      #
      legend.position = "bottom",
      legend.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend bg
      legend.box.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend panel
      legend.text = element_text(size=12),
      strip.text = element_text(size=12),
      text = element_text(size = 12),
      strip.background = element_rect( fill=NA, color="gray" ),
      plot.title = element_text(hjust = 0.5, size = 11.5)
    )+
    scale_y_continuous(breaks = c(0.0, 0.1, 0.2),  limits = c(0,0.28)) +
    scale_x_continuous(breaks = 25+c(0,5,10,15,20),  labels=c("","","","",""),
                       limits = c(25,48)) +
    xlab("  ")  + ylab("Density") +
    ggtitle(paste0("Hospital ", g))
  
  
  pp2[[g]] = ggplot(data[data$Group==g,], aes(x = y, y = jit, col=cl)) + 
    geom_point(size = 1.6, alpha = 0.5) +
    scale_color_manual(values = c("1" = "turquoise4","3" = "orange2","2" = "#0520ca"))+
    theme_minimal() +
    theme(
      plot.margin = unit(c(-20,1,0,0.1), "pt"),
      panel.background = element_rect(fill='transparent', color=NA), #transparent panel bg
      plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
      panel.grid.major.y = element_blank(), #remove major gridlines
      panel.grid.minor = element_blank(), #remove minor gridlines
      panel.border = element_rect(color = "darkgray ", fill=NA),
      # axis.line.y.left = element_line(color="gray"),
      axis.line.x.bottom = element_line(color="gray"),
      #
      legend.position = "none",
      legend.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend bg
      legend.box.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend panel
      legend.text = element_text(size=12),
      strip.text = element_text(size=12),
      text = element_text(size = 12),
      strip.background = element_rect( fill=NA, color="gray" ),
      plot.title = element_text(hjust = 0.5)
    )+
    xlab("  ")  + ylab("  ") +
    scale_x_continuous(breaks = 25+c(0,5,10,15,20),  limits = c(25,48)) +
    scale_y_continuous(labels = c(), limits = c(-0.35,0.35))  
  
}



for(g in c(2,3,5,6,8,9,11,12)){
  pp[[g]] = ggplot(data = CPPdensity_CI_thinnedDDP[CPPdensity_CI_thinnedDDP$Group==g,], 
                   aes(x = Seq, y = mean)) + 
    geom_histogram(data = data[data$Group==g,], aes(x = y, y=..density..),
                   binwidth = 1.5, colour="gray30", fill="gray80",lwd=0.3, alpha = 0.3) +
    geom_ribbon(aes(ymin = lower, ymax = upper), 
                alpha=0.2, col = "maroon", fill = "maroon", lwd= 0.3) +
    geom_line(aes(col = "Thinned DDP"), lwd = 1, col = "deeppink4") +
    theme_minimal() +
    theme(
      plot.margin = unit(c(1,1,-0.5,0.1), "pt"),
      panel.background = element_rect(fill='transparent', color=NA), #transparent panel bg
      plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
      # panel.grid.major = element_blank(), #remove major gridlines
      panel.grid.minor = element_blank(), #remove minor gridlines
      panel.border = element_rect(color = "darkgray ", fill=NA),
      # axis.line.y.left = element_line(color="gray"),
      axis.line.x.bottom = element_line(color="gray"),
      #
      legend.position = "bottom",
      legend.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend bg
      legend.box.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend panel
      legend.text = element_text(size=12),
      strip.text = element_text(size=12),
      text = element_text(size = 12),
      strip.background = element_rect( fill=NA, color="gray" ),
      plot.title = element_text(hjust = 0.5, size = 11.5)
    )+
    scale_y_continuous(breaks = c(0.0, 0.1, 0.2),  limits = c(0,0.28)) +
    scale_x_continuous(breaks = 25+c(0,5,10,15,20),  labels=c("","","","",""),
                       limits = c(25,48)) +
    xlab("  ")  + ylab("") +
    ggtitle(paste0("Hospital ", g))
  
  
  pp2[[g]] = ggplot(data[data$Group==g,], aes(x = y, y = jit, col=cl)) + 
    geom_point(size = 1.6, alpha = 0.5) +
    scale_color_manual(values = c("1" = "turquoise4","3" = "orange2","2" = "#0520ca"))+
    theme_minimal() +
    theme(
      plot.margin = unit(c(-20,1,0,0.1), "pt"),
      panel.background = element_rect(fill='transparent', color=NA), #transparent panel bg
      plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
      panel.grid.major.y = element_blank(), #remove major gridlines
      panel.grid.minor = element_blank(), #remove minor gridlines
      panel.border = element_rect(color = "darkgray ", fill=NA),
      # axis.line.y.left = element_line(color="gray"),
      axis.line.x.bottom = element_line(color="gray"),
      #
      legend.position = "none",
      legend.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend bg
      legend.box.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend panel
      legend.text = element_text(size=12),
      strip.text = element_text(size=12),
      text = element_text(size = 12),
      strip.background = element_rect( fill=NA, color="gray" ),
      plot.title = element_text(hjust = 0.5)
    )+
    xlab("  ")  + ylab("  ") +
    scale_x_continuous(breaks = 25+c(0,5,10,15,20), 
                       limits = c(25,48)) +
    scale_y_continuous(labels = c(), limits = c(-0.35,0.35))  
  
}

for(g in c(10,11,12)){

  pp2[[g]] = ggplot(data[data$Group==g,], aes(x = y, y = jit, col=cl)) +
    geom_point(size = 1.6, alpha = 0.5) +
    scale_color_manual(values = c("1" = "turquoise4","3" = "orange2","2" = "#0520ca"))+
    theme_minimal() +
    theme(
      plot.margin = unit(c(-20,1,0,0.1), "pt"),
      panel.background = element_rect(fill='transparent', color=NA), #transparent panel bg
      plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
      panel.grid.major.y = element_blank(), #remove major gridlines
      panel.grid.minor = element_blank(), #remove minor gridlines
      panel.border = element_rect(color = "darkgray ", fill=NA),
      # axis.line.y.left = element_line(color="gray"),
      axis.line.x.bottom = element_line(color="gray"),
      #
      legend.position = "none",
      legend.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend bg
      legend.box.background = element_rect(fill='transparent', color = 'transparent'), #transparent legend panel
      legend.text = element_text(size=12),
      strip.text = element_text(size=12),
      text = element_text(size = 12),
      strip.background = element_rect( fill=NA, color="gray" ),
      plot.title = element_text(hjust = 0.5)
    )+
    xlab("Gestational age")  + ylab("  ") +
    scale_x_continuous(breaks = 25+c(0,5,10,15,20),
                       limits = c(25,48)) +
    scale_y_continuous(labels = c(), limits = c(-0.35,0.35))

}



plot_grid(pp[[1]], pp[[2]], pp[[3]], 
          pp2[[1]], pp2[[2]], pp2[[3]],
          pp[[4]], pp[[5]],pp[[6]],
          pp2[[4]], pp2[[5]],pp2[[6]],
          pp[[7]], pp[[8]], pp[[9]],
          pp2[[7]], pp2[[8]], pp2[[9]],
          pp[[10]], pp[[11]], pp[[12]],
          pp2[[10]], pp2[[11]], pp2[[12]],
          ncol  = 3,
          align = 'v',
          axis  = 'tb',
          labels = c(),
          rel_heights = c(2, 0.5)
          )


ggsave("02_Application/output_images/CPP_density_clusters.pdf", width = 8, height = 9)




