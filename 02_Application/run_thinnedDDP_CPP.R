library(devtools)
# library(mclust)
library(salso)
# library(mcclust.ext)
library(ggplot2)
library(cowplot)

library(TeachingDemos)
document()

load("../02_Application_CPP/CPP.Rdata")
str(CPP)

summary(CPP$gest)
var(CPP$gest)

plot(density(CPP$gest))
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

nrep = 10000
burnin = (floor(nrep/2))


# run_gibbs_CPP = Bthin::blocked_gibbs_BernThin_burnin_prior(nrep = nrep, 
#                                                            burnin = burnin, 
#                                                            y = CPP$gest,
#                                                            group = CPP$hosp-1, 
#                                                            trunc = trunc,
#                                                            m0 = mean(CPP$gest), tau0 = tau0, 
#                                                            gamma0 = gam0, lambda0 = lam0, 
#                                                            alpha = 1, 
#                                                            a_beta = 3, b_beta = 3,
#                                                            mu_start = mu_start,
#                                                            sigma2_start = sigma2_start,
#                                                            cl_start = cl_start-1,
#                                                            progressbar = T)
namesave = paste0("../02_Application_CPP/run_CPP.Rdata")
# save(run_gibbs_CPP, file = namesave)
load(namesave)
str(run_gibbs_CPP)

time_CPP = run_gibbs_CPP$time
namesave = paste0("../02_Application_CPP/time_CPP.Rdata")
# save(time_CPP, file = namesave)
load(namesave)

#-----# Check convergence with traceplots #-----# 
plot(run_gibbs_CPP$mu[,1], type="l" , ylim=c(30,45))
for(j in 2:4) { lines(run_gibbs_CPP$mu[,j], col=j) }

plot(run_gibbs_CPP$sigma2[,1], type="l")

nrep = nrow(run_gibbs_CPP$mu)

i=1
g=2
plot(unlist(sapply(1:nrep, function(x) run_gibbs_CPP$mu[x,run_gibbs_CPP$cl[x,i]+1])), type="l")
lines(1:nrep,
      cumsum(unlist(sapply(1:nrep, function(x) run_gibbs_CPP$mu[x,run_gibbs_CPP$cl[x,i]+1])))/1:nrep, col=2)

plot(unlist(sapply(1:nrep, function(x) run_gibbs_CPP$sigma2[x,run_gibbs_CPP$cl[x,i]+1])), type="l")
lines(1:nrep,
      cumsum(unlist(sapply(1:nrep, function(x) run_gibbs_CPP$sigma2[x,run_gibbs_CPP$cl[x,i]+1])))/1:nrep, col=2)

plot(1:nrep, round(sapply(1:nrep, function(iter) sum(run_gibbs_CPP$pi[,g,iter])),6), type="l")



## label of the first occupied cluster
plot(1:nrep, apply(run_gibbs_CPP$cl, 1, min), type="l")
## label of the last occupied cluser
plot(1:nrep, apply(run_gibbs_CPP$cl, 1, max), type="l")




#-----# Estimate partition of observations #-----# 
# clus_mat = salso::psm((run_gibbs_CPP$cl+1), nCores = 3 )


# cl_point_est_CPP = salso::salso(run_gibbs_CPP$cl+1, maxNClusters = 10)
# length(unique(cl_point_est_CPP))
# plot(CPP$gest, CPP$hosp, col = cl_point_est_CPP)


namesave = paste0("../02_Application_CPP/est_cl_CPP.Rdata")
# save(cl_point_est_CPP, file = namesave)
load(namesave)


tapply(CPP$gest, cl_point_est_CPP, mean)


#-----# Estimate group-specific partition of observations #-----# 

# cluster_est = c(cl_point_est_CPP)
# 
# cluster_est = numeric(length(CPP$hosp))
# for(g in 1:length(unique(CPP$hosp))) {
#   submat_cl = run_gibbs_CPP$cl[,CPP$hosp==g]+1
#   cluster_est_g = salso(submat_cl, loss = VI(), maxNClusters = 7)
#   mean_gest = tapply(CPP$gest[CPP$hosp==g], cluster_est_g, mean)
#   sort_gest = sort(mean_gest, index.return = T)
# 
#   for(id in 1:length(unique(sort_gest$ix))) {
#     cluster_est_g[cluster_est_g == sort_gest$ix[id]] = id-5
#   }
#   cluster_est[CPP$hosp==g] = cluster_est_g
# }
# cluster_est = cluster_est +5
plot(CPP$gest, CPP$hosp, col = cluster_est)

tapply(cluster_est, CPP$hosp, function(x) length(unique(x)))


#-----# #-----# Compute density estimate #-----# #-----#
seqq = seq(range(CPP$gest)[1]-2, range(CPP$gest)[2]+2, length.out = 300)
# density_est_CPP = Bthin::compute_density(seqq,
#                                          weight = run_gibbs_CPP$pi,
#                                          means = t(run_gibbs_CPP$mu),
#                                          variances = t(run_gibbs_CPP$sigma2)
# )
namesave = paste0("../02_Application_CPP/density_est_CPP.Rdata")
# save(density_est_CPP, file = namesave)
load(namesave)

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

namesave = paste0("../02_Application_CPP/density_CI_CPP.Rdata")
# save(CPPdensity_CI_thinnedDDP, file = namesave)
load(namesave)






#-----# #-----# Cluster of hospitals #-----# #-----#

iter_hosp = matrix(rep(1:12,5000), 5000, 12, byrow = T)
for(iter in 1:5000) {
  for(g1 in 1:12) {
    for(g2 in 1:g1) {
      cl1 = sort(unique(run_gibbs_CPP$cl[iter, run_gibbs_CPP$group==(g1-1)]))
      cl2 =  sort(unique(run_gibbs_CPP$cl[iter, run_gibbs_CPP$group==(g2-1)]))
      if(length(cl1) == length(cl2)) {
        if(prod(cl1==cl2)) { iter_hosp[iter,g1] = iter_hosp[iter,g2] }}
    }
  }
}

namesave = paste0("../02_Application_CPP/iter_hospitals_CPP.Rdata")
# save(iter_hosp, file = namesave)
load(namesave)

psm_hosp = salso::psm(iter_hosp, nCores = 3 )

hosp_point_est_CPP = salso::salso(iter_hosp, loss = VI(a=1.97), maxNClusters = 12)
length(unique(hosp_point_est_CPP))
hosp_point_est_CPP
which(hosp_point_est_CPP == 1)
which(hosp_point_est_CPP == 2)
which(hosp_point_est_CPP == 3)

# mat_hosp = matrix(0, 12, 12)
# for(iter in 1:5000) {
#   for(g1 in 1:12) {
#     for(g2 in 1:g1) {
#       cl1 = sort(unique(run_gibbs_CPP$cl[iter, run_gibbs_CPP$group==(g1-1)]))
#       cl2 =  sort(unique(run_gibbs_CPP$cl[iter, run_gibbs_CPP$group==(g2-1)]))
#       if(length(cl1) == length(cl2)) { 
#         mat_hosp[g1,g2] = mat_hosp[g1,g2] + prod(cl1==cl2)
#         mat_hosp[g2,g1] = mat_hosp[g1,g2] }
#     }
#   }
# }
# mat_hosp = mat_hosp/5000

namesave = paste0("../02_Application_CPP/mat_clus_hospitals_CPP.Rdata")
# save(mat_hosp, file = namesave)
load(namesave)

minVI(mat_hosp)$cl
which(minVI(mat_hosp)$cl == 2)
idclh = sort(minVI(mat_hosp)$cl, index.return=T)


# matrix TV distance
TV_distance = function(seq, dens1, dens2) {
  sum(.5 * abs(dens1-dens2))
}

mat_hospTV = matrix(0, 12, 12)
for(iter in 1:5000) {
  for(g1 in 1:12) {
    for(g2 in 1:g1) {
      mat_hospTV[g1,g2] =  mat_hospTV[g1,g2] + TV_distance(seqq, density_est_CPP$density_mcmc[,g1,iter],
                                                           density_est_CPP$density_mcmc[,g2,iter])
      mat_hospTV[g2,g1] = mat_hospTV[g1,g2]
    }
  }
}

namesave = paste0("../02_Application_CPP/mat_TV_hosp_CPP.Rdata")
mat_hospTV = mat_hospTV / 5000
# save(mat_hospTV, file = namesave)
load(namesave)

isSymmetric(mat_hospTV)


row.names(mat_hospTV) = factor(1:12)
colnames(mat_hospTV) = factor(1:12)

library(pheatmap)
library("RColorBrewer")

pheatmap(mat_hospTV, treeheight_row = 0, treeheight_col = 0, angle_col = 0, 
         color = brewer.pal(n = 9, name = "YlGnBu")[9:1] )




#-----# #-----# Plots #-----# #-----#
# cl_point_est_CPP = cluster_est
# 
# data = data.frame("Group" = CPP$hosp, "y" = CPP$gest, 
#                   "cl" = as.factor(cl_point_est_CPP), "jit" = runif(length(cl_point_est_CPP), -0.3,0.3))
# 
# saveRDS(data, file = "../02_Application_CPP/dataset_plot_densities_cpp.RDS")

data = readRDS("../02_Application_CPP/dataset_plot_densities_cpp.RDS")
namesave = paste0("../02_Application_CPP/density_CI_CPP.Rdata")
load(namesave)


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
    scale_color_manual(values = c("1" = "turquoise4","2" = "orange2","3" = "#0520ca"))+
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
    scale_color_manual(values = c("1" = "turquoise4","2" = "orange2","3" = "#0520ca"))+
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
    scale_color_manual(values = c("1" = "turquoise4","2" = "orange2","3" = "#0520ca"))+
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


ggsave("../02_Application_CPP/CPP_density_clusters2.pdf", width = 8, height = 9)









# 
# pp = list()
# pp2 = list()

for(g in c(7,9,11)){
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
      legend.background = element_rect(fill='transparent', color = 'white'), #transparent legend bg
      legend.box.background = element_rect(fill='transparent'), #transparent legend panel
      legend.text = element_text(size=10),
      strip.text = element_text(size=10),
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
    scale_color_manual(values = c("1" = "turquoise4","2" = "orange2","3" = "#0520ca"))+
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
      legend.background = element_rect(fill='transparent', color = 'white'), #transparent legend bg
      legend.box.background = element_rect(fill='transparent'), #transparent legend panel
      legend.text = element_text(size=10),
      strip.text = element_text(size=10),
      strip.background = element_rect( fill=NA, color="gray" ),
      plot.title = element_text(hjust = 0.5)
    )+
    xlab("  ")  + ylab("  ") +
    scale_x_continuous(breaks = 25+c(0,5,10,15,20),  limits = c(25,48)) +
    scale_y_continuous(labels = c(), limits = c(-0.35,0.35))  
  
}

for(g in c(11)){
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
      legend.background = element_rect(fill='transparent', color = 'white'), #transparent legend bg
      legend.box.background = element_rect(fill='transparent'), #transparent legend panel
      legend.text = element_text(size=10),
      strip.text = element_text(size=10),
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
    scale_color_manual(values = c("1" = "turquoise4","2" = "orange2","3" = "#0520ca"))+
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
      legend.background = element_rect(fill='transparent', color = 'white'), #transparent legend bg
      legend.box.background = element_rect(fill='transparent'), #transparent legend panel
      legend.text = element_text(size=10),
      strip.text = element_text(size=10),
      strip.background = element_rect( fill=NA, color="gray" ),
      plot.title = element_text(hjust = 0.5)
    )+
    xlab("Gestational age")  + ylab("  ") +
    scale_x_continuous(breaks = 25+c(0,5,10,15,20),  limits = c(25,48)) +
    scale_y_continuous(labels = c(), limits = c(-0.35,0.35))  
  
}

for(g in c(8,10,12)){
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
      legend.background = element_rect(fill='transparent', color = 'white'), #transparent legend bg
      legend.box.background = element_rect(fill='transparent'), #transparent legend panel
      legend.text = element_text(size=10),
      strip.text = element_text(size=10),
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
    scale_color_manual(values = c("1" = "turquoise4","2" = "orange2","3" = "#0520ca"))+
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
      legend.background = element_rect(fill='transparent', color = 'white'), #transparent legend bg
      legend.box.background = element_rect(fill='transparent'), #transparent legend panel
      legend.text = element_text(size=10),
      strip.text = element_text(size=10),
      strip.background = element_rect( fill=NA, color="gray" ),
      plot.title = element_text(hjust = 0.5)
    )+
    xlab("  ")  + ylab("  ") +
    scale_x_continuous(breaks = 25+c(0,5,10,15,20), 
                       limits = c(25,48)) +
    scale_y_continuous(labels = c(), limits = c(-0.35,0.35))  
  
}

for(g in c(12)){
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
      legend.background = element_rect(fill='transparent', color = 'white'), #transparent legend bg
      legend.box.background = element_rect(fill='transparent'), #transparent legend panel
      legend.text = element_text(size=10),
      strip.text = element_text(size=10),
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
    scale_color_manual(values = c("1" = "turquoise4","2" = "orange2","3" = "#0520ca"))+
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
      legend.background = element_rect(fill='transparent', color = 'white'), #transparent legend bg
      legend.box.background = element_rect(fill='transparent'), #transparent legend panel
      legend.text = element_text(size=10),
      strip.text = element_text(size=10),
      strip.background = element_rect( fill=NA, color="gray" ),
      plot.title = element_text(hjust = 0.5)
    )+
    xlab("Gestational age")  + ylab("  ") +
    scale_x_continuous(breaks = 25+c(0,5,10,15,20), 
                       limits = c(25,48)) +
    scale_y_continuous(labels = c(), limits = c(-0.35,0.35))  
  
}



plot_grid(pp[[7]], pp[[8]], 
          pp2[[7]], pp2[[8]], 
          pp[[9]], pp[[10]],
          pp2[[9]], pp2[[10]],
          pp[[11]], pp[[12]],
          pp2[[11]], pp2[[12]],
          ncol  = 2,
          align = 'v',
          axis  = 'tb',
          labels = c(),
          rel_heights = c(2, 0.5)
)


ggsave("../02_Application_CPP/CPP_density_clusters2.pdf", width = 8, height = 9)






plot_grid(pp[[1]], pp[[2]], 
          pp2[[1]], pp2[[2]], 
          pp[[3]], pp[[4]],
          pp2[[3]], pp2[[4]],
          pp[[5]], pp[[6]],
          pp2[[5]], pp2[[6]],
  pp[[7]], pp[[8]], 
          pp2[[7]], pp2[[8]], 
          pp[[9]], pp[[10]],
          pp2[[9]], pp2[[10]],
          pp[[11]], pp[[12]],
          pp2[[11]], pp2[[12]],
          ncol  = 2,
          align = 'v',
          axis  = 'tb',
          labels = c(),
          rel_heights = c(2, 0.5)
)

ggsave("../02_Application_CPP/CPP_density_clusters3.pdf", width = 8, height = 12)
