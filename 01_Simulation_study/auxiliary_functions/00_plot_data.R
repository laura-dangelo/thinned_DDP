library(ggplot2)
library(viridis)
library(ggpubr)



plotdata = function(data) {
  df = data.frame(x = c(data$y),
                  y = data$group/70,
                  group = as.factor(c(data$group)),
                  obs_cl = as.factor(c(data$cl)))
  p1 = ggplot(df, aes(x=x)) + 
    geom_density(aes( color=group, fill = group), linewidth = 0.9, alpha=.4) + 
    scale_color_manual(values = plasma(length(unique(data$g)))) +
    scale_fill_manual(values = plasma(length(unique(data$g)))) +
    xlim(range(df$x) + c(-5,5)) +
    theme_bw() + 
    theme(panel.grid.minor = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) 
  
  p2 = ggplot(df, aes(x=x, y=y)) + 
    geom_point(col="black", cex=1,  ) +
    geom_point(aes( color=group), cex=0.8  ) +
    scale_color_manual(values = plasma(length(unique(data$g)))) +
    scale_fill_manual(values = plasma(length(unique(data$g)))) +
    xlim(range(df$x) + c(-5,5)) +
    ylim(range(df$y) ) +
    theme_bw() + 
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank())+
    labs(y = "")+
    guides(colour = guide_legend(nrow = 1))
  
  ggarrange(p1, p2, ncol=1, nrow=2, 
            align = "v",
            heights = c(2.5, 1),
            common.legend = TRUE, legend="bottom")
}
