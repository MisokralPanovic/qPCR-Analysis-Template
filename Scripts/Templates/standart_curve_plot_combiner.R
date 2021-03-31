# .libPaths('C:/r_packages')
library(tidyverse)
library(ggpubr)

###########
# change 'standard_curves', '_x', plot names
###########

##########
dpi <- 600
width <- 35
height <- 20
file_name <- 'standard_curves'
##########

combined_plot_x <- ggarrange(plot_standard1, plot_standard2, plot_standard3, plot_standard5, # what to plot
                             labels = c('A', 'B', 'C', 'D'), # labels
                             ncol = 2, nrow = 2, # dimension of new plot
                             common.legend = T, legend = "bottom") # shared legend and legend position

combined_plot_x

ggsave(filename = paste(file_name, '.svg', sep = ''),
       plot = combined_plot, 
       device = 'svg', 
       path = 'Figures/Combined-Plots', 
       dpi = dpi, 
       height = height, 
       width = width, 
       units = 'cm')
ggsave(filename = paste(file_name, '.png', sep = ''),
       plot = combined_plot, 
       device = 'png', 
       path = 'Figures/Combined-Plots', 
       dpi = dpi, 
       height = height, 
       width = width, 
       units = 'cm')
