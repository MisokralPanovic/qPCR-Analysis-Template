.libPaths('C:/r_packages')
library(tidyverse)
library(ggpubr)

###########
# change 'human_plots', '_x', plot names
###########

##########
dpi <- 600
width <- 45
height <- 30
file_name <- 'human_plots'
##########


combined_plot_x <- ggarrange(plot_a549_1, plot_a549_2, plot_a549_3, plot_a549_5, # what to plot
                           labels = c('A', 'B', 'C', 'D'), # labels
                           ncol = 4, nrow = 1, # dimension of new plot
                           common.legend = T, legend = "bottom") # shared legend and legend position

combined_plot_x

ggsave(filename = paste(file_name, '.png', sep = ''),
       plot = combined_plot_x, 
       device = 'png', 
       path = 'Figures/Combined Plots', 
       dpi = dpi, 
       height = height, 
       width = width, 
       units = 'cm')
