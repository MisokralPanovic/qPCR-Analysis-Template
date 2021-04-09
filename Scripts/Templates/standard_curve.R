# .libPaths('C:/r_packages')
library(tidyverse)
library(ggplot2)
library(scales)

###########
# change 'bIFIT2', 'PS22' and 'scb2_2'
###########

# data prep ----
standards_data <- read.csv(
  paste('Data/', 
        
        ############
        'standard_curves_data', 
        ############
        
        '.csv', 
        sep = '')
  )
      
scb2_2 <- standards_data %>%
  
  ###########
  filter(Target == 'bIFIT2',
         Primer_set == 'PS22'
         )
  ###########

scb2_2

model_lmscb2_2 <- lm(
  log10(Copy_number)~Ct, 
  data = scb2_2
  )
plot(model_lmscb2_2)

# put the -# into efficency equation
lm(Ct~log10(Copy_number), 
   data = scb2_2)

###########
rate_scb2_2 <- -3.3
###########

efficiency_scb2_2 <- paste(
  round(
    (10^(-1/ rate_scb2_2 ) -1)*100, 
    digits = 2), 
  '% Amplification Efficiency', 
  sep = ''
  )

model.predict_scb2_2 <- 10^predict(
  model_lmscb2_2, 
  interval = 'prediction'
  )
data.combined.predict_scb2_2 <- cbind(
  scb2_2,
  model.predict_scb2_2
  )

# data plotting ----

###########
x_annotation_position <- 1000000
y_annotation_position <- 37.5
top_range <- 40
plot_title <- 'bIFIT2 - PS22 Standard Curve'
y_axis_title <- 'Cycle Threshold'
x_axis_title <- 'Copy Number'
###########

plot_scb2_2 <- ggplot(
  data = data.combined.predict_scb2_2,
  mapping = aes(x = Copy_number,
                y = Ct)) +
  geom_point() +
  stat_smooth(method = lm) +
  labs(
    title = plot_title,
    y = y_axis_title,
    x = x_axis_title
  ) +
  scale_y_continuous(
    breaks = seq(from = 0, 
                 to = 40, 
                 by = 10),
    limits = c(0,top_range)
  ) +
  scale_x_log10(
    labels = trans_format("log10", 
                          math_format(10^.x)),
    breaks = trans_breaks("log10", 
                          function(x) 10^x)
  ) +
  annotation_logticks(sides='b') +
  theme(
    plot.title = element_text(
      size = 20, 
      face = 'bold', 
      margin = margin(10, 0, 10, 0), 
      hjust = 0.5
    ),
    axis.text.y = element_text(
      angle = 0, 
      size = 12, 
      vjust = 0.5),
    axis.text.x.bottom = element_text(
      angle = 0, 
      size = 12, 
      vjust = 0.5),
    axis.title.x = element_text(
      size = 15, 
      face = 'bold', 
      vjust = -0.5, 
      margin = margin(0, 10, 0, 0)),
    axis.title.y = element_text(
      size = 15, 
      face = 'bold', 
      vjust = -0.5, 
      margin = margin(0, 10, 0, 0)),
    aspect.ratio = 1/2
  ) +
  annotate('text',
           y = y_annotation_position, 
           x = x_annotation_position, 
           label = efficiency_scb2_2, 
           size = 5)

plot_scb2_2

# saving plot ------------------------------

###########
dpi <- 600
width <- 16
height <- 10
file_name <- 'scb2_2'
###########

ggsave(filename = paste(
  file_name, '.svg', sep = ''), 
       plot = plot_scb2_2, 
       device = 'svg', 
       path = 'Figures', 
       dpi = dpi, 
       height = height, 
       width = width, 
       units = 'cm')
ggsave(filename = paste(
  file_name, '.png', sep = ''), 
       plot = plot_scb2_2, 
       device = 'png', 
       path = 'Figures', 
       dpi = dpi, 
       height = height, 
       width = width, 
       units = 'cm')
