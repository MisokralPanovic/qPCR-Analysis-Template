# .libPaths('C:/r_packages')
library(tidyverse)
library(scales)
library(data.table)

figure_theme <- theme(
  plot.title = element_text(
    size = 15, 
    face = 'bold', 
    margin = margin(8, 0, 8, 0), 
    hjust = 0.5
  ),
  axis.text.y = element_text(
    angle = 0, 
    size = 9, 
    vjust = 0.5),
  axis.text.x.bottom = element_text(
    angle = 0, 
    size = 9, 
    vjust = 0.5),
  axis.title.x = element_text(
    size = 12, 
    face = 'bold', 
    vjust = -0.5, 
    margin = margin(0, 8, 0, 0)),
  axis.title.y = element_text(
    size = 12, 
    face='bold', 
    vjust=-0.5, 
    margin = margin(0, 8, 0, 0)),
  aspect.ratio = 1/2
)

###########
# change 'bIFIT1', 'PS11' and 'scb1'
###########

# data prep ----
scb1 <- fread(
  paste('Data/', 
        
        ############
        'standard_curves_data', 
        ############
        
        '.csv', 
        sep = '')
  ) %>%
    
  ###########
  filter(Target == 'bIFIT1',
         Primer_set == 'PS11'
         )
  ###########

scb1

model_lmscb1 <- lm(
  log10(Copy_number)~Ct, 
  data = scb1
  )
plot(model_lmscb1)

efficiency_scb1 <- paste(
  round(
    (10^(-1/ lm(Ct~log10(Copy_number), data = scb1)[[1]][2]) -1)*100, 
    digits = 2), 
  '% Amplification Efficiency', 
  sep = ''
  )

model.predict_scb1 <- 10^predict(
  model_lmscb1, 
  interval = 'prediction'
  )
data.combined.predict_scb1 <- cbind(
  scb1,
  model.predict_scb1
  )

# data plotting ----

###########
plot_title <- 'bIFIT1 - PS11 Standard Curve'
x_annotation_position <- 1000000
y_annotation_position <- 37.5
top_range <- 40
y_axis_title <- 'Cycle Threshold'
x_axis_title <- 'Copy Number'
###########

plot_scb1 <- ggplot(
  data = data.combined.predict_scb1,
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
                 to = top_range, 
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
  theme_bw() +
  figure_theme +
  annotate('text',
           y = y_annotation_position, 
           x = x_annotation_position, 
           label = efficiency_scb1, 
           size = 5)

plot_scb1

# saving plot ------------------------------

###########
dpi <- 600
width <- 16
height <- 20
file_name <- 'scb1'
###########

ggsave(filename = paste(
  file_name, '.svg', sep = ''), 
       plot = plot_scb1, 
       device = 'svg', 
       path = 'Figures', 
       dpi = dpi, 
       height = height, 
       width = width, 
       units = 'cm')
ggsave(filename = paste(
  file_name, '.png', sep = ''), 
       plot = plot_scb1, 
       device = 'png', 
       path = 'Figures', 
       dpi = dpi, 
       height = height, 
       width = width, 
       units = 'cm')
