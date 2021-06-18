# .libPaths('C:/r_packages')
library(tidyverse)
library(ggsignif)
library(data.table)
library(ggsignif)
library(scales)
library(ggplot2)
library(ggpubr)

###########
# change 'hRSV N', ' 'mdbk_brsvn_48h'
###########

###########
list_of_conditions <- c(
  'Mock',
  'cbRSV 0.1-48',
  'cbRSV 1-48',
  'cbRSV 2-48'
)
list_of_colours <- c(
  '#999999', 
  '#e200e2', 
  '#a800a8', 
  '#6d006d'
)
file_emmit <- c(
  T,T,T,
  T,T,T,
  T,T,T,
  T,T,T
)
replicates <- 3
###########

# data prep ----
gene_of_interest_data <- read.csv(
  paste('Data/', 
        
        ############
        'ct_data', 
        ############
        
        '.csv', 
        sep = '')
)
housekeeping_control <- gene_of_interest_data %>% 
  
  ###########
filter(Target == 'bGAPDH',
       Time_point == '48h',
       Additional_info != 'bMx1',
       Condition %in% list_of_conditions
)
###########

housekeeping_control <- aggregate(
  housekeeping_control[1],
  list(housekeeping_control$Condition),
  mean) %>% 
  arrange(match(
    Group.1, 
    list_of_conditions))
housekeeping_control <- rep(
  housekeeping_control$Ct, 
  each=replicates)
housekeeping_control_vector <- housekeeping_control[file_emmit]

mdbk_brsvn_48h <- gene_of_interest_data %>%
  
  ###########
filter(
  Target == 'bRSV_N',
  Time_point == '48h',
  Condition %in% list_of_conditions
) %>% 
  ###########

arrange(match(
  Condition, 
  list_of_conditions)) %>%
  mutate(
    Control_ct_mean = housekeeping_control_vector,
    dCt = Ct - Control_ct_mean,
    Control_dct_mean = mean(
      dCt[Condition == list_of_conditions[1]],
      na.rm = T),
    ddCt = dCt - Control_dct_mean,
    log2_ddct = 2^ (- ddCt),
    mock_mean_log = mean(
      log2_ddct[Condition == list_of_conditions[1]],
      na.rm = T),
    Value_norm = log2_ddct / mock_mean_log
  )

mdbk_brsvn_48h

# data analysis ----
boxplot(Value_norm~Condition, mdbk_brsvn_48h)
plot(lm(Value_norm~Condition, mdbk_brsvn_48h))

# test normality
shapiro.test(mdbk_brsvn_48h$Value_norm[1:3])
shapiro.test(mdbk_brsvn_48h$Value_norm[4:6])
shapiro.test(mdbk_brsvn_48h$Value_norm[7:9])
shapiro.test(mdbk_brsvn_48h$Value_norm[10:12])

plot(residuals(
  lm(Value_norm~Condition, 
     mdbk_brsvn_48h)))
shapiro.test(residuals(
  lm(Value_norm~Condition, 
     mdbk_brsvn_48h)))

### normal distribution

# if distribution normal
bartlett.test(Value_norm~Condition, mdbk_brsvn_48h)

### non equal variance

# for normal distribution but non equal variance - multiple comparison
library(rstatix)
library(dplyr)
oneway.test(Value_norm~Condition, 
            mdbk_brsvn_48h, 
            var.equal = F)
mdbk_brsvn_48h %>% games_howell_test(Value_norm~Condition)

# plot data ----
boxplot(Value_norm~Condition, mdbk_brsvn_48h)
mdbk_brsvn_48h
list_of_conditions

###########
p_val <- c(
  0.076,
  0.137,
  0.032
)
range_y <- 9
plot_title <- 'MDBK - 48h - bRSV N'
y_axis_title <- 'Relative mRNA Levels'
###########

textsize_values <- c()

for (value in p_val) {
  if (value > 0.05) {
    textsize_values <- append(
      textsize_values, 4)
  } else {
    textsize_values <- append(
      textsize_values, 5)
  }
}

plot_mdbk_brsvn_48h <- ggplot(
  mdbk_brsvn_48h, 
  aes(
    x = Condition, 
    y = Value_norm,
    fill = Condition)) +
  geom_violin(
    trim = FALSE,
    alpha = 0.5,
    scale = 'width',
    adjust = 0.7) +
  stat_summary(
    fun.data = mean_se, 
    fun.args = list(mult=1), 
    geom = "pointrange", 
    color = "black",
    show.legend = F) +
  scale_x_discrete(
    limits = list_of_conditions) +
  scale_fill_manual(
    breaks = list_of_conditions,
    values = list_of_colours) +
  theme(
    plot.title = element_text(
      size = 20, 
      face = 'bold', 
      margin = margin(10, 0, 10, 0), 
      hjust = 0.5
    ),
    legend.text = element_text(
      size=15),  
    legend.title = element_blank(),
    axis.text.y = element_text(
      angle=0, 
      size=12, 
      vjust=0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(
      size = 15, 
      face='bold', 
      vjust=-0.5, 
      margin = margin(0, 10, 0, 0)),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    aspect.ratio = 2/1
  ) +
  labs(
    title = plot_title,
    y = y_axis_title,
    x = NULL
  ) +
  scale_y_log10(labels = trans_format("log10", 
                                      math_format(10^.x)),
                breaks = trans_breaks("log10", 
                                      function(x) 10^x)
  ) +
  annotation_logticks(sides='l') +
  geom_signif(
    comparisons = list(c(
      list_of_conditions[1], 
      list_of_conditions[2])), 
    annotation = p_val[1], 
    y_position = 0.93*range_y - 2*(range_y*0.075), 
    tip_length = 0, 
    vjust= -0.2, 
    size = 0.7,
    textsize = textsize_values[1]) +
  
  geom_signif(
    comparisons = list(c(
      list_of_conditions[1], 
      list_of_conditions[3])), 
    annotation = p_val[2], 
    y_position = 0.93*range_y - 1*(range_y*0.075), 
    tip_length = 0, 
    vjust= -0.2, 
    size = 0.7, 
    textsize = textsize_values[2]) +
  
  geom_signif(
    comparisons = list(c(list_of_conditions[1], 
                         list_of_conditions[4])), 
    annotation = p_val[3], 
    y_position = 0.93*range_y - 0*(range_y*0.075), 
    tip_length = 0, 
    vjust= -0.2, 
    size = 0.7, 
    textsize = textsize_values[3])

plot_mdbk_brsvn_48h

# saving data and plot ------------------------------

###########
dpi <- 600
width <- 16
height <- 20
file_name <- 'mdbk_brsvn_48h'
###########

ggsave(filename = paste(
  file_name, '.svg', sep = ''), 
  plot = plot_mdbk_brsvn_48h, 
  device = 'svg', 
  path = 'Figures', 
  dpi = dpi, 
  height = height, 
  width = width, 
  units = 'cm')
ggsave(filename = paste(
  file_name, '.png', sep = ''), 
  plot = plot_mdbk_brsvn_48h, 
  device = 'png', 
  path = 'Figures', 
  dpi = dpi, 
  height = height, 
  width = width, 
  units = 'cm')

fwrite(mdbk_brsvn_48h, 
       paste('Adjusted-Data/', 
             file_name, 
             '.csv', 
             sep = ''))
