# .libPaths('C:/r_packages')
library(tidyverse)
library(ggsignif)
library(data.table)
library(scales)

figure_theme <-  theme(
  plot.title = element_text(
    size = 15, 
    face = 'bold', 
    margin = margin(8, 0, 8, 0), 
    hjust = 0.5
  ),
  legend.text = element_text(
    size=10),  
  legend.title = element_blank(),
  axis.text.y = element_text(
    angle=0, 
    size=9, 
    vjust=0.5),
  axis.title.x = element_blank(),
  axis.title.y = element_text(
    size = 12, 
    face='bold', 
    vjust=-0.5, 
    margin = margin(0, 8, 0, 0)),
  axis.text.x=element_blank(),
  axis.ticks.x=element_blank(),
  aspect.ratio = 2/1
)
###########
# change 'MDBK', 'bRSV_N', 'mdbk_bRSVN', 'bRSV N',
###########

###########
list_of_conditions <- c(
  'Mock',
  'bRSV dSH 1-24'
)
list_of_colours <- c(
  '#999999', 
  '#6d006d'
)
###########

# data prep -------------------------
mdbk_bRSVN <- fread(
  paste('Data/', 
        
        ############
        'ct_data', 
        ############
        
        '.csv', 
        sep = '')
  ) %>%

# all dataset normalisation - reverse values so mock low and infection high
  
  ###########
  filter(Target == 'bRSV_N',
         Cell_line == "MDBK",
         Condition %in% list_of_conditions
         ) %>%
  ###########

  arrange(match(
    Condition, 
    list_of_conditions)) %>%
  mutate(
    control_mean_ct = mean(
      Ct[Condition == list_of_conditions[1]], 
      na.rm = T),
    log2_dCt = 2^ (- (Ct - control_mean_ct)),
    control_mean_log = mean(
      log2_dCt[Condition == list_of_conditions[1]],
      na.rm = T),
    Value_norm = log2_dCt / control_mean_log
    )
                         
mdbk_bRSVN

# data analysis ----
boxplot(Value_norm~Condition, mdbk_bRSVN)
plot(lm(Value_norm~Condition, mdbk_bRSVN))

# test normality
shapiro.test(mdbk_bRSVN$Value_norm[1:3])
shapiro.test(mdbk_bRSVN$Value_norm[4:6])

plot(residuals(lm(Value_norm~Condition, mdbk_bRSVN)))
shapiro.test(residuals(lm(Value_norm~Condition, mdbk_bRSVN)))

### non normal distribution

# if non normal distribution
library(car)
leveneTest(Value_norm~Condition, 
           Value_norm~Condition, 
           mdbk_bRSVN)

### equal variance
# for non normal distribution but equal variance
t.test(Value_norm~Condition, 
       data=mdbk_bRSVN, 
       alternative='two.sided',
       var.equal=T)

# plot data ----
boxplot(Value_norm~Condition, 
        mdbk_bRSVN)
mdbk_bRSVN
list_of_conditions

###########
p_val <- c(
  0.0464754
  )
plot_title <- 'MDBK - bRSV N'
y_axis_title <- 'Relative mRNA Levels'
###########

textsize_values <- c()

for (value in p_val) {
  if (value > 0.05) {p
    textsize_values <- append(
      textsize_values, 3)
  } else {
    textsize_values <- append(
      textsize_values, 4)
  }
}

plot_mdbk_bRSVN <- ggplot(
  mdbk_bRSVN, aes(
    x = Condition, 
    y = Value_norm, 
    fill = Condition)) +
  geom_violin(
    trim=FALSE,
    alpha = 0.5,
    scale = 'width',
    adjust = 0.7) +
  stat_summary(
    fun.data=mean_se, 
    fun.args = list(mult=1), 
    geom="pointrange", 
    color="black",
    show.legend = F) +
  scale_x_discrete(
    limits = list_of_conditions) +
  scale_fill_manual(
    breaks = list_of_conditions,
    values = list_of_colours) +
  theme_bw() +
  figure_theme +
  labs(
    title = plot_title,
    y = y_axis_title,
    x = NULL
  ) +
  scale_y_log10(labels = trans_format("log10", 
                                      math_format(10^.x)),
                breaks = trans_breaks("log10", 
                                      function(x) 10^x, n = 7),
                limits = c(10^-1, 10^9)
  ) +
  annotation_logticks(sides='l') +
 
  geom_signif(
    comparisons = list(c(
      list_of_conditions[1], 
      list_of_conditions[2])), 
    annotation = p_val[1], 
    y_position = 7, 
    tip_length = 0, 
    vjust= -0.2, 
    size = 0.7, 
    textsize = textsize_values[1])

plot_mdbk_bRSVN

# saving data and plot ------------------------------

###########
dpi <- 600
width <- 16
height <- 20
file_name <- 'mdbk_bRSVN'
###########

ggsave(filename = paste(
  file_name, '.svg', sep = ''), 
       plot = plot_mdbk_bRSVN, 
       device = 'svg', 
       path = 'Figures', 
       dpi = dpi, 
       height = height, 
       width = width, 
       units = 'cm')
ggsave(filename = paste(
  file_name, '.png', sep = ''), 
       plot = plot_mdbk_bRSVN, 
       device = 'png', 
       path = 'Figures', 
       dpi = dpi, 
       height = height, 
       width = width, 
       units = 'cm')
                
fwrite(mdbk_bRSVN, 
       paste('Adjusted-Data/', 
             file_name, 
             '.csv', 
             sep = ''))
