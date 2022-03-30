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
# change 'bifit1', 'MDBK', 'bIFIT1', 'model_lmscb1',
###########

###########
list_of_conditions <- c(
  'Mock',
  'bRSV dSH 1-24'
  )
list_of_colours <- c(
  '#999999', 
  '#ffc034'
  )
###########

# data prep ----
housekeeping_gene_data <- fread(
  paste('Data/', 
        ############
        'ct_data', 
        ############
        
        '.csv', 
        sep = ''
        )
  ) %>% 
                                       
  ########################
  filter(Target == 'bGAPDH',
         Cell_line == "MDBK",
         Condition %in% list_of_conditions
         ) %>%
  ########################
  
  mutate(
    control_mean = mean(
       Ct[Condition == list_of_conditions[1]],
       na.rm = T),
    log2_dCt = 2^ (- (Ct - control_mean)),
    control_mean_log = mean(
      log2_dCt[Condition == list_of_conditions[1]],
      na.rm = T),
    Value_norm = log2_dCt / control_mean_log
    ) %>%
  group_by(Condition) %>%
  mutate(mean_Ct = mean(Ct)) %>%
  ungroup() %>%
  arrange(match(
    Condition, 
    list_of_conditions))

bifit1 <- fread(
  paste('Data/', 
        
        ############
        'copy_number_extrapolation_data', 
        ############
        
        '.csv', 
        sep = ''
        )
  )  %>% 
  ###########
  filter(Target == 'bIFIT1',
       Cell_line == "MDBK",
       Condition %in% list_of_conditions
  ) %>%
  ###########
  arrange(
    match(Condition, 
          list_of_conditions))
  
bifit1 <- bifit1 %>% 
  mutate(
    Copy_number = 10^predict(model_lmscb1, 
                             newdata = bifit1),
    Factor = housekeeping_gene_data$mean_Ct,
    Copy_number_mod = Copy_number / Factor,
    Control_mean = mean(
      Copy_number_mod[Condition == list_of_conditions[1]], 
      na.rm = T
      ),
    Value_norm = Copy_number_mod / Control_mean,
    Value_norm_old = Copy_number / Control_mean
    )       

bifit1

# data analysis ----
boxplot(Value_norm~Condition, bifit1)
plot(lm(Value_norm~Condition, bifit1))

# test normality
shapiro.test(bifit1$Value_norm[1:3])
shapiro.test(bifit1$Value_norm[4:6])

plot(residuals(lm(Value_norm~Condition, bifit1)))
shapiro.test(residuals(lm(Value_norm~Condition, bifit1)))

### normal distribution

# if distribution normal
bartlett.test(Value_norm~Condition, bifit1)

### equal variance

# for normal distribution and variance
library(DTK)
anova(lm(Value_norm~Condition, bifit1))
TukeyHSD(aov(Value_norm~Condition, bifit1))

# plot data ----
boxplot(Value_norm~Condition, bifit1)
bifit1
list_of_conditions

###########
p_val <- c(
  0.8031417
  )
plot_title <- 'bIFIT1'
y_axis_title <- 'Relative mRNA Levels'
###########

textsize_values <- c()

for (value in p_val) {
  if (value > 0.05) {
    textsize_values <- append(
      textsize_values, 3)
  } else {
    textsize_values <- append(
      textsize_values, 4)
  }
}

plot_bifit1 <- ggplot(
  bifit1, 
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
  theme_bw() +
  figure_theme +
  labs(
    title = plot_title,
    y = y_axis_title,
    x = NULL
  ) +
  scale_y_continuous(trans = log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x, n = 8),
                     labels = trans_format("log2", math_format(2^.x)),
                     limits = c(2^-7,2^9),
                     sec.axis = sec_axis(trans = identity,
                                         breaks = c(2^-6, 2^-4, 2^-2, 2^0, 2^2, 2^4, 2^6, 2^8),
                                         labels = c(0.016, 0.062, 0.25, 1, 4, 16, 64, 256)
                     )) +
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

plot_bifit1

# saving data and plot ------------------------------

###########
dpi <- 600
width <- 16
height <- 20
file_name <- 'bifit1'
###########

ggsave(filename = paste(
  file_name, '.svg', sep = ''), 
       plot = plot_bifit1, 
       device = 'svg', 
       path = 'Figures', 
       dpi = dpi, 
       height = height, 
       width = width, 
       units = 'cm')
ggsave(filename = paste(
  file_name, '.png', sep = ''), 
       plot = plot_bifit1, 
       device = 'png', 
       path = 'Figures', 
       dpi = dpi, 
       height = height, 
       width = width, 
       units = 'cm')

fwrite(bifit1, 
       paste('Adjusted-Data/', 
             file_name, 
             '.csv', 
             sep = ''))
