.libPaths('C:/r_packages')
library(tidyverse)
library(ggsignif)
library(data.table)

###########
# change 'b24_1', 'MDBK', 'bIFIT1', 'model_lmsc1', 'gapdh_ratios_24'
###########

###########
list_of_conditions <- c('Mock',
                        'cbRSV dNS1 0.001-24',
                        'cbRSV dNS2 0.01-24',
                        'cbRSV dNS1/2 0.001-24')
list_of_colours <- c('#999999', 
                     '#ffc034', 
                     '#9a6a00', 
                     '#009e73')
###########

# data prep ----
initial_data <- read.csv(paste('Data/', 
      
                               ############
                               'bRSV_data', 
                               ############

                               '.csv', 
                               sep = ''))

b24_1 <- initial_data %>%
  
  ###########
  filter(TimePoint == 24) %>%
  filter(Target == 'bIFIT1')
  ###########

b24_1 <- b24_1 %>%
  mutate(Copy_number = 10^predict(model_lmsc1, newdata = b24_1)) %>%
  arrange(match(Condition, list_of_conditions)) %>%

    mutate(Factor = gapdh_ratios_24) %>%
  mutate(Copy_number_mod = Copy_number / Factor) %>%
  mutate(Mock_mean = mean(Copy_number_mod[Condition == list_of_conditions[1]],
                          na.rm = T),
         Value_norm = Copy_number_mod / Mock_mean) %>%
  
  mutate(Value_norm_old = Copy_number / Mock_mean)
  
b24_1

# data analysis ----
boxplot(Value_norm~Condition, b24_1)
plot(lm(Value_norm~Condition, b24_1))

# test normality
shapiro.test(b24_1$Value_norm[1:3])
shapiro.test(b24_1$Value_norm[4:6])
shapiro.test(b24_1$Value_norm[7:9])
shapiro.test(b24_1$Value_norm[10:12])

plot(residuals(lm(Value_norm~Condition, b24_1)))
shapiro.test(residuals(lm(Value_norm~Condition, b24_1)))

### normal distribution

# if distribution normal
bartlett.test(Value_norm~Condition, b24_1)

### equal variance

# for normal distribution and variance
library(DTK)
anova(lm(Value_norm~Condition, b24_1))
TukeyHSD(aov(Value_norm~Condition, b24_1))

# plot data ----
boxplot(Value_norm~Condition, b24_1)
b24_1

###########
p_val <- c(0.8031417,
           0.9655795,
           0.0564754)
range_y <- 2
breaks_y <- 0.5
plot_title <- 'bIFIT1 - 24h'
y_axis_title <- 'mRNA Fold Change'
###########

textsize_values <- c()

for (value in p_val) {
  if (value > 0.05) {
    textsize_values <- append(textsize_values, 4)
  } else {
    textsize_values <- append(textsize_values, 5)
  }
}

plot_b24_1 <- ggplot(b24_1, aes(x = Condition, 
                                y = Value_norm, 
                                fill = Condition)) +
  geom_boxplot(varwidth = T) +
  scale_x_discrete(limits = list_of_conditions) +
  scale_fill_manual(breaks = list_of_conditions,
                    values = list_of_colours) +
  theme(
    plot.title = element_text(
      size=20, 
      face='bold', 
      margin = margin(10, 0, 10, 0), 
      hjust = 0.5
    ),
    legend.text = element_text(size=15),  
    legend.title=element_blank(),
    axis.text.y=element_text(angle=0, 
                             size=12, 
                             vjust=0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 15, 
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
  scale_y_continuous(breaks= seq(0, range_y, breaks_y), 
                     limits = c(0, range_y)) +
  
  geom_signif(comparisons = list(c(list_of_conditions[1], list_of_conditions[2])), 
              annotation = p_val[1], 
              y_position = 0.93*range_y - 2*(range_y*0.075), 
              tip_length = 0, 
              vjust= -0.2, 
              size = 0.7, 
              textsize = textsize_values[1]) +
  geom_signif(comparisons = list(c(list_of_conditions[1], list_of_conditions[3])), 
              annotation = p_val[2], 
              y_position = 0.93*range_y - 1*(range_y*0.075), 
              tip_length = 0, 
              vjust= -0.2, 
              size = 0.7, 
              textsize = textsize_values[2]) +
  geom_signif(comparisons = list(c(list_of_conditions[1], list_of_conditions[4])), 
              annotation = p_val[3], 
              y_position = 0.93*range_y - 0*(range_y*0.075), 
              tip_length = 0, 
              vjust= -0.2, 
              size = 0.7, 
              textsize = textsize_values[3])


plot_b24_1

# saving data and plot ------------------------------

###########
dpi <- 600
width <- 16
height <- 20
file_name <- 'b24_1'
###########

ggsave(filename = paste(file_name, '.png', sep = ''), 
       plot = plot_b24_1, 
       device = 'png', 
       path = 'Figures', 
       dpi = dpi, 
       height = height, 
       width = width, 
       units = 'cm')

fwrite(b24_1, 
       paste('Adjusted Data/', 
             file_name, 
             '.csv', 
             sep = ''))
