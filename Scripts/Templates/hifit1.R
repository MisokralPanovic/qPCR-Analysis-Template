.libPaths('C:/r_packages')
library(tidyverse)
library(ggsignif)
library(data.table)

###########
# change 'hIFIT1', 'A549', 'hifit1'
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

hifit1 <- initial_data %>%
  
  ###########
  filter(Target == 'hIFIT1') %>%
  filter(CellLine == 'A549') %>%
  ###########

  arrange(match(Condition, list_of_conditions)) %>%
  mutate(Mock_mean = mean(Value[Condition == list_of_conditions[1]],
                          na.rm = T),
         Value_norm = Value / Mock_mean)

hifit1

# data analysis ----
boxplot(Value_norm~Condition, hifit1)
plot(lm(Value_norm~Condition, hifit1))

# test normality
shapiro.test(hifit1$Value_norm[1:3])
shapiro.test(hifit1$Value_norm[4:6])
shapiro.test(hifit1$Value_norm[7:9])
shapiro.test(hifit1$Value_norm[10:12])

plot(residuals(lm(Value_norm~Condition, hifit1)))
shapiro.test(residuals(lm(Value_norm~Condition, hifit1)))

### normal distribution

# if distribution normal
bartlett.test(Value_norm~Condition, hifit1)

### non equal variance

# for normal distribution but non equal variance
library(userfriendlyscience)
oneway.test(Value_norm~Condition, hifit1, var.equal = FALSE)
oneway(hifit1$Value_norm, as.factor(hifit1$Condition), posthoc = 'games-howell')


# plot data ----
boxplot(Value_norm~Condition, hifit1)
hifit1

###########
p_val <- c(0.8031417,
           0.9655795,
           0.0564754)
range_y <- 20
breaks_y <- 5
plot_title <- 'A549 - hIFIT1'
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

plot_mdbk_bi1 <- ggplot(hifit1, aes(x = Condition, 
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


plot_hifit1


# saving data and plot ------------------------------

###########
dpi <- 600
width <- 16
height <- 20
file_name <- 'A549_hIFIT1'
###########

ggsave(filename = paste(file_name, '.png', sep = ''), 
       plot = plot_hifit1, 
       device = 'png', 
       path = 'Figures', 
       dpi = dpi, 
       height = height, 
       width = width, 
       units = 'cm')

fwrite(hifit1, 
       paste('Adjusted Data/', 
             file_name, 
             '.csv', 
             sep = ''))
