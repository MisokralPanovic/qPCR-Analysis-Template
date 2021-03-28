.libPaths('C:/r_packages')
library(tidyverse)
library(ggsignif)
library(data.table)
library(ggsignif)
library(scales)
library(ggplot2)
library(ggpubr)

###########
# change 'MDBK', 'bRSV_N', 'mdbk_bRSVN', 'bRSV N',
###########

###########
list_of_conditions <- c('Mock',
                        'bRSV dNS1 0.001-24',
                        'bRSV dNS2 0.01-24',
                        'bRSV dNS1/2 0.001-24')
list_of_colours <- c('#999999', 
                     '#ffc034', 
                     '#9a6a00', 
                     '#009e73')
###########

# data prep -------------------------
all_RSV_N_data <- read.csv(paste('Data/', 
      
                                 ############
                                 'bRSV_data', 
                                 ############

                                 '.csv', 
                                 sep = ''))

# all dataset normalisation - reverse values so mock low and infection high
all_RSV_N_data <- all_RSV_N_data %>%
  mutate(normalisator = 45) %>%
  mutate(Value_adju = 2^abs(Ct - normalisator))

mdbk_bRSVN <- all_RSV_N_data %>%
  
  ###########
  filter(CellLine == 'MDBK') %>%
  filter(Target == 'bRSV_N') %>%
  ###########
  
  mutate(Mock_mean = mean(Value_adju[Condition == list_of_conditions[1]],
                          na.rm = T),
         Value_norm = Value_adju / Mock_mean) %>% 
  
  arrange(match(Condition, list_of_conditions))

mdbk_bRSVN

# data analysis ----
boxplot(Value_norm~Condition, mdbk_bRSVN)
plot(lm(Value_norm~Condition, mdbk_bRSVN))

# test normality
shapiro.test(mdbk_bRSVN$Value_norm[1:3])
shapiro.test(mdbk_bRSVN$Value_norm[4:6])
shapiro.test(mdbk_bRSVN$Value_norm[7:9])
shapiro.test(mdbk_bRSVN$Value_norm[10:12])

plot(residuals(lm(Value_norm~Condition, mdbk_bRSVN)))
shapiro.test(residuals(lm(Value_norm~Condition, mdbk_bRSVN)))

### non normal distribution

# if non normal distribution
library(car)
leveneTest(Value_norm~Condition, Value_norm~Condition, mdbk_bRSVN)

### equal variance
# for non normal distribution but equal variance
t.test(Value_norm~Condition, data=mdbk_bRSVN, alternative='two.sided',
       var.equal=T)

# plot data ----
boxplot(Value_norm~Condition, mdbk_bRSVN)
mdbk_bRSVN

###########
p_val <- c(0.8031417,
           0.9655795,
           0.0564754)
range_y <- 8
plot_title <- 'MDBK - bRSV N'
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

plot_mdbk_bRSVN <- ggplot(mdbk_bRSVN, aes(x = Condition, 
                                          y =Value_norm, 
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
  scale_y_log10(labels = trans_format("log10", 
                                      math_format(10^.x)),
                breaks = trans_breaks("log10", 
                                      function(x) 10^x)
  ) +
  annotation_logticks(sides='l') +
 
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


plot_mdbk_bRSVN

# saving data and plot ------------------------------

###########
dpi <- 600
width <- 20
height <- 20
file_name <- 'mdbk_bRSVN'
###########

ggsave(filename = paste(file_name, '.png', sep = ''), 
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