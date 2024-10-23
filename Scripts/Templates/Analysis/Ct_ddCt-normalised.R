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
# change 'hIFIT1', 'A549', 'hifit1'
###########

###########
list_of_conditions <- c(
  'Mock',
  'bRSV dSH 1-24',
  'hRSV 1-24 ROXO 5-24',
  "hIFNa 1000-24"
  )
list_of_colours <- c(
  '#999999', 
  '#9a6a00', 
  '#009e73',
  "#0484a5"
  )
###########

# data prep ----
housekeeping_control <- fread(
  paste('Data/', 
        
        ############
        'ct_data', 
        ############
        
        '.csv', 
        sep = '')
  ) %>% 
  
  ###########
  filter(Target == 'hGAPDH',
         Cell_line == "A549",
         Condition %in% list_of_conditions
         ) %>%
  ###########

  group_by(Condition) %>%
  mutate(mean_Ct = mean(Ct)) %>%
  ungroup() %>%
  arrange(match(
    Condition, 
    list_of_conditions))

hifit1 <- gene_of_interest_data %>%
  
  ###########
  filter(
    Target == 'hIFIT1',
    Cell_line == "A549",
    Condition %in% list_of_conditions
         ) %>%
  ###########

  arrange(match(
    Condition, 
    list_of_conditions)) %>%
  mutate(
    Control_ct_mean = housekeeping_control$mean_Ct,
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

hifit1

# data analysis ----
boxplot(Value_norm~Condition, hifit1)
plot(lm(Value_norm~Condition, hifit1))

# test normality
shapiro.test(hifit1$Value_norm[1:3])
shapiro.test(hifit1$Value_norm[4:6])
shapiro.test(hifit1$Value_norm[7:9])
shapiro.test(hifit1$Value_norm[10:12])

plot(residuals(
  lm(Value_norm~Condition, 
     hifit1)))
shapiro.test(residuals(
  lm(Value_norm~Condition, 
     hifit1)))

### normal distribution

# if distribution normal
bartlett.test(Value_norm~Condition, 
              hifit1)

### non equal variance

# for normal distribution but non equal variance
library(userfriendlyscience)
oneway.test(Value_norm~Condition, 
            hifit1, 
            var.equal = FALSE)
oneway(hifit1$Value_norm, 
       as.factor(hifit1$Condition), 
       posthoc = 'games-howell')

# plot data ----
boxplot(Value_norm~Condition, hifit1)
hifit1
list_of_conditions

###########
p_val <- c(
  0.8031417,
  0.0464754,
  0.3186545
  )
plot_title <- 'A549 - hIFIT1'
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

plot_hifit1 <- ggplot(
  hifit1, 
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
    y_position = 6, 
    tip_length = 0, 
    vjust= -0.2, 
    size = 0.7,
    textsize = textsize_values[1]) +
  
  geom_signif(
    comparisons = list(c(
      list_of_conditions[1], 
      list_of_conditions[3])), 
    annotation = p_val[2], 
    y_position = 7, 
    tip_length = 0, 
    vjust= -0.2, 
    size = 0.7, 
    textsize = textsize_values[2]) +
  
  geom_signif(
    comparisons = list(c(
      list_of_conditions[1], 
      list_of_conditions[4])), 
    annotation = p_val[3], 
    y_position = 8, 
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
file_name <- 'hifit1'
###########

ggsave(filename = paste(
  file_name, '.svg', sep = ''), 
       plot = plot_hifit1, 
       device = 'svg', 
       path = 'Figures', 
       dpi = dpi, 
       height = height, 
       width = width, 
       units = 'cm')
ggsave(filename = paste(
  file_name, '.png', sep = ''), 
       plot = plot_hifit1, 
       device = 'png', 
       path = 'Figures', 
       dpi = dpi, 
       height = height, 
       width = width, 
       units = 'cm')

fwrite(hifit1, 
       paste('Adjusted-Data/', 
             file_name, 
             '.csv', 
             sep = ''))
