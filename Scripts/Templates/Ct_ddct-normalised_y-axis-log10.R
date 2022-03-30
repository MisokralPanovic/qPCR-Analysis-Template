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
# change 'bRSV N', ' 'mdbk_brsv'
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

# data prep ----
gene_of_interest_data <- fread(
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
         Cell_line == "MDBK",
        Condition %in% list_of_conditions
) %>%
  ###########

  group_by(Condition) %>%
  mutate(mean_Ct = mean(Ct)) %>%
  ungroup() %>%
  arrange(match(
    Condition, 
    list_of_conditions))

mdbk_brsv <- gene_of_interest_data %>%
  
  ###########
filter(
  Target == 'bRSV_N',
  Cell_line == "MDBK",
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

mdbk_brsv

# data analysis ----
boxplot(Value_norm~Condition, mdbk_brsv)
plot(lm(Value_norm~Condition, mdbk_brsv))

# test normality
shapiro.test(mdbk_brsv$Value_norm[1:3])
shapiro.test(mdbk_brsv$Value_norm[4:6])

plot(residuals(
  lm(Value_norm~Condition, 
     mdbk_brsv)))
shapiro.test(residuals(
  lm(Value_norm~Condition, 
     mdbk_brsv)))

### normal distribution

# if distribution normal
bartlett.test(Value_norm~Condition, mdbk_brsv)

### non equal variance

# for normal distribution but non equal variance - multiple comparison
library(rstatix)
library(dplyr)
oneway.test(Value_norm~Condition, 
            mdbk_brsv, 
            var.equal = F)
mdbk_brsv %>% games_howell_test(Value_norm~Condition)

# plot data ----
boxplot(Value_norm~Condition, mdbk_brsv)
mdbk_brsv
list_of_conditions

###########
p_val <- c(
  0.036,
  0.232
)
plot_title <- 'MDBK - bRSV N'
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

plot_mdbk_brsv <- ggplot(
  mdbk_brsv, 
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

plot_mdbk_brsv

# saving data and plot ------------------------------

###########
dpi <- 600
width <- 16
height <- 20
file_name <- 'mdbk_brsv'
###########

ggsave(filename = paste(
  file_name, '.svg', sep = ''), 
  plot = plot_mdbk_brsv, 
  device = 'svg', 
  path = 'Figures', 
  dpi = dpi, 
  height = height, 
  width = width, 
  units = 'cm')
ggsave(filename = paste(
  file_name, '.png', sep = ''), 
  plot = plot_mdbk_brsv, 
  device = 'png', 
  path = 'Figures', 
  dpi = dpi, 
  height = height, 
  width = width, 
  units = 'cm')

fwrite(mdbk_brsv, 
       paste('Adjusted-Data/', 
             file_name, 
             '.csv', 
             sep = ''))
