## ---- plot_theme ----
plot_theme <- theme(
  plot.title = element_text(
    size = 20, 
    face = 'bold', 
    margin = margin(10, 0, 10, 0), 
    hjust = 0.5
  ),
  legend.text = element_text(size=15),
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
    margin = margin(0, 10, 0, 0)))
## ---- lm_figure_loading ----
plot_scdata <- ggplot(
  data = cbind(scdata, 10^predict(model_lmscdata, interval = 'prediction')),
  mapping = aes(x = Copy_number, y = Ct)) +
  geom_point() +
  stat_smooth(method = lm) +
  labs(title = paste(params$target_gene, 'Standard Curve', sep = " "),
       y = 'Cycle Threshold',
       x = 'Copy Number'
  ) +
  scale_y_continuous(
    breaks = seq(from = 0, 
                 to = 45, 
                 by = 10),
    limits = c(0,45)
  ) +
  theme_minimal()+
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)), breaks = trans_breaks("log10", function(x) 10^x)) +
  annotation_logticks(sides='b') +
  plot_theme +
  annotate('text',
           y = 37.5, 
           x = 1000000, 
           label = paste(round((10^(-1/ lm(Ct~log10(Copy_number),  data = scdata)[[1]][2] ) -1)*100,
                               digits = 2), '% Amplification Efficiency', sep = ''), 
           size = 5)
## ---- final_plotting ----
textsize_values <- c()

for (value in params$pvals) {
  if (value > 0.05) {
    textsize_values <- append(
      textsize_values, 4)
  } else {
    textsize_values <- append(
      textsize_values, 5)
  }
}

plot_target_data <- ggplot(
  target_data, 
  aes(x = Condition, y = Value_norm,fill = Condition)) +
  geom_violin(trim = T, alpha = 0.5, scale = 'width', adjust = 0.7) +
  stat_summary(fun.data = mean_se, fun.args = list(mult=1), geom = "pointrange", color = "black", show.legend = F) +
  scale_x_discrete(limits = params$conditions) +
  scale_fill_manual(breaks = params$conditions, values = params$colors) +
  theme_minimal()+
  plot_theme +
  labs(
    title = params$target_gene,
    y = 'Relative mRNA Levels',
    x = NULL
  ) +
  coord_cartesian(ylim = c(0, params$range)) +
  scale_y_continuous(breaks= seq(0, params$range, params$breaks)) +
  theme(aspect.ratio = 2/1, axis.text.x=element_blank(), axis.ticks.x=element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())