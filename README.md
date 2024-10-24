Example of qPCR Pipeline Usage
================
Michal Varga
Last edited: 2024-10-24

# Introduction

Purpose:

based on scripts:

- standard curve
- control gene Ct based
- Copy number GAPDH normalised dataset

parameters incuded in script:

- `target_gene:` bIFIT1
- `housekeeping_gene:` bGAPDH
- `conditions:` !r c(‘Mock’, ‘phRSV 1-24’,‘bIFNa 5-6’)
- `colors:` !r c(‘\#999999’, ‘\#ffc034’, ‘\#9a6a00’)
- `statistics:` “Normal distribution and equal variance”
- `pvals:` !r c(0.9935156, 0.1817141)
- `range:` 14
- `breaks:` 2

Link to stats pipeline?

## Libraries (maybe put them when they are used?)

``` r
library(tidyverse)
library(data.table)
library(ggsignif)
library(scales)
```

# 1. Standard Curve

## Introduction

Some text about what the fuck

## Data Loading

Load `standard_curves_data.csv` from `/Data` folder, and filter the
target gene and primer set used based on the declared parameters in
YAML.

``` r
standard_curve_data <- fread('../Data/standard_curves_data.csv') %>%
  filter(Target == params$target_gene,
         Primer_set == params$primer_set)
```

|        Ct | Copy_number | Primer_set | Target | Experiment | Additional_info |
|----------:|------------:|:-----------|:-------|:-----------|:----------------|
|  7.664552 |       1e+08 | PS11       | bIFIT1 | NA         | NA              |
| 11.265419 |       1e+07 | PS11       | bIFIT1 | NA         | NA              |
| 14.448029 |       1e+06 | PS11       | bIFIT1 | NA         | NA              |
| 17.917984 |       1e+05 | PS11       | bIFIT1 | NA         | NA              |
| 21.056202 |       1e+04 | PS11       | bIFIT1 | NA         | NA              |
| 22.973955 |       1e+03 | PS11       | bIFIT1 | NA         | NA              |
| 26.622187 |       1e+02 | PS11       | bIFIT1 | NA         | NA              |

## Linear Model Production

### Construct Linear Model

``` r
model_standard_curve <- lm(log10(Copy_number)~Ct, data = standard_curve_data)
```

    ## 
    ## Call:
    ## lm(formula = log10(Copy_number) ~ Ct, data = standard_curve_data)
    ## 
    ## Coefficients:
    ## (Intercept)           Ct  
    ##     10.5888      -0.3208

### Construct a Prediction Based on the Linear Model

``` r
prediction_model_standard_curve <- 10^predict(
  model_standard_curve, 
  interval = 'prediction'
  )
```

    ## Warning in predict.lm(model_standard_curve, interval = "prediction"): predictions on current data refer to _future_ responses

|          fit |          lwr |          upr |
|-------------:|-------------:|-------------:|
| 1.348840e+08 | 4.340240e+07 | 4.191860e+08 |
| 9.436060e+06 | 3.300113e+06 | 2.698066e+07 |
| 8.990821e+05 | 3.288738e+05 | 2.457929e+06 |
| 6.928302e+04 | 2.569156e+04 | 1.868371e+05 |
| 6.821448e+03 | 2.478090e+03 | 1.877743e+04 |
| 1.654456e+03 | 5.848611e+02 | 4.680130e+03 |
| 1.117612e+02 | 3.650013e+01 | 3.422060e+02 |

``` r
data_combined_standard_curve <- cbind(
  standard_curve_data,
  prediction_model_standard_curve
  ) # save it to a variabel
```

| Ct | Copy_number | Primer_set | Target | Experiment | Additional_info | fit | lwr | upr |
|---:|---:|:---|:---|:---|:---|---:|---:|---:|
| 7.664552 | 1e+08 | PS11 | bIFIT1 | NA | NA | 1.348840e+08 | 4.340240e+07 | 4.191860e+08 |
| 11.265419 | 1e+07 | PS11 | bIFIT1 | NA | NA | 9.436060e+06 | 3.300113e+06 | 2.698066e+07 |
| 14.448029 | 1e+06 | PS11 | bIFIT1 | NA | NA | 8.990821e+05 | 3.288738e+05 | 2.457929e+06 |
| 17.917984 | 1e+05 | PS11 | bIFIT1 | NA | NA | 6.928302e+04 | 2.569156e+04 | 1.868371e+05 |
| 21.056202 | 1e+04 | PS11 | bIFIT1 | NA | NA | 6.821448e+03 | 2.478090e+03 | 1.877743e+04 |
| 22.973955 | 1e+03 | PS11 | bIFIT1 | NA | NA | 1.654456e+03 | 5.848611e+02 | 4.680130e+03 |
| 26.622187 | 1e+02 | PS11 | bIFIT1 | NA | NA | 1.117612e+02 | 3.650013e+01 | 3.422060e+02 |

### Calculate Amplification Efficiency

Some text about what the fuck are we doing

Put amplification efficiency equation

``` math
\mbox{Amplification Efficiency} = 10^{-1/\mbox{slope}}-1
```

``` r
efficiency_standard_curve <- paste(
  round(
    (10^(-1/ lm(Ct~log10(Copy_number), 
                data = standard_curve_data)[[1]][2]) -1)*100, 
    digits = 2), 
  '% Amplification Efficiency', 
  sep = ''
  ) # save it as a variable
```

    ## [1] "110% Amplification Efficiency"

## Figure

### Constructing Figure

``` r
plot_title_standardCurve <- paste(params$target_gene, 
                    '-', 
                    params$primer_set,
                    'Standard Curve', 
                    sep = " ")
x_annotation_position_standardCurve <- 1000000
y_annotation_position_standardCurve <- 37.5
top_range_standardCurve <- 40
y_axis_title_standardCurve <- 'Cycle Threshold'
x_axis_title_standardCurve <- 'Copy Number'
```

``` r
plot_standard_curve <- ggplot(
  data = data_combined_standard_curve,
  mapping = aes(x = Copy_number,
                y = Ct)) +
  geom_point() +
  stat_smooth(method = lm) +
  labs(
    title = plot_title_standardCurve,
    y = x_axis_title_standardCurve,
    x = y_axis_title_standardCurve
  )
```

``` r
plot_standard_curve <- plot_standard_curve +
  scale_y_continuous(
    breaks = seq(from = 0, 
                 to = top_range_standardCurve, 
                 by = 10),
    limits = c(0,top_range_standardCurve)
  )
```

``` r
plot_standard_curve <- plot_standard_curve +
  scale_x_log10(
    labels = trans_format("log10", 
                          math_format(10^.x)),
    breaks = trans_breaks("log10", 
                          function(x) 10^x)
  ) +
  annotation_logticks(sides='b')
```

``` r
plot_standard_curve <- plot_standard_curve +
  theme_bw() +
  theme(
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
```

``` r
plot_standard_curve <- plot_standard_curve +
  annotate('text',
           y = y_annotation_position_standardCurve, 
           x = x_annotation_position_standardCurve, 
           label = efficiency_standard_curve, 
           size = 5)
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](Reports/copy_number_gapdh_normalised_Example_files/figure-gfm/standard_curve_plot_display,%20-1.png)<!-- -->

### Saving Figure

# 2. Housekeeping Gene Control

## Introduction

Some text about what the fuck

## Data Loading

``` r
housekeeping_gene_data <- fread('../Data/ct_data.csv') %>% 
  filter(Target == params$housekeeping_gene, Condition %in% params$conditions)
```

## Processing Data

## Make Factors

# 3. Factorised Copy Number Extrapolation

## Introduction

Some text about what the fuck

## Data loading and filtering

``` r
target_data <- fread('../Data/copy_number_extrapolation_data.csv') %>%
  filter(Target == params$target_gene, Condition %in% params$conditions)
```

## Copy number extrapolation

## Factorisation based on housekeeping gene levels

## Statistics

## Plotting

# Libraries and Data Loading

the libraries needed and why

# Data Manipulation

## Linear Model Creation

``` r
model_lmscdata <- lm(log10(Copy_number)~Ct,  data = scdata)
```

## Linear Model Figure

``` r
plot_scdata
```

## Housekeeping Log2 Data Generation

``` r
housekeeping_gene_data <-  housekeeping_gene_data %>% 
  mutate(control_mean = mean(Ct[Condition == params$conditions[1]], na.rm = T),
        log2_dCt = 2^ (- (Ct - control_mean)),
        control_mean_log = mean(log2_dCt[Condition == params$conditions[1]], na.rm = T),
        Value_norm = log2_dCt / control_mean_log) %>%
  group_by(Condition) %>%
  mutate(mean_Ct = mean(Ct)) %>%
  ungroup() %>%
  arrange(match(Condition, params$conditions))
```

## Extrapolation From Standard Curve and Factorisation by Housekeeping Gene Data

``` r
target_data <- target_data %>% 
  arrange(match(Condition, params$conditions)) %>%
  mutate(Copy_number = 10^predict(model_lmscdata, newdata = target_data),
         Factor = housekeeping_gene_data$mean_Ct,
         Copy_number_mod = Copy_number / Factor,
         Control_mean = mean(Copy_number_mod[Condition == params$conditions[1]], na.rm = T),
         Value_norm = Copy_number_mod / Control_mean)
```

``` r
kable(target_data, caption = "Final Table")
```

## Dataset Save

``` r
print("SAVE THE DATAAA")
```

Data was saved to **HELL**

# Statistics

## Normal distribution and equal variance

``` r
library(DTK, quietly = T)
anova(lm(Value_norm~Condition, target_data))
TukeyHSD(aov(Value_norm~Condition, target_data))
```

# Data Visualisation

``` r
plot_target_data +
  geom_signif(
    comparisons = list(c(
      params$conditions[1], 
      params$conditions[2])), 
    annotation = params$pvals[1], 
    y_position = 0.93*params$range - 1*(params$range*0.075), 
    tip_length = 0, vjust= -0.2, size = 0.7,
    textsize = textsize_values[1]) +
              
  geom_signif(
    comparisons = list(c(
      params$conditions[1], 
      params$conditions[3])), 
    annotation = params$pvals[2], 
    y_position = 0.93*params$range - 0*(params$range*0.075), 
    tip_length = 0, vjust= -0.2, size = 0.7, 
    textsize = textsize_values[2])
```

## Figure Save

``` r
print("SAVE THE FIGUREEEEE")
```

Figure was saved to **HELL**
