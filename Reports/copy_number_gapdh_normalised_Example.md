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
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.3     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.2     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter()     masks stats::filter()
    ## ✖ dplyr::group_rows() masks kableExtra::group_rows()
    ## ✖ dplyr::lag()        masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(data.table)
```

    ## 
    ## Attaching package: 'data.table'
    ## 
    ## The following objects are masked from 'package:lubridate':
    ## 
    ##     hour, isoweek, mday, minute, month, quarter, second, wday, week,
    ##     yday, year
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     between, first, last
    ## 
    ## The following object is masked from 'package:purrr':
    ## 
    ##     transpose

``` r
library(ggsignif)
library(scales)
```

    ## 
    ## Attaching package: 'scales'
    ## 
    ## The following object is masked from 'package:purrr':
    ## 
    ##     discard
    ## 
    ## The following object is masked from 'package:readr':
    ## 
    ##     col_factor

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

## Model Production

``` r
lm(log10(Copy_number)~Ct, data = standard_curve_data) # create a model
```

    ## 
    ## Call:
    ## lm(formula = log10(Copy_number) ~ Ct, data = standard_curve_data)
    ## 
    ## Coefficients:
    ## (Intercept)           Ct  
    ##     10.5888      -0.3208

``` r
model_standard_curve <- lm(log10(Copy_number)~Ct, data = standard_curve_data) # save it to a variabel
```

### Calculate Amplification Efficiency

Some text about what the fuck are we doing

Put amplification efficiency equation

``` math
\mbox{Amplification Efficiency} = 10^{-1/\mbox{slope}}-1
```

``` r
paste(
  round(
    (10^(-1/ lm(Ct~log10(Copy_number), 
                data = standard_curve_data)[[1]][2]) -1)*100, 
    digits = 2), 
  '% Amplification Efficiency', 
  sep = ''
  ) # calculate the amplification efficiency
```

    ## [1] "110% Amplification Efficiency"

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

## Figure

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
