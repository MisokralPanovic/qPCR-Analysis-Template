.libPaths('C:/r_packages')
library(tidyverse)

###########
# change '_condition1', 'c('Mock', 'cbRSV 0.1-24', 'cbRSV 1-24', 'cbRSV 2-24')'
###########

# start with control
########################
conditions <- c('Mock',
                'cbRSV 0.1-24',
                'cbRSV 1-24',
                'cbRSV 2-24')
file_emmit <- c(T,T,T,
                T,T,T,
                T,T,T,
                T,T,T)
replicates <- 3
########################

gapdh_data <- read.csv(paste('Data/', 
      
                             ############
                             'bRSV_data', 
                             ############

                             '.csv', 
                             sep = ''))

# --------------------------------------------------------------------------------------------------

gapdh_data_condition1 <- gapdh_data %>%
  
  ########################
  filter(Time == '24h') %>%
  ########################
  
  # run the rest to get list of factor values
  mutate(normalisator = 45) %>%
  mutate(Value_adju = 2^abs(Value - normalisator))

gapdh_data_condition1 <- aggregate(gapdh_data_condition1[-1],
                           list(gapdh_data_condition1$Condition),
                           mean)
gapdh_data_condition1 <- gapdh_data_condition1 %>%
  mutate(Mock_mean = Value_adju[Group.1 == conditions[1]],
         Value_norm = Value_adju / Mock_mean) %>%
  arrange(match(Group.1, conditions))


gapdh_ratios_condition1 <- rep(gapdh_data_condition1$Value_nor, each=replicates)
gapdh_ratios_mdbk <- gapdh_ratios_mdbk[file_emmit]

# final vector
gapdh_ratios_condition1

