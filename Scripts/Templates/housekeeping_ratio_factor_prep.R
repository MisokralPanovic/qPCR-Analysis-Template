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
  
  mutate(control_mean = mean([Condition == list_of_conditions[1]],
                             na.rm = T),
         log2_dCt = 2^ (- (Ct - control_mean)),
         control_mean_log = mean(log2_dCt[Condition == list_of_conditions[1]],
                              na.rm = T),
         Value_norm = log2_ddct / mock_mean_log)
         

gapdh_data_condition1 <- aggregate(gapdh_data_condition1[-1],
                           list(gapdh_data_condition1$Condition),
                           mean)
gapdh_data_condition1 <- gapdh_data_condition1 %>%
  arrange(match(Group.1, conditions))

gapdh_ratios_condition1 <- rep(gapdh_data_condition1$Value_norm, each=replicates)
gapdh_ratios_condition1 <- gapdh_ratios_condition1[file_emmit]

# final vector
gapdh_ratios_condition1
