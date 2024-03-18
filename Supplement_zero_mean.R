rm(list = ls())

setwd("[Enter main path]")

library(scoringRules)
library(dplyr)
library(tidyr)

# load functions
source("gdp_procs23.R")

flnm <- "inf_us" # either "inf_us", "gdp_us" or "gdp_de"
dat <- read.csv(paste0("data/", flnm, ".csv")) %>% na.omit
yrs <- dat %>% pull(target_year) %>% unique %>% sort

# initialize empty data frame
res_all <- data.frame()

# Loop over years (cross-validation)
for (yy in yrs){
  # test sample
  test <- dat %>% filter(target_year == yy)
  # training sample
  train <- dat %>% filter(target_year != yy) %>% 
    select(e, h) %>% na.omit
  # start data frame for this year
  res_tmp <- test %>% select(rlz, forecast, target_year, e, h) %>%
    mutate(pred1_l = NA, pred1_u = NA, score1 = NA, cov1 = NA,
           pred2_l = NA, pred2_u = NA, score2 = NA, cov2 = NA)
 
  # Gaussian regression
  fit1 <- gaussreg_pp(y_training = train$e, h_training = train$h, 
                      h_test = test$h, tau_gauss = c(.1, .9))
  res_tmp[, c("pred1_l", "pred1_u")] <- fit1$pred[,1:2]
  res_tmp$score1 <- ints(y = test$e, 
                         x_lower = fit1$pred[, 1],
                         x_upper = fit1$pred[, 2], target_coverage = .8)
  res_tmp$cov1 <- (test$e >= fit1$pred[,1]) & 
    (test$e <= fit1$pred[,2])
  
  # Gaussian regression (zero mean)
  fit2 <- gaussreg_pp(y_training = train$e, h_training = train$h, 
                      h_test = test$h, tau_gauss = c(.1, .9), 
                      zero_mean = TRUE)
  res_tmp[, c("pred2_l", "pred2_u")] <- fit2$pred[,1:2]
  res_tmp$score2 <- ints(y = test$e, 
                         x_lower = fit2$pred[, 1],
                         x_upper = fit2$pred[, 2], target_coverage = .8)
  res_tmp$cov2 <- (test$e >= fit2$pred[,1]) & 
    (test$e <= fit2$pred[,2])
  
  # expand data frame
  res_all <- rbind(res_all, res_tmp)
}

# Compute evaluation statistics
qs <- colMeans(res_all[, paste0("score", 1:2)], na.rm = TRUE)
qs_print <- print_helper(qs)
covg_print <- print_helper(100*
                             colMeans(res_all[, paste0("cov", 1:2)],
                                      na.rm = TRUE))
len_print <- res_all %>% transmute(len1 = pred1_u - pred1_l,
                                   len2 = pred2_u - pred2_l,
                                   score1) %>%
  na.omit %>%
  summarise(len1 = mean(len1), len2 = mean(len2)) %>% print_helper

meth_ind <- c(1:2) # indexes of methods (as presented in table)
# names of models 1-4
meth_names <- c("Gaussian", "Gaussian (zero mean)") 
# Make tex output for table
paste(meth_names[meth_ind], "&", 
      paste0(covg_print[meth_ind], "\\%"), "&", 
      len_print[meth_ind], "&", 
      qs_print[meth_ind], "\\\\") %>%
  writeLines
