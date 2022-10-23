rm(list = ls())
setwd("[ENTER PATH HERE]")

library(quantreg)
library(dplyr)
library(crch)
library(murphydiagram)
library(scoringRules)
library(ggplot2)
library(tidyr)
library(colorspace)
library(sandwich)

# Settings

# Choice of data set
data_choice <- "US_inf" # "DE" for German data, "US_gdp" and "US_inf" for US data
# Quantile levels 
quantile_levels <- c(.1, .9) # Equal to c(.1, .9) throughout the paper
# Whether or not to use ceiling of forecast horizon variable. 
# Equal to FALSE in baseline setting
ceil_horizon <- FALSE

# Load data
if (data_choice == "DE"){
  dat <- read.csv("gdp_forecasts_DE.csv") %>%
    mutate(e = truth_preliminary - forecast, h = horizon, 
           rlz = truth_preliminary) %>%
    select(-horizon)
} else if (data_choice == "US_gdp"){
  dat <- read.csv("gdp_forecasts_US.csv") %>% na.omit
} else if (data_choice == "US_inf"){
  dat <- read.csv("inf_forecasts_US.csv") %>% na.omit
}

# Optionally: Use ceiling of horizon (so that there are integer horizons only)
if (ceil_horizon){
  dat <- dat %>% mutate(h = ceiling(h))
}

# Initialize stuff
yrs <- unique(dat$target_year)
eval_df <- data.frame()
theta_choices <- data.frame(yy = yrs, model2 = NA, model4 = NA)

# Load R functions
source("quantile_procs.R")

# Loop over years (cross-validation; different test-sample year yy in each iteration)
for (yy in yrs){
  
  # "training" and "test" set for this iteration
  dat_train <- dat %>% filter(target_year != yy)
  dat_test <- dat %>% filter(target_year == yy)
  
  # model fitting
  
  # alphasigma model
  model0 <- fit_alphasigma(dat_train, dat_test, quantile_levels = quantile_levels)
  
  # crch, fixed theta
  model1 <- fit_fixed_theta(dat_train, dat_test, theta = 12, 
                            quantile_levels = quantile_levels, 
                            model = "crch")
  
  # crch, optimized theta
  crch_theta <- find_best_theta(dat_train, quantile_levels = quantile_levels, 
                                model = "crch")
  model2 <- fit_fixed_theta(dat_train, dat_test, theta = crch_theta$best_theta, 
                            quantile_levels = quantile_levels, 
                            model = "crch")
  theta_choices$model2[which(yrs == yy)] <- crch_theta$best_theta
  
  # rq, fixed theta
  model3 <- fit_fixed_theta(dat_train, dat_test, theta = 12, 
                            quantile_levels = quantile_levels, 
                            model = "rq")
  
  # rq, optimized theta
  rq_theta <- find_best_theta(dat_train, quantile_levels = quantile_levels, 
                              model = "rq")
  model4 <- fit_fixed_theta(dat_train, dat_test, theta = rq_theta$best_theta, 
                            quantile_levels = quantile_levels, 
                            model = "rq")
  theta_choices$model4[which(yrs == yy)] <- rq_theta$best_theta
  
  if (yy == 2020){
    model0_2020 <- model0
    model2_2020 <- model2
    model4_2020 <- model4
  }
  
  # model evaluation
  
  # put point forecasts into matrix format
  aux_mat <- matrix(dat_test$forecast, nrow = nrow(dat_test), 
                    ncol = length(quantile_levels))
  
  # add point forecasts to predicted quantiles of forecast error
  pred0 <- aux_mat + model0$pred
  pred1 <- aux_mat + model1$pred
  pred2 <- aux_mat + model2$pred
  pred3 <- aux_mat + model3$pred
  pred4 <- aux_mat + model4$pred

  # expand evaluation data frame
  eval_df <- rbind(eval_df, 
                   data.frame(dat_test[, c("target_year", "h", 
                                           "forecast", "rlz")], 
                              score0 = model0$loss, score1 = model1$loss, 
                              score2 = model2$loss, score3 = model3$loss, 
                              score4 = model4$loss,
                              pred0_lower = pred0[,1],
                              pred0_upper = pred0[,length(quantile_levels)],
                              pred1_lower = pred1[,1],
                              pred1_upper = pred1[,length(quantile_levels)],
                              pred2_lower = pred2[,1], 
                              pred2_upper = pred2[,length(quantile_levels)],
                              pred3_lower = pred3[,1],
                              pred3_upper = pred3[,length(quantile_levels)],
                              pred4_lower = pred4[,1],
                              pred4_upper = pred4[,length(quantile_levels)]))
}

# Define helper function for printing
print_helper <- function(x, nd = 2){
  format(round(x, nd), nsmall = nd)
}

# Plot quantile scores (in the order appearing in the paper)
# Multiply all scores by 2/(interval coverage level) to get the standard scaling of the 
# interval score
qs <- colMeans(eval_df[, paste0("score", c(0, 2, 1, 4, 3))])*(2/(2*quantile_levels[1]))
qs_print <- print_helper(qs, nd = 2)

# Compute coverage rates
covg_df <- eval_df %>% mutate(truth = rlz,
                              covg0 = (truth >= pred0_lower) & (truth <= pred0_upper),
                              covg1 = (truth >= pred1_lower) & (truth <= pred1_upper),
                              covg2 = (truth >= pred2_lower) & (truth <= pred2_upper), 
                              covg3 = (truth >= pred3_lower) & (truth <= pred3_upper),
                              covg4 = (truth >= pred4_lower) & (truth <= pred4_upper)) 
covg_print <- covg_df %>%
  summarise(covg0 = mean(covg0), covg2 = mean(covg2), covg1 = mean(covg1), 
            covg4 = mean(covg4), covg3 = mean(covg3), ) %>%
  (function(x) print_helper(100*x))

# Compute mean length of prediction intervals
len_print <- eval_df %>%
  summarise(l0 = mean(pred0_upper-pred0_lower), 
            l2 = mean(pred2_upper-pred2_lower), 
            l1 = mean(pred1_upper-pred1_lower), 
            l4 = mean(pred4_upper-pred4_lower),
            l3 = mean(pred3_upper-pred3_lower)) %>%
  print_helper

# Make tex output for table
nms <- c("AR1", "Gauss ($\\theta$ estimated)", "Gauss ($\\theta = 12$)", 
         "QR ($\\theta$ estimated)", "QR ($\\theta = 12$)")
paste(nms, "&", paste0(covg_print, "\\%"), "&", len_print, "&", qs_print, "\\\\") %>%
  writeLines

# Compute Diebold-Mariano tests (AR1 as benchmark; different variants of standard errors)
# `cluster' and `cluster2' should yield the same results (different implementations of the same estimator)
dms <- data.frame(standard = rep(NA, 4), cluster = rep(NA, 4), 
                  cluster2 = rep(NA, 4))
loss1 <- eval_df$score0
cluster_var <- eval_df$target_year
inds <- c(2, 1, 4, 3)
for (jj in 1:4){
  loss2 <- eval_df[, paste0("score", inds[jj])]
  dms$standard[jj] <- t.test(loss2, loss1)$statistic
  dms$cluster[jj] <- dm_cluster(loss2, loss1, cluster_var)$t
  dms$cluster2[jj] <- dm_cluster(loss2, loss1, cluster_var, type = "sandwich")$t
}
# print Diebold-Mariano p-values
dms[,1:2] %>% print_helper %>% 
  (function(x) apply(x, 1, function(z) paste(z, collapse = " & ")))

# print range of estimates for theta
apply(theta_choices[,-1],2, range) %>% unname %>% 
  (function(x) apply(x, 2, function(z) paste(z, collapse = "-")))

# Plot forecasts of two illustrative models
# First get predicted quantiles of forecast error distribution
plot_df <- eval_df %>% mutate(e = rlz - forecast, 
                              pred2_lower_e = pred2_lower - forecast, 
                              pred2_upper_e = pred2_upper - forecast,
                              pred4_lower_e = pred4_lower - forecast,
                              pred4_upper_e = pred4_upper - forecast) %>%
  pivot_longer(cols = c(contains("pred"))) %>%
  filter(grepl("_e", name)) %>% mutate(pred_type = substr(name, 1, 5)) %>%
  mutate(target_year == as.character(target_year)) %>%
  mutate(target_year = if_else(target_year == "2020", 
                               "2020", "other"))
plot_df_2 <- plot_df %>% filter(target_year == "2020") 

# Plot fitted quantiles (use models estimated without 2020 data)
cl <- rainbow_hcl(n = 2)
ggplot(plot_df, aes(x = h, y = e, shape = target_year)) + 
  geom_point(col = grey(.1, .1)) + theme_minimal(base_size = 11) + 
  geom_line(data = plot_df_2, mapping = aes(x = h, y = value, color = name)) +
  theme(legend.position = "top") + xlab("Horizon") + 
  scale_shape_manual(name = "Target Year", breaks = c("2020", "other"), 
                     values = c(12, 20)) + ylab("Forecast Error") + 
  scale_color_manual(name = "", breaks = c("pred2_lower_e", "pred2_upper_e", 
                                           "pred4_lower_e", "pred4_upper_e"),
                     values = rep(cl, each = 2)) + 
  guides(color = "none") +
  annotate("text", x = 20, y = 2.5, 
           label = as.character(expression("Gauss (" * theta * " est.)")),
           parse = TRUE, color = cl[1], size = I(2.7)) + 
  annotate("text", x = 22, y = -4, 
           label = as.character(expression("QR (" * theta * " est.)")),
           parse = TRUE,
           color = cl[2], size = I(2.7)) 

# If in baseline setting: Save results to disk
if (!ceil_horizon){
  # Save plot in pdf format
  ggsave(filename = paste0("2020errors_", data_choice, ".pdf"))
  
  # Save evaluation data frame as csv file
  write.table(eval_df, paste0("evaluation_", data_choice, ".csv"), sep = ",",
              col.names = TRUE, row.names = FALSE)  
} else {
  print("Not in baseline setting - results not printed to disk")
}
