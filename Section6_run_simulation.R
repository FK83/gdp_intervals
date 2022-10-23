rm(list = ls())

setwd("[ENTER PATH HERE]")

library(quantreg)
library(dplyr)
library(crch)
library(murphydiagram)
library(scoringRules)
library(foreach)
library(doParallel)

source("quantile_procs.R")

# get counter for current job
job_id <- 1 # when using a cluster, can be set to something like as.integer(as.character(Sys.getenv("SLURM_ARRAY_TASK_ID")))

# nr of cores
n_cores <- 4 # use more cores if available (e.g. on cluster)

# settings for simulation
h_max <- 24 # maximal forecast horizon
a <- .9 # ar parameter
s2 <- .1 # error term variance
quantile_levels <- c(.1, .9) # quantile levels (set to c(.1, .9) throughout paper)
include_kinks <- TRUE # whether to include models with estimated theta

# nr of months to be simulated
nr_y <- 5 + 20 # include "pre-sample" period of five years. Additional nr of years set to either 20 or 40 in paper.
nr_m <- nr_y*12
aux_ind <- 1:nr_m

# nr of mc iterations
n_mc <- 20

# sample size (= nr of forecast errors) in each iteration 
n <- 300 # either 300 or 600 in paper

# set random seed
set.seed(job_id)

# register Cores
registerDoParallel(n_cores)

# loop over MC iterations
simulation_results <- foreach(jj = 1:n_mc, .combine = rbind) %dopar% {
  
  sel <- data.frame(year = sample(6:nr_y, n, replace = TRUE), 
                    h = sample.int(h_max, size = n, replace = TRUE), 
                    e = NA)
  
  # simulated monthly history
  sim <- data.frame(value = arima.sim(nr_m, model = list(ar = a), sd = sqrt(s2)), 
                    year = rep(1:nr_y, each = 12), 
                    month = rep(1:12, nr_y))
  
  # simulated yearly observations
  y_year <- sim %>% group_by(year) %>% summarise(value = sum(value)) %>%
    ungroup

  # loop to make forecast error observations for years and horizons in 'sel'
  for (ll in 1:n){
    
    # index of forecast origin date
    t0 <- sel$year[ll]*12-sel$h[ll]
    
    # nr of known months from year
    known_m <- max(c(12-sel$h[ll], 0))
    
    # index of known months
    t_known <- (t0 - known_m + 1):t0
    
    # latest available observation
    latest_y <- sim$value[t0]
    
    # get forecast
    if ((t0 > 0) & (min(t_known) > 0)){
      if (sel$h[ll] >= 12){
        # all constituents are unknown
        fc <- (latest_y * (a^(1:sel$h[ll]))) %>% tail(12) %>% sum
      } else {
        # some constituents are known
        fc_aux <- c(sim$value[t_known], (latest_y * (a^(1:sel$h[ll]))))
        if (length(fc_aux) != 12){
          stop()
        } else {
          fc <- sum(fc_aux)
        }
      }
      # get outcome and forecast error
      rlz <- y_year$value[y_year$year == sel$year[ll]]
      sel$e[ll] <- rlz - fc
    }
  }
  
  # Drop NA's, change format of h
  sel <- sel %>% na.omit %>% mutate(h = as.numeric(h))
  
  # Run cross-validation procedure 
  yrs <- unique(sel$year) %>% sort
  
  # Loop over years (cross-validation; different year yy is left out in each iteration)
  # initialize data frame (to be appended year by year)
  eval_df <- data.frame()
  
  for (yy in yrs){
    
    # "training" and "test" set for this iteration
    dat_train <- sel %>% filter(year != yy) 
    dat_test <- sel %>% filter(year == yy)
    
    if (nrow(dat_test) == 0) stop("No test sample available")
    
    # model fitting
    
    # alphasigma model
    bench_fixed <- 
      fit_alphasigma(dat_train, dat_test, quantile_levels = quantile_levels, 
                     params = c(a, sqrt(s2)))
    
    bench_est <- fit_alphasigma(dat_train, dat_test, quantile_levels = quantile_levels, 
                                params = NULL)
    
    # crch, fixed theta
    model1 <- fit_fixed_theta(dat_train, dat_test, theta = 12, 
                              quantile_levels = quantile_levels, 
                              model = "crch")
    
    # initialize empty data frames for kink models
    crch_theta <- rq_theta <- data.frame(best_theta = NA)
    model2 <- model4 <- data.frame(loss = NA)
    if (include_kinks){
      # crch, optimized theta
      crch_theta <- find_best_theta(dat_train, quantile_levels = quantile_levels, 
                                    model = "crch") %>% try(silent = TRUE)
      if (class(crch_theta) == "try-error"){
        crch_theta <- data.frame(best_theta = NA)
      } else {
        # fit crch with selected theta
        model2 <- fit_fixed_theta(dat_train, dat_test, theta = crch_theta$best_theta, 
                                  quantile_levels = quantile_levels, 
                                  model = "crch")
      }
      # rq, optimized theta
      rq_theta <- find_best_theta(dat_train, quantile_levels = quantile_levels, 
                                  model = "rq") %>% try(silent = TRUE)
      if (class(crch_theta) == "try-error"){
        rq_theta <- data.frame(best_theta = NA)
      } else {
        # fit rq with selected theta
        model4 <- fit_fixed_theta(dat_train, dat_test, theta = rq_theta$best_theta, 
                                  quantile_levels = quantile_levels, 
                                  model = "rq")
      }
    }

    # rq, fixed theta
    model3 <- fit_fixed_theta(dat_train, dat_test, theta = 12, 
                              quantile_levels = quantile_levels, 
                              model = "rq")
    
    # expand evaluation data frame
    eval_df <- rbind(eval_df, 
                     data.frame(iter = jj, dat_test[, c("year", "h")],
                                alpha_est = bench_est$theta[1], 
                                s2_est = bench_est$theta[2]^2,
                                scoreb_fixed = bench_fixed$loss, 
                                scoreb_est = bench_est$loss, 
                                score1 = model1$loss, 
                                score2 = model2$loss, score3 = model3$loss, 
                                score4 = model4$loss,
                                theta2 = crch_theta$best_theta, 
                                theta4 = rq_theta$best_theta))
  }
  eval_df
}

# save results
save.image(file = paste0("fe_sim/mc_n", n, "_id", job_id, "_a", 10*a, ".RData"))
