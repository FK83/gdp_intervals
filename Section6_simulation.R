rm(list = ls())

setwd("[Enter main path]")

library(scoringRules)
library(isodistrreg)

# load functions
source("gdp_procs23.R")

# load Alexander Henzi's code for estimation under convex order constraint
source("estimation.R")

# helper function used below
cov_ind <- function(y, x_lower, x_upper){
  ((y >= x_lower) & (y <= x_upper))
}

# data frame holding all settings to be run
all_settings <- data.frame(n = c(500, 1000, 500), 
                           n_years_eval = c(30, 30, 60))

# get iteration (= relevant row number of all_settings; either 1, 2 or 3)
iter <- 3

# set random seed
set.seed(iter)

sim_start <- Sys.time()
n <- all_settings$n[iter]
n_years_eval <- all_settings$n_years_eval[iter] # nr of years used for evaluation
n_years_init <- 30 # initial burn-in period (to remove impact of prior)
n_years <- n_years_eval + n_years_init
n_mc <- 2e3

dat_all <- ar1_model(n_years = n_years, n_mc = n_mc, 
                     type_agg = "annual_average")
sel <- apply(dat_all$e, 1, function(z) !any(is.na(z))) &
  (dat_all$e[,1] > n_years_init)
e_all <- dat_all$e[sel,]
s_true_all <- dat_all$s[sel,]
res_mc_all <- matrix(NA, n_mc, 15)

for (jj in 1:n_mc){
  # print progress
  if ((jj %% 100) == 0){
    print(paste("Now running iteration", jj))
  }
  ind <- sample.int(nrow(e_all), size = n)
  dat <- e_all[ind, c(1:2, jj+2)] 
  s_true <- s_true_all[ind, ]
  yrs <- unique(dat[,1])
  scores1 <- scores2 <- scores3 <- scores4 <- scores5 <- 
    cov1 <- cov2 <- cov3 <- cov4 <- cov5 <- 
    len1 <- len2 <- len3 <- len4 <- len5 <- NULL
  # cross-validation for this MC iteration
  for (yy in yrs){
    # test sample
    test <- dat[dat[,1] == yy,,drop = FALSE]
    # training sample
    train <- dat[dat[,1] != yy,,drop = FALSE]
    # simple IDR
    fit1 <- isodistrreg_pp(y_training = train[,3], 
                           h_training = train[,2], simple = TRUE, 
                           h_test = test[,2], tau_idr = c(.1, .9))
    pred1 <- fit1$pred
    scores1 <- c(scores1, 
                 ints(y = test[, 3], x_lower = pred1[, 1],
                      x_upper = pred1[, 2], target_coverage = .8))
    cov1 <- c(cov1,
              cov_ind(y = test[,3], x_lower = pred1[, 1],
                      x_upper = pred1[, 2]))
    len1 <- c(len1, pred1[, 2] - pred1[, 1])
    # convex order IDR
    fit2 <- iso_icv_icx(y = train[,3], X = data.frame(x = train[,2]), 
                        concave = FALSE)
    aux2 <- predict(fit2, data = data.frame(x = test[,2]))
    pred2 <- qpred(aux2, quantiles = c(.1, .9))
    scores2 <- c(scores2, 
                 ints(y = test[, 3], x_lower = pred2[, 1],
                      x_upper = pred2[, 2], target_coverage = .8))
    cov2 <- c(cov2, 
              cov_ind(y = test[, 3], x_lower = pred2[, 1],
                      x_upper = pred2[, 2]))
    len2 <- c(len2, pred2[, 2] - pred2[, 1])
    # gaussian regression
    fit3 <- gaussreg_pp(y_training = train[,3], h_training = train[,2], 
                        h_test = test[,2], tau_gauss = c(.1, .9))
    pred3 <- fit3$pred
    scores3 <- c(scores3, 
                 ints(y = test[, 3], x_lower = pred3[, 1],
                      x_upper = pred3[, 2], target_coverage = .8))
    cov3 <- c(cov3, 
              cov_ind(y = test[,3], x_lower = pred3[, 1],
                      x_upper = pred3[, 2]))
    len3 <- c(len3, pred3[, 2] - pred3[, 1])
    # true model
    tmp4 <- s_true[dat[,1] == yy,3]
    pred4 <- cbind(qnorm(.1, sd = tmp4), qnorm(.9, sd = tmp4))
    scores4 <- c(scores4, 
                 ints(y = test[, 3], x_lower = pred4[, 1], 
                      x_upper = pred4[, 2], target_coverage = .8))
    cov4 <- c(cov4, 
              cov_ind(y = test[, 3], x_lower = pred4[, 1],
                      x_upper = pred4[, 2]))
    len4 <- c(len4, pred4[, 2] - pred4[, 1])
    # combination of models 1-3
    pred5 <- 1/3*(pred1[,1:2, drop = FALSE] + pred2[,1:2, drop = FALSE] + 
                    pred3[,1:2, drop = FALSE])
    scores5 <- c(scores5, 
                 ints(y = test[,3], x_lower = pred5[, 1], 
                      x_upper = pred5[, 2], target_coverage = .8))
    cov5 <- c(cov5, 
              cov_ind(y = test[,3], x_lower = pred5[, 1], 
                      x_upper = pred5[, 2]))
    len5 <- c(len5, pred5[, 2] - pred5[, 1])
  }
  res_mc_all[jj, 1:5] <- cbind(scores1, scores2, scores3, scores4, scores5) |> 
    colMeans()
  res_mc_all[jj, 6:10] <- cbind(cov1, cov2, cov3, cov4, cov5) |> 
    colMeans()
  res_mc_all[jj, 11:15] <- cbind(len1, len2, len3, len4, len5) |> 
    colMeans()
}

boxplot(res_mc_all)
sim_end <- Sys.time()

# make filename (based on time at which simulation was finished)
sv_nm <- sim_end |> as.character() |>
  (function(x) substr(x, 1, 16))() |>
  (function(x) gsub(" ", "_", x))() |>
  (function(x) paste0("sim_", x, ".RData"))()

save.image(file = paste0("fe_sim/", sv_nm))
