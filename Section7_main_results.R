rm(list = ls())

setwd("[Enter main path]")

library(scoringRules)
library(dplyr)
library(tidyr)
library(colorspace)
library(ggplot2)
source("gdp_procs23.R")
library(isodistrreg)
# use Alexander Henzi's code for estimation under convex order constraint
source("estimation.R")

flnm <- "inf_us"
dat <- read.csv(paste0("data/", flnm, ".csv")) %>% na.omit
yrs <- dat %>% pull(target_year) %>% unique %>% sort

res_all <- data.frame()

for (yy in yrs){
  # test sample
  test <- dat %>% filter(target_year == yy)
  # training sample
  train <- dat %>% filter(target_year != yy) %>% 
    select(e, h) %>% na.omit
  # start data frame for this year
  res_tmp <- test %>% select(rlz, forecast, target_year, e, h) %>%
    mutate(pred1_l = NA, pred1_u = NA, score1 = NA, cov1 = NA,
           pred2_l = NA, pred2_u = NA, score2 = NA, cov2 = NA, 
           pred3_l = NA, pred3_u = NA, score3 = NA, cov3 = NA,
           pred4_l = NA, pred4_u = NA, score4 = NA, cov4 = NA)
  # decomposition method
  fit1 <- isodistrreg_pp(y_training = train$e, 
                         h_training = train$h, simple = TRUE, 
                         h_test = test$h, tau_idr = c(.1, .9))
  res_tmp[, c("pred1_l", "pred1_u")] <- fit1$pred[,1:2]
  res_tmp$score1 <- ints(y = test$e, x_lower = fit1$pred[, 1],
                         x_upper = fit1$pred[, 2], 
                         target_coverage = .8)
  res_tmp$cov1 <- (test$e >= fit1$pred[,1]) & 
    (test$e <= fit1$pred[,2])

  # flexible method
  fit2 <- iso_icv_icx(y = train$e, X = data.frame(x = train$h), 
                      concave = FALSE)
  aux2 <- predict(fit2, data = data.frame(x = test$h))
  pred2 <- res_tmp[, c("pred2_l", "pred2_u")] <- 
    qpred(aux2, quantiles = c(.1, .9))
  res_tmp$score2 <- ints(y = test$e, x_lower = pred2[, 1],
                         x_upper = pred2[, 2], target_coverage = .8)
  res_tmp$cov2 <- (test$e >= pred2[,1]) & 
    (test$e <= pred2[,2])
  
  # Gaussian regression
  fit3 <- gaussreg_pp(y_training = train$e, h_training = train$h, 
                      h_test = test$h, tau_gauss = c(.1, .9))
  res_tmp[, c("pred3_l", "pred3_u")] <- fit3$pred[,1:2]
  res_tmp$score3 <- ints(y = test$e, 
                         x_lower = fit3$pred[, 1],
                         x_upper = fit3$pred[, 2], target_coverage = .8)
  res_tmp$cov3 <- (test$e >= fit3$pred[,1]) & 
    (test$e <= fit3$pred[,2])
  
  # quantile combination
  pred_ens <- (fit1$pred[,1:2] + pred2 + fit3$pred[,1:2])/3
  res_tmp[, c("pred4_l", "pred4_u")] <- pred_ens
  res_tmp$score4 <- ints(y = test$e, 
                         x_lower = pred_ens[, 1],
                         x_upper = pred_ens[, 2], target_coverage = .8)
  res_tmp$cov4 <- (test$e >= pred_ens[,1]) & 
    (test$e <= pred_ens[,2])
  
  # expand data frame
  res_all <- rbind(res_all, res_tmp)
}

# compute evaluation statistics
qs <- colMeans(res_all[, paste0("score", 1:4)], na.rm = TRUE)
qs_print <- print_helper(qs)
covg_print <- print_helper(100*
                             colMeans(res_all[, paste0("cov", 1:4)], 
                                      na.rm = TRUE))
len_print <- res_all %>% transmute(len1 = pred1_u - pred1_l,
                                   len2 = pred2_u - pred2_l, 
                                   len3 = pred3_u - pred3_l, 
                                   len4 = pred4_u - pred4_l,
                                   score1) %>%
  na.omit %>%
  summarise(len1 = mean(len1), len2 = mean(len2), 
            len3 = mean(len3), len4 = mean(len4)) %>% print_helper

meth_ind <- c(3, 1, 2, 4) # indexes of methods (as presented in table)
# names of models 1-4
meth_names <- c("Decomposition", "Flexible", "Gaussian", "Combination") 
# make tex output for Table 3 in paper
paste(meth_names[meth_ind], "&", 
      paste0(covg_print[meth_ind], "\\%"), "&", 
      len_print[meth_ind], "&", 
      qs_print[meth_ind], "\\\\") %>%
  writeLines

# plot coverage against horizon (Figure S2 in online supplement)
for (mm in 1:4){
  fit <- loess(formula(paste0("cov", mm, "~ h")), data = res_all, 
               degree = 1)
  pred_df <- data.frame(h = seq(from = 0, to = 104, by = .5))
  pred_df$pred <- predict(fit, newdata = pred_df)
  pdf(paste0("plots/", flnm, "_calib", mm, ".pdf"))
  par(mar = c(5.1, 6.1, 4.1, 2.1))
  plot(pred_df$h, pred_df$pred, type = "l", ylim = c(0, 1), bty = "n", 
       xlab = "Horizon (weeks)", ylab = "Coverage", lwd = 4, 
       cex.axis = 2, cex.lab = 2)
  points(jitter(res_all$h), res_all[,paste0("cov", mm)], col = grey(.1, .1), 
         pch = 20, cex = 2)
  abline(h = .8, lty = 2, lwd = 2)
  dev.off()
}

# illustration plot (model excluding 2020)
# (Figure 5 in paper)
meth_sel <- 3:4 # choose methods for quantile predictions shown in plot
lab_sel <- meth_names[meth_sel]
pch_2020 <- 17 # symbol type for 2020 and other years
pch_other <- 20
col_2020 <- "black" # color for 2020 and other years
col_other <- grey(.6, .6)
plot_df <- res_all %>% select(target_year, e, h, 
                              contains(c(paste0("pred", meth_sel[1]),
                                         paste0("pred", meth_sel[2])))) %>%
  mutate(target_year_2020 = (target_year == 2020), 
         aux_pch = if_else(target_year_2020, pch_2020, pch_other), 
         aux_col = if_else(target_year_2020, col_2020, col_other))
# Print plot
pdf(paste0("plots/2020errors_", flnm, ".pdf"), 
    height = 6, width = 16)
par(mar=c(5,5,4,1)+.1)
plot(plot_df$h, plot_df$e, pch = plot_df$aux_pch, bty = "n", 
     col = plot_df$aux_col, xlab = "Horizon (weeks)", 
     ylab = "Forecast Error", cex.lab = 1.6, cex.axis = 1.6, cex = 1.6, 
     ylim = c(-7, 7), xlim = c(0, 104))
plot_df_2 <- plot_df %>% filter(target_year_2020) %>% arrange(h)
matplot(plot_df_2$h, plot_df_2[,grepl("pred", names(plot_df))], add = TRUE, 
        col = 1, type = "l", lty = c(1, 1, 2, 2), lwd = 3)
dev.off()

# construct legend as separate plot
# see https://stackoverflow.com/questions/48966645/how-can-i-create-a-legend-without-a-plot-in-r
pdf(paste0("plots/2020errors_legend.pdf"), 
    height = 6, width = 16)
plot(NULL, xaxt = "n", yaxt = "n", bty = "n", ylab = "", xlab = "", 
     ylim = c(0, 1), xlim = c(0, 1))
legend("center", c("2020", lab_sel[1], "other years", lab_sel[2]), 
       lty = c(NA, 1, NA, 2), col = c(col_2020, 1, col_other, 1), 
       lwd = c(NA, 3, NA, 3), ncol = 2, bty = "n", 
       pch = c(pch_2020, NA, pch_other, NA), 
       cex = 4)
dev.off()

# save results as csv file
write.csv(res_all, paste0("data/", flnm, "_eval.csv"), row.names = FALSE)