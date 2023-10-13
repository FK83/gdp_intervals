rm(list = ls())

setwd("[Enter main path]")

library(scoringRules)
library(dplyr)
library(tidyr)
library(colorspace)
library(ggplot2)
library(isodistrreg)

# load functions
source("gdp_procs23.R")

# load Alexander Henzi's code for estimation under convex order constraint
source("estimation.R")

flnm <- "inf_us" # either "inf_us", "gdp_us" or "gdp_de"
dat <- read.csv(paste0("data/", flnm, ".csv")) %>% na.omit
yrs <- dat %>% pull(target_year) %>% unique %>% sort

# initialize empty data frame
res_all <- data.frame()

# loop over years (cross-validation)
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
  # Decomposition method
  fit1 <- isodistrreg_pp(y_training = train$e, 
                         h_training = train$h, simple = TRUE, 
                         h_test = test$h, tau_idr = c(.1, .9))
  res_tmp[, c("pred1_l", "pred1_u")] <- fit1$pred[,1:2]
  res_tmp$score1 <- ints(y = test$e, x_lower = fit1$pred[, 1],
                         x_upper = fit1$pred[, 2], 
                         target_coverage = .8)
  res_tmp$cov1 <- (test$e >= fit1$pred[,1]) & 
    (test$e <= fit1$pred[,2])

  # Flexible method
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
  
  # Quantile combination
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

# Compute evaluation statistics
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
# Make tex output for table
paste(meth_names[meth_ind], "&", 
      paste0(covg_print[meth_ind], "\\%"), "&", 
      len_print[meth_ind], "&", 
      qs_print[meth_ind], "\\\\") %>%
  writeLines

# plot coverage against horizon
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
meth_sel <- 3:4
lab_sel <- meth_names[meth_sel]
plot_df <- res_all %>% select(target_year, e, h, 
                              contains(c(paste0("pred", meth_sel[1]),
                                         paste0("pred", meth_sel[2])))) %>%
  pivot_longer(cols = contains("pred")) %>%
  mutate(target_year = if_else(target_year == 2020, "2020", "other"))
plot_df_2 <- plot_df %>% filter(target_year == "2020") 
cl <- c("blue4", "brown4")
ggplot(filter(plot_df, target_year != "2020"), 
       aes(x = h, y = e)) + 
  geom_point(col = grey(.1, .1), size = I(3)) + 
  theme_minimal(base_size = 20) + 
  geom_line(data = plot_df_2, mapping = aes(x = h, y = value, 
                                            color = name, linetype = name), 
            lwd = I(1.5)) +
  geom_point(plot_df_2, mapping = aes(x = h, y = e), 
             color = "green4", size = I(6), pch = 20) +
  theme(legend.position = "top") + xlab("Horizon") + 
  scale_color_manual(name = "", 
                     breaks = c(paste0("pred", meth_sel[1], c("_l", "_u")), 
                                paste0("pred", meth_sel[2], c("_l", "_u"))),
                     values = rep(cl, each = 2)) + 
  scale_linetype_manual(name = "", 
                        breaks = c(paste0("pred", meth_sel[1], c("_l", "_u")), 
                                   paste0("pred", meth_sel[2], c("_l", "_u"))),
                        values = rep(1:2, each = 2)) + 
  guides(color = "none", linetype = "none")
ggsave(filename = paste0("plots/2020errors_", flnm, ".pdf"), 
       height = 6, width = 16)

# print legend as separate file
pdf("plots/2020errors_legend.pdf", height = 4, width = 16)
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("top", c("2020", lab_sel[1], "other years", lab_sel[2]), 
       lty = c(NA, 1, NA, 2), lwd = c(NA, 7, NA, 7), 
       pch = c(16, NA, 16, NA), bty = "n", cex = 4, pt.cex = c(7, NA, 7, NA),
       col = c("green4", cl[1], "grey", cl[2]), ncol = 2)
dev.off()

# write evaluation results to csv file
write.csv(res_all, paste0(flnm, "_eval.csv"), row.names = FALSE)
