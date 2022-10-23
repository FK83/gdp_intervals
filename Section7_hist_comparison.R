rm(list = ls())

setwd("[ENTER PATH HERE]")

library(dplyr)
library(tidyr)
library(colorspace)
library(ggplot2)
library(murphydiagram)

source("quantile_procs.R")

var <- "gdp" 	# choose variable (either "gdp" or "inf")
model_nr <- "0" # choose model to be compared to histograms ("0" yields AR1 baseline as used in the paper)
if (model_nr == "4"){
  lb <- expression("QR (" * theta * " estimated)")
} else if (model_nr == "3"){
  lb <- expression("QR (" * theta * " fixed)")
} else if (model_nr == "2"){
  lb <- expression("Gauss (" * theta * " estimated)")
} else if (model_nr == "1"){
  lb <- expression("Gauss (" * theta * " fixed)")
} else if (model_nr == "0"){
  lb <- "AR1"
}

# load quantile regression and histogram forecasts
dat1 <- read.csv(paste0("evaluation_US_", var, ".csv")) %>%
  select(target_year, h, forecast, rlz, contains(model_nr))
# Remove model nr from variable names (easier to access below)
names(dat1) <- gsub(model_nr, "", names(dat1))
dat2 <- read.csv(paste0("forecast_histograms_", var, ".csv"))

# merge data sets (to get set of common forecast cases)
dat <- merge(dat1, dat2)
dat$hist_score <- quantile_eval(as.matrix(dat[, c("hist_lower", "hist_upper")]), 
                                dat$rlz, c(.1, .9))

# coverage rates of both forecasts
hist_cvg <- mean( (dat$rlz >= dat$hist_lower) & (dat$rlz <= dat$hist_upper) )
qr_cvg <- mean( (dat$rlz >= dat$pred_lower) & (dat$rlz <= dat$pred_upper) )
round(100*hist_cvg, 2)
round(100*qr_cvg, 2)

# correlation
cor(dat$pred_lower, dat$hist_lower)
cor(dat$pred_upper, dat$hist_upper)
cor(dat$forecast, dat$hist_fcst)

# scores
# multiply scores by 10 to get common scaling of interval score 
colMeans(10*dat[, c("hist_score", "score")]) %>% (function(x) round(x, 2))
t1 <- t.test(dat$hist_score-dat$score)
t2 <- dm_cluster(dat$hist_score, dat$score, dat$target_year)
round(c(t1$statistic, t2$t), 2)

# sharpness
dat <- dat %>% mutate(qr_length = pred_upper - pred_lower, 
                      hist_length = hist_upper - hist_lower)
dat %>% summarise(qr_length = mean(qr_length), hist_length = mean(hist_length)) %>%
  (function(x) round(x, 2))

# plot length of prediction interval against forecast horizon
dat %>% select(target_year, h, hist = hist_length, qr = qr_length) %>%
  pivot_longer(c("hist", "qr")) %>% 
  ggplot(aes(x = h, y = value, color = name)) + 
  geom_point(size = I(2), alpha = .6) + 
  theme_minimal(base_size = 14) + xlab("Horizon") + ylab("Length of PI") + 
  theme(legend.position = "top") + 
  scale_color_discrete_qualitative(name = "", palette = "Dark 3", 
                                   breaks = c("hist", "qr"), 
                                   labels = c("SPF Histogram", lb)) 
# save plot as pdf file
ggsave(paste0("hist_dots_", var, ".pdf"))
