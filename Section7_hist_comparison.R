rm(list = ls())

setwd("[Enter main path]")

library(dplyr)
library(tidyr)
library(colorspace)
library(ggplot2)
library(murphydiagram)
library(forecastcalibration)

source("gdp_procs23.R")

var <- "gdp" 	# choose variable (either "gdp" or "inf")
model_nr <- "4" # choose model to be compared to histograms ("4" yields ensemble baseline as used in the paper)
if (model_nr == "1"){
  lb <- expression("Decomposition")
} else if (model_nr == "2"){
  lb <- expression("Flexible")
} else if (model_nr == "3"){
  lb <- expression("Gaussian")
} else if (model_nr == "4"){
  lb <- "Comb."
}
cov_ind <- function(y, x_lower, x_upper){
  ((y >= x_lower) & (y <= x_upper))
}

# load postprocessing and histogram forecasts
dat1 <- read.csv(paste0("data/", var, "_us_eval.csv")) %>%
  select(target_year, h, forecast, rlz, contains(model_nr))
# Remove model nr from variable names (easier to access below)
names(dat1) <- gsub(model_nr, "", names(dat1))
dat2 <- read.csv(paste0("data/histograms_", var, ".csv")) %>%
  transmute(target_year, h, hist1_lower = hist_lower, hist1_upper = hist_upper)
dat3 <- read.csv(paste0("data/individual_histograms_", var, ".csv")) %>%
  na.omit %>%
  group_by(target_year, h) %>% 
  summarise(hist2_lower = mean(hist_lower), hist2_upper = mean(hist_upper), 
            hist3_lower = median(hist_lower), hist3_upper = median(hist_upper))

# merge data sets (to get set of common forecast cases)
dat <- merge(dat1, merge(dat2, dat3)) %>% na.omit %>%
  mutate(pp_lower = forecast + pred_l, 
         pp_upper = forecast + pred_u) %>%
  mutate(hist1_score = ints(y = rlz, x_lower = hist1_lower, 
                            x_upper = hist1_upper,
                            target_coverage = .8),
         hist2_score = ints(y = rlz, x_lower = hist2_lower, 
                            x_upper = hist2_upper, 
                            target_coverage = .8), 
         hist3_score = ints(y = rlz, x_lower = hist3_lower, 
                            x_upper = hist3_upper, 
                            target_coverage = .8),
         hist1_cov = cov_ind(y = rlz, x_lower = hist1_lower, 
                             x_upper = hist1_upper), 
         hist2_cov = cov_ind(y = rlz, x_lower = hist2_lower, 
                             x_upper = hist2_upper), 
         hist3_cov = cov_ind(y = rlz, x_lower = hist3_lower, 
                             x_upper = hist3_upper))

# coverage rates of both forecasts
round(100*colMeans(dat[,c("cov", paste0("hist", 1:3, "_cov"))]))

# correlation
cor(dat[, grepl("lower", names(dat))])
summary(dat[, grepl("lower", names(dat))])
cor(dat[, grepl("upper", names(dat))])
summary(dat[, grepl("upper", names(dat))])

# scores
summary(dat[, c("score", paste0("hist", 1:3, "_score"))])
scores <- colMeans(dat[, c("score", paste0("hist", 1:3, "_score"))]) %>%
  unname

# double-check scores for postprocessing method
# (now based on outcomes, rather than forecast errors; should be identical)
dat$score_check <- ints(y = dat$rlz, x_lower = dat$pp_lower, 
                        x_upper = dat$pp_upper, 
                        target_coverage = .8)
all.equal(dat$score, dat$score_check)

# sharpness
dat <- dat %>% mutate(pred_length = pred_u - pred_l, 
                      hist1_length = hist1_upper - hist1_lower, 
                      hist2_length = hist2_upper - hist2_lower, 
                      hist3_length = hist3_upper - hist3_lower)
len <- dat %>% summarise(pred_length = mean(pred_length),
                         hist1_length = mean(hist1_length),
                         hist2_length = mean(hist2_length),
                         hist3_length = mean(hist3_length)) %>%
  unname

# tex code for table
paste(c("SPF Histogram", "SPF Histogram (quantile mean)", 
        "SPF Histogram (quantile median)"), "&",
      paste0(print_helper(unname(100*colMeans(dat[,paste0("hist", 1:3, "_cov")]))), 
                          "\\%"), "&", 
      print_helper(len[2:4]), "&", 
      print_helper(scores[2:4]), "\\\\") %>%
  writeLines

# plot length of prediction interval against forecast horizon
# (Figure 4 in paper)
hist_sel <- "hist1"
len_name <- paste0(hist_sel,"_length")
dat_plot <- dat[, c("target_year", "h", len_name, "pred_length")]
pdf(paste0("plots/", hist_sel, "_dots_", var, ".pdf"))
par(mar=c(5,5,4,1)+.1)
plot(dat_plot$h, dat_plot$hist1_length, bty = "n", xlab = "Horizon (weeks)", 
     ylab = "PI Length", col = grey(.55, .55), pch = 20, ylim = c(0, 16), 
     cex.lab = 1.6, cex.axis = 1.6, cex = 1.6)
legend("topleft", c("Combination", "SPF Histogram"), pch = c(2, 20), 
       bty = "n", cex = 1.6, col = c(1, grey(.55, .55)))
points(dat_plot$h, dat_plot$pred_length, pch = 2, cex = 1.6)
dev.off()

# separate DM tests for each horizon
hs <- dat$h %>% unique %>% sort
df_dm <- data.frame(h = hs, t = NA, p = NA, n = NA)
for (hh in hs){
  ind <- which(hs == hh)
  # score difference for this horizon
  d <- dat %>% filter(h == hh) %>% 
    mutate(d = score-.[,paste0(hist_sel,"_score")]) %>% 
    arrange(target_year) %>% pull(d)
  t_test <- t_hac(d, implementation = "sandwich")
  df_dm$n[ind] <- length(d)
  df_dm$t[ind] <- t_test$t
  df_dm$p[ind] <- t_test$pval
}
df_dm$p_adj <- p.adjust(df_dm$p, method = "bonferroni")
# Plot for Figure 3 in paper
pdf(paste0("plots/pvals_", var, "_", hist_sel, ".pdf"))
plot(df_dm$h, df_dm$p, pch = 3, xlab = "Horizon (weeks)", ylab = "P-value", 
     bty = "n", ylim = c(0, 1), xlim = c(0, 104), 
     cex.lab = 1.4, cex.axis = 1.4, cex = 2)
points(df_dm$h, df_dm$p_adj, pch = 2, col = 1, cex = 2)
abline(h = .05, lty = 2)
legend("topleft", c("Raw", "Adj. (Bonferroni)"), col = c(1, 1), 
       bty = "n", pch = c(3, 2), cex = 1.5)
dev.off()
# print indexes of score differences that are in favor of postprocessing method
which(df_dm$t < 0)