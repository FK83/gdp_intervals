rm(list = ls())

setwd("[Enter main path]")

# load packages
library(dplyr)
library(forecastcalibration)
library(tidyr)

# load functions
source("gdp_procs23.R")
source("dm_procs24.R")

# helper function that collects four types of p-values
dm_pvals <- function(d, h){
  res <- rep(NA, 4)
  # implementation in sandwich package
  res[1] <- t_hac(d, implementation = "sandwich")$pval
  # implementation used by Clark and McCracken (2013)
  res[2] <- dmtest(d = d, h = h)$p1
  # EWC estimator recommended by Lazarus et al (JBES, 2018)
  res[3] <- ewc(lm(d~1))$p
  # IID t-test
  res[4] <- t.test(d)$p.value 
  res
}

var <- "gdp" 	# choose variable (either "gdp" or "inf")
model_nr <- "4" # choose model to be compared to histograms ("4" yields combination as used in the paper)
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

# plot length of prediction interval against forecast horizon
hist_sel <- "hist1"

# separate DM tests for each horizon
hs <- dat$h %>% unique %>% sort
df_dm <- data.frame(h = hs, p1 = NA, p2 = NA, 
                    p3 = NA, p4 = NA, rho_d = NA, 
                    n = NA)
for (hh in hs){
  ind <- which(hs == hh)
  # score difference for this horizon
  d <- dat %>% filter(h == hh) %>% 
    mutate(d = score-.[,paste0(hist_sel,"_score")]) %>% 
    arrange(target_year) %>% pull(d)
  # horizon (in years)
  h_y <- ifelse(sum(hs <= hh) <= 4, 1, 2)
  df_dm[ind,paste0("p", 1:4)] <- dm_pvals(d = d, h = h_y)
  df_dm[ind, "rho_d"] <- acf(d, plot = FALSE)$acf[2]
  df_dm$n[ind] <- length(d)
}

# Plots for top row of Figure S3 in online supplement
pdf(paste0("plots/DM_pvals_", var, ".pdf"))
plot(df_dm$h, df_dm$p1, bty = "n", xlab = "Horizon (weeks)",
     ylab = "DM p-value", pch = 1, ylim = c(0, 1), 
     cex.lab = 1.4, cex.axis = 1.6, cex = 2)
for (jj in 2:4){
  points(jitter(df_dm$h), df_dm[, paste0("p", jj)], 
         col = 1, pch = jj, cex = 2)
}
legend("topleft", c("sandwich", "CM13", "EWC", "IID"), 
       col = 1, pch = 1:4, bty = "n", cex = 1.6)
abline(h = .05, lty = 2)
dev.off()

# Plots for bottom row of Figure S3 in online supplement
pdf(paste0("plots/DM_persistence_", var, ".pdf"))
plot(df_dm$h, df_dm$rho_d, bty = "n", xlab = "Horizon (weeks)", 
     ylab = "Autocorr. of Loss Diff.", pch = 1, ylim = c(-.5, .5), 
     cex.axis = 1.6, cex.lab = 1.4, cex = 2)
abline(h = 0, lty = 2)
points(df_dm$h, qnorm(.975)/sqrt(df_dm$n), pch = 2)
lines(df_dm$h, qnorm(.975)/sqrt(df_dm$n), lwd = .5)
points(df_dm$h, -qnorm(.975)/sqrt(df_dm$n), pch = 2)
lines(df_dm$h, -qnorm(.975)/sqrt(df_dm$n), lwd = .5)
dev.off()

# Tests that average across horizons
# Get set of common forecast cases
dat_all <- dat %>% 
  transmute(target_year, h, d = score-.[,paste0(hist_sel,"_score")]) %>% 
  arrange(target_year) %>% 
  pivot_wider(names_from = "h", values_from = "d") %>%
  na.omit

# Compute average loss diff. across horizons (for each year)
d_av <- cbind(dat_all$target_year, rowMeans(dat_all[,-1]))

# Make Figure S4 in online supplement
pdf(paste0("plots/avg_diff_", var, ".pdf"))
plot(d_av[,1], d_av[,2], xlab = "Year", 
     ylab = "Avg. Loss. Diff.", bty = "n", 
     cex = 2, cex.axis = 1.6, cex.lab = 1.4, 
     ylim = c(-5, 5))
abline(h = 0, lty = 2)
abline(h = mean(d_av[,2]), lty = 1)
dev.off()
round(dm_pvals(d = d_av[,2], h = 1), 3)
