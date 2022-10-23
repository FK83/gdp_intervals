rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)

setwd("[ENTER PATH HERE]/fe_sim/")
fn <- "simulation_results.RData"

if (file.exists(fn)){ # load summary data file if already available
  load(fn)
} else {  # if summary file is not available, loop through individual files to create it
  l <- list.files() 
  l <- l[l != "_old"]
  res_all <- matrix(NA, 1.8e6, 11)
  n_rows_all <- 0
  
  for (jj in l){
    load(jj)
    n <- substr(jj, 5,7) %>% as.numeric
    if (grepl("a5", jj)){
      a <- .5
    } else {
      a <- .9
    }
    n_rows <- nrow(simulation_results)
    inds <- (n_rows_all + 1):(n_rows_all + n_rows)
    # Collect results. NOTE: Switch order of models with estimated vs fixed theta (score2, score1, score4, score3)
    res_all[inds, 1:9] <- 
      simulation_results %>% select(scoreb_fixed, scoreb_est, score2,
                                    score1, score4, score3, iter, year, h) %>% 
      as.matrix %>% unname
    res_all[inds, 10] <- a
    res_all[inds, 11] <- n
    n_rows_all <- n_rows_all + n_rows
    rm(list = setdiff(ls(), c("l", "jj", "res_all", "n_rows_all")))
  }
  
  res_all <- data.frame(Truth = res_all[,1], Truth_est = res_all[,2], 
                        Gauss1 = res_all[,3], Gauss2 = res_all[,4], 
                        QR1 = res_all[,5], QR2 = res_all[,6], 
                        iter = res_all[,7], year = res_all[,8], h = res_all[,9],
                        a = res_all[,10], n_smpl = as.factor(res_all[,11])) %>% 
    na.omit
  save(res_all, file = "simulation_results.RData")
}
  
# Make plot for one choice of the AR parameter (called a_sel here)
pd <- position_dodge(.1)
a_sel <- .5
res_all %>% filter(a == a_sel) %>% select(-iter, -year, -h) %>%
  pivot_longer(-c("n_smpl", "a")) %>%
  mutate(name = factor(name, levels = c("Truth", "Truth_est", "Gauss1", "Gauss2", 
                                        "QR1", "QR2"), 
                             labels = c("Truth", "Truth (est.)", "Gauss1", "Gauss2",
                                         "QR1", "QR2"))) %>%
  group_by(name, n_smpl) %>%
  # Compute mean and multiply by 10 (to get standard scaling of interval score)
  summarise(mean = 10*mean(value)) %>%
  ungroup %>%
  ggplot(aes(x = name, y = mean, color = n_smpl, group = n_smpl)) + 
  geom_line(position = pd, size = .2) + 
  geom_point(position = pd, size = 1.2) + 
  theme_minimal(base_size = 9) + theme(legend.position = "top") + 
  scale_color_discrete(name = "Sample Size") + ylab("Mean Score") + 
  xlab("Method") 
# save plot as pdf file
ggsave(paste0("simulation_results_a", 10*a_sel, ".pdf"))
