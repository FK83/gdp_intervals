rm(list = ls())

library(tidyr)
library(ggplot2)
library(dplyr)

setwd("[Enter main path]")

# load simulation results (three files, one for each setting)

# note: to avoid large files, the saved simulation files in fe_sim (sim1.RData, sim2.RData, sim3.RData) 
# contain only the simulation results that directly enter Table 1. When running the simulation anew
# (file Section6_simulation.R), much larger files containing additional outputs (and named based on 
# the simulation's finishing time, see variable sv_nm in Section6_simulation.R) are produced. 

# first setup (n = 500, T = 30)
load("fe_sim/sim1.RData")
res500 <- as.data.frame(res_mc_all) 
rm(list = setdiff(ls(), "res500"))

# second setup (n = 1000, T = 30)
load("fe_sim/sim2.RData")
res1000 <- as.data.frame(res_mc_all)
rm(list = setdiff(ls(), c("res500", "res1000")))

# third setup (n = 500, T = 60)
load("fe_sim/sim3.RData")
res60 <- as.data.frame(res_mc_all)
rm(list = setdiff(ls(), c("res500", "res1000", "res60")))

source("gdp_procs23.R")

# method labels
nms <- c("Decomposition", "Flexible", "Gaussian", "True", "Combination")
method_order <- c(3, 1, 2, 5, 4) # ordering of methods in table

# put results into common data frame
df1 <- data.frame(method = nms, 
                  score = colMeans(res500)[1:5],
                  covg = colMeans(res500)[6:10], 
                  len = colMeans(res500)[11:15])
df2 <- data.frame(method = nms, 
                  score = colMeans(res1000)[1:5],
                  covg = colMeans(res1000)[6:10], 
                  len = colMeans(res1000)[11:15])
df3 <- data.frame(method = nms, 
                  score = colMeans(res60)[1:5],
                  covg = colMeans(res60)[6:10], 
                  len = colMeans(res60)[11:15])
res <- rbind(df1, df2, df3) %>% 
  mutate(method = factor(method, levels = nms[method_order]))
rownames(res) <- NULL

# format numbers for printing
covg_print <- print_helper(100*res$covg)
len_print <- print_helper(res$len)
qs_print <- print_helper(res$score)

# print results for first MC setup
paste(nms[method_order], "&", 
      paste0(covg_print[1:5][method_order], "\\%"), "&", 
      len_print[1:5][method_order], "&", 
      qs_print[1:5][method_order], "\\\\") %>%
  writeLines

# print results for second MC setup
paste(nms[method_order], "&", 
      paste0(covg_print[6:10][method_order], "\\%"), "&", 
      len_print[6:10][method_order], "&", 
      qs_print[6:10][method_order], "\\\\") %>%
  writeLines

# print results for third MC setup
paste(nms[method_order], "&", 
      paste0(covg_print[11:15][method_order], "\\%"), "&", 
      len_print[11:15][method_order], "&", 
      qs_print[11:15][method_order], "\\\\") %>%
  writeLines

