rm(list = ls())
library(dplyr)
library(ggplot2)
library(knitr)

setwd("[ENTER PATH HERE]")

dat1 <- read.csv("gdp_forecasts_DE.csv") %>% 
  transmute(target_year, h = horizon, type = "DE")
dat2 <- read.csv("gdp_forecasts_US.csv") %>% na.omit %>%
  transmute(target_year, h, type = "US_gdp")
dat3 <- read.csv("inf_forecasts_US.csv") %>% na.omit %>%
  transmute(target_year, h, type = "US_inf")
dat <- rbind(dat1, dat2, dat3) 

dat %>% group_by(type) %>% 
  summarise(from = min(target_year), to = max(target_year),
            n = n(), min_h = min(h), q25 = quantile(h, .25), 
            q50 = quantile(h, .5), q75 = quantile(h, .75), 
            max_h = max(h)) %>% 
  ungroup %>% kable(digits = 1, format = "latex", booktabs = TRUE)

