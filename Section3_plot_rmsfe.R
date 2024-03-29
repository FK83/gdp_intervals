rm(list = ls())

setwd("[Enter main path]")

# set random seed
set.seed(2023)

# load empirical RMSEs for Germany
rmse_gdp <- read.csv("data/rmse_gdp_de.csv")

# load functions
source("gdp_procs23.R")

# PT GDP: s2e = .04, s2w = .004, phi = .2
# PT Inflation: s2e = .012, s2w = 2e-4, phi = .2
# DE GDP: s2e = .09, s2w = .003, phi = .3
model1 <- ar1_model(n_mc = 2, n_years = 10, s2e = .09, 
                    s2w = .003, 
                    phi = .3)
# Plot for Figure 1 in paper
pdf("plots/plot_rmse.pdf", 
    width = 16, height = 9)
par(mar = c(5.1, 6.1, 4.1, 2.1))
plot(model1$rmse, ylim = c(0, max(c(model1$rmse, rmse_gdp$rmse))), 
     bty = "n", 
     ylab = "Root Mean Squared Forecast Error", 
     xlab = "Forecast Horizon (in Weeks)", 
     type = "l" ,
     cex.axis = 1.6, cex.lab = 1.6, lwd = 3)
abline(v = 52, lty = 2)
points(rmse_gdp$h, rmse_gdp$rmse, 
       cex = 2, pch = 20)
dev.off()

# higher persistence
model2 <- ar1_model(n_mc = 2, n_years = 10, s2e = .09, 
                    s2w = .003, 
                    phi = .6)
# higher AR(1) error term variance
model3 <- ar1_model(n_mc = 2, n_years = 10, s2e = .18, 
                    s2w = .003, 
                    phi = .3)
# higher observation error term 
model4 <- ar1_model(n_mc = 2, n_years = 10, s2e = .09, 
                    s2w = .006, 
                    phi = .3)

rmse_all <- cbind(model1$rmse, model2$rmse, model3$rmse, model4$rmse)

lb1 <- expression(rho * " = 0.3," * sigma[epsilon]^2 * " = 0.09," *
                      sigma[eta]^2 * " = 0.003")
lb2 <- expression("Doubling " * rho)
lb3 <- expression("Doubling " * sigma[epsilon]^2)
lb4 <- expression("Doubling " * sigma[eta]^2)

# Plot for Figure S1 in online supplement
pdf("plots/plot_rmse_all.pdf", 
    width = 16, height = 9)
par(mar = c(5.1, 6.1, 4.1, 2.1))
matplot(rmse_all, type = "l", bty = "n", 
        ylab = "Root Mean Squared Forecast Error", 
        xlab = "Forecast Horizon (in Weeks)", 
        pch = 20, 
        cex.axis = 1.6, cex.lab = 1.6, cex = 2, lty = 1:4, lwd = 3, col = 1)
legend("topleft", c(lb1, lb2, lb3, lb4), lty = 1:4, lwd = 3, bty = "n", 
       cex = 1.5, col = 1)
dev.off()