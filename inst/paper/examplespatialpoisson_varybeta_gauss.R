


rm(list = ls())

library(midasINLA)
library(ggplot2)
library(midasr)
library(INLA)
library(ggplot2)
library(tidyr)

data("data_spatialpoisson_varybeta_gauss", package = "midasINLA")
str(data_spatialpoisson_varybeta_gauss)

for_plot <- data_spatialpoisson_varybeta_gauss$data_y
for_plot$loc <- as.factor(for_plot$loc)

ggplot(for_plot, aes(y = y, x= Time, group = loc, color = loc)) +
  geom_line() +
  geom_point() +
  theme_minimal()


data_y <- data_spatialpoisson_varybeta_gauss[["data_y"]]
data_y[which(data_y[["Time"]] %in% 50:60),"y"] <- NA




#### Model fitting ####

Midas_objects <- fit_Minla_spatial_varybeta(xdata = data_spatialpoisson_varybeta_gauss[["data_x"]][["x"]],
                                            ydata = data_y[["y"]],
                                            loc_x = data_spatialpoisson_varybeta_gauss[["data_x"]][["loc"]],
                                            loc_y = data_y[["loc"]],
                                            constraint = "gaussian",
                                            K = 0:29,
                                            m = 30,
                                            family = "poisson")


res = inla(y ~ 1 +
             f(idx, model = Midas_objects[["rgen"]],
               n = nrow(Midas_objects[["data"]])) +
             f(iid_loc, model = "iid"),
           data = data.frame(y = Midas_objects[["data"]][["y"]],
                             idx = 1:nrow(Midas_objects[["data"]]),
                             iid_loc = Midas_objects[["idx_loc"]]),
           verbose = TRUE,
           family = "poisson",
           control.compute=list(config = TRUE))

summary(res)



png("inst/paper/figures/Example_spatialpoisson_varybeta_gauss_paramestimates.png", width=30, height=10, units = 'cm', res = 300)

par(mfrow=c(2,4))

plot(inla.smarginal(res$marginals.fixed[["(Intercept)"]]),
     type="l", lwd=3, col="red", xlab=expression(beta[0]), ylab="",
     cex.lab = 2.2, cex.axis=1.5)
abline(v = data_spatialpoisson_varybeta_gauss$beta0, col = 'blue', lty=1, lwd = 2)
abline(v = quantile(inla.rmarginal(200, res$marginals.fixed$`(Intercept)`), prob = 0.025), lty = 2)
abline(v = quantile(inla.rmarginal(200, res$marginals.fixed$`(Intercept)`), prob = 0.975), lty = 2)

plot(inla.smarginal(res$marginals.hyperpar[['Theta3 for idx']]),
     type="l", lwd=3, col="red", xlab=expression(beta[{1(1)}]), ylab="",
     cex.lab = 2.2, cex.axis=1.5)
abline(v = data_spatialpoisson_varybeta_gauss$beta[1], col = 'blue', lty = 1, lwd = 2)
abline(v = quantile(inla.rmarginal(200, res$marginals.hyperpar$`Theta3 for idx`), prob = 0.025), lty = 2)
abline(v = quantile(inla.rmarginal(200, res$marginals.hyperpar$`Theta3 for idx`), prob = 0.975), lty = 2)

plot(inla.smarginal(res$marginals.hyperpar[['Theta4 for idx']]),
     type="l", lwd=3, col="red", xlab=expression(beta[{1(2)}]), ylab="",
     cex.lab = 2.2, cex.axis=1.5)
abline(v = data_spatialpoisson_varybeta_gauss$beta[2], col = 'blue', lty = 1, lwd = 2)
abline(v = quantile(inla.rmarginal(200, res$marginals.hyperpar$`Theta4 for idx`), prob = 0.025), lty = 2)
abline(v = quantile(inla.rmarginal(200, res$marginals.hyperpar$`Theta4 for idx`), prob = 0.975), lty = 2)

plot(inla.smarginal(res$marginals.hyperpar[['Theta5 for idx']]),
     type="l", lwd=3, col="red", xlab=expression(beta[{1(3)}]), ylab="",
     cex.lab = 2.2, cex.axis=1.5)
abline(v = data_spatialpoisson_varybeta_gauss$beta[3], col = 'blue', lty = 1, lwd = 2)
abline(v = quantile(inla.rmarginal(200, res$marginals.hyperpar$`Theta5 for idx`), prob = 0.025), lty = 2)
abline(v = quantile(inla.rmarginal(200, res$marginals.hyperpar$`Theta5 for idx`), prob = 0.975), lty = 2)

plot(inla.smarginal(res$marginals.hyperpar[['Theta6 for idx']]),
     type="l", lwd=3, col="red", xlab=expression(beta[{1(4)}]), ylab="",
     cex.lab = 2.2, cex.axis=1.5)
abline(v = data_spatialpoisson_varybeta_gauss$beta[4], col = 'blue', lty = 1, lwd = 2)
abline(v = quantile(inla.rmarginal(200, res$marginals.hyperpar$`Theta6 for idx`), prob = 0.025), lty = 2)
abline(v = quantile(inla.rmarginal(200, res$marginals.hyperpar$`Theta6 for idx`), prob = 0.975), lty = 2)

plot(inla.smarginal(res$marginals.hyperpar[['Theta7 for idx']]),
     type="l", lwd=3, col="red", xlab=expression(beta[{1(5)}]), ylab="",
     cex.lab = 2.2, cex.axis=1.5)
abline(v = data_spatialpoisson_varybeta_gauss$beta[5], col = 'blue', lty = 1, lwd = 2)
abline(v = quantile(inla.rmarginal(200, res$marginals.hyperpar$`Theta7 for idx`), prob = 0.025), lty = 2)
abline(v = quantile(inla.rmarginal(200, res$marginals.hyperpar$`Theta7 for idx`), prob = 0.975), lty = 2)

plot(inla.smarginal(res$marginals.hyperpar[['Precision for iid_loc']]),
     type="l", lwd=3, col="red", xlab=expression(tau), ylab="",
     cex.lab = 2.2, cex.axis=1.5)
abline(v = 1/data_spatialpoisson_gauss$phi_sd^2, col = 'blue', lty = 1, lwd = 2)
abline(v = quantile(inla.rmarginal(200, res$marginals.hyperpar$`Precision for iid_loc`), prob = 0.025), lty = 2)
abline(v = quantile(inla.rmarginal(200, res$marginals.hyperpar$`Precision for iid_loc`), prob = 0.975), lty = 2)

dev.off()



#### Compute weights ####

res_weights <- compute_weights(model = res,
                               constraint = "gaussian",
                               lag_k = 29)
str(res_weights)

png("inst/paper/figures/Example_spatialpoisson_varybeta_gauss_weights.png", width=35, height=14, units = 'cm', res = 300)
res_weights$lag <- factor(res_weights$lag, levels=unique(res_weights$lag))
ggplot(res_weights, aes(x=lag, y=mean)) +
  geom_point(aes(col="Posterior mean")) +
  geom_errorbar(aes(ymin=q2.5, ymax=q97.5), width=.2,
                position=position_dodge(0.05), col = "black") +
  geom_point(aes(x=lag, y=data_spatialpoisson_gauss$weights, color = "True value")) +
  scale_color_manual(name='',
                     breaks=c('Posterior mean', 'True value'),
                     values=c('Posterior mean'='red', 'True value'='blue')) +
  theme_bw() +
  theme(axis.text=element_text(size=17),
        axis.title=element_text(size=20,face="bold"),
        legend.position = "bottom",
        legend.text=element_text(size=24, face = "plain"),
        legend.title=element_blank()) +
  ylab("")
dev.off()




