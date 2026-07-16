
rm(list = ls())

library(midasINLA)
library(ggplot2)
library(midasr)
library(INLA)

data("data_spatialpoisson_example", package = "midasINLA")
str(data_spatialpoisson_example)

for_plot <- data_spatialpoisson_example$data_y
for_plot <- for_plot[which(for_plot$loc %in% c(1:6)),]
for_plot$loc <- as.factor(for_plot$loc)
ggplot(for_plot, aes(y = y, x= Time, group = loc, color = loc)) +
  geom_line() +
  geom_point() +
  theme_minimal()

for_plot <- data_spatialpoisson_example$data_x1
for_plot <- for_plot[which(for_plot$loc %in% c(1:4)),]
for_plot$loc <- as.factor(for_plot$loc)
ggplot(for_plot, aes(y = x1, x= Time, group = loc, color = loc)) +
  geom_line() +
  theme_minimal()


data_y <- data_spatialpoisson_example[["data_y"]]
data_y$y_all <- data_y$y
data_y[which(data_y[["Time"]] %in% 182:192),"y"] <- NA






#### Model fitting ####

response_data <- data_y

g <- inla.read.graph(filename = "inst/map.adj")
Midas_x1 <- prepare_Minla_spatial(x = data_spatialpoisson_example$data_x1$x1,
                                  loc_x = data_spatialpoisson_example$data_x1$loc,
                                  constraint = "hyperbolic",
                                  K = 0:29,
                                  m = 30,
                                  svc = TRUE,
                                  svc_prior = "icar",
                                  g = g)
Midas_x2 <- prepare_Minla_spatial(x = data_spatialpoisson_example$data_x2$x2,
                                  loc_x = data_spatialpoisson_example$data_x2$loc,
                                  constraint = "hyperbolic",
                                  K = 0:45,
                                  m = 30,
                                  svc = FALSE)


fit_res <- fit_Minla_spatial(formula = y ~ 1,
                             data = response_data,
                             loc_var = "loc",
                             time_var = "Time",
                             family = "poisson",
                             hf_input = list(Midas_x1,Midas_x2),
                             inla_options = list(verbose = T))




#png("inst/paper/figures/Example_binomial_hyperbolic_paramestimates.png", width=30, height=10, units = 'cm', res = 300)
par(mfrow=c(1,4))

plot(inla.smarginal(fit_res$res$marginals.fixed[["(Intercept)"]]),
     type="l", lwd=3, col="red", xlab=expression(beta[0]), ylab="",
     cex.lab = 2.2, cex.axis=1.5)
abline(v = data_binomial_example$beta0, col = 'blue', lty=1, lwd = 2)
abline(v = quantile(inla.rmarginal(200, fit_res$res$marginals.fixed$`(Intercept)`), prob = 0.025), lty = 2)
abline(v = quantile(inla.rmarginal(200, fit_res$res$marginals.fixed$`(Intercept)`), prob = 0.975), lty = 2)

plot(inla.smarginal(fit_res$res$marginals.fixed[['trend']]),
     type="l", lwd=3, col="red", xlab=expression(beta[1]), ylab="",
     cex.lab = 2.2, cex.axis=1.5)
abline(v = data_binomial_example$beta1, col = 'blue', lty = 1, lwd = 2)
abline(v = quantile(inla.rmarginal(200, fit_res$res$marginals.fixed[['trend']]), prob = 0.025), lty = 2)
abline(v = quantile(inla.rmarginal(200, fit_res$res$marginals.fixed[['trend']]), prob = 0.975), lty = 2)

plot(inla.smarginal(fit_res$res$marginals.hyperpar[['Theta2 for hf_idx_1']]),
     type="l", lwd=3, col="red", xlab=expression(beta[2]), ylab="",
     cex.lab = 2.2, cex.axis=1.5)
abline(v = data_binomial_example$beta2, col = 'blue', lty = 1, lwd = 2)
abline(v = quantile(inla.rmarginal(200, fit_res$res$marginals.hyperpar$`Theta2 for hf_idx_1`), prob = 0.025), lty = 2)
abline(v = quantile(inla.rmarginal(200, fit_res$res$marginals.hyperpar$`Theta2 for hf_idx_1`), prob = 0.975), lty = 2)

plot(inla.smarginal(fit_res$res$marginals.hyperpar[['Theta3 for hf_idx_2']]),
     type="l", lwd=3, col="red", xlab=expression(beta[3]), ylab="",
     cex.lab = 2.2, cex.axis=1.5)
abline(v =data_binomial_example$beta3, col = 'blue', lty = 1, lwd = 2)
abline(v = quantile(inla.rmarginal(200, fit_res$res$marginals.hyperpar$`Theta3 for hf_idx_2`), prob = 0.025), lty = 2)
abline(v = quantile(inla.rmarginal(200, fit_res$res$marginals.hyperpar$`Theta3 for hf_idx_2`), prob = 0.975), lty = 2)
#dev.off()



#### Compute weights ####

res_weights <- compute_weights(fit_res)
str(res_weights)

#png("inst/paper/figures/data_binomial_example_weights_covariate1.png", width=27, height=12, units = 'cm', res = 300)
res_weights$hf_1$lag <- factor(res_weights$hf_1$lag, levels=unique(res_weights$hf_1$lag))
ggplot(res_weights$hf_1, aes(x=lag, y=mean)) +
  geom_point(aes(col="Posterior mean")) +
  geom_errorbar(aes(ymin=q2.5, ymax=q97.5), width=.2,
                position=position_dodge(0.05), col = "black") +
  geom_point(aes(x=lag, y=data_spatialpoisson_example$weights1, color = "True value")) +
  scale_color_manual(name='',
                     breaks=c('Posterior mean', 'True value'),
                     values=c('Posterior mean'='red', 'True value'='blue')) +
  theme_bw() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        legend.position = "bottom",
        legend.text=element_text(size=24, face = "plain"),
        legend.title=element_blank()) +
  ylab("")
#dev.off()

#png("inst/paper/figures/data_binomial_example_weights_covariate2.png", width=27, height=12, units = 'cm', res = 300)
res_weights$hf_2$lag <- factor(res_weights$hf_2$lag, levels=unique(res_weights$hf_2$lag))
ggplot(res_weights$hf_2, aes(x=lag, y=mean)) +
  geom_point(aes(col="Posterior mean")) +
  geom_errorbar(aes(ymin=q2.5, ymax=q97.5), width=.2,
                position=position_dodge(0.05), col = "black") +
  geom_point(aes(x=lag, y=data_spatialpoisson_example$weights2, color = "True value")) +
  scale_color_manual(name='',
                     breaks=c('Posterior mean', 'True value'),
                     values=c('Posterior mean'='red', 'True value'='blue')) +
  theme_bw() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        legend.position = "bottom",
        legend.text=element_text(size=24, face = "plain"),
        legend.title=element_blank()) +
  ylab("")
#dev.off()


#### Compare observed versus predicted values ####

pred_res <- predict_midas(model = fit_res,
                          family = "poisson",
                          Ntrials = NULL,
                          nsamples = 1000)
both_ts <- data.frame(observed = fit_res$data_final$y_all,
                      predicted = pred_res$computed_y$mean,
                      loc = fit_res$data_final$loc,
                      Time = fit_res$data_final$Time)
PI <- data.frame(lower = pred_res$computed_y$q2.5,
                 upper = pred_res$computed_y$q97.5,
                 loc = fit_res$data_final$loc,
                 Time = fit_res$data_final$Time)

library(tidyr)
both_long <- both_ts |>
  pivot_longer(
    cols = c(observed, predicted),
    names_to = "series",
    values_to = "value"
  )
first_na_idx <- which(diff(c(FALSE, is.na(fit_res$data_final$y))) == 1)
idx <- seq_along(fit_res$data_final$y)
non_na <- !is.na(fit_res$data_final$y)
segment <- cumsum(non_na != dplyr::lag(non_na, default = TRUE))
segment[!non_na] <- NA
rel_idx <- ave(idx, segment, FUN = seq_along)
first_na <- which(diff(c(FALSE, is.na(fit_res$data_final$y))) == 1)

vlines_df <- data.frame(cut = rel_idx[first_na - 1] + 1,
                        loc = 1:16)

both_long_sub <- both_long[which(both_long$loc %in% 1:4),]
vlines_df_sub <- vlines_df[which(vlines_df$loc %in% 1:4),]
PI_sub <- PI[which(PI$loc %in% 1:4),]



#png("inst/paper/figures/Example_spatialpoisson_gauss_predsVSobs.png", width=35, height=20, units = 'cm', res = 300)
ggplot(both_long_sub, aes(x = Time, y = value, colour = series)) +
  geom_line() +
  geom_ribbon(data = PI_sub,
              aes(x = Time, ymin = lower, ymax = upper),
              inherit.aes = FALSE,
              alpha = 0.2,
              fill  = "red") +
  theme_minimal() +
  facet_wrap(~loc, ncol = 2,
             labeller = labeller(loc = function(x) paste("Loc =", x))) +
  scale_colour_manual(
    values = c("blue","red")
  ) +
  ylab("y") +
  geom_vline(
    data = vlines_df_sub,
    aes(xintercept = cut),
    colour = "black",
    size = 1
  ) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"),
        legend.position = "bottom",
        legend.text=element_text(size=20, face = "plain"),
        legend.title=element_blank(),
        strip.text = element_text(size=20))
#dev.off()



