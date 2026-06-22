
rm(list = ls())

library(midasINLA)
library(ggplot2)
library(midasr)
library(INLA)

data("data_poiss_gaussianc", package = "midasINLA")
str(data_poiss_gaussianc)

y_ts <- ts(data_poiss_gaussianc$y, frequency = 1)
x_ts <- ts(data_poiss_gaussianc$x1, frequency = 1)
x2_ts <- ts(data_poiss_gaussianc$x2, frequency = 1)

#png("paper/figures/Example_binomial_hyperbolic_covariatedata.png", width=17, height=10, units = 'cm', res = 300)
autoplot(x_ts, color = "brown") + theme_bw() + geom_point(color = "red", size = .01) +
  theme(axis.text.x=element_text(size=16),
        axis.title.x=element_text(size=16,face="bold"),
        axis.text.y=element_text(size=16),
        axis.title.y=element_blank(),
        legend.position = "bottom",
        legend.text=element_text(size=20, face = "plain"),
        legend.title=element_blank())
#dev.off()
autoplot(x2_ts, color = "brown") + theme_bw() + geom_point(color = "red", size = .01) +
  theme(axis.text.x=element_text(size=16),
        axis.title.x=element_text(size=16,face="bold"),
        axis.text.y=element_text(size=16),
        axis.title.y=element_blank(),
        legend.position = "bottom",
        legend.text=element_text(size=20, face = "plain"),
        legend.title=element_blank())
#png("paper/figures/Example_binomial_hyperbolic_response.png", width=17, height=10, units = 'cm', res = 300)
autoplot(y_ts, color = "black") + theme_bw() + geom_point(color = "gray", size = .01) +
  theme(axis.text.x=element_text(size=16),
        axis.title.x=element_text(size=16,face="bold"),
        axis.text.y=element_text(size=16),
        axis.title.y=element_blank(),
        legend.position = "bottom",
        legend.text=element_text(size=20, face = "plain"),
        legend.title=element_blank())
#dev.off()


x <- data_poiss_gaussianc[["x1"]]
y <- data_poiss_gaussianc[["y"]]
y[(length(y)-15):length(y)] <- NA # remove last 16 time points as forecast test set

#### Model fitting ####

Midas_objects <- fit_Minla_c(xdata = x,
                           ydata = y,
                           constraint = "gaussian",
                           K = 0:29,
                           m = 30,
                           lagY = 0,
                           family = "poisson")

res = inla(y ~ 1 +f(idx,
                    model = Midas_objects[["cmodel"]]) +
             f(x2,
               model = "linear"),
           data = data.frame(y = Midas_objects[["data"]][["y"]],
                             idx = 1:nrow(Midas_objects[["data"]]),
                             x2 = data_poiss_gaussianc$x2),
           verbose = FALSE,
           family = "poisson",
           control.compute=list(config = TRUE))

summary(res)


res = inla(y ~ 1 +f(x2,
               model = "linear"),
           data = data.frame(y = Midas_objects[["data"]][["y"]],
                             x2 = data_poiss_gaussianc$x2),
           verbose = TRUE,
           family = "poisson",
           control.compute=list(config = TRUE))
summary(res)


res = inla(y ~ 1 +f(idx,
                    model = Midas_objects[["cmodel"]]),
           data = data.frame(y = Midas_objects[["data"]][["y"]],
                             idx = 1:nrow(Midas_objects[["data"]])),
           verbose = FALSE,
           family = "poisson",
           control.compute=list(config = TRUE))
summary(res)


res = inla(y ~ 1 +f(idx,
                    model = Midas_objects[["cmodel"]]) +
             f(Time_iid, model = "iid"),
           data = data.frame(y = Midas_objects[["data"]][["y"]],
                             idx = 1:nrow(Midas_objects[["data"]]),
                             Time_iid = 1:nrow(Midas_objects[["data"]])),
           verbose = FALSE,
           family = "poisson",
           control.compute=list(config = TRUE))
summary(res)


#png("inst/paper/figures/Example_binomial_hyperbolic_paramestimates.png", width=30, height=10, units = 'cm', res = 300)
par(mfrow=c(1,3))

plot(inla.smarginal(res$marginals.fixed[["(Intercept)"]]),
     type="l", lwd=3, col="red", xlab=expression(beta[0]), ylab="",
     cex.lab = 2.2, cex.axis=1.5)
abline(v = data_poiss_gaussianc$beta0, col = 'blue', lty=1, lwd = 2)
abline(v = quantile(inla.rmarginal(200, res$marginals.fixed$`(Intercept)`), prob = 0.025), lty = 2)
abline(v = quantile(inla.rmarginal(200, res$marginals.fixed$`(Intercept)`), prob = 0.975), lty = 2)

plot(inla.smarginal(res$marginals.hyperpar[['Theta3 for idx']]),
     type="l", lwd=3, col="red", xlab=expression(beta[2]), ylab="",
     cex.lab = 2.2, cex.axis=1.5)
abline(v = data_poiss_gaussianc$beta1, col = 'blue', lty = 1, lwd = 2)
abline(v = quantile(inla.rmarginal(200, res$marginals.hyperpar$`Theta3 for idx`), prob = 0.025), lty = 2)
abline(v = quantile(inla.rmarginal(200, res$marginals.hyperpar$`Theta3 for idx`), prob = 0.975), lty = 2)

plot(inla.smarginal(res$marginals.fixed[['x2']]),
     type="l", lwd=3, col="red", xlab=expression(beta[2]), ylab="",
     cex.lab = 2.2, cex.axis=1.5)
abline(v = data_poiss_gaussianc$beta2, col = 'blue', lty = 1, lwd = 2)
abline(v = quantile(inla.rmarginal(200, res$marginals.fixed[['x2']]), prob = 0.025), lty = 2)
abline(v = quantile(inla.rmarginal(200, res$marginals.fixed[['x2']]), prob = 0.975), lty = 2)


#dev.off()



#### Compute weights ####

res_weights <- compute_weights(model = res,
                               constraint = "gaussian",
                               lag_k = 29)
str(res_weights)

#png("inst/paper/figures/Example_binomial_hyperbolic_weights.png", width=27, height=12, units = 'cm', res = 300)
res_weights$lag <- factor(res_weights$lag, levels=unique(res_weights$lag))
ggplot(res_weights, aes(x=lag, y=mean)) +
  geom_point(aes(col="Posterior mean")) +
  geom_errorbar(aes(ymin=q2.5, ymax=q97.5), width=.2,
                position=position_dodge(0.05), col = "black") +
  geom_point(aes(x=lag, y=data_poiss_gaussianc$weights, color = "True value")) +
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

pred_data <- Midas_objects[["data"]]
pred_res <- predict_midas(model = res,
                          data = pred_data,
                          family = "poisson",
                          nsamples = 30)
both_ts <- ts(data.frame(observed = as.vector(data_poiss_gaussianc$y)[(length(data_poiss_gaussianc$y)-length(Midas_objects$data$y)+1):length(data_poiss_gaussianc$y)],
                         predicted = pred_res$computed_y$mean),
              start = 1, end = length(pred_res$computed_y$mean))
PI <- data.frame(lower = pred_res$computed_y$q2.5,
                 upper = pred_res$computed_y$q97.5,
                 Time = 1:length(pred_res$computed_y$mean))
PI$Time <- time(both_ts)


#png("inst/paper/figures/Example_binomial_hyperbolic_predsVSobs.png", width=25, height=12, units = 'cm', res = 300)
autoplot(both_ts) +
  geom_ribbon(data = PI,
              aes(x = Time, ymin = lower, ymax = upper),
              inherit.aes = FALSE,
              alpha = 0.2,
              fill  = "red") +
  theme_bw() +
  scale_colour_manual(
    values = c("blue","red")
  ) +
  ylab("y") +
  geom_vline(xintercept = max(which(!is.na(Midas_objects$data$y))) + 0.5, linetype = "solid", colour = "black", linewidth = 1) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"),
        legend.position = "bottom",
        legend.text=element_text(size=20, face = "plain"),
        legend.title=element_blank())
#dev.off()




#### Fit a baseline approach ####

df <- data.frame(y = Midas_objects[["data"]][["y"]],
                 x = Midas_objects[["data"]][["lag0"]],
                 trend = 1:length(Midas_objects[["data"]][["y"]]),
                 Ntrials = Midas_objects[["Ntrials"]])

res_baseline = inla(y ~ 1 + f(x, model = "linear") + f(trend, model = "linear"),
                    data = df,
                    family = "binomial",
                    Ntrials = df$Ntrials,
                    verbose = TRUE,
                    control.compute=list(config = TRUE))

summary(res_baseline)

pred_res_baseline <- predict_midas(model = res_baseline,
                                   data = df,
                                   family = "binomial",
                                   Ntrials = Midas_objects[["Ntrials"]],
                                   nsamples = 30)

both_ts_baseline <- ts(data.frame(observed = as.vector(data_binom_hyperbolic$y)[(length(data_binom_hyperbolic$y)-length(Midas_objects$data$y)+1):length(data_binom_hyperbolic$y)],
                                  predicted = pred_res_baseline$computed_y$mean),
                       start = 1, end = length(pred_res_baseline$computed_y$mean))


png("inst/paper/figures/PredsVsObs_hyperbolic_baseline_binomial.png", width=25, height=12, units = 'cm', res = 300)
autoplot(both_ts_baseline) +
  geom_ribbon(data = cbind(as.data.frame(pred_res_baseline$computed_y),
                           Time = 1:length(pred_res_baseline$computed_y$mean)),
              aes(x = Time, ymin = q2.5, ymax = q97.5),
              inherit.aes = FALSE,
              alpha = 0.2,
              fill  = "red") +
  theme_bw() +
  scale_colour_manual(
    values = c("blue","red")
  ) +
  geom_vline(xintercept = max(which(!is.na(Midas_objects$data$y))) + 0.5, linetype = "solid", colour = "black", linewidth = 1) +
  #scale_x_continuous(limits = c(180, NA)) +
  ylab("y") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"),
        legend.position = "bottom",
        legend.text=element_text(size=20, face = "plain"),
        legend.title=element_blank())
dev.off()




#### Compute scores ####


res_scores <- compute_forecast_scores(y = data_binom_hyperbolic[["y"]],
                                      Midas_object = Midas_objects,
                                      pred_res = pred_res,
                                      family = "binomial")
str(res_scores)

res_scores$scores$Approach <- "Lag-distributed model"
res_scores$scores


res_scores_baseline <- compute_forecast_scores(y = data_binom_hyperbolic$y,
                                               Midas_object = Midas_objects,
                                               pred_res = pred_res_baseline,
                                               family = "binomial")
res_scores_baseline$scores$Approach <- "Model comparator"
res_scores_baseline$scores

summary(res_scores$scores)
summary(res_scores_baseline$scores)


for_plot <- rbind(res_scores$scores,
                  res_scores_baseline$scores)
for_plot$Approach <- factor(for_plot$Approach)


library(dplyr)
means_DS <- for_plot %>%
  group_by(Approach) %>%
  summarise(mean_DS = mean(DS, na.rm = TRUE))
means_NLS <- for_plot %>%
  group_by(Approach) %>%
  summarise(mean_NLS = mean(NLS, na.rm = TRUE))

png("C:/Users/sv20/Downloads/Comparison_DS_binomial.png", width=18, height=12, units = 'cm', res = 300)
ggplot(for_plot, aes(x = DS, fill = Approach)) +
  geom_histogram(alpha = 0.5, position = "identity") +
  geom_vline(data = means_DS,
             aes(xintercept = mean_DS, color = Approach, linetype = "Mean"),
             linewidth = 1) +
  scale_linetype_manual(name = "", values = c(Mean = "dashed")) +
  theme_bw() +
  theme(legend.text = element_text(size = 15),
        legend.title = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 15),
        legend.position = "bottom")
dev.off()

png("C:/Users/sv20/Downloads/Comparison_NLS_binomial.png", width=18, height=12, units = 'cm', res = 300)
ggplot(for_plot, aes(x = NLS, fill = Approach)) +
  geom_histogram(alpha = 0.5, position = "identity") +
  geom_vline(data = means_NLS,
             aes(xintercept = mean_NLS, color = Approach, linetype = "Mean"),
             linewidth = 1) +
  scale_linetype_manual(name = "", values = c(Mean = "dashed")) +
  theme_bw() +
  theme(legend.text = element_text(size = 15),
        legend.title = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 15),
        legend.position = "bottom")
dev.off()


