
rm(list = ls())

library(midasINLA)
library(ggplot2)
library(midasr)
library(INLA)

data("data_binomial_example", package = "midasINLA")
str(data_binomial_example)

y_ts <- ts(data_binomial_example$y, frequency = 1)
x1_ts <- ts(data_binomial_example$x1, frequency = 1)
x2_ts <- ts(data_binomial_example$x2, frequency = 1)


png("inst/paper/figures/data_binomial_example_covariate1.png", width=17, height=10, units = 'cm', res = 300)
autoplot(x1_ts, color = "brown") + theme_bw() + geom_point(color = "red", size = .01) +
  theme(axis.text.x=element_text(size=16),
        axis.title.x=element_text(size=16,face="bold"),
        axis.text.y=element_text(size=16),
        axis.title.y=element_blank(),
        legend.position = "bottom",
        legend.text=element_text(size=20, face = "plain"),
        legend.title=element_blank())
dev.off()
png("inst/paper/figures/data_binomial_example_covariate2.png", width=17, height=10, units = 'cm', res = 300)
autoplot(x2_ts, color = "brown") + theme_bw() + geom_point(color = "red", size = .01) +
  theme(axis.text.x=element_text(size=16),
        axis.title.x=element_text(size=16,face="bold"),
        axis.text.y=element_text(size=16),
        axis.title.y=element_blank(),
        legend.position = "bottom",
        legend.text=element_text(size=20, face = "plain"),
        legend.title=element_blank())
dev.off()
png("inst/paper/figures/data_binomial_example_response.png", width=17, height=10, units = 'cm', res = 300)
autoplot(y_ts, color = "black") + theme_bw() + geom_point(color = "gray", size = .01) +
  theme(axis.text.x=element_text(size=16),
        axis.title.x=element_text(size=16,face="bold"),
        axis.text.y=element_text(size=16),
        axis.title.y=element_blank(),
        legend.position = "bottom",
        legend.text=element_text(size=20, face = "plain"),
        legend.title=element_blank())
dev.off()



data = data.frame(y = data_binomial_example$y,
                  Ntrials = data_binomial_example$Ntrials,
                  trend = 1:length(data_binomial_example$y))
data$y_all <- data$y
data$y[(length(data$y)-15):length(data$y)] <- NA # remove last 16 time points as forecast test set


#### Model fitting ####


Midas_x1 <- prepare_Minla(x = data_binomial_example$x1,
                          constraint = "hyperbolic",
                          K = 0:13,
                          m = 7)
Midas_x2 <- prepare_Minla(x = data_binomial_example$x2,
                          constraint = "gaussian",
                          K = 0:20,
                          m = 7)


fit_res <- fit_Minla(formula = y ~ 1 + f(trend, model = "linear"),
                     data = data,
                     family = "binomial",
                     Ntrials = data$Ntrials,
                     hf_input = list(Midas_x1,Midas_x2),
                     inla_options = list(verbose = T))

summary(fit_res$res)


png("inst/paper/figures/Example_binomial_hyperbolic_paramestimates.png", width=30, height=10, units = 'cm', res = 300)
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

dev.off()



#### Compute weights ####

res_weights <- compute_weights(fit_res)
str(res_weights)

png("inst/paper/figures/data_binomial_example_weights_covariate1.png", width=27, height=12, units = 'cm', res = 300)
res_weights$hf_1$lag <- factor(res_weights$hf_1$lag, levels=unique(res_weights$hf_1$lag))
ggplot(res_weights$hf_1, aes(x=lag, y=mean)) +
  geom_point(aes(col="Posterior mean")) +
  geom_errorbar(aes(ymin=q2.5, ymax=q97.5), width=.2,
                position=position_dodge(0.05), col = "black") +
  geom_point(aes(x=lag, y=data_binomial_example$weights1, color = "True value")) +
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
dev.off()

png("inst/paper/figures/data_binomial_example_weights_covariate2.png", width=27, height=12, units = 'cm', res = 300)
res_weights$hf_2$lag <- factor(res_weights$hf_2$lag, levels=unique(res_weights$hf_2$lag))
ggplot(res_weights$hf_2, aes(x=lag, y=mean)) +
  geom_point(aes(col="Posterior mean")) +
  geom_errorbar(aes(ymin=q2.5, ymax=q97.5), width=.2,
                position=position_dodge(0.05), col = "black") +
  geom_point(aes(x=lag, y=data_binomial_example$weights2, color = "True value")) +
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
dev.off()


#### Compare observed versus predicted values ####

pred_res <- predict_midas(model = fit_res,
                          family = "binomial",
                          Ntrials = NULL,
                          nsamples = 1000)
both_ts <- ts(data.frame(observed = fit_res$data_final$y_all,
                         predicted = pred_res$computed_y$mean),
              start = 1, end = length(pred_res$computed_y$mean))
PI <- data.frame(lower = pred_res$computed_y$q2.5,
                 upper = pred_res$computed_y$q97.5,
                 Time = 1:length(pred_res$computed_y$mean))
PI$Time <- time(both_ts)




png("inst/paper/figures/data_binomial_example_predsVSobs.png", width=25, height=12, units = 'cm', res = 300)
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
  geom_vline(xintercept = max(which(!is.na(data$y))) + 0.5, linetype = "solid", colour = "black", linewidth = 1) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"),
        legend.position = "bottom",
        legend.text=element_text(size=20, face = "plain"),
        legend.title=element_blank())
dev.off()




#### Fit a baseline approach ####

data_baseline <- data.frame(y = fit_res$data_final$y,
                            x1 = fit_res$hf_input[[1]]$X_matrix[-c(1:2),1],
                            x2 = fit_res$hf_input[[2]]$X_matrix[-c(1:2),1],
                            trend = fit_res$data_final$trend,
                            Ntrials = fit_res$data_final$Ntrials)

res_baseline = inla(y ~ 1 +
                      f(x1, model = "linear") +
                      f(x2, model = "linear") +
                      f(trend, model = "linear"),
                    data = data_baseline,
                    family = "binomial",
                    Ntrials = data_baseline$Ntrials,
                    verbose = TRUE,
                    control.compute=list(config = TRUE))
summary(res_baseline)

nsamples <- 1000
temp <- INLA::inla.posterior.sample(n = nsamples, res_baseline)
latent_names <- rownames(temp[[1]]$latent)
predictor_idx <- grep("^Predictor", latent_names)
predictor_idx <- predictor_idx[order(as.numeric(
  sub("Predictor:", "", latent_names[predictor_idx])
))]
latent_predictor <- sapply(seq_len(nsamples), function(i) {
  temp[[i]]$latent[predictor_idx]
})
samples <- list(
  latent_predictor = latent_predictor
)
sample_y <- sapply(seq_len(nsamples), function(i) {
  stats::rbinom(
    n = nrow(data_baseline),
    size = data_baseline$Ntrials,
    prob = exp(latent_predictor[, i]) / (1 + exp(latent_predictor[, i]))
  )
})

predicted_y <- list(
  mean = rowMeans(sample_y),
  sd = apply(sample_y, 1, stats::sd),
  q2.5 = matrixStats::rowQuantiles(sample_y, probs = 0.025),
  q97.5 = matrixStats::rowQuantiles(sample_y, probs = 0.975)
)


both_ts <- ts(data.frame(observed = fit_res$data_final$y_all,
                         predicted = predicted_y$mean),
              start = 1, end = length(predicted_y$mean))
PI <- data.frame(lower = predicted_y$q2.5,
                 upper = predicted_y$q97.5,
                 Time = 1:length(predicted_y$mean))
PI$Time <- time(both_ts)

png("inst/paper/figures/data_binomial_example_predsVSobs_baseline.png", width=25, height=12, units = 'cm', res = 300)
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
  geom_vline(xintercept = max(which(!is.na(data$y))) + 0.5, linetype = "solid", colour = "black", linewidth = 1) +
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


