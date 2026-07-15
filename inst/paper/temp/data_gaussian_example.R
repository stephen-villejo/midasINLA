rm(list=ls())

library(midasr)

#
##### Hyperbolic Scheme Polynomial #####
#

# model parameters

n <- 260
m <- 7
beta0 = 10
beta1 = 3
beta2 = -1

gamma_param <- 0.9
mu_f <- 8
sigma_f <- 7
lag_k_1 <- 13 # 0:lag_k, this implies (lag_k + 1) high freq covariate associated with each response y
lag_k_2 <- 20


# compute weights
xi <- c()
for(lag in 0:lag_k_1){
  temp <- gamma(lag+gamma_param)/(gamma(lag+1)*gamma(gamma_param))
  xi <- c(xi, temp)
}
weights1 <- xi/sum(xi)
#plot(weights1, type = "l")

xi <- c()
for(lag in 0:lag_k_2){
  temp <- exp( -(lag-mu_f)^2 / (2 * (sigma_f^2)) )
  xi <- c(xi, temp)
}
weights2 <- xi/sum(xi)
#plot(weights2, type = "l")


# compute covariate, Ntrials, and linear predictor

set.seed(9912)
x1 <- rnorm(m*n,5,3) # covariate 1
x2 <- rnorm(m*n,2,4) # covariate 2
mu <- beta0 + beta1*mls(x1,0:lag_k_1,m)%*%weights1 + beta2*mls(x2,0:lag_k_2,m)%*%weights2

y <- matrix(NA, nrow = length(mu), ncol = 1)
y[!is.na(mu)] <- rnorm(
  sum(!is.na(mu)),
  mean = mu[!is.na(mu)],
  sd = 1
)

data_gaussian_example <- list(x1 = x1[-seq_len(14)],
                             x2 = x2[-seq_len(14)],
                             y = y[!is.na(lambda)],
                             mu = mu[!is.na(mu)],
                             weights1 = weights1,
                             weights2 = weights2,
                             beta0 = beta0,
                             beta1 = beta1,
                             beta2 = beta2)






y_ts <- ts(data_gaussian_example$y, frequency = 1)
x1_ts <- ts(data_gaussian_example$x1, frequency = 1)
x2_ts <- ts(data_gaussian_example$x2, frequency = 1)


#png("inst/paper/figures/data_binomial_example_covariate1.png", width=17, height=10, units = 'cm', res = 300)
autoplot(x1_ts, color = "brown") + theme_bw() + geom_point(color = "red", size = .01) +
  theme(axis.text.x=element_text(size=16),
        axis.title.x=element_text(size=16,face="bold"),
        axis.text.y=element_text(size=16),
        axis.title.y=element_blank(),
        legend.position = "bottom",
        legend.text=element_text(size=20, face = "plain"),
        legend.title=element_blank())
#dev.off()
#png("inst/paper/figures/data_binomial_example_covariate2.png", width=17, height=10, units = 'cm', res = 300)
autoplot(x2_ts, color = "brown") + theme_bw() + geom_point(color = "red", size = .01) +
  theme(axis.text.x=element_text(size=16),
        axis.title.x=element_text(size=16,face="bold"),
        axis.text.y=element_text(size=16),
        axis.title.y=element_blank(),
        legend.position = "bottom",
        legend.text=element_text(size=20, face = "plain"),
        legend.title=element_blank())
#dev.off()
#png("inst/paper/figures/data_binomial_example_response.png", width=17, height=10, units = 'cm', res = 300)
autoplot(y_ts, color = "black") + theme_bw() + geom_point(color = "gray", size = .01) +
  theme(axis.text.x=element_text(size=16),
        axis.title.x=element_text(size=16,face="bold"),
        axis.text.y=element_text(size=16),
        axis.title.y=element_blank(),
        legend.position = "bottom",
        legend.text=element_text(size=20, face = "plain"),
        legend.title=element_blank())
#dev.off()


data = data.frame(y = data_gaussian_example$y)
data$y_all <- data$y
data$y[(length(data$y)-15):length(data$y)] <- NA # remove last 16 time points as forecast test set


#### Model fitting ####


Midas_x1 <- prepare_Minla(x = data_gaussian_example$x1,
                          constraint = "hyperbolic",
                          K = 0:13,
                          m = 7)
Midas_x2 <- prepare_Minla(x = data_gaussian_example$x2,
                          constraint = "gaussian",
                          K = 0:20,
                          m = 7)


fit_res <- fit_Minla(formula = y ~ 1,
                     data = data,
                     family = "gaussian",
                     hf_input = list(Midas_x1,Midas_x2),
                     inla_options = list(verbose = T))

summary(fit_res$res)


#png("inst/paper/figures/Example_binomial_hyperbolic_paramestimates.png", width=30, height=10, units = 'cm', res = 300)
par(mfrow=c(1,3))

plot(inla.smarginal(fit_res$res$marginals.fixed[["(Intercept)"]]),
     type="l", lwd=3, col="red", xlab=expression(beta[0]), ylab="",
     cex.lab = 2.2, cex.axis=1.5)
abline(v = data_gaussian_example$beta0, col = 'blue', lty=1, lwd = 2)
abline(v = quantile(inla.rmarginal(200, fit_res$res$marginals.fixed$`(Intercept)`), prob = 0.025), lty = 2)
abline(v = quantile(inla.rmarginal(200, fit_res$res$marginals.fixed$`(Intercept)`), prob = 0.975), lty = 2)

plot(inla.smarginal(fit_res$res$marginals.hyperpar[['Theta2 for hf_idx_1']]),
     type="l", lwd=3, col="red", xlab=expression(beta[2]), ylab="",
     cex.lab = 2.2, cex.axis=1.5)
abline(v = data_gaussian_example$beta1, col = 'blue', lty = 1, lwd = 2)
abline(v = quantile(inla.rmarginal(200, fit_res$res$marginals.hyperpar$`Theta2 for hf_idx_1`), prob = 0.025), lty = 2)
abline(v = quantile(inla.rmarginal(200, fit_res$res$marginals.hyperpar$`Theta2 for hf_idx_1`), prob = 0.975), lty = 2)

plot(inla.smarginal(fit_res$res$marginals.hyperpar[['Theta3 for hf_idx_2']]),
     type="l", lwd=3, col="red", xlab=expression(beta[3]), ylab="",
     cex.lab = 2.2, cex.axis=1.5)
abline(v =data_gaussian_example$beta2, col = 'blue', lty = 1, lwd = 2)
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
#dev.off()

#png("inst/paper/figures/data_binomial_example_weights_covariate2.png", width=27, height=12, units = 'cm', res = 300)
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
#dev.off()


#### Compare observed versus predicted values ####

pred_res <- predict_midas(model = fit_res,
                          family = "gaussian",
                          Ntrials = NULL,
                          nsamples = 1000)
both_ts <- ts(data.frame(observed = fit_res$data_final$y_all,
                         predicted = pred_res$computed_y$mean),
              start = 1, end = length(pred_res$computed_y$mean))
PI <- data.frame(lower = pred_res$computed_y$q2.5,
                 upper = pred_res$computed_y$q97.5,
                 Time = 1:length(pred_res$computed_y$mean))
PI$Time <- time(both_ts)




#png("inst/paper/figures/data_binomial_example_predsVSobs.png", width=25, height=12, units = 'cm', res = 300)
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
#dev.off()




#### Fit a baseline approach ####

data_baseline <- data.frame(y = fit_res$data_final$y,
                            x1 = fit_res$hf_input[[1]]$X_matrix[-c(1:2),1],
                            x2 = fit_res$hf_input[[2]]$X_matrix[-c(1:2),1])

res_baseline = inla(y ~ 1 +
                      f(x1, model = "linear") +
                      f(x2, model = "linear"),
                    data = data_baseline,
                    family = "gaussian",
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
  stats::rnorm(
    n = nrow(data_baseline),
    mean = latent_predictor[, i],
    sd = sqrt(1/res_baseline$summary.hyperpar$mean))
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

#png("inst/paper/figures/data_binomial_example_predsVSobs_baseline.png", width=25, height=12, units = 'cm', res = 300)
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
#dev.off()



