
rm(list = ls())

library(midasINLA)
library(ggplot2)
library(midasr)
library(INLA)
library(ltc)


data("data_spatialpoisson_example", package = "midasINLA")
str(data_spatialpoisson_example)

pal=ltc("paloma",5,"continuous")
a <- ggplot(data_spatialpoisson_example$grid_sf) +
  geom_sf(aes(fill = icar)) +
  scale_fill_gradientn(colours = pal,
                       name = expression(b[i])) +
  geom_sf_text(
    aes(label = id),
    size = 2.5
  ) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 11, face = "bold"),
        #legend.key.height = unit(1, "cm"),
        #legend.key.width = unit(1.5, "cm"),
        legend.position = "left")
a

df_weights <- data.frame(
  lag = seq_along(data_spatialpoisson_example$weights1),
  weight = data_spatialpoisson_example$weights1
)
b <- ggplot(df_weights, aes(x = lag, y = weight)) +
  geom_line(linewidth = 0.3) +
  geom_point(size = .3) +
  labs(
    x = "Lag",
    y = "Weight",
    title = expression(w[1*k])
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))

df_weights <- data.frame(
  lag = seq_along(data_spatialpoisson_example$weights2),
  weight = data_spatialpoisson_example$weights2
)
c <- ggplot(df_weights, aes(x = lag, y = weight)) +
  geom_line(linewidth = 0.3) +
  geom_point(size = .3) +
  labs(
    x = "Lag",
    y = "Weight",
    title = expression(w[2*k])
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))


library(patchwork)
final_plot <- a + (b / c) +
  plot_layout(widths = c(.6, 1))
png("inst/paper/figures/data_spatialpoisson_example.png", width=17, height=10, units = 'cm', res = 300)
final_plot
dev.off()





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


response_data <- data_spatialpoisson_example[["data_y"]]
response_data$y_all <- response_data$y
response_data[which(response_data[["Time"]] %in% 183:192),"y"] <- NA






#### Model fitting ####


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
                                  constraint = "gaussian",
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
summary(fit_res$res)


beta_results <- compute_beta_spatial(model = fit_res,
                                     n_loc = 16)

beta_results$hf_index_1$summary.global.beta
beta_results$hf_index_1$summary.icar.beta
beta_results$hf_index_2$summary.beta

png("inst/paper/figures/Example_spatialpoisson_betaestimates.png", width=20, height=9, units = 'cm', res = 300)

par(mgp = c(4, 1, 0))
par(mar = c(5, 4, 2, 1),
    mgp = c(4, 1, 0))
layout(matrix(c(1, 2, 3, 4, 4, 4),
              nrow = 2,
              byrow = TRUE),
       heights = c(1, 0.2))

# Panel 1
plot(inla.smarginal(fit_res$res$marginals.fixed[["(Intercept)"]]),
     type = "l", lwd = 3, col = "red",
     xlab = expression(beta[0]), ylab = "",
     cex.lab = 2.2, cex.axis = 1.5)

abline(v = data_spatialpoisson_example$beta0,
       col = "blue", lty = 1, lwd = 2)
abline(v = fit_res$res$summary.fixed["(Intercept)", "mean"],
       col = "black", lty = 1, lwd = 2)
abline(v = fit_res$res$summary.fixed["(Intercept)", "0.025quant"],
       lty = 2)
abline(v = fit_res$res$summary.fixed["(Intercept)", "0.975quant"],
       lty = 2)

# Panel 2
plot(inla.smarginal(beta_results$hf_index_1$marginal.global.beta),
     type = "l", lwd = 3, col = "red",
     xlab = expression(beta[1]^"*"), ylab = "",
     cex.lab = 2.2, cex.axis = 1.5)
abline(v = data_spatialpoisson_example$beta1,
       col = "blue", lty = 1, lwd = 2)
abline(v = beta_results$hf_index_1$summary.global.beta$Mean,
       col = "black", lty = 1, lwd = 2)
abline(v = beta_results$hf_index_1$summary.global.beta$`2.5%`,
       lty = 2)
abline(v = beta_results$hf_index_1$summary.global.beta$`97.5%`,
       lty = 2)

# Panel 3
plot(inla.smarginal(beta_results$hf_index_2$marginal.beta),
     type = "l", lwd = 3, col = "red",
     xlab = expression(beta[2]), ylab = "",
     cex.lab = 2.2, cex.axis = 1.5)

abline(v = data_spatialpoisson_example$beta2,
       col = "blue", lty = 1, lwd = 2)
abline(v = beta_results$hf_index_2$summary.beta$Mean,
       col = "black", lty = 1, lwd = 2)
abline(v = beta_results$hf_index_2$summary.beta$`2.5%`,
       lty = 2)
abline(v = beta_results$hf_index_2$summary.beta$`97.5%`,
       lty = 2)

# Legend panel
par(mar = c(0, 0, 0, 0))
plot.new()

legend("center",
       legend = c("True value",
                  "Posterior mean",
                  "95% credible interval"),
       col = c("blue", "black", "black"),
       lty = c(1, 1, 2),
       lwd = c(2, 2, 1),
       bty = "n",
       horiz = TRUE,
       cex = 1.7)

dev.off()


#### Compute weights ####

res_weights <- compute_weights(fit_res)
str(res_weights)


png("inst/paper/figures/Example_spatialpoisson_weights_covariate1.png", width=30, height=15, units = 'cm', res = 300)
res_weights$hf_1$lag <- factor(res_weights$hf_1$lag, levels=unique(res_weights$hf_1$lag))
ggplot(res_weights$hf_1, aes(x = lag, y = mean)) +
  geom_errorbar(
    aes(
      ymin = q2.5,
      ymax = q97.5,
      color = "95% credible interval"
    ),
    width = .7,
    position = position_dodge(0.05)
  ) +
  geom_point(
    aes(color = "Posterior mean"),
    size = 1
  ) +
  geom_point(
    aes(
      x = lag,
      y = data_spatialpoisson_example$weights1,
      color = "True value"
    ),
    size = 1
  ) +
  scale_color_manual(
    name = "",
    breaks = c(
      "Posterior mean",
      "95% credible interval",
      "True value"
    ),
    values = c(
      "Posterior mean" = "red",
      "95% credible interval" = "black",
      "True value" = "blue"
    )
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        size = c(3, 1, 3)
      )
    )
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 16, face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(size = 16, face = "plain"),
    legend.title = element_blank()
  ) +
  ylab("")
dev.off()


png("inst/paper/figures/Example_spatialpoisson_weights_covariate2.png", width=30, height=15, units = 'cm', res = 300)
res_weights$hf_2$lag <- factor(res_weights$hf_2$lag, levels=unique(res_weights$hf_2$lag))
ggplot(res_weights$hf_2, aes(x = lag, y = mean)) +
  geom_errorbar(
    aes(
      ymin = q2.5,
      ymax = q97.5,
      color = "95% credible interval"
    ),
    width = .7,
    position = position_dodge(0.05)
  ) +
  geom_point(
    aes(color = "Posterior mean"),
    size = 1
  ) +
  geom_point(
    aes(
      x = lag,
      y = data_spatialpoisson_example$weights2,
      color = "True value"
    ),
    size = 1
  ) +
  scale_color_manual(
    name = "",
    breaks = c(
      "Posterior mean",
      "95% credible interval",
      "True value"
    ),
    values = c(
      "Posterior mean" = "red",
      "95% credible interval" = "black",
      "True value" = "blue"
    )
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        size = c(3, 1, 3)
      )
    )
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 16, face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(size = 16, face = "plain"),
    legend.title = element_blank()
  ) +
  ylab("")
dev.off()



#### Compare observed versus predicted values ####

pred_res <- predict_midas(model = fit_res,
                          family = "poisson",
                          Ntrials = NULL,
                          nsamples = 1000)
str(pred_res)
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



png("inst/paper/figures/Example_spatialpoisson_forecast.png", width=35, height=20, units = 'cm', res = 300)
ggplot(both_long_sub, aes(x = Time, y = value, colour = series)) +
  geom_line() +
  geom_ribbon(
    data = PI_sub,
    aes(
      x = Time,
      ymin = lower,
      ymax = upper,
      fill = "95% credible interval"
    ),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  theme_minimal() +
  facet_wrap(
    ~loc,
    ncol = 2,
    labeller = labeller(loc = function(x) paste("Loc =", x))
  ) +
  scale_colour_manual(
    values = c("blue", "red")
  ) +
  scale_fill_manual(
    name = "",
    values = c("95% credible interval" = "red")
  ) +
  ylab("y") +
  geom_vline(
    data = vlines_df_sub,
    aes(xintercept = cut),
    colour = "black",
    size = 1
  ) +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16, face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(size = 20, face = "plain"),
    legend.title = element_blank(),
    strip.text = element_text(size = 20)
  )
dev.off()




#### Model comparator ####



