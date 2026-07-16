rm(list = ls())

library(midasr)
library(sf)
library(spdep)
library(MASS)
library(INLA)
library(midasINLA)

#
##### Hyperbolic Scheme Polynomial #####
#

#Create a square boundary

boundary <- st_as_sf(
  st_sfc(
    st_polygon(list(rbind(
      c(0, 0),
      c(4, 0),
      c(4, 4),
      c(0, 4),
      c(0, 0)
    )))
  )
)

# Generate a 4 x 4 grid
grid <- st_make_grid(boundary, n = c(4, 4))

# Convert to an sf object
grid_sf <- st_sf(
  id = 1:length(grid),
  geometry = grid
)

# Plot
plot(grid_sf["id"])

# model parameters

locs <- nrow(grid_sf)
n <- 12 * locs
m <- 30
lag_k1 <- 29
lag_k2 <- 45
beta0 = 1
beta1 = -1.1
beta2 = 2

# simulate my spatially varying betas from an iCAR model
# Construct queen adjacency
nb <- poly2nb(grid_sf, queen = TRUE)
# Adjacency matrix
W <- nb2mat(nb, style = "B", zero.policy = TRUE)
# ICAR precision matrix (tau = 1)
D <- Diagonal(x = rowSums(W))
Q <- D - W

tau <- 10
# eigen decomposition
eig <- eigen(as.matrix(Q))
# remove null space (constant vector)
keep <- eig$values > 1e-10
V <- eig$vectors[, keep]
Lambda <- eig$values[keep]
# ICAR covariance on constrained space
Sigma <- V %*% diag(1 / (tau * Lambda)) %*% t(V)

set.seed(123)
icar <- MASS::mvrnorm(1, rep(0, nrow(Q)), Sigma)
# enforce sum-to-zero (numerical safety)
icar <- icar - mean(icar)
grid_sf$icar <- as.numeric(icar)
plot(grid_sf["icar"])

beta_i <- icar


# compute weights

gamma_param <- 0.9
mu_val <- 10
sigma_val <- 12

xi <- c()
for(lag in 0:lag_k1){
  temp <- gamma(lag+gamma_param)/(gamma(lag+1)*gamma(gamma_param))
  xi <- c(xi, temp)
}
weights1 <- xi/sum(xi)
plot(weights1, type = "l")

xi <- c()
for(lag in 0:lag_k2){
  temp <- exp( -(lag-mu_val)^2 / (2 * (sigma_val^2)) )
  xi <- c(xi, temp)
}
weights2 <- xi/sum(xi)
plot(weights2, type = "l")

# compute covariate, Ntrials, and linear predictor

x1_list <- vector("list", locs)
x2_list <- vector("list", locs)
y_list <- vector("list", locs)

for(i in 1:locs){

  set.seed(9912 + i)
  x1 <- rnorm(m*n,3,2)
  x2 <- rnorm(m*n,1,1)
  eta <- beta0 + (beta1 + beta_i[i])*mls(x1,0:lag_k1,m)%*%weights1 +
    beta2*mls(x2,0:lag_k2,m)%*%weights2
  mu <- exp(eta)

  y <- matrix(NA, nrow = length(mu), ncol = 1)
  set.seed(9912 + i)
  y[!is.na(mu)] <- rpois(
    sum(!is.na(mu)),
    lambda = mu[!is.na(mu)]
  )

  #y <- y[-which(is.na(y))]
  #x <- x[-c(1:30)]

  x1_list[[i]] <- data.frame(
    x1 = x1,
    loc = i,
    Time = seq_along(x1)
  )
  x2_list[[i]] <- data.frame(
    x2 = x2,
    loc = i,
    Time = seq_along(x2)
  )

  y_list[[i]] <- data.frame(
    y = y,
    loc = i,
    Time = seq_along(y)
  )

}


data_x1 <- do.call(rbind, x1_list)
data_x2 <- do.call(rbind, x2_list)
data_y <- do.call(rbind, y_list)


data_spatialpoisson_example <- list(data_x1 = data_x1,
                                    data_x2 = data_x2,
                                    data_y = data_y,
                                    weights1 = weights1,
                                    weights2 = weights2,
                                    eta = eta,
                                    beta0 = beta0,
                                    beta1 = beta1,
                                    beta2 = beta2,
                                    icar = icar,
                                    tau = tau)





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
                                  constraint = "gaussian",
                                  K = 0:45,
                                  m = 30,
                                  svc = FALSE)



fit_res <- fit_Minla_spatial(formula = y ~ 1,
                             data = response_data,
                             loc_var = "loc",
                             time_var = "Time",
                             family = "poisson",
                             hf_input = list(Midas_x1, Midas_x2),
                             inla_options = list(verbose = T))

summary(fit_res$res)


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


