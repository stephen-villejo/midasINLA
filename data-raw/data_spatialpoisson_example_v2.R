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
beta1 = 1.1
beta2 = -2

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


library(ggplot2)

grid_sf$icar <- as.numeric(icar)

ggplot(grid_sf) +
  geom_sf(aes(fill = icar)) +
  geom_sf_text(aes(label = id)) +
  theme_minimal()

beta_i <- icar


# compute weights

gamma_param <- 0.9
mu_val <- 10
sigma_val <- 12

par(mfrow=c(1,1))

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

  set.seed(423 + i)
  x1 <- rnorm(m*n,2,2)
  x2 <- rnorm(m*n,0,1)
  eta <- beta0 + (beta1 + beta_i[i])*mls(x1,0:lag_k1,m)%*%weights1 +
    beta2*mls(x2,0:lag_k2,m)%*%weights2
  mu <- exp(eta)

  y <- matrix(NA, nrow = length(mu), ncol = 1)
  set.seed(753 + i)
  y[!is.na(mu)] <- rpois(
    sum(!is.na(mu)),
    lambda = mu[!is.na(mu)]
  )


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


data_spatialpoisson_example_v2 <- list(data_x1 = data_x1,
                                       data_x2 = data_x2,
                                       data_y = data_y,
                                       weights1 = weights1,
                                       weights2 = weights2,
                                       eta = eta,
                                       beta0 = beta0,
                                       beta1 = beta1,
                                       beta2 = beta2,
                                       icar = icar,
                                       tau = tau,
                                       grid_sf = grid_sf)


usethis::use_data(data_spatialpoisson_example_v2, overwrite = TRUE)


#source("data-raw/data_spatialpoisson_example_v2.R")

