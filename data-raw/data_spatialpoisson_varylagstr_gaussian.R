

library(midasr)

#
##### Hyperbolic Scheme Polynomial #####
#

# model parameters

locs <- 5
n <- 12 * locs
m <- 30
beta0 = 1
beta = c(0.9, 1, 1.2, 1.3,1.5)

set.seed(5109)
phi_sd <- 0.5
phi <- rnorm(locs,0,phi_sd)

lag_k <- rep(29, times = 5)

# compute weights

mu_f <- c(2,4,6,8,10)
sigma_f <- rep(7, times = 5)

weights <- list()
for(i in 1:locs){
  xi <- c()
  for(lag in 0:lag_k[i]){
    temp <- exp( -(lag-mu_f[i])^2 / (2 * (sigma_f[i]^2)) )
    xi <- c(xi, temp)
  }
  weights[[i]] <- xi/sum(xi)
}

par(mfrow = c(2,2))
plot(weights[[1]], type = "l")
plot(weights[[2]], type = "l")
plot(weights[[3]], type = "l")
plot(weights[[5]], type = "l")

# compute covariate, Ntrials, and linear predictor

x_list <- vector("list", locs)
y_list <- vector("list", locs)


for(i in 1:locs){

  set.seed(9912 + i)
  x <- rnorm(m*n,3,2)
  eta <- beta0 + beta[i]*mls(x,0:lag_k[i],m)%*%weights[[i]] + phi[i]
  mu <- exp(eta)

  y <- matrix(NA, nrow = length(mu), ncol = 1)
  set.seed(9912 + i)
  y[!is.na(mu)] <- rpois(
    sum(!is.na(mu)),
    lambda = mu[!is.na(mu)]
  )

  #y <- y[-which(is.na(y))]
  #x <- x[-c(1:30)]

  x_list[[i]] <- data.frame(
    x = x,
    loc = i,
    Time = seq_along(x)
  )

  y_list[[i]] <- data.frame(
    y = y,
    loc = i,
    Time = seq_along(y)
  )

}


data_x <- do.call(rbind, x_list)
data_y <- do.call(rbind, y_list)


data_spatialpoisson_varylagstr_gauss <- list(data_x = data_x,
                                             data_y = data_y,
                                             weights = weights,
                                             eta = eta,
                                             beta0 = beta0,
                                             beta = beta,
                                             phi_sd = phi_sd,
                                             phi = phi)


usethis::use_data(data_spatialpoisson_varylagstr_gauss, overwrite = TRUE)


#source("data-raw/data_spatialpoisson_varylagstr_gaussian.R")

