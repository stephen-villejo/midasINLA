rm(list=ls())

library(midasr)

#
##### Hyperbolic Scheme Polynomial #####
#

# model parameters

n <- 260
m <- 7
beta0 = -2.5
beta1 = 0.01 # trend coefficient
beta2 = 1.2
beta3 = -1

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
x1 <- rnorm(m*n,2,3) # covariate 1
x2 <- rnorm(m*n,0,3) # covariate 2
eta <- beta0 + beta1*c(1:n) + beta2*mls(x1,0:lag_k_1,m)%*%weights1 + beta3*mls(x2,0:lag_k_2,m)%*%weights2
p <- as.vector(exp(eta)/(1+exp(eta)))
Ntrials <- sample(50:150, size = n, replace = TRUE)

y <- matrix(NA, nrow = length(p), ncol = 1)
y[!is.na(p)] <- rbinom(
  sum(!is.na(p)),
  size = Ntrials[!is.na(p)],
  prob = p[!is.na(p)]
)

data_binomial_example <- list(x1 = x1[-seq_len(14)],
                              x2 = x2[-seq_len(14)],
                              y = y[!is.na(p)],
                              p = p[!is.na(p)],
                              Ntrials = Ntrials[!is.na(p)],
                              weights1 = weights1,
                              weights2 = weights2,
                              beta0 = beta0,
                              beta1 = beta1,
                              beta2 = beta2,
                              beta3 = beta3)


usethis::use_data(data_binomial_example, overwrite = TRUE)


#source("data-raw/data_binomial_example.R")
