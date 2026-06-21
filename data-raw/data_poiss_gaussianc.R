
library(midasr)

#
##### Gaussian Polynomial #####
#

# model parameters

n <- 60
m <- 30
beta0 = 1
beta1 = 1.2
beta2 = -0.5

mu_f <- 8
sigma_f <- 7
lag_k <- 29 # 0:lag_k, this implies (lag_k + 1) high freq covariate associated with each response y

# compute weights

xi <- c()
for(lag in 0:lag_k){
  temp <- exp( -(lag-mu_f)^2 / (2 * (sigma_f^2)) )
  xi <- c(xi, temp)
}
weights <- xi/sum(xi)
#plot(weights, type = "l")

# compute covariate, Ntrials, and linear predictor

set.seed(9912)
x1 <- rnorm(m*n,3,2) # covariate data
x2 <- rnorm(n,0,1) # low frequency covariate
eta <- beta0 + beta1*mls(x1,0:lag_k,m)%*%weights + beta2*x2
mu <- exp(eta)

y <- matrix(NA, nrow = length(mu), ncol = 1)
set.seed(9912)
y[!is.na(mu)] <- rpois(
  sum(!is.na(mu)),
  lambda = mu[!is.na(mu)]
)

data_poiss_gaussianc <- list(x1 = x1,
                             x2 = x2,
                             y = y,
                             mu = mu,
                             weights = weights,
                             beta0 = beta0,
                             beta1 = beta1,
                             beta2 = beta2)


usethis::use_data(data_poiss_gaussianc, overwrite = TRUE)


#source("data-raw/data_poiss_gaussianc.R")



