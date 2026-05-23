
library(midasr)

#
##### Hyperbolic Scheme Polynomial #####
#

# model parameters

n <- 260
m <- 7
beta0 = -2.5
trend = 0.01
beta1 = 1.2

gamma_param <- 0.9
lag_k <- 13 # 0:lag_k, this implies (lag_k + 1) high freq covariate associated with each response y

# compute weights

xi <- c()
for(lag in 0:lag_k){
  temp <- gamma(lag+gamma_param)/(gamma(lag+1)*gamma(gamma_param))
  xi <- c(xi, temp)
}
weights <- xi/sum(xi)
plot(weights, type = "l")

# compute covariate, Ntrials, and linear predictor

set.seed(9912)
x <- rnorm(m*n,2,3) # covariate data
eta <- beta0 + trend*c(1:n) + beta1*mls(x,0:lag_k,m)%*%weights
p <- as.vector(exp(eta)/(1+exp(eta)))
Ntrials <- sample(50:150, size = n, replace = TRUE)

y <- matrix(NA, nrow = length(p), ncol = 1)
y[!is.na(p)] <- rbinom(
  sum(!is.na(p)),
  size = Ntrials[!is.na(p)],
  prob = p[!is.na(p)]
)

data_binom_hyperbolic <- list(x = x,
                              y = y,
                              p = p,
                              Ntrials = Ntrials,
                              weights = weights,
                              beta0 = beta0,
                              beta1 = beta1,
                              trend = trend)

usethis::use_data(data_binom_hyperbolic, overwrite = TRUE)


#source("data-raw/data_binom_hyperbolic.R")

