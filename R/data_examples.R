#' Example MIDAS dataset
#'
#' Simulated dataset generated using a binomial sampling model and hyperbolic weighting scheme.
#'
#' @description
#' This dataset contains simulated high-frequency covariates and a
#' binomial response constructed using MIDAS lag weights.
#'
#' @format A list with the following components:
#' \describe{
#'   \item{x}{Numeric vector of high-frequency covariates}
#'   \item{y}{Integer vector of binomial response counts}
#'   \item{p}{Underlying success probabilities}
#'   \item{Ntrials}{Number of trials for each observation}
#'   \item{weights}{Lag weights used in the MIDAS structure}
#'   \item{beta0}{True value of intercept \eqn{\beta_0}}
#'   \item{beta1}{True value of \eqn{\beta_1}}
#'   \item{beta2}{True value of \eqn{\beta_2}}
#' }
#'
#' @details
#' The dataset is generated using a hyperbolic weighting scheme
#' with parameter gamma = 0.9 and lag length 13.
#'
#' @source Simulated data
"data_binom_hyperbolic"

#devtools::document()





#' Example MIDAS dataset
#'
#' Simulated dataset generated using a spatial poisson likelihood and gaussian weighting scheme.
#'
#' @description
#' This dataset contains simulated high-frequency covariates and a
#' poisson response constructed using MIDAS lag weights.
#'
#' @format A list with the following components:
#' \describe{
#'   \item{data_x}{Data frame containing covariate data}
#'   \item{data_y}{Data frame containing response data}
#'   \item{weights}{True values of weights}
#'   \item{eta}{True values of linear predictor}
#'   \item{beta0}{True value of intercept \eqn{\beta_0}}
#'   \item{beta1}{True value of \eqn{\beta_1}}
#'   \item{phi_sd}{True value of \eqn{\phi} standard deviation}
#'   \item{phi}{True value of \eqn{\phi} effects}
#' }
#'
#' @details
#' The dataset is generated using a gaussian weighting scheme
#' with parameters mu = 8 and sigma = 7 and lag length 30
#'
#' @source Simulated data
"data_spatialpoisson_gauss"

#devtools::document()





#' Example MIDAS dataset
#'
#' Simulated dataset generated using a spatial poisson likelihood,
#' varying beta, and gaussian weighting scheme.
#'
#' @description
#' This dataset contains simulated high-frequency covariates and a
#' poisson response constructed using MIDAS lag weights. The betas
#' vary across locations.
#'
#' @format A list with the following components:
#' \describe{
#'   \item{data_x}{Data frame containing covariate data}
#'   \item{data_y}{Data frame containing response data}
#'   \item{weights}{True values of weights}
#'   \item{eta}{True values of linear predictor}
#'   \item{beta0}{True value of intercept \eqn{\beta_0}}
#'   \item{beta}{True values of \eqn{\beta_i}}
#'   \item{phi_sd}{True value of \eqn{\phi} standard deviation}
#'   \item{phi}{True value of \eqn{\phi} effects}
#' }
#'
#' @details
#' The dataset is generated using a gaussian weighting scheme
#' with parameters mu = 8 and sigma = 7 and lag length 30
#'
#' @source Simulated data
"data_spatialpoisson_varybeta_gauss"

#devtools::document()
