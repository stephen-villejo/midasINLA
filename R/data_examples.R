#' Example MIDAS dataset
#'
#' Simulated dataset generated using a binomial sampling model
#'
#' @description
#' This dataset contains simulated high-frequency covariates and a
#' binomial response constructed using MIDAS lag weights.
#'
#' @format A list with the following components:
#' \describe{
#'   \item{x1}{Numeric vector of high-frequency covariate 1}
#'   \item{x2}{Numeric vector of high-frequency covariate 2}
#'   \item{y}{Integer vector of binomial response counts}
#'   \item{p}{Underlying success probabilities}
#'   \item{Ntrials}{Number of trials for each observation}
#'   \item{weights1}{Lag weights used in the MIDAS structure for covariate 1}
#'   \item{weights2}{Lag weights used in the MIDAS structure for covariate 2}
#'   \item{beta0}{True value of intercept \eqn{\beta_0}}
#'   \item{beta1}{True value of \eqn{\beta_1}}
#'   \item{beta2}{True value of \eqn{\beta_2}}
#'   \item{beta3}{True value of \eqn{\beta_3}}
#' }
#'
#' @details
#' The dataset is generated using two covariates. The first one has a hyperbolic
#' weighting scheme with parameter gamma = 0.9 and lag length 13. The second one
#' has gaussian weighting scheme with parameters mu = 8, sigma = 7, and
#' lag length of 20.
#'
#' @source Simulated data
"data_binomial_example"

#devtools::document()

