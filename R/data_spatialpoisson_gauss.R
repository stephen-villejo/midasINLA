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
