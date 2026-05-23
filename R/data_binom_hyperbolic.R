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
#' }
#'
#' @details
#' The dataset is generated using a hyperbolic weighting scheme
#' with parameter gamma = 0.9 and lag length 13.
#'
#' @source Simulated data
"data_binom_hyperbolic"

#devtools::document()
