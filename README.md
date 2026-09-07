<img width="468" height="67" alt="image" src="https://github.com/user-attachments/assets/2fda6668-6b65-4f1e-b670-e23e2784a5ff" />
<!-- README.md is generated from README.Rmd. Please edit that file -->

# midasINLA

<!-- badges: start -->

<!-- badges: end -->

`midasINLA` provides tools for fitting Mixed Data Sampling (MIDAS) regression models with the Integrated Nested Laplace Approximation (INLA). The package is designed for settings where a response is observed at a lower frequency than one or more explanatory variables, and supports both global and spatially varying MIDAS coefficients.

## Overview

MIDAS models allow high-frequency covariates to be incorporated into models for lower-frequency responses through weighted distributed lags. midasINLA combines this framework with INLA, allowing MIDAS regression models to be fitted efficiently within a latent Gaussian modelling framework.


## Installation

The development version can be installed from GitHub using `remotes`:

``` r
install.packages("remotes")
remotes::install_github("stephen-villejo/midasINLA")
```

## Getting started

The following example uses the spatial Poisson dataset included with the package.

``` r
library(midasINLA)
library(INLA)

data("data_spatialpoisson_example")

# Read the spatial adjacency graph
g <- INLA::inla.read.graph(
  filename = system.file("map.adj", package = "midasINLA")
)

# Prepare the first high-frequency covariate
Midas_x1 <- prepare_Minla_spatial(
  x = data_spatialpoisson_example$data_x1$x1,
  loc_x = data_spatialpoisson_example$data_x1$loc,
  constraint = "hyperbolic",
  K = 0:29,
  m = 30,
  svc = TRUE,
  svc_prior = "icar",
  g = g
)

# Prepare the second high-frequency covariate
Midas_x2 <- prepare_Minla_spatial(
  x = data_spatialpoisson_example$data_x2$x2,
  loc_x = data_spatialpoisson_example$data_x2$loc,
  constraint = "gaussian",
  K = 0:45,
  m = 30,
  svc = FALSE
)

# Create a response with the final observations held out
response_data <- data_spatialpoisson_example$data_y
response_data$y_all <- response_data$y
response_data[response_data$Time %in% 183:192,"y"] <- NA

# Fit the model
fit_res <- fit_Minla_spatial(
  formula = y ~ 1,
  data = response_data,
  loc_var = "loc",
  time_var = "Time",
  family = "poisson",
  hf_input = list(Midas_x1, Midas_x2),
  inla_options = list(
    verbose = FALSE,
    num.threads = 1
  )
)

# Posterior summaries
summary(fit_res$res)

# MIDAS coefficient summaries
beta_results <- compute_beta_spatial(
  model = fit_res,
  n_loc = 16
)

# MIDAS lag weights
weights <- compute_weights(fit_res)

# Predictions
pred_res <- predict_midas(
  model = fit_res,
  family = "poisson",
  Ntrials = NULL,
  nsamples = 1000
)

```

For a complete walkthrough, including model specification, interpretation of MIDAS weights, spatially varying coefficients, and prediction, see the package vignette:

``` r
vignette("midasINLA", package = "midasINLA")
```
