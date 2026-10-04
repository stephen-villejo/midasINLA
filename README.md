# midasINLA

<!-- badges: start -->

<!-- badges: end -->

`midasINLA` provides tools for fitting mixed-data sampling (MIDAS) regression models using Integrated Nested Laplace Approximation (INLA). The package is designed for settings where the response is observed at a lower frequency than one or more explanatory variables, and supports both constant and spatially varying regression coefficients.

## Overview

MIDAS models allow high-frequency covariates to be incorporated into models for lower-frequency responses through weighted distributed lags. midasINLA combines this framework with INLA, allowing MIDAS regression models to be fitted efficiently within a latent Gaussian modelling framework.

**Conceptual framework for spatial distributed-lag MIDAS modelling.**  
  High-frequency covariates are linked to lower-frequency responses through
weighted distributed lags, with MIDAS coefficients potentially varying across
spatial locations.

<p align="center">
  <img src="man/figures/schematic_paper.png" width="70%">
</p>

## Installation

The latest released version of midasINLA can be installed from CRAN:

```r
install.packages("midasINLA")
```

The development version can be installed from GitHub using `remotes`:

``` r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

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
    num.threads = 1,
    control.predictor = list(
          compute = TRUE,
          link = 1
        )
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

For a complete walkthrough, see the package vignette:

``` r
vignette("midasINLA", package = "midasINLA")
```

## Model specification

A MIDAS model incorporates high-frequency covariates through weighted lagged terms. For example, for a response $y_{it}$, a high-frequency covariate $x$ can enter the model through
```math
\sum_{k=0}^{K} w_k x_{i,s(t)-k}
```
where $w_k$ denotes the MIDAS weight associated with lag $k$, and $s(t)$ denotes the high-frequency time index corresponding to $y_{it}$.

The weights are parameterised using constraint functions, allowing the lag-response pattern to be estimated with substantially fewer parameters than an unconstrained distributed-lag model.

The MIDAS effect can additionally vary across spatial locations, with spatial variation modelled using either an iid or an intrinsic conditional autoregressive (iCAR) prior.

## Documentation

The package vignette provides a detailed introduction to the modelling workflow:
``` r
vignette("midasINLA", package = "midasINLA")
```

Individual functions can also be explored using R’s help system:
`?prepare_Minla_spatial`
`?fit_Minla_spatial`
`?compute_beta_spatial`
`?compute_weights`
`?predict_midas`

## Applications

midasINLA was developed to facilitate Bayesian modelling of relationships between variables observed at different temporal resolutions, with particular emphasis on applications involving spatial data.

The package is particularly useful in environmental epidemiology, climate–health modelling, and other settings where high-frequency exposures need to be linked to lower-frequency outcomes.

## Status

`midasINLA` is currently under development.

## Citation

If you use `midasINLA` in your research, please cite the package and the associated methodological work once the relevant publication is available.





