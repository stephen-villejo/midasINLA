# Fit a spatial MIDAS model using INLA

Fits a Mixed Data Sampling (MIDAS) regression model with optional
spatially varying coefficients using Integrated Nested Laplace
Approximation (INLA). High-frequency covariates are supplied as MIDAS
objects created by
[`prepare_Minla_spatial()`](https://stephen-villejo.github.io/midasINLA/reference/prepare_Minla_spatial.md).

## Usage

``` r
fit_Minla_spatial(
  formula,
  data,
  loc_var,
  time_var,
  family,
  hf_input = NULL,
  Ntrials = NULL,
  E = NULL,
  inla_options = list()
)
```

## Arguments

- formula:

  A model formula specifying the response and other covariates. The
  response variable must be a column in `data`.

- data:

  A data frame containing the response and any additional model
  covariates. It must contain the variables specified by `loc_var` and
  `time_var`.

- loc_var:

  Character string specifying the name of the location variable in
  `data`.

- time_var:

  Character string specifying the name of the time variable in `data`.

- family:

  Character string specifying the likelihood family for the response.
  For example, `"poisson"` or `"binomial"`.

- hf_input:

  A list of MIDAS objects returned by
  [`prepare_Minla_spatial()`](https://stephen-villejo.github.io/midasINLA/reference/prepare_Minla_spatial.md).
  Each object specifies a high-frequency covariate and its MIDAS lag
  structure. Multiple MIDAS objects can be supplied.

- Ntrials:

  Optional vector specifying the number of trials for a binomial
  response. Used only when `family = "binomial"`.

- E:

  Optional vector of expected counts or exposure values for a Poisson
  model. Its length must match the number of rows in `data`. Used only
  when `family = "poisson"`.

- inla_options:

  A named list of additional arguments passed to
  [`INLA::inla()`](https://rdrr.io/pkg/INLA/man/inla.html). By default,
  `control.compute$config` is set to `TRUE` if it is not already
  specified.

## Value

A list containing:

- formula_final:

  The final INLA model formula, including the MIDAS components.

- data_final:

  The response data used for model fitting after ordering by location
  and time and removing rows with incomplete lagged covariate
  information.

- rm_max:

  The maximum number of initial observations removed across locations
  because of incomplete MIDAS lagged covariates.

- res:

  The fitted INLA model returned by
  [`INLA::inla()`](https://rdrr.io/pkg/INLA/man/inla.html).

- hf_input:

  The list of MIDAS objects supplied through `hf_input`.

## Details

The function constructs the INLA model by incorporating the MIDAS
components specified in `hf_input`. Multiple high-frequency covariates
can be included by supplying multiple MIDAS objects in `hf_input`.

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("INLA", quietly = TRUE)) {
  data(data_spatialpoisson_example)

  # Read the spatial adjacency graph
  g <- INLA::inla.read.graph(
    filename = system.file("map.adj", package = "midasINLA")
  )

  # Prepare a spatial MIDAS predictor
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

  # Fit the spatial Poisson MIDAS model
  fit <- fit_Minla_spatial(
    formula = y ~ 1,
    data = data_spatialpoisson_example$data_y,
    loc_var = "loc",
    time_var = "Time",
    family = "poisson",
    hf_input = list(Midas_x1),
    inla_options = list(
      verbose = FALSE,
      control.predictor = list(
        compute = TRUE,
        link = 1
      )
    )
  )

  # Inspect the fitted INLA model
  fit$res
}
} # }
```
