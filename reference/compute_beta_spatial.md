# Compute posterior summaries of MIDAS coefficients

Computes posterior summaries of the MIDAS regression coefficients from a
fitted spatial MIDAS model returned by
[`fit_Minla_spatial()`](https://stephen-villejo.github.io/midasINLA/reference/fit_Minla_spatial.md).
The output depends on whether the MIDAS coefficient is spatially varying
and on the specified spatial prior.

## Usage

``` r
compute_beta_spatial(model, n_loc)
```

## Arguments

- model:

  A fitted spatial MIDAS model returned by
  [`fit_Minla_spatial()`](https://stephen-villejo.github.io/midasINLA/reference/fit_Minla_spatial.md).

- n_loc:

  Integer specifying the number of spatial locations for which
  coefficient summaries should be computed.

## Value

A list containing one element for each high-frequency MIDAS covariate in
`model$hf_input`. The elements are named `"hf_index_1"`, `"hf_index_2"`,
and so on. The contents of each element depend on the spatial structure
of the corresponding MIDAS covariate:

- Spatially varying coefficient with ICAR prior:

  A list containing `summary.global.beta`, `marginal.global.beta`,
  `summary.icar.beta`, `marginal.icar.beta`, `summary.total.beta`, and
  `sample.total.beta`.

- Spatially varying coefficient with IID prior:

  A list containing `summary.beta` and `marginal.beta`.

- Non-spatially varying coefficient:

  A list containing `summary.beta` and `marginal.beta`.

The summary data frames contain the posterior mean, standard deviation,
and 2.5%, 50%, and 97.5% posterior quantiles.

## Details

For spatially varying coefficients with an ICAR prior, the function
returns summaries and marginal distributions for the global coefficient,
the location-specific spatial deviations, and the resulting
location-specific total coefficients. Posterior samples of the total
coefficients are also returned.

For spatially varying coefficients with an IID prior, the function
returns posterior summaries and marginal distributions for the
location-specific coefficients.

For non-spatially varying coefficients, the function returns the
posterior marginal distribution and summary of the MIDAS coefficient
associated with each high-frequency covariate.

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("INLA", quietly = TRUE)) {
  INLA::inla.setOption(num.threads = 1)

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

  # Compute posterior summaries of the MIDAS coefficients
  beta_summary <- compute_beta_spatial(
    model = fit,
    n_loc = 16
  )

  # Inspect summaries of the location-specific total coefficients
  beta_summary$hf_index_1$summary.total.beta
}
} # }
```
