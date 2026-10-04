# Compute posterior estimates of MIDAS lag weights

Computes posterior estimates of the normalized MIDAS lag weights from a
fitted MIDAS model returned by
[`fit_Minla_spatial()`](https://stephen-villejo.github.io/midasINLA/reference/fit_Minla_spatial.md).
The function draws samples from the posterior marginal distributions of
the MIDAS hyperparameters and uses these samples to obtain the
corresponding lag-weight functions.

## Usage

``` r
compute_weights(model, n.samples = 200)
```

## Arguments

- model:

  A fitted MIDAS model returned by
  [`fit_Minla_spatial()`](https://stephen-villejo.github.io/midasINLA/reference/fit_Minla_spatial.md).

- n.samples:

  Positive integer specifying the number of posterior samples drawn from
  the MIDAS hyperparameter marginal distributions to estimate the
  lag-weight distribution. Defaults to `200`.

## Value

A list containing one data frame for each high-frequency covariate in
`model$hf_input`. The elements are named `"hf_1"`, `"hf_2"`, and so on.
Each data frame contains:

- lag:

  The MIDAS lag index.

- mean:

  The posterior mean of the normalized lag weight.

- q2.5:

  The 2.5% posterior quantile of the lag weight.

- q97.5:

  The 97.5% posterior quantile of the lag weight.

## Details

The supported MIDAS lag constraints are `"hyperbolic"`, `"gaussian"`,
`"beta1"`, `"beta2"`, and `"almon2"`. For each high-frequency covariate,
the resulting weights are normalized to sum to one across all included
lags.

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("INLA", quietly = TRUE)) {
  data(data_spatialpoisson_example)

  g <- INLA::inla.read.graph(
    filename = system.file("map.adj", package = "midasINLA")
  )

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

  weights <- compute_weights(
    model = fit,
    n.samples = 200
  )

  head(weights$hf_1)
}
} # }
```
