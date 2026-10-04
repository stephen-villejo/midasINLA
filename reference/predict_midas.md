# Generate posterior predictive samples from a MIDAS model

Generates posterior predictive samples from a fitted MIDAS model using
posterior samples of the latent linear predictor. The function supports
Gaussian, Poisson, and binomial response distributions.

## Usage

``` r
predict_midas(model, family = "gaussian", Ntrials = NULL, nsamples = 1000)
```

## Arguments

- model:

  A fitted MIDAS model returned by
  [`fit_Minla_spatial()`](https://stephen-villejo.github.io/midasINLA/reference/fit_Minla_spatial.md).

- family:

  Character string specifying the likelihood family. Supported values
  are `"gaussian"`, `"poisson"`, and `"binomial"`. Defaults to
  `"gaussian"`.

- Ntrials:

  Optional vector specifying the number of trials for each observation
  when `family = "binomial"`. If omitted, the function attempts to use
  the `Ntrials` column in `model$data_final`.

- nsamples:

  Positive integer specifying the number of posterior samples used to
  generate the predictive distribution. Defaults to `1000`.

## Value

A list with two components:

- computed_y:

  A list containing posterior predictive summaries:

  mean

  : Posterior predictive mean for each observation.

  sd

  : Posterior predictive standard deviation for each observation.

  q2.5

  : 2.5% posterior predictive quantile for each observation.

  q97.5

  : 97.5% posterior predictive quantile for each observation.

- samples:

  A list containing posterior samples of the latent predictor and, for
  Gaussian responses, the posterior samples of the observation variance.

## Details

For a Gaussian response, posterior samples of the observation variance
are obtained from the posterior samples of the Gaussian precision
parameter and used to generate predictive observations. For Poisson and
binomial responses, observations are generated using the appropriate
inverse-link function.

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("INLA", quietly = TRUE)) {
  data(data_spatialpoisson_example)

  Midas_x1 <- prepare_Minla_spatial(
    x = data_spatialpoisson_example$data_x1$x1,
    loc_x = data_spatialpoisson_example$data_x1$loc,
    constraint = "hyperbolic",
    K = 0:29,
    m = 30,
    svc = FALSE
  )

  fit <- fit_Minla_spatial(
    formula = y ~ 1,
    data = data_spatialpoisson_example$data_y,
    loc_var = "loc",
    time_var = "Time",
    family = "poisson",
    hf_input = list(Midas_x1),
    inla_options = list(verbose = FALSE)
  )

  predictions <- predict_midas(
    model = fit,
    family = "poisson",
    nsamples = 1000
  )

  # Inspect posterior predictive summaries
  head(predictions$computed_y$mean)
  head(predictions$computed_y$q2.5)
  head(predictions$computed_y$q97.5)
}
} # }
```
