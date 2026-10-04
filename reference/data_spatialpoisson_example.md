# Example MIDAS dataset

This dataset contains simulated high-frequency covariates and a poisson
response constructed using MIDAS lag weights.

## Usage

``` r
data_spatialpoisson_example
```

## Format

A list with the following components:

- data_x1:

  High-frequency covariate 1 data frame

- data_x2:

  High-frequency covariate 2 data frame

- data_y:

  Response data frame

- weights1:

  Lag weights used in the MIDAS structure for covariate 1

- weights2:

  Lag weights used in the MIDAS structure for covariate 2

- eta:

  Linear predictor log lambda

- beta0:

  True value of intercept \\\beta_0\\

- beta1:

  True value of \\\beta_1\\

- beta2:

  True value of \\\beta_2\\

- icar:

  Vector of region-specific deviations from \\\beta_1\\

- tau:

  Precision parameter of the ICAR model

## Source

Simulated data

## Details

Simulated dataset generated using a spatial poisson sampling model

The dataset is generated using two covariates. The first one has a
hyperbolic weighting scheme with parameter gamma = 0.9 and lag length
29. The second one has Gaussian constraint with parameters mu = 10 and
sigma = 12.
