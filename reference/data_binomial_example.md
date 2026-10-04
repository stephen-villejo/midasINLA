# Example MIDAS dataset

This dataset contains simulated high-frequency covariates and a binomial
response constructed using MIDAS lag weights.

## Usage

``` r
data_binomial_example
```

## Format

A list with the following components:

- x1:

  Numeric vector of high-frequency covariate 1

- x2:

  Numeric vector of high-frequency covariate 2

- y:

  Integer vector of binomial response counts

- p:

  Underlying success probabilities

- Ntrials:

  Number of trials for each observation

- weights1:

  Lag weights used in the MIDAS structure for covariate 1

- weights2:

  Lag weights used in the MIDAS structure for covariate 2

- beta0:

  True value of intercept \\\beta_0\\

- beta1:

  True value of \\\beta_1\\

- beta2:

  True value of \\\beta_2\\

- beta3:

  True value of \\\beta_3\\

## Source

Simulated data

## Details

Simulated dataset generated using a binomial sampling model

The dataset is generated using two covariates. The first one has a
hyperbolic weighting scheme with parameter gamma = 0.9 and lag length
13. The second one has gaussian weighting scheme with parameters mu = 8,
sigma = 7, and lag length of 20.
