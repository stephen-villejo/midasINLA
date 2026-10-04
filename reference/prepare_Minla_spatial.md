# Prepare a spatial MIDAS object for INLA estimation

Constructs the MIDAS design matrix and associated model specifications
for use with
[`fit_Minla_spatial()`](https://stephen-villejo.github.io/midasINLA/reference/fit_Minla_spatial.md).
The function creates lagged high-frequency covariate values for each
location and optionally specifies a spatially varying coefficient (SVC)
component.

## Usage

``` r
prepare_Minla_spatial(
  x,
  loc_x,
  constraint,
  K,
  m,
  svc = FALSE,
  svc_prior = "iid",
  g = NULL
)
```

## Arguments

- x:

  Numeric vector of high-frequency covariate observations.

- loc_x:

  Vector identifying the location associated with each observation in
  `x`. The length of `loc_x` must equal the length of `x`.

- constraint:

  Character string specifying the constraint used for the MIDAS
  lag-response association. Supported constraints include
  `"hyperbolic"`, `"gaussian"`, `"beta1"`, `"beta2"`, and `"almon2"`.

- K:

  Numeric vector specifying the lags to be included in the MIDAS
  representation.

- m:

  Numeric vector specifying the number of high-frequency covariate
  observations associated with each response observation. A single value
  can be supplied when the number of observations is constant over time,
  or a vector can be supplied when this number varies across response
  times.

- svc:

  Logical; if `TRUE`, specifies a spatially varying coefficient model.
  Defaults to `FALSE`.

- svc_prior:

  Character string specifying the prior for the spatially varying
  coefficient component. Must be either `"icar"` or `"iid"`. Defaults to
  `"iid"`.

- g:

  An INLA graph object used for the `"icar"` prior. Required when
  `svc_prior = "icar"`.

## Value

A list containing the MIDAS design matrix and model specifications. The
returned object includes:

- X_matrix:

  The MIDAS design matrix, including a location index.

- constraint:

  The MIDAS constraint used for the lag-response association.

- lag_k:

  The maximum lag specified in `K`.

- K:

  The vector of lags used to construct the design matrix.

- m:

  The number of high-frequency observations associated with each
  response observation.

- rm.row:

  The number of initial rows removed from each location because of
  incomplete lagged observations.

- svc:

  Whether a spatially varying coefficient component is specified.

- svc_prior:

  The prior specified for the spatially varying coefficient component.

- g:

  The INLA graph object, included when `svc_prior = "icar"`.

## Examples

``` r
if (requireNamespace("INLA", quietly = TRUE)) {
  data(data_spatialpoisson_example)

  # Prepare a MIDAS object using a hyperbolic lag constraint
  # and a spatially varying coefficient with an ICAR prior.
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

  # Inspect the resulting MIDAS design matrix
  head(Midas_x1$X_matrix)
}
#>            lag0        lag1       lag2       lag3       lag4     lag5     lag6
#> [1,]  3.3476186  0.01817311  0.3681445  2.3016830  0.5067161 4.408850 1.909036
#> [2,]  2.4864753  4.06072089  5.1092570  2.7308815  3.7797569 1.103207 2.333270
#> [3,]  1.9115982 -1.33000836  0.9097436  5.5770802 -0.2284686 3.655859 1.085481
#> [4,]  3.2369918  0.67073035 -0.3463683  1.5294842  3.1199743 1.160617 3.414520
#> [5,] -0.6901928  6.46649267  1.2216607 -0.9564444  1.6109662 1.859666 4.117308
#> [6,]  3.4602856  1.85359448  4.5318530  0.8639463  4.0706130 2.605277 3.418157
#>            lag7       lag8      lag9     lag10       lag11     lag12    lag13
#> [1,] -0.5614291  2.3128501 0.7916824 1.5002006 -1.18624790 0.8533801 4.508304
#> [2,]  4.1363923  2.6400443 1.8004433 0.6138136  3.13127227 5.0012469 3.851161
#> [3,]  3.0006452  3.7296628 3.4812855 3.7257405  1.49294182 4.7230649 3.237976
#> [4,]  1.8423261  2.5860938 4.2605953 3.4196465  0.76643148 1.1544946 1.542278
#> [5,]  1.3754550  6.7709315 0.5750773 3.9950856 -0.07849181 4.0278543 2.197294
#> [6,]  0.2346903 -0.4553831 1.4213427 1.6819284 -2.44594168 5.5570183 4.084163
#>           lag14     lag15      lag16       lag17     lag18    lag19      lag20
#> [1,]  3.8724709 7.0140009 -1.9837898 -1.02440548 1.7938687 2.257050 -1.1134601
#> [2,]  1.9162171 3.7875422  3.8491306  2.97534561 2.9917505 2.022454  4.8039730
#> [3,]  1.8599673 2.4572730 -1.5404110  1.03967652 2.6341253 1.000411  1.3400006
#> [4,]  0.5440659 2.7045564  3.6376599  2.75351264 0.7318604 1.621922  0.5517963
#> [5,]  7.6717373 4.5855835 -0.1195613 -0.06711773 0.9373837 3.486460  2.4655293
#> [6,] -0.2521852 0.5551746  1.2628123  0.34729417 0.9935365 4.112260  1.0543692
#>           lag21     lag22       lag23      lag24     lag25     lag26
#> [1,]  6.5242674 -1.522246 -0.06789194  0.1984434 2.3612870 0.3694894
#> [2,] -0.7897482  1.419727  1.74709201  1.0057677 0.4601563 2.8873097
#> [3,]  1.4696951  2.749143  3.60456228 -0.6737688 3.2214203 1.3291661
#> [4,]  4.7079697  3.864386  3.10587044 -0.2165249 1.6472653 4.5915974
#> [5,]  1.4005922  4.007637 -0.88910391  5.3616978 2.1111839 2.2114814
#> [6,]  3.4279198  6.526910  2.16798238  0.4524545 2.1833393 4.9359112
#>            lag27        lag28     lag29 loc
#> [1,]  0.56899585  0.022331169 5.4273026   1
#> [2,]  2.55395439 -0.573808337 1.8285856   1
#> [3,]  0.03075021  2.511625857 0.6536545   1
#> [4,]  2.08972004  1.949383782 3.6419254   1
#> [5,] -0.07553904  0.008450714 1.3235851   1
#> [6,]  2.14586541 -2.816771919 1.7349923   1
```
