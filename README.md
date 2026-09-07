
<!-- README.md is generated from README.Rmd. Please edit that file -->

# midasINLA

<!-- badges: start -->

<!-- badges: end -->

midasINLA provides tools for fitting Mixed Data Sampling (MIDAS) regression models with the Integrated Nested Laplace Approximation (INLA). The package is designed for settings where a response is observed at a lower frequency than one or more explanatory variables, and supports both global and spatially varying MIDAS coefficients.<img width="468" height="95" alt="image" src="https://github.com/user-attachments/assets/a53ebc07-7466-4446-b076-7b443a8bf213" />


## Installation

You can install the development version of midasINLA from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("stephen-villejo/midasINLA")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(midasINLA)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" alt="" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
