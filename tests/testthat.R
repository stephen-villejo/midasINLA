library(testthat)
library(midasINLA)

INLA::inla.setOption(num.threads = 1)

test_check("midasINLA")
