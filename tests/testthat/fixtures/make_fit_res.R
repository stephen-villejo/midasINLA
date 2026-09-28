

library(midasINLA)
library(INLA)

INLA::inla.setOption(num.threads = 1)

data("data_spatialpoisson_example", package = "midasINLA")

response_data <- data_spatialpoisson_example[["data_y"]]
response_data$y_all <- response_data$y
response_data[
  response_data[["Time"]] %in% 183:192,
  "y"
] <- NA

g <- INLA::inla.read.graph(
  filename = system.file("map.adj", package = "midasINLA")
)

K1 <- 0:29

Midas_x1 <- prepare_Minla_spatial(
  x = data_spatialpoisson_example$data_x1$x1,
  loc_x = data_spatialpoisson_example$data_x1$loc,
  constraint = "hyperbolic",
  K = K1,
  m = 30,
  svc = TRUE,
  svc_prior = "icar",
  g = g
)

K2 <- 0:45

Midas_x2 <- prepare_Minla_spatial(
  x = data_spatialpoisson_example$data_x2$x2,
  loc_x = data_spatialpoisson_example$data_x2$loc,
  constraint = "gaussian",
  K = K2,
  m = 30,
  svc = FALSE
)

hf_input <- list(Midas_x1, Midas_x2)

fit_res <- fit_Minla_spatial(
  formula = y ~ 1,
  data = response_data,
  loc_var = "loc",
  time_var = "Time",
  family = "poisson",
  hf_input = hf_input,
  inla_options = list(verbose = FALSE,
                      control.predictor = list(
                        compute = TRUE,link = 1))
)


saveRDS(
  hf_input,
  file = file.path(
    "tests", "testthat", "fixtures", "hf_input.rds"
  )
)

saveRDS(
  fit_res,
  file = file.path(
    "tests", "testthat", "fixtures", "fit_res.rds"
  )
)

