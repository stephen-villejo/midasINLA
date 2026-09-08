


test_that("package functions are available", {
  expect_true(is.function(prepare_Minla_spatial))
  expect_true(is.function(fit_Minla_spatial))
  expect_true(is.function(compute_beta_spatial))
  expect_true(is.function(compute_weights))
  expect_true(is.function(predict_midas))
})


test_that("compute_weights returns valid normalized weights", {
  skip_if_not_installed("INLA")

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

  expect_equal(ncol(Midas_x1$X_matrix), length(K1) + 1)

  K2 <- 0:45

  Midas_x2 <- prepare_Minla_spatial(
    x = data_spatialpoisson_example$data_x2$x2,
    loc_x = data_spatialpoisson_example$data_x2$loc,
    constraint = "gaussian",
    K = K2,
    m = 30,
    svc = FALSE
  )

  expect_equal(ncol(Midas_x2$X_matrix), length(K2) + 1)

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

  hyperpar <- summary(fit_res$res)$hyperpar
  hyperpar_names <- rownames(hyperpar)

  expected_hf <- paste0(
    "hf_idx_",
    seq_along(hf_input)
  )

  for (hf in expected_hf) {
    expect_true(
      any(grepl(hf, hyperpar_names, fixed = TRUE)),
      info = paste("No hyperparameters found for", hf)
    )
  }

  actual_hf <- unique(
    sub(
      ".*(hf_idx_[0-9]+).*",
      "\\1",
      hyperpar_names[grepl("hf_idx_[0-9]+", hyperpar_names)]
    )
  )

  expect_setequal(actual_hf, expected_hf)


  # compute beta estimates


  beta_results <- compute_beta_spatial(model = fit_res,
                                       n_loc = 16)

  expected_beta_names <- paste0(
    "hf_index_",
    seq_along(hf_input)
  )

  expect_setequal(
    names(beta_results),
    expected_beta_names
  )


  # compute weights

  w <- compute_weights(fit_res)

  expect_type(w, "list")
  expect_true("hf_1" %in% names(w))
  expect_true("hf_2" %in% names(w))

  expect_true(is.data.frame(w$hf_1))
  expect_true(is.data.frame(w$hf_2))

  expect_true(all(is.finite(w$hf_1$mean)))
  expect_true(all(is.finite(w$hf_2$mean)))

  expect_equal(
    sum(w$hf_1$mean),
    1,
    tolerance = 1e-8
  )

  expect_equal(
    sum(w$hf_2$mean),
    1,
    tolerance = 1e-8
  )

  expect_true(all(c("lag", "mean", "q2.5", "q97.5") %in%
                    names(w$hf_1)))

  expect_true(all(c("lag", "mean", "q2.5", "q97.5") %in%
                    names(w$hf_2)))

  expect_true(all(w$hf_1$q2.5 <= w$hf_1$mean))
  expect_true(all(w$hf_1$mean <= w$hf_1$q97.5))

  expect_true(all(w$hf_2$q2.5 <= w$hf_2$mean))
  expect_true(all(w$hf_2$mean <= w$hf_2$q97.5))

})
