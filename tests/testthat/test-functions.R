


test_that("package functions are available", {
  expect_true(is.function(prepare_Minla_spatial))
  expect_true(is.function(fit_Minla_spatial))
  expect_true(is.function(compute_beta_spatial))
  expect_true(is.function(compute_weights))
  expect_true(is.function(predict_midas))
})


test_that("prepare_Minla_spatial works for a non-spatial MIDAS model", {

  data(
    data_spatialpoisson_example,
    package = "midasINLA"
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

  expect_type(Midas_x2, "list")

  expect_true(
    all(c(
      "X_matrix",
      "constraint",
      "lag_k",
      "K",
      "m",
      "rm.row",
      "svc",
      "svc_prior"
    ) %in% names(Midas_x2))
  )

  expect_equal(Midas_x2$constraint, "gaussian")
  expect_equal(Midas_x2$lag_k, 45)
  expect_equal(Midas_x2$K, K2)
  expect_equal(Midas_x2$m, 30)
  expect_false(Midas_x2$svc)

  expect_true(is.matrix(Midas_x2$X_matrix))
  expect_true("loc" %in% colnames(Midas_x2$X_matrix))
})


test_that("precomputed vignette results have expected structure", {

  result_file <- system.file(
    "extdata",
    "vignette_results.rds",
    package = "midasINLA"
  )

  skip_if(
    !nzchar(result_file),
    "Precomputed vignette results not available"
  )

  results <- readRDS(result_file)

  beta_results <- results$beta_results
  w <- results$res_weights

  # beta results
  expect_type(beta_results, "list")

  expected_beta_names <- c(
    "hf_index_1",
    "hf_index_2"
  )

  expect_setequal(
    names(beta_results),
    expected_beta_names
  )

  # weights
  expect_type(w, "list")

  expect_true(all(c("hf_1", "hf_2") %in% names(w)))

  expect_true(is.data.frame(w$hf_1))
  expect_true(is.data.frame(w$hf_2))

  expect_true(
    all(c("lag", "mean", "q2.5", "q97.5") %in% names(w$hf_1))
  )

  expect_true(
    all(c("lag", "mean", "q2.5", "q97.5") %in% names(w$hf_2))
  )

  expect_true(all(is.finite(w$hf_1$mean)))
  expect_true(all(is.finite(w$hf_2$mean)))

  expect_equal(sum(w$hf_1$mean), 1, tolerance = 1e-8)
  expect_equal(sum(w$hf_2$mean), 1, tolerance = 1e-8)

  expect_true(all(w$hf_1$q2.5 <= w$hf_1$mean))
  expect_true(all(w$hf_1$mean <= w$hf_1$q97.5))

  expect_true(all(w$hf_2$q2.5 <= w$hf_2$mean))
  expect_true(all(w$hf_2$mean <= w$hf_2$q97.5))
})
