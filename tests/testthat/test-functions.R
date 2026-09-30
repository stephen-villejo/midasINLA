


test_that("package functions are available", {
  expect_true(is.function(prepare_Minla_spatial))
  expect_true(is.function(fit_Minla_spatial))
  expect_true(is.function(compute_beta_spatial))
  expect_true(is.function(compute_weights))
  expect_true(is.function(predict_midas))
})


test_that("downstream functions work", {


  skip_if_not_installed("INLA")

  INLA::inla.setOption(num.threads = 1)

  hf_input <- readRDS(testthat::test_path("fixtures", "hf_input.rds"))

  fit_file <- system.file(
    "extdata",
    "fit_res.rds",
    package = "midasINLA"
  )

  skip_if(!nzchar(fit_file), "Precomputed fit not available")

  fit_res <- readRDS(fit_file)

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
