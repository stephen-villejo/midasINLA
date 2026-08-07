


#' Create lag matrix
#'
#' @param tsdata High-frequency covariate data
#' @param lags Lags to be considered
#' @param frequency Number of high-frequency covariate values associated
#'   with each response. Can be a single value or a vector.
#' @return Matrix of lagged values
#' @export
create_lag_Xmatrix <- function(tsdata,
                               lags,
                               frequency) {

  tsdata <- as.numeric(tsdata)

  if (any(frequency <= 0) || any(frequency != as.integer(frequency))) {
    stop("'frequency' must contain positive integers.")
  }

  # Scalar frequency
  if (length(frequency) == 1) {

    if (length(tsdata) %% frequency != 0) {
      stop(
        "When 'frequency' is a scalar, its value must divide ",
        "the length of 'tsdata'."
      )
    }

    frequency <- rep(
      frequency,
      length(tsdata) / frequency
    )

    # Vector frequency
  } else {

    if (sum(frequency) != length(tsdata)) {
      stop(
        "When 'frequency' is a vector, its sum must equal ",
        "the length of 'tsdata'."
      )
    }
  }

  end_idx <- cumsum(frequency)

  X_matrix <- matrix(
    NA_real_,
    nrow = length(frequency),
    ncol = length(lags)
  )

  for (i in seq_along(end_idx)) {
    for (j in seq_along(lags)) {

      idx <- end_idx[i] - lags[j]

      if (idx >= 1 && idx <= length(tsdata)) {
        X_matrix[i, j] <- tsdata[idx]
      }
    }
  }

  colnames(X_matrix) <- paste0("lag", lags)

  X_matrix
}


#' Prepare MIDAS objects for INLA estimation
#' @param x High-frequency covariate data
#' @param constraint constraint function for lag-response association
#' @param K Lags to be considered
#' @param m Number of covariates values associated with each response
#' @return A list of objects that will be used for model fitting
#' @export
prepare_Minla <- function(x,
                          constraint,
                          K,
                          m) {

  X_matrix <- create_lag_Xmatrix(tsdata = x,
                                 lags = K,
                                 frequency = m)

  temp_data <- as.data.frame(X_matrix)

  bad_rows <- which(!stats::complete.cases(temp_data))

  rm.row <- if (length(bad_rows) == 0) integer(0) else max(bad_rows)

  return(list(
    X_matrix = X_matrix,
    constraint = constraint,
    lag_k = max(K),
    K = K,
    m = m,
    rm.row = rm.row
  ))
}


#' Fit MIDAS model using INLA
#' @param formula model formula
#' @param data dataframe
#' @param family Likelihood family for response data
#' @param hf_input list of output objects from prepare_Minla
#' @param Ntrials Number of trials for a binomial family response
#' @param inla_options arguments in inla function
#' @return MIDAS output
#' @export
fit_Minla <- function(formula,
                      data,
                      family,
                      hf_input = NULL,
                      Ntrials = data$Ntrials,
                      inla_options = list()) {


  build_rgen <- function(x, constraint) {
    if (constraint == "beta") {
      INLA::inla.rgeneric.define(model = rgeneric.Beta.midas, x = x)
    } else if (constraint == "beta2") {
      INLA::inla.rgeneric.define(model = rgeneric.Beta2.midas, x = x)
    } else if (constraint == "almon2") {
      INLA::inla.rgeneric.define(model = rgeneric.Almon2.midas, x = x)
    } else if (constraint == "almon3") {
      INLA::inla.rgeneric.define(model = rgeneric.Almon3.midas, x = x)
    } else if (constraint == "hyperbolic") {
      INLA::inla.rgeneric.define(model = rgeneric.Hyperbolic.midas, x = x)
    } else if (constraint == "gaussian") {
      INLA::inla.rgeneric.define(model = rgeneric.Gaussian.midas, x = x)
    } else {
      stop("Unknown constraint.")
    }
  }

  response_name <- all.vars(formula[[2]])

  if (family == "binomial") {
    if ("Ntrials" %in% names(data)) {
      data$Ntrials <- NULL
      data$Ntrials <- Ntrials
    }else{
      data$Ntrials <- Ntrials
    }
  }

  get_rm_row <- function(obj) {
    if (length(obj$rm.row) == 0) {
      0
    } else {
      obj$rm.row
    }
  }

  rm_max <- if (is.null(hf_input) || length(hf_input) == 0) {
    0
  } else {
    max(vapply(hf_input, get_rm_row, numeric(1)))
  }

  data_final <- if (rm_max > 0) {
    data[-seq_len(rm_max), , drop = FALSE]
  } else {
    data
  }



  formula_final <- formula
  final_n <- nrow(data_final)

  formula_env <- environment(formula_final)
  if (is.null(formula_env)) {
    formula_env <- parent.frame()
  }

  if (!is.null(hf_input) && length(hf_input) > 0) {
    stopifnot(is.list(hf_input))

    for (i in seq_along(hf_input)) {
      obj <- hf_input[[i]]

      idx_name <- paste0("hf_idx_", i)
      model_name <- paste0("hf_model_", i)

      temp_data <- if (rm_max > 0) {
        cbind(obj$X_matrix[-seq_len(rm_max), , drop = FALSE],data_final[,which(names(data_final) == response_name)])
      } else {
        cbind(obj$X_matrix,data_final[,which(names(data_final) == response_name)])
      }
      temp_data <- as.data.frame(temp_data)
      rgen = build_rgen(x = temp_data, constraint = obj$constraint)

      data_final[[idx_name]] <- seq_len(final_n)
      assign(model_name, rgen, envir = formula_env)

      term_txt <- sprintf(
        "f(%s, model = %s, n = %d)",
        idx_name, model_name, final_n
      )

      formula_final <- stats::update(formula_final, paste(". ~ . +", term_txt))
    }
  }

  environment(formula_final) <- formula_env

  args_inla <- list(
    formula = formula_final,
    data = data_final,
    family = family,
    control.compute = list(config = TRUE)
  )

  if (family == "binomial") {
    args_inla$Ntrials <- data_final$Ntrials
  }

  args_inla <- c(args_inla, inla_options)

  res <- do.call(INLA::inla, args_inla)


  return(list(formula_final = formula_final,
              data_final = data_final,
              rm_max = rm_max,
              res = res,
              hf_input = hf_input))

}



#' Make predictions form a MIDAS model output
#' @param model MIDAS model output
#' @param family  Likelihood family for response data
#' @param Ntrials Number of trials for a binomial family response
#' @param nsamples Number of posterior samples
#' @return A list containing the predictions and the samples
#' @export
predict_midas <- function(model,
                          family = "gaussian",
                          Ntrials = NULL,
                          nsamples = 1000) {

  fit <- if (!is.null(model$res)) model$res else model

  if (is.null(data)) {
    if (!is.null(model$data_final)) {
      data <- model$data_final
    } else {
      stop("`data` must be supplied if not available in `model$data_final`.")
    }
  }

  data <- model$data_final
  n_obs <- nrow(data)

  temp <- INLA::inla.posterior.sample(n = nsamples, fit)

  latent_names <- rownames(temp[[1]]$latent)

  predictor_idx <- grep("^Predictor", latent_names)

  if (length(predictor_idx) < n_obs) {
    stop("Could not identify enough Predictor entries in posterior samples.")
  }

  predictor_idx <- predictor_idx[seq_len(n_obs)]

  latent_predictor <- sapply(seq_len(nsamples), function(i) {
    temp[[i]]$latent[predictor_idx]
  })

  if (family == "gaussian") {

    sigma2 <- sapply(seq_len(nsamples), function(i) {
      1 / unname(temp[[i]]$hyperpar[
        which(names(temp[[i]]$hyperpar) == "Precision for the Gaussian observations")
      ])
    })

    samples <- list(
      latent_predictor = latent_predictor,
      sigma2 = sigma2
    )

    sample_y <- sapply(seq_len(nsamples), function(i) {
      latent_predictor[, i] +
        stats::rnorm(n_obs, mean = 0, sd = sqrt(sigma2[i]))
    })

  } else if (family == "poisson") {

    samples <- list(
      latent_predictor = latent_predictor
    )

    sample_y <- sapply(seq_len(nsamples), function(i) {
      stats::rpois(n = n_obs, lambda = exp(latent_predictor[, i]))
    })

  } else if (family == "binomial") {

    if (is.null(Ntrials)) {
      if ("Ntrials" %in% names(data)) {
        Ntrials <- data$Ntrials
      } else {
        stop("`Ntrials` must be supplied for binomial predictions.")
      }
    }

    samples <- list(
      latent_predictor = latent_predictor
    )

    sample_y <- sapply(seq_len(nsamples), function(i) {
      stats::rbinom(
        n = n_obs,
        size = Ntrials,
        prob = exp(latent_predictor[, i]) / (1 + exp(latent_predictor[, i]))
      )
    })

  } else {
    stop("Unsupported family.")
  }

  computed_y <- list(
    mean = rowMeans(sample_y),
    sd = apply(sample_y, 1, stats::sd),
    q2.5 = matrixStats::rowQuantiles(sample_y, probs = 0.025),
    q97.5 = matrixStats::rowQuantiles(sample_y, probs = 0.975)
  )

  return(list(
    computed_y = computed_y,
    samples = samples
  ))
}


#' Compute weights estimates
#' @param model midas model result
#' @param n.samples number of posterior samples
#' @return estimates of the weights
#' @export
compute_weights <- function(model, n.samples = 200) {

  if (is.null(model$res)) {
    stop("`model$res` not found. `model` must be the output of `fit_Minla()`.")
  }

  if (is.null(model$hf_input) || length(model$hf_input) == 0) {
    stop("`model$hf_input` is empty or missing.")
  }

  fit <- model$res
  hf_input <- model$hf_input

  out <- vector("list", length(hf_input))

  for (i in seq_along(hf_input)) {

    obj <- hf_input[[i]]

    if (is.null(obj$constraint)) {
      stop(sprintf("`constraint` missing in hf_input[[%d]].", i))
    }

    if (is.null(obj$lag_k)) {
      stop(sprintf("`lag_k` missing in hf_input[[%d]].", i))
    }

    constraint <- obj$constraint
    lag_k <- obj$lag_k

    theta1_name <- sprintf("Theta1 for hf_idx_%d", i)
    theta2_name <- sprintf("Theta2 for hf_idx_%d", i)
    theta3_name <- sprintf("Theta3 for hf_idx_%d", i)

    if (!(theta1_name %in% names(fit$marginals.hyperpar))) {
      stop(sprintf("Hyperparameter `%s` not found.", theta1_name))
    }

    if (constraint == "hyperbolic") {

      theta1_samples <- INLA::inla.rmarginal(
        n.samples,
        fit$marginals.hyperpar[[theta1_name]]
      )

      gamma_samples <- exp(theta1_samples) / (1 + exp(theta1_samples))

      psi_mat <- sapply(0:lag_k, function(lag) {
        gamma(lag + gamma_samples) /
          (gamma(lag + 1) * gamma(gamma_samples))
      })

      w_mat <- psi_mat / rowSums(psi_mat)

    } else if (constraint == "gaussian") {

      if (!(theta2_name %in% names(fit$marginals.hyperpar))) {
        stop(sprintf("Hyperparameter `%s` not found.", theta2_name))
      }

      theta1_samples <- INLA::inla.rmarginal(
        n.samples,
        fit$marginals.hyperpar[[theta1_name]]
      )

      theta2_samples <- INLA::inla.rmarginal(
        n.samples,
        fit$marginals.hyperpar[[theta2_name]]
      )

      mu_val_samples <- lag_k * (1 / (1 + exp(-theta1_samples)))
      sigma_val_samples <- exp(theta2_samples)

      psi_mat <- sapply(0:lag_k, function(lag) {
        exp(-(lag - mu_val_samples)^2 / (2 * sigma_val_samples^2))
      })

      w_mat <- psi_mat / rowSums(psi_mat)

    } else if (constraint == "beta1"){

      theta1_samples <- INLA::inla.rmarginal(
        n.samples,
        fit$marginals.hyperpar[[theta1_name]]
      )

      gamma2_val_samples <- exp(theta1_samples) + 1
      xi <- 1e-4
      psi_mat <- sapply(0:lag_k, function(lag) {
        1 + (1-(xi + (1 - 2 * xi) * (lag / lag_k)))^(gamma2_val_samples-1)
      })

      w_mat <- psi_mat / rowSums(psi_mat)

    } else {
      stop(sprintf("Constraint `%s` not yet implemented.", constraint))
    }

    out[[i]] <- data.frame(
      lag = 0:lag_k,
      mean = colMeans(w_mat),
      q2.5 = apply(w_mat, 2, stats::quantile, probs = 0.025),
      q97.5 = apply(w_mat, 2, stats::quantile, probs = 0.975)
    )
  }

  names(out) <- paste0("hf_", seq_along(hf_input))

  return(out)
}





#' Prepare MIDAS objects (spatial case) for INLA estimation
#' @param x High-frequency covariate data
#' @param loc_x index for locations
#' @param constraint constraint function for lag-response association
#' @param K Lags to be considered
#' @param m Number of covariates values associated with each response
#' @param svc TRUE/FALSE, whether a spatially varying coefficent model or not
#' @param svc_prior Prior for the svc component, either "icar" or "iid"
#' @param g graph for the icar prior
#' @return A list of objects that will be used for model fitting
#' @export
prepare_Minla_spatial <- function(x,
                                  loc_x,
                                  constraint,
                                  K,
                                  m,
                                  svc = FALSE,
                                  svc_prior = "iid",
                                  g = NULL) {

  if (length(x) != length(loc_x)) {
    stop("`x` and `loc_x` must have the same length.")
  }

  if (svc_prior == "icar" && is.null(g)) {
    stop("`g` must be supplied when `svc_prior = 'icar'`.")
  }

  compile_X_matrix <- NULL
  unique_loc_x <- unique(loc_x)
  rm.row.loc <- c()

  for(i in unique_loc_x){

    temp_xdata <- x[which(loc_x == i)]

    X_matrix <- create_lag_Xmatrix(tsdata = temp_xdata,
                                   lags = K,
                                   frequency = m)

    bad_rows_i <- which(!stats::complete.cases(X_matrix))
    rm_i <- if (length(bad_rows_i) == 0) 0 else max(bad_rows_i)
    rm.row.loc <- c(rm.row.loc, rm_i)

    X_matrix <- cbind(X_matrix,loc = rep(i, nrow(X_matrix)))

    compile_X_matrix <- rbind(compile_X_matrix,
                              X_matrix)

  }

  if (length(unique(rm.row.loc)) != 1) {
    stop("`rm.row` is not the same across locations for this covariate.")
  }
  rm.row <- unique(rm.row.loc)

  out <- list(
    X_matrix = compile_X_matrix,
    constraint = constraint,
    lag_k = max(K),
    K = K,
    m = m,
    rm.row = rm.row,
    svc = svc,
    svc_prior = svc_prior
  )

  if (svc_prior == "icar") {
    if (is.null(g)) {
      stop("`g` must be supplied when `svc_prior = 'icar'`.")
    }
    out$g <- g
  }

  return(out)

}



#' Fit MIDAS model (spatial case) using INLA
#' @param formula model formula
#' @param data response data frame
#' @param loc_var variable name for location in response data frame
#' @param time_var variable name for time in response data frame
#' @param family Likelihood family for response data
#' @param hf_input list of output objects from prepare_Minla_spatial
#' @param Ntrials Number of trials for a binomial family response
#' @param E expected cases for Poisson model
#' @param inla_options arguments in inla function
#' @return MIDAS output
#' @export
fit_Minla_spatial <- function(formula,
                              data,
                              loc_var,
                              time_var,
                              family,
                              hf_input = NULL,
                              Ntrials = NULL,
                              E = NULL,
                              inla_options = list()) {


  build_rgen <- function(temp_data, hf_info) {

    if (is.null(temp_data)) {
      stop("`hf_info$temp_data` is missing.")
    }
    if (is.null(hf_info$constraint)) {
      stop("`hf_info$constraint` is missing.")
    }
    if (is.null(hf_info$svc)) {
      stop("`hf_info$svc` is missing.")
    }

    if(!hf_info$svc)
    {
      if (hf_info$constraint == "beta1") {
        return(INLA::inla.rgeneric.define(model = rgeneric.globalbeta.Beta1.midas,
                                          x = temp_data))
      } else if (hf_info$constraint == "beta2") {
        return(INLA::inla.rgeneric.define(model = rgeneric.globalbeta.Beta2.midas,
                                          x = temp_data))
      } else if (hf_info$constraint == "almon2") {
        return(INLA::inla.rgeneric.define(model = rgeneric.globalbeta.Almon2.midas,
                                          x = temp_data))
      } else if (hf_info$constraint == "almon3") {
        return(INLA::inla.rgeneric.define(model = rgeneric.globalbeta.Almon3.midas,
                                          x = temp_data))
      } else if (hf_info$constraint == "hyperbolic") {
        return(INLA::inla.rgeneric.define(model = rgeneric.globalbeta.Hyperbolic.midas,
                                          x = temp_data))
      } else if (hf_info$constraint == "gaussian") {
        return(INLA::inla.rgeneric.define(model = rgeneric.globalbeta.Gaussian.midas,
                                          x = temp_data))
      } else {
        stop("Unknown constraint.")
      }

    } else {

      if (is.null(hf_info$svc_prior)) {
        stop("`hf_info$svc_prior` is missing.")
      }

      if(hf_info$svc_prior == "icar"){

        if (is.null(hf_info$g)) {
          stop("`hf_info$g` is required when `svc_prior = 'icar'`.")
        }

        if (hf_info$constraint == "beta1") {
          return(INLA::inla.rgeneric.define(model = rgeneric.svc.Beta1.midas.icar,
                                            x = temp_data,
                                            g = hf_info$g))
        } else if (hf_info$constraint == "beta2") {
          return(INLA::inla.rgeneric.define(model = rgeneric.svc.Beta2.midas.icar,
                                            x = temp_data,
                                            g = hf_info$g))
        } else if (hf_info$constraint == "almon2") {
          return(INLA::inla.rgeneric.define(model = rgeneric.svc.Almon2.midas.icar,
                                            x = temp_data,
                                            g = hf_info$g))
        } else if (hf_info$constraint == "almon3") {
          return(INLA::inla.rgeneric.define(model = rgeneric.svc.Almon3.midas.icar,
                                            x = temp_data,
                                            g = hf_info$g))
        } else if (hf_info$constraint == "hyperbolic") {
          return(INLA::inla.rgeneric.define(model = rgeneric.svc.Hyperbolic.midas.icar,
                                            x = temp_data,
                                            g = hf_info$g))
        } else if (hf_info$constraint == "gaussian") {
          return(INLA::inla.rgeneric.define(model = rgeneric.svc.Gaussian.midas.icar,
                                            x = temp_data,
                                            g = hf_info$g))
        } else {
          stop("Unknown constraint.")
        }

      } else if(hf_info$svc_prior == "iid") {

        if (hf_info$constraint == "beta1") {
          return(INLA::inla.rgeneric.define(model = rgeneric.svc.Beta1.midas.iid,
                                            x = temp_data))
        } else if (hf_info$constraint == "beta2") {
          return(INLA::inla.rgeneric.define(model = rgeneric.svc.Beta2.midas.iid,
                                            x = temp_data))
        } else if (hf_info$constraint == "almon2") {
          return(INLA::inla.rgeneric.define(model = rgeneric.svc.Almon2.midas.iid,
                                            x = temp_data))
        } else if (hf_info$constraint == "almon3") {
          return(INLA::inla.rgeneric.define(model = rgeneric.svc.Almon3.midas.iid,
                                            x = temp_data))
        } else if (hf_info$constraint == "hyperbolic") {
          return(INLA::inla.rgeneric.define(model = rgeneric.svc.Hyperbolic.midas.iid,
                                            x = temp_data))
        } else if (hf_info$constraint == "gaussian") {
          return(INLA::inla.rgeneric.define(model = rgeneric.svc.Gaussian.midas.iid,
                                            x = temp_data))
        } else {
          stop("Unknown constraint.")
        }
      }else {
        stop("Unknown `svc_prior`.")
      }
    }
  }

  response_name <- all.vars(formula[[2]])[1]

  if (family == "binomial") {
    data$Ntrials <- Ntrials
  }else if (family == "poisson" & !is.null(E)){
    if (length(E) != nrow(data)) {
      stop("Length of `E` must match the number of rows in `data`.")
    }
    data$E <- E
  }


  rm_max <- max(vapply(hf_input, function(obj) obj$rm.row, numeric(1)))


  if (!loc_var %in% names(data)) {
    stop(sprintf("`loc_var` = '%s' not found in `data`.", loc_var))
  }

  if (!time_var %in% names(data)) {
    stop(sprintf("`time_var` = '%s' not found in `data`.", time_var))
  }

  data <- data[order(data[[loc_var]], data[[time_var]]), , drop = FALSE]

  data_final <- if (rm_max > 0) {
    do.call(
      rbind,
      lapply(split(data, data[[loc_var]]), function(d) {
        d[-seq_len(rm_max), , drop = FALSE]
      })
    )
  } else {
    data
  }

  rownames(data_final) <- NULL


  formula_final <- formula
  final_n <- nrow(data_final)

  formula_env <- environment(formula_final)
  if (is.null(formula_env)) {
    formula_env <- parent.frame()
  }

  if (!is.null(hf_input) && length(hf_input) > 0) {
    stopifnot(is.list(hf_input))

    for (i in seq_along(hf_input)) {
      obj <- hf_input[[i]]

      idx_name <- paste0("hf_idx_", i)
      model_name <- paste0("hf_model_", i)


      X_df <- as.data.frame(obj$X_matrix)

      X_trimmed <- if (rm_max > 0) {
        do.call(
          rbind,
          lapply(split(X_df, X_df$loc), function(d) {
            d[-seq_len(rm_max), , drop = FALSE]
          })
        )
      } else {
        X_df
      }

      rownames(X_trimmed) <- NULL

      temp_data <- data.frame(
        X_trimmed,
        data_final[, response_name, drop = FALSE]
      )

      names(temp_data)[names(temp_data) == "loc"] <- "region"
      rgen <- build_rgen(temp_data = temp_data, hf_info = obj)

      n_region <- max(temp_data$region)
      if(!obj$svc){
        data_final[[idx_name]] <- seq_len(final_n)
      }else{
        if(obj$svc_prior == "iid"){
          data_final[[idx_name]] <- n_region + seq_len(final_n)
        }else if(obj$svc_prior == "icar"){
          data_final[[idx_name]] <- (1 + n_region) + seq_len(final_n)
        }
      }

      assign(model_name, rgen, envir = formula_env)

      if(!obj$svc){
        term_txt <- sprintf(
          "f(%s, model = %s, n = %d)",
          idx_name, model_name, final_n
        )
      }else{
        if(obj$svc_prior == "iid"){
          term_txt <- sprintf(
            "f(%s, model = %s, n = %d)",
            idx_name, model_name, n_region + final_n
          )
        }else if(obj$svc_prior == "icar"){
          A <- Matrix::Matrix(
            0,
            nrow = 1,
            ncol = 1 + n_region + final_n,
            sparse = TRUE
          )
          A[1, 2:(n_region + 1)] <- 1
          extraconstr <- list(
            A = A,
            e = 0
          )
          constr_name <- paste0("hf_extraconstr_", i)
          assign(constr_name, extraconstr, envir = formula_env)
          term_txt <- sprintf(
            "f(%s, model = %s, n = %d, extraconstr = %s)",
            idx_name, model_name, 1 + n_region + final_n, constr_name
          )
        }
      }

      formula_final <- stats::update(formula_final, paste(". ~ . +", term_txt))
    }
  }

  environment(formula_final) <- formula_env

  if (is.null(inla_options$control.compute)) {
    inla_options$control.compute <- list(config = TRUE)
  } else if (is.null(inla_options$control.compute$config)) {
    inla_options$control.compute$config <- TRUE
  }
  args_inla <- list(
    formula = formula_final,
    data = data_final,
    family = family
  )

  if (family == "binomial") {
    args_inla$Ntrials <- data_final$Ntrials
  }

  if (family == "poisson" && !is.null(E)) {
    args_inla$E <- data_final$E
  }

  args_inla <- c(args_inla, inla_options)

  res <- do.call(INLA::inla, args_inla)


  return(list(formula_final = formula_final,
              data_final = data_final,
              rm_max = rm_max,
              res = res,
              hf_input = hf_input))

}





#' Compute beta summaries for the spatial model
#' @param model output from fit_Minla_spatial
#' @param n_loc number of locations
#' @export
compute_beta_spatial <- function(model,
                                 n_loc){

  length_hf <- length(model$hf_input)

  hf_summary_output <- vector(mode = "list", length = length_hf)

  for(i in seq_len(length_hf)){

    temp <- model$hf_input[[i]]

    idx_name <- paste0("hf_idx_", i)

    if(temp$svc){

      if(temp$svc_prior == "icar"){

        res <- vector(mode = "list", length = 4)
        names(res) <- c("summary.global.beta",
                        "marginal.global.beta",
                        "summary.icar.beta",
                        "marginal.icar.beta")
        res[["marginal.icar.beta"]] <- vector(mode = "list", length = n_loc)


        marg <- INLA::inla.rmarginal(1000,marginal = model$res$marginals.random[[idx_name]][[1]])
        summary.global.beta <- data.frame(
          Mean = mean(marg),
          SD = stats::sd(marg),
          `2.5%` = stats::quantile(marg, probs = 0.025),
          `50%` = stats::quantile(marg, probs = 0.5),
          `97.5%` = stats::quantile(marg, probs = 0.975),
          row.names = "beta",
          check.names = FALSE
        )

        res[["summary.global.beta"]] <- summary.global.beta
        res[["marginal.global.beta"]] <- model$res$marginals.random[[idx_name]][[1]]


        for(loc_id in seq_len(n_loc)){

          res[["marginal.icar.beta"]][[loc_id]] <- model$res$marginals.random[[idx_name]][[1+loc_id]]

        }


        temp <- lapply(seq_len(n_loc), function(x){
          marg <- INLA::inla.rmarginal(1000,marginal = res[["marginal.icar.beta"]][[x]])
          data.frame(
            Mean = mean(marg),
            SD = stats::sd(marg),
            `2.5%` = stats::quantile(marg, probs = 0.025),
            `50%` = stats::quantile(marg, probs = 0.5),
            `97.5%` = stats::quantile(marg, probs = 0.975),
            row.names = paste0("b",x),
            check.names = FALSE
          )
        })

        res[["summary.icar.beta"]]<- do.call(rbind,temp)


        hf_summary_output[[i]] <- res


      }else if(temp$svc_prior == "iid"){


        res <- vector(mode = "list", length = 2)
        names(res) <- c("summary.beta",
                        "marginal.beta")
        res[["marginal.beta"]] <- vector(mode = "list", length = n_loc)

        for(loc_id in seq_len(n_loc)){

          res[["marginal.beta"]][[loc_id]] <- model$res$marginals.random[[idx_name]][[loc_id]]

        }

        temp <- lapply(seq_len(n_loc), function(x){
          marg <- INLA::inla.rmarginal(1000,marginal = res[["marginal.beta"]][[x]])
          data.frame(
            Mean = mean(marg),
            SD = stats::sd(marg),
            `2.5%` = stats::quantile(marg, probs = 0.025),
            `50%` = stats::quantile(marg, probs = 0.5),
            `97.5%` = stats::quantile(marg, probs = 0.975),
            row.names = paste0("beta",x),
            check.names = FALSE
          )
        })

        res[["summary.beta"]] <- do.call(rbind,temp)


        hf_summary_output[[i]] <- res

      }

    }else{

      if(temp$constraint == "gaussian"){

        res <- vector(mode = "list", length = 2)
        names(res) <- c("summary.beta",
                        "marginal.beta")

        res[["marginal.beta"]] <- model$res$marginals.hyperpar[[paste0("Theta3 for ", idx_name)]]

        marg <- INLA::inla.rmarginal(1000,marginal = res[["marginal.beta"]])
        summary.beta <- data.frame(
          Mean = mean(marg),
          SD = stats::sd(marg),
          `2.5%` = stats::quantile(marg, probs = 0.025),
          `50%` = stats::quantile(marg, probs = 0.5),
          `97.5%` = stats::quantile(marg, probs = 0.975),
          row.names = "beta",
          check.names = FALSE
        )

        res[["summary.beta"]] <- summary.beta

      }else if(temp$constraint == "hyperbolic"){

        res <- vector(mode = "list", length = 2)
        names(res) <- c("summary.beta",
                        "marginal.beta")

        res[["marginal.beta"]] <- model$res$marginals.hyperpar[[paste0("Theta2 for ", idx_name)]]

        marg <- INLA::inla.rmarginal(1000,marginal = res[["marginal.beta"]])
        summary.beta <- data.frame(
          Mean = mean(marg),
          SD = stats::sd(marg),
          `2.5%` = stats::quantile(marg, probs = 0.025),
          `50%` = stats::quantile(marg, probs = 0.5),
          `97.5%` = stats::quantile(marg, probs = 0.975),
          row.names = "beta",
          check.names = FALSE
        )

        res[["summary.beta"]] <- summary.beta

      }else if(temp$consrtaint == "beta1"){

        res <- vector(mode = "list", length = 2)
        names(res) <- c("summary.beta",
                        "marginal.beta")

        res[["marginal.beta"]] <- model$res$marginals.hyperpar[[paste0("Theta2 for ", idx_name)]]

        marg <- INLA::inla.rmarginal(1000,marginal = res[["marginal.beta"]])
        summary.beta <- data.frame(
          Mean = mean(marg),
          SD = stats::sd(marg),
          `2.5%` = stats::quantile(marg, probs = 0.025),
          `50%` = stats::quantile(marg, probs = 0.5),
          `97.5%` = stats::quantile(marg, probs = 0.975),
          row.names = "beta",
          check.names = FALSE
        )

        res[["summary.beta"]] <- summary.beta

      }

      hf_summary_output[[i]] <- res

    }

  }

  names(hf_summary_output) <- paste0("hf_index_", seq_len(length_hf))

  return(hf_summary_output)

}

