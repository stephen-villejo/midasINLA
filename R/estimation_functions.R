





#' Create lag matrix
#' @param tsdata High-frequency covariate data
#' @param lags Lags to be considered
#' @param frequency Number of covariates values associated with each response
#' @return Matrix of lagged values
#' @export
create_lag_Xmatrix <- function(tsdata,
                               lags,
                               frequency){
  X_matrix <- midasr::mls(tsdata, lags, frequency)
  compile_names <- c()
  for(i in 0:(ncol(X_matrix)-1)){
    compile_names <- c(compile_names,paste0("lag",i))
  }
  colnames(X_matrix) <- compile_names

  return(X_matrix)
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


#' Prepare MIDAS objects for INLA estimation
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
#' @param data Data object to fit MIDAS
#' @param family  Likelihood family for response data
#' @param Ntrials Number of trials for a binomial family response
#' @param nsamples Number of posterior samples
#' @return A list containing the predictions and the samples
#' @export
predict_midas <- function(model,
                          data = NULL,
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



#' Compute scores for model assessment
#' @param y Response vector
#' @param Midas_object Midas object
#' @param pred_res Predictions from the MIDAS model
#' @param family  Likelihood family for response data
#' @return Computed scores
#' @export
compute_forecast_scores <- function(y,
                                    Midas_object,
                                    pred_res,
                                    family = "binomial"){


  if(family == "gaussian"){

    y_forecast <- y[-c(1:(length(y)-length(Midas_object$data$y)))][-c(1:max(which(!is.na(Midas_object$data$y))))]

    SE = (y_forecast - pred_res$computed_y$mean[-c(1:max(which(!is.na(Midas_object$data$y))))])^2

    compile_logscore <- matrix(NA,nrow=length(y_forecast),ncol=ncol(pred_res$samples$latent_predictor))
    for(i in 1:ncol(pred_res$samples$latent_predictor)){
      score <- stats::dnorm(y_forecast,
                            mean = pred_res$samples$latent_predictor[-c(1:max(which(!is.na(Midas_object$data$y)))),i],
                            sd = sqrt(pred_res$samples$sigma2[i]))
      compile_logscore[,i] <- log(score)
    }

    post_E <- rowMeans(pred_res$samples$latent_predictor[-c(1:max(which(!is.na(Midas_object$data$y)))),])
    post_Var <- mean(pred_res$samples$sigma2) + apply(pred_res$samples$latent_predictor, 1, stats::var)[-c(1:max(which(!is.na(Midas_object$data$y))))]

    scores <- data.frame(
      SE = SE,
      DS = (y_forecast - post_E)^2 / post_Var + log(post_Var),
      NLS = -rowMeans(compile_logscore)
    )

  }else if(family == "poisson"){

    y_forecast <- y[-c(1:(length(y)-length(Midas_object$data$y)))][-c(1:max(which(!is.na(Midas_object$data$y))))]

    SE = (y_forecast - pred_res$computed_y$mean[-c(1:max(which(!is.na(Midas_object$data$y))))])^2

    compile_logscore <- matrix(NA,nrow=length(y_forecast),ncol=ncol(pred_res$samples$latent_predictor))
    for(i in 1:ncol(pred_res$samples$latent_predictor)){
      score <- stats::dpois(y_forecast,
                            lambda = exp(pred_res$samples$latent_predictor[-c(1:max(which(!is.na(Midas_object$data$y)))),i]))
      compile_logscore[,i] <- log(score)
    }

    post_E <- rowMeans(exp(pred_res$samples$latent_predictor[-c(1:max(which(!is.na(Midas_object$data$y)))),]))
    post_Var <- rowMeans(exp(pred_res$samples$latent_predictor[-c(1:max(which(!is.na(Midas_object$data$y)))),])) +
      apply(exp(pred_res$samples$latent_predictor), 1, stats::var)[-c(1:max(which(!is.na(Midas_object$data$y))))]

    scores <- data.frame(
      SE = SE,
      DS = (y_forecast - post_E)^2 / post_Var + log(post_Var),
      NLS = -rowMeans(compile_logscore)
    )

  }else if(family == "binomial"){

    y_forecast <- y[-c(1:(length(y)-length(Midas_object$data$y)))][-c(1:max(which(!is.na(Midas_object$data$y))))]

    SE = (y_forecast - pred_res$computed_y$mean[-c(1:max(which(!is.na(Midas_object$data$y))))])^2

    compile_logscore <- matrix(NA,nrow=length(y_forecast),ncol=ncol(pred_res$samples$latent_predictor))
    for(i in 1:ncol(pred_res$samples$latent_predictor)){
      score <- stats::dbinom(y_forecast,
                             size = Midas_object$Ntrials[-c(1:max(which(!is.na(Midas_object$data$y))))],
                             prob = exp(pred_res$samples$latent_predictor[-c(1:max(which(!is.na(Midas_object$data$y)))),i])/(1+exp(pred_res$samples$latent_predictor[-c(1:max(which(!is.na(Midas_object$data$y)))),i])))
      compile_logscore[,i] <- log(score)
    }

    p_temp <- exp(pred_res$samples$latent_predictor[-c(1:max(which(!is.na(Midas_object$data$y)))),])/(1+exp(pred_res$samples$latent_predictor[-c(1:max(which(!is.na(Midas_object$data$y)))),]))
    post_E <- rowMeans(p_temp) * Midas_object$Ntrials[-c(1:max(which(!is.na(Midas_object$data$y))))]
    post_Var <- (Midas_object$Ntrials[-c(1:max(which(!is.na(Midas_object$data$y))))] * rowMeans(p_temp * (1 - p_temp))) +
      ((Midas_object$Ntrials[-c(1:max(which(!is.na(Midas_object$data$y))))]^2) * apply(p_temp, 1, stats::var))

    scores <- data.frame(
      SE = SE,
      DS = (y_forecast - post_E)^2 / post_Var + log(post_Var),
      NLS = -rowMeans(compile_logscore)
    )

  }

  return(list(scores = scores))

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





