


# Create lag matrix
#
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




#' Generate posterior predictive samples from a MIDAS model
#'
#' Generates posterior predictive samples from a fitted MIDAS model
#' using posterior samples of the latent linear predictor. The function
#' supports Gaussian, Poisson, and binomial response distributions.
#'
#' For a Gaussian response, posterior samples of the observation
#' variance are obtained from the posterior samples of the Gaussian
#' precision parameter and used to generate predictive observations.
#' For Poisson and binomial responses, observations are generated using
#' the appropriate inverse-link function.
#'
#' @param model A fitted MIDAS model returned by
#'   [fit_Minla_spatial()].
#' @param family Character string specifying the likelihood family.
#'   Supported values are `"gaussian"`, `"poisson"`, and `"binomial"`.
#'   Defaults to `"gaussian"`.
#' @param Ntrials Optional vector specifying the number of trials for
#'   each observation when `family = "binomial"`. If omitted, the
#'   function attempts to use the `Ntrials` column in
#'   `model$data_final`.
#' @param nsamples Positive integer specifying the number of posterior
#'   samples used to generate the predictive distribution. Defaults to
#'   `1000`.
#'
#' @return A list with two components:
#'   \describe{
#'     \item{computed_y}{
#'       A list containing posterior predictive summaries:
#'       \describe{
#'         \item{mean}{Posterior predictive mean for each observation.}
#'         \item{sd}{Posterior predictive standard deviation for each
#'           observation.}
#'         \item{q2.5}{2.5% posterior predictive quantile for each
#'           observation.}
#'         \item{q97.5}{97.5% posterior predictive quantile for each
#'           observation.}
#'       }
#'     }
#'     \item{samples}{
#'       A list containing posterior samples of the latent predictor and,
#'       for Gaussian responses, the posterior samples of the observation
#'       variance.
#'     }
#'   }
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("INLA", quietly = TRUE)) {
#'   data(data_spatialpoisson_example)
#'
#'   Midas_x1 <- prepare_Minla_spatial(
#'     x = data_spatialpoisson_example$data_x1$x1,
#'     loc_x = data_spatialpoisson_example$data_x1$loc,
#'     constraint = "hyperbolic",
#'     K = 0:29,
#'     m = 30,
#'     svc = FALSE
#'   )
#'
#'   fit <- fit_Minla_spatial(
#'     formula = y ~ 1,
#'     data = data_spatialpoisson_example$data_y,
#'     loc_var = "loc",
#'     time_var = "Time",
#'     family = "poisson",
#'     hf_input = list(Midas_x1),
#'     inla_options = list(verbose = FALSE)
#'   )
#'
#'   predictions <- predict_midas(
#'     model = fit,
#'     family = "poisson",
#'     nsamples = 1000
#'   )
#'
#'   # Inspect posterior predictive summaries
#'   head(predictions$computed_y$mean)
#'   head(predictions$computed_y$q2.5)
#'   head(predictions$computed_y$q97.5)
#' }
#' }
#' @export
predict_midas <- function(model,
                          family = "gaussian",
                          Ntrials = NULL,
                          nsamples = 1000) {

  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop(
      "Package 'INLA' is required for `predict_midas()`. ",
      "Please install INLA from the R-INLA repository."
    )
  }

  data <- model$data_final

  if (is.null(data)) {
    stop("`model$data_final` is missing.")
  }

  fit <- if (!is.null(model$res)) model$res else model

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



#' Compute posterior estimates of MIDAS lag weights
#'
#' Computes posterior estimates of the normalized MIDAS lag weights from
#' a fitted MIDAS model returned by [fit_Minla_spatial()]. The function
#' draws samples from the posterior marginal distributions of the
#' MIDAS hyperparameters and uses these samples to obtain the
#' corresponding lag-weight functions.
#'
#' The supported MIDAS lag constraints are `"hyperbolic"`, `"gaussian"`,
#' `"beta1"`, `"beta2"`, and `"almon2"`. For each high-frequency
#' covariate, the resulting weights are normalized to sum to one across
#' all included lags.
#'
#' @param model A fitted MIDAS model returned by
#'   [fit_Minla_spatial()].
#' @param n.samples Positive integer specifying the number of posterior
#'   samples drawn from the MIDAS hyperparameter marginal distributions
#'   to estimate the lag-weight distribution. Defaults to `200`.
#'
#' @return A list containing one data frame for each high-frequency
#'   covariate in `model$hf_input`. The elements are named `"hf_1"`,
#'   `"hf_2"`, and so on. Each data frame contains:
#'   \describe{
#'     \item{lag}{The MIDAS lag index.}
#'     \item{mean}{The posterior mean of the normalized lag weight.}
#'     \item{q2.5}{The 2.5% posterior quantile of the lag weight.}
#'     \item{q97.5}{The 97.5% posterior quantile of the lag weight.}
#'   }
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("INLA", quietly = TRUE)) {
#'   data(data_spatialpoisson_example)
#'
#'   g <- INLA::inla.read.graph(
#'     filename = system.file("map.adj", package = "midasINLA")
#'   )
#'
#'   Midas_x1 <- prepare_Minla_spatial(
#'     x = data_spatialpoisson_example$data_x1$x1,
#'     loc_x = data_spatialpoisson_example$data_x1$loc,
#'     constraint = "hyperbolic",
#'     K = 0:29,
#'     m = 30,
#'     svc = TRUE,
#'     svc_prior = "icar",
#'     g = g
#'   )
#'
#'   fit <- fit_Minla_spatial(
#'     formula = y ~ 1,
#'     data = data_spatialpoisson_example$data_y,
#'     loc_var = "loc",
#'     time_var = "Time",
#'     family = "poisson",
#'     hf_input = list(Midas_x1),
#'     inla_options = list(
#'       verbose = FALSE,
#'       control.predictor = list(
#'         compute = TRUE,
#'         link = 1
#'       )
#'     )
#'   )
#'
#'   weights <- compute_weights(
#'     model = fit,
#'     n.samples = 200
#'   )
#'
#'   head(weights$hf_1)
#' }
#' }
#' @export
compute_weights <- function(model, n.samples = 200) {

  if (!requireNamespace("INLA", quietly = TRUE)){
    stop(
      "Package 'INLA' is required for `compute_weights()`. ",
      "Please install INLA from the R-INLA repository."
    )
  }

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

    } else if (constraint == "beta2"){

      theta1_samples <- INLA::inla.rmarginal(
        n.samples,
        fit$marginals.hyperpar[[theta1_name]]
      )

      theta2_samples <- INLA::inla.rmarginal(
        n.samples,
        fit$marginals.hyperpar[[theta2_name]]
      )

      gamma1_val_samples <- exp(theta1_samples) + 1
      gamma2_val_samples <- exp(theta2_samples) + 1
      xi <- 1e-4
      psi_mat <- sapply(0:lag_k, function(lag) {
        ((xi + (1 - 2 * xi) * (lag / lag_k))^(gamma1_val_samples-1))*((1-(xi + (1 - 2 * xi) * (lag / lag_k)))^(gamma2_val_samples-1))
      })

      w_mat <- psi_mat / rowSums(psi_mat)


      } else if(constraint == "almon2"){

        theta1_samples <- INLA::inla.rmarginal(
          n.samples,
          fit$marginals.hyperpar[[theta1_name]]
        )

        theta2_samples <- INLA::inla.rmarginal(
          n.samples,
          fit$marginals.hyperpar[[theta2_name]]
        )

        gamma1_val_samples <- 0.01 * tanh(theta1_samples)
        gamma2_val_samples <- 0.01 * tanh(theta2_samples)
        psi_mat <- sapply(0:lag_k, function(lag) {
          exp(gamma1_val_samples*(lag^1) + gamma2_val_samples*(lag^2))
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





#' Prepare a spatial MIDAS object for INLA estimation
#'
#' Constructs the MIDAS design matrix and associated model specifications
#' for use with [fit_Minla_spatial()]. The function creates lagged
#' high-frequency covariate values for each location and optionally
#' specifies a spatially varying coefficient (SVC) component.
#'
#' @param x Numeric vector of high-frequency covariate observations.
#' @param loc_x Vector identifying the location associated with each
#'   observation in `x`. The length of `loc_x` must equal the length of `x`.
#' @param constraint Character string specifying the constraint used for
#'   the MIDAS lag-response association. Supported constraints include
#'   `"hyperbolic"`, `"gaussian"`, `"beta1"`, `"beta2"`, and `"almon2"`.
#' @param K Numeric vector specifying the lags to be included in the
#'   MIDAS representation.
#' @param m Numeric vector specifying the number of high-frequency
#'   covariate observations associated with each response observation.
#'   A single value can be supplied when the number of observations is
#'   constant over time, or a vector can be supplied when this number
#'   varies across response times.
#' @param svc Logical; if `TRUE`, specifies a spatially varying
#'   coefficient model. Defaults to `FALSE`.
#' @param svc_prior Character string specifying the prior for the spatially
#'   varying coefficient component. Must be either `"icar"` or `"iid"`.
#'   Defaults to `"iid"`.
#' @param g An INLA graph object used for the `"icar"` prior. Required
#'   when `svc_prior = "icar"`.
#'
#' @return A list containing the MIDAS design matrix and model
#'   specifications. The returned object includes:
#'   \describe{
#'     \item{X_matrix}{The MIDAS design matrix, including a location
#'       index.}
#'     \item{constraint}{The MIDAS constraint used for the lag-response
#'       association.}
#'     \item{lag_k}{The maximum lag specified in `K`.}
#'     \item{K}{The vector of lags used to construct the design matrix.}
#'     \item{m}{The number of high-frequency observations associated with
#'       each response observation.}
#'     \item{rm.row}{The number of initial rows removed from each location
#'       because of incomplete lagged observations.}
#'     \item{svc}{Whether a spatially varying coefficient component is
#'       specified.}
#'     \item{svc_prior}{The prior specified for the spatially varying
#'       coefficient component.}
#'     \item{g}{The INLA graph object, included when `svc_prior = "icar"`.}
#'   }
#'
#' @examples
#' if (requireNamespace("INLA", quietly = TRUE)) {
#'   data(data_spatialpoisson_example)
#'
#'   # Prepare a MIDAS object using a hyperbolic lag constraint
#'   # and a spatially varying coefficient with an ICAR prior.
#'   g <- INLA::inla.read.graph(
#'     filename = system.file("map.adj", package = "midasINLA")
#'   )
#'
#'   Midas_x1 <- prepare_Minla_spatial(
#'     x = data_spatialpoisson_example$data_x1$x1,
#'     loc_x = data_spatialpoisson_example$data_x1$loc,
#'     constraint = "hyperbolic",
#'     K = 0:29,
#'     m = 30,
#'     svc = TRUE,
#'     svc_prior = "icar",
#'     g = g
#'   )
#'
#'   # Inspect the resulting MIDAS design matrix
#'   head(Midas_x1$X_matrix)
#' }
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



#' Fit a spatial MIDAS model using INLA
#'
#' Fits a Mixed Data Sampling (MIDAS) regression model with optional
#' spatially varying coefficients using Integrated Nested Laplace
#' Approximation (INLA). High-frequency covariates are supplied as
#' MIDAS objects created by [prepare_Minla_spatial()].
#'
#' The function constructs the INLA model by incorporating the
#' MIDAS components specified in `hf_input`. Multiple high-frequency
#' covariates can be included by supplying multiple MIDAS objects in
#' `hf_input`.
#'
#' @param formula A model formula specifying the response and other
#'   covariates. The response variable must be a column in `data`.
#' @param data A data frame containing the response and any additional
#'   model covariates. It must contain the variables specified by
#'   `loc_var` and `time_var`.
#' @param loc_var Character string specifying the name of the location
#'   variable in `data`.
#' @param time_var Character string specifying the name of the time
#'   variable in `data`.
#' @param family Character string specifying the likelihood family for
#'   the response. For example, `"poisson"` or `"binomial"`.
#' @param hf_input A list of MIDAS objects returned by
#'   [prepare_Minla_spatial()]. Each object specifies a high-frequency
#'   covariate and its MIDAS lag structure. Multiple MIDAS objects can
#'   be supplied.
#' @param Ntrials Optional vector specifying the number of trials for a
#'   binomial response. Used only when `family = "binomial"`.
#' @param E Optional vector of expected counts or exposure values for a
#'   Poisson model. Its length must match the number of rows in `data`.
#'   Used only when `family = "poisson"`.
#' @param inla_options A named list of additional arguments passed to
#'   [INLA::inla()]. By default, `control.compute$config` is set to
#'   `TRUE` if it is not already specified.
#'
#' @return A list containing:
#'   \describe{
#'     \item{formula_final}{The final INLA model formula, including the
#'       MIDAS components.}
#'     \item{data_final}{The response data used for model fitting after
#'       ordering by location and time and removing rows with incomplete
#'       lagged covariate information.}
#'     \item{rm_max}{The maximum number of initial observations removed
#'       across locations because of incomplete MIDAS lagged covariates.}
#'     \item{res}{The fitted INLA model returned by [INLA::inla()].}
#'     \item{hf_input}{The list of MIDAS objects supplied through
#'       `hf_input`.}
#'   }
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("INLA", quietly = TRUE)) {
#'   data(data_spatialpoisson_example)
#'
#'   # Read the spatial adjacency graph
#'   g <- INLA::inla.read.graph(
#'     filename = system.file("map.adj", package = "midasINLA")
#'   )
#'
#'   # Prepare a spatial MIDAS predictor
#'   Midas_x1 <- prepare_Minla_spatial(
#'     x = data_spatialpoisson_example$data_x1$x1,
#'     loc_x = data_spatialpoisson_example$data_x1$loc,
#'     constraint = "hyperbolic",
#'     K = 0:29,
#'     m = 30,
#'     svc = TRUE,
#'     svc_prior = "icar",
#'     g = g
#'   )
#'
#'   # Fit the spatial Poisson MIDAS model
#'   fit <- fit_Minla_spatial(
#'     formula = y ~ 1,
#'     data = data_spatialpoisson_example$data_y,
#'     loc_var = "loc",
#'     time_var = "Time",
#'     family = "poisson",
#'     hf_input = list(Midas_x1),
#'     inla_options = list(
#'       verbose = FALSE,
#'       control.predictor = list(
#'         compute = TRUE,
#'         link = 1
#'       )
#'     )
#'   )
#'
#'   # Inspect the fitted INLA model
#'   fit$res
#' }
#' }
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

  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop(
      "Package 'INLA' is required for this function. ",
      "Please install INLA from the R-INLA repository."
    )
  }

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





#' Compute posterior summaries of MIDAS coefficients
#'
#' Computes posterior summaries of the MIDAS regression coefficients from
#' a fitted spatial MIDAS model returned by [fit_Minla_spatial()]. The
#' output depends on whether the MIDAS coefficient is spatially varying
#' and on the specified spatial prior.
#'
#' For spatially varying coefficients with an ICAR prior, the function
#' returns summaries and marginal distributions for the global coefficient,
#' the location-specific spatial deviations, and the resulting
#' location-specific total coefficients. Posterior samples of the total
#' coefficients are also returned.
#'
#' For spatially varying coefficients with an IID prior, the function
#' returns posterior summaries and marginal distributions for the
#' location-specific coefficients.
#'
#' For non-spatially varying coefficients, the function returns the
#' posterior marginal distribution and summary of the MIDAS coefficient
#' associated with each high-frequency covariate.
#'
#' @param model A fitted spatial MIDAS model returned by
#'   [fit_Minla_spatial()].
#' @param n_loc Integer specifying the number of spatial locations for
#'   which coefficient summaries should be computed.
#'
#' @return A list containing one element for each high-frequency MIDAS
#'   covariate in `model$hf_input`. The elements are named
#'   `"hf_index_1"`, `"hf_index_2"`, and so on. The contents of each
#'   element depend on the spatial structure of the corresponding MIDAS
#'   covariate:
#'   \describe{
#'     \item{Spatially varying coefficient with ICAR prior}{
#'       A list containing `summary.global.beta`,
#'       `marginal.global.beta`, `summary.icar.beta`,
#'       `marginal.icar.beta`, `summary.total.beta`, and
#'       `sample.total.beta`.
#'     }
#'     \item{Spatially varying coefficient with IID prior}{
#'       A list containing `summary.beta` and `marginal.beta`.
#'     }
#'     \item{Non-spatially varying coefficient}{
#'       A list containing `summary.beta` and `marginal.beta`.
#'     }
#'   }
#'
#'   The summary data frames contain the posterior mean, standard
#'   deviation, and 2.5%, 50%, and 97.5% posterior quantiles.
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("INLA", quietly = TRUE)) {
#'   INLA::inla.setOption(num.threads = 1)
#'
#'   data(data_spatialpoisson_example)
#'
#'   # Read the spatial adjacency graph
#'   g <- INLA::inla.read.graph(
#'     filename = system.file("map.adj", package = "midasINLA")
#'   )
#'
#'   # Prepare a spatial MIDAS predictor
#'   Midas_x1 <- prepare_Minla_spatial(
#'     x = data_spatialpoisson_example$data_x1$x1,
#'     loc_x = data_spatialpoisson_example$data_x1$loc,
#'     constraint = "hyperbolic",
#'     K = 0:29,
#'     m = 30,
#'     svc = TRUE,
#'     svc_prior = "icar",
#'     g = g
#'   )
#'
#'   # Fit the spatial Poisson MIDAS model
#'   fit <- fit_Minla_spatial(
#'     formula = y ~ 1,
#'     data = data_spatialpoisson_example$data_y,
#'     loc_var = "loc",
#'     time_var = "Time",
#'     family = "poisson",
#'     hf_input = list(Midas_x1),
#'     inla_options = list(
#'       verbose = FALSE,
#'       control.predictor = list(
#'         compute = TRUE,
#'         link = 1
#'       )
#'     )
#'   )
#'
#'   # Compute posterior summaries of the MIDAS coefficients
#'   beta_summary <- compute_beta_spatial(
#'     model = fit,
#'     n_loc = 16
#'   )
#'
#'   # Inspect summaries of the location-specific total coefficients
#'   beta_summary$hf_index_1$summary.total.beta
#' }
#' }
#' @export
compute_beta_spatial <- function(model,
                                 n_loc){

  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop(
      "Package 'INLA' is required for this function. ",
      "Please install INLA from the R-INLA repository."
    )
  }

  length_hf <- length(model$hf_input)

  hf_summary_output <- vector(mode = "list", length = length_hf)

  for(i in seq_len(length_hf)){

    temp <- model$hf_input[[i]]

    idx_name <- paste0("hf_idx_", i)

    if(temp$svc){

      if(temp$svc_prior == "icar"){

        res <- vector(mode = "list", length = 6)
        names(res) <- c("summary.global.beta",
                        "marginal.global.beta",
                        "summary.icar.beta",
                        "marginal.icar.beta",
                        "summary.total.beta",
                        "sample.total.beta")

        res[["marginal.icar.beta"]] <- vector(mode = "list", length = n_loc)
        res[["sample.total.beta"]]  <- vector(mode = "list", length = n_loc)

        # Global beta summary
        global_marg <- model$res$marginals.random[[idx_name]][[1]]
        global_draws <- INLA::inla.rmarginal(1000, marginal = global_marg)

        res[["marginal.global.beta"]] <- global_marg
        res[["summary.global.beta"]] <- data.frame(
          Mean = mean(global_draws),
          SD = stats::sd(global_draws),
          `2.5%` = stats::quantile(global_draws, probs = 0.025),
          `50%` = stats::quantile(global_draws, probs = 0.5),
          `97.5%` = stats::quantile(global_draws, probs = 0.975),
          row.names = "beta",
          check.names = FALSE
        )

        # Local iCAR marginals and summaries
        local_summaries <- lapply(seq_len(n_loc), function(loc_id){

          local_marg <- model$res$marginals.random[[idx_name]][[1 + loc_id]]
          res[["marginal.icar.beta"]][[loc_id]] <<- local_marg

          local_draws <- INLA::inla.rmarginal(1000, marginal = local_marg)

          data.frame(
            Mean = mean(local_draws),
            SD = stats::sd(local_draws),
            `2.5%` = stats::quantile(local_draws, probs = 0.025),
            `50%` = stats::quantile(local_draws, probs = 0.5),
            `97.5%` = stats::quantile(local_draws, probs = 0.975),
            row.names = paste0("b", loc_id),
            check.names = FALSE
          )
        })

        res[["summary.icar.beta"]] <- do.call(rbind, local_summaries)

        # Joint posterior samples for beta_i = beta* + b_i
        post_samp <- INLA::inla.posterior.sample(
          n = 1000,
          result = model$res
        )

        latent_names <- rownames(post_samp[[1]]$latent)

        # Identify latent entries corresponding to this MIDAS covariate
        relevant_idx <- grep(paste0("^", idx_name), latent_names)

        # Assumption: first is global beta, remaining are location-specific deviations
        global_idx <- relevant_idx[1]
        local_idx  <- relevant_idx[-1]

        total_beta_draws <- lapply(seq_len(n_loc), function(loc_id){
          sapply(post_samp, function(s){
            s$latent[global_idx, 1] + s$latent[local_idx[loc_id], 1]
          })
        })

        res[["sample.total.beta"]] <- total_beta_draws

        total_summaries <- lapply(seq_len(n_loc), function(loc_id){
          draws <- total_beta_draws[[loc_id]]
          data.frame(
            Mean = mean(draws),
            SD = stats::sd(draws),
            `2.5%` = stats::quantile(draws, probs = 0.025),
            `50%` = stats::quantile(draws, probs = 0.5),
            `97.5%` = stats::quantile(draws, probs = 0.975),
            row.names = paste0("beta", loc_id),
            check.names = FALSE
          )
        })

        res[["summary.total.beta"]] <- do.call(rbind, total_summaries)

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

      }else if(temp$constraint == "beta1"){

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

      }else if(temp$constraint == "beta2"){

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

      }else if(temp$constraint == "almon2"){

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

      hf_summary_output[[i]] <- res

    }

  }

  names(hf_summary_output) <- paste0("hf_index_", seq_len(length_hf))

  return(hf_summary_output)

}


