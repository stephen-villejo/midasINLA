
#' Rgeneric MIDAS spatially-varying coefficient model with hyperbolic
#' scheme constraint
#'
#' Defines a custom \code{rgeneric} model for use with the \code{INLA} framework,
#' implementing MIDAS-type lag weights using a hyperbolic scheme structure,
#' and where \code{beta[i]} varies for each location
#' following an iCAR model
#'

#' @param cmd Character string indicating the INLA command.
#'   One of \code{"graph"}, \code{"Q"}, \code{"mu"}, \code{"initial"},
#'   \code{"log.norm.const"}, \code{"log.prior"}, or \code{"quit"}.
#' @param theta Numeric vector of hyperparameters controlling the MIDAS weights.
#'
#' @return Depends on \code{cmd}:
#' \itemize{
#'   \item \code{graph}: Sparse precision structure
#'   \item \code{Q}: Precision matrix
#'   \item \code{mu}: Mean vector
#'   \item \code{initial}: Initial values for \code{theta}
#'   \item \code{log.norm.const}: Normalizing constant (numeric(0))
#'   \item \code{log.prior}: Log prior density
#' }
#'
#' @details
#' The MIDAS lag weights are constructed using a normalized Hyperbolic scheme
#' polynomial transformation of lag indices. The parameter \code{theta[1]}
#' control the shape of the lag weighting function.
#'
#' @importFrom Matrix Diagonal
#'
rgeneric.svc.Hyperbolic.midas.icar = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                                                      "log.prior", "quit"),
                                              theta = NULL){

  envir <- parent.env(environment())

  x <- envir$x
  g <- envir$g

  regions <- sort(unique(x$region))
  n_regions <- length(regions)

  region_id <- match(x$region, regions)

  if (any(is.na(region_id))) {
    stop("Some observations have region labels not matched to 'regions'.")
  }
  if (!all(region_id %in% seq_len(n_regions))) {
    stop("region_id values are outside 1:n_regions.")
  }

  build_Q_icar_from_nbs <- function(g) {

    nbs <- g$nbs
    N <- length(nbs)

    ii <- integer(0)
    jj <- integer(0)

    for (i in seq_len(N)) {
      nei <- nbs[[i]]
      if (length(nei) > 0L) {
        ii <- c(ii, rep(i, length(nei)))
        jj <- c(jj, nei)
      }
    }

    W <- Matrix::sparseMatrix(
      i = ii,
      j = jj,
      x = 1,
      dims = c(N, N)
    )

    W <- (W + Matrix::t(W)) > 0
    W <- Matrix::Matrix(W * 1, sparse = TRUE)

    D <- Matrix::Diagonal(x = Matrix::rowSums(W))
    Q_icar <- D - W

    return(Q_icar)
  }

  Q_icar <- build_Q_icar_from_nbs(g)
  stopifnot(n_regions == nrow(Q_icar))
  stopifnot(all(region_id %in% seq_len(n_regions)))

  prec.high <- exp(15)
  prec_beta <- 1 # this is the precision for the global beta parameter

  interpret.theta <- function() {

    lag_cols <- grep("^lag[0-9]+$", names(x), value = TRUE)
    lag_k <- length(lag_cols) - 1L

    gamma_val <- stats::plogis(theta[1L])
    tau_icar  <- exp(theta[2L])

    psi <- numeric(lag_k + 1L)
    for (lag in 0:lag_k) {
      psi[lag + 1L] <- exp(lgamma(lag + gamma_val) -
                             lgamma(lag + 1L) -
                             lgamma(gamma_val))
    }

    w <- psi / sum(psi)

    out_list <- list()

    out_list[["w"]] <- w
    out_list[["gamma"]] <- gamma_val
    out_list[["tau_icar"]] <- tau_icar

    return(out_list)
  }

  graph <- function() {

    n_obs <- nrow(x)
    N <- nrow(Q_icar)

    # global beta node
    Sgg <- Matrix::Matrix(1, nrow = 1, ncol = 1, sparse = TRUE)

    # global-to-region links: must be present because Qgb is nonzero
    Sgb <- Matrix::Matrix(1, nrow = 1, ncol = N, sparse = TRUE)

    # region random effects block
    Sbb <- Q_icar
    Sbb@x[] <- 1
    Sbb <- Sbb + Matrix::Diagonal(n = N, x = 1)
    Sbb@x[] <- 1

    # global beta to obs links
    Sge <- Matrix::Matrix(1, nrow = 1, ncol = n_obs, sparse = TRUE)

    # region beta_i to obs links
    Sbe <- Matrix::sparseMatrix(
      i = region_id,
      j = seq_len(n_obs),
      x = 1,
      dims = c(N, n_obs)
    )

    # observation block
    See <- Matrix::Diagonal(n = n_obs, x = 1)

    G <- rbind(
      cbind(Sgg, Sgb, Sge),
      cbind(Matrix::t(Sgb), Sbb, Sbe),
      cbind(Matrix::t(Sge), Matrix::t(Sbe), See)
    )

    G@x[] <- 1
    return(G)
  }


  Q <- function() {

    par <- interpret.theta()

    n_obs <- nrow(x)
    N <- nrow(Q_icar)

    lag_cols <- grep("^lag[0-9]+$", names(x), value = TRUE)

    Zmat <- as.matrix(x[, lag_cols, drop = FALSE])
    z <- as.vector(Zmat %*% par$w)

    prec_beta <- 1

    # global beta prior
    Qgg <- Matrix::Matrix(prec_beta + prec.high * sum(z^2),
                          nrow = 1, ncol = 1, sparse = TRUE)

    # global-region block
    zb_sum_by_region <- rowsum(z^2, group = region_id, reorder = FALSE)
    zb_sum_by_region <- as.numeric(zb_sum_by_region)

    if (length(zb_sum_by_region) < N) {
      tmp <- numeric(N)
      tmp[sort(unique(region_id))] <- zb_sum_by_region
      zb_sum_by_region <- tmp
    }

    Qgb <- Matrix::Matrix(prec.high * zb_sum_by_region,
                          nrow = 1, ncol = N, sparse = TRUE)

    # global-observation block
    Qge <- Matrix::sparseMatrix(
      i = rep(1L, n_obs),
      j = seq_len(n_obs),
      x = -prec.high * z,
      dims = c(1, n_obs)
    )

    # region-region block
    beta_diag <- rowsum(z^2, group = region_id, reorder = FALSE)
    beta_diag <- as.numeric(beta_diag)

    if (length(beta_diag) < N) {
      tmp <- numeric(N)
      tmp[sort(unique(region_id))] <- beta_diag
      beta_diag <- tmp
    }

    Qbb <- par$tau_icar * Q_icar + Matrix::Diagonal(n = N, x = prec.high * beta_diag)

    # region-observation block
    Qbe <- Matrix::sparseMatrix(
      i = region_id,
      j = seq_len(n_obs),
      x = -prec.high * z,
      dims = c(N, n_obs)
    )

    # observation block
    Qee <- Matrix::Diagonal(n = n_obs, x = prec.high)

    # assemble
    Qfull <- rbind(
      cbind(Qgg, Qgb, Qge),
      cbind(Matrix::t(Qgb), Qbb, Qbe),
      cbind(Matrix::t(Qge), Matrix::t(Qbe), Qee)
    )

    return(Qfull)
  }

  mu <- function() {
    n_obs <- nrow(x)
    N <- nrow(Q_icar)
    return(numeric(1 + N + n_obs))
  }

  log.norm.const <- function() {
    return(numeric(0))
  }

  log.prior <- function() {
    val <- (stats::dnorm(theta[1L], mean=0, sd=1, log=TRUE) +
              stats::dnorm(theta[2L], mean=0, sd=1, log=TRUE))
    return(val)
  }

  initial <- function() {
    return(c(0, 0))
  }

  quit <- function() {
    return(invisible())
  }

  val <- do.call(match.arg(cmd), args = list())

  return(val)

}




#' Rgeneric MIDAS spatially-varying coefficient model with hyperbolic
#' scheme constraint
#'
#' Defines a custom \code{rgeneric} model for use with the \code{INLA} framework,
#' implementing MIDAS-type lag weights using a hyperbolic scheme structure,
#' and where \code{beta} is constant for all locations
#'

#' @param cmd Character string indicating the INLA command.
#'   One of \code{"graph"}, \code{"Q"}, \code{"mu"}, \code{"initial"},
#'   \code{"log.norm.const"}, \code{"log.prior"}, or \code{"quit"}.
#' @param theta Numeric vector of hyperparameters controlling the MIDAS weights.
#'
#' @return Depends on \code{cmd}:
#' \itemize{
#'   \item \code{graph}: Sparse precision structure
#'   \item \code{Q}: Precision matrix
#'   \item \code{mu}: Mean vector
#'   \item \code{initial}: Initial values for \code{theta}
#'   \item \code{log.norm.const}: Normalizing constant (numeric(0))
#'   \item \code{log.prior}: Log prior density
#' }
#'
#' @details
#' The MIDAS lag weights are constructed using a normalized Hyperbolic scheme
#' polynomial transformation of lag indices. The parameter \code{theta[1]}
#' control the shape of the lag weighting function.
#'
#' @importFrom Matrix Diagonal
#'
rgeneric.globalbeta.Hyperbolic.midas = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                                                        "log.prior", "quit"),
                                                theta = NULL){

  envir = parent.env(environment())
  x <- envir$x

  ## artificial high precision to be added to the mean-model
  prec.high = exp(15)

  interpret.theta = function() {

    lag_cols <- grep("^lag[0-9]+$", names(x), value = TRUE)
    lag_k <- length(lag_cols) - 1L

    gamma_val = exp(theta[1L])/(1+exp(theta[1L]))
    beta1 <- theta[2L]

    psi <- numeric(lag_k + 1)
    for(lag in 0:lag_k){
      psi[lag + 1] <- gamma(lag+gamma_val) / (gamma(lag+1)*gamma(gamma_val))
    }

    w <- psi / sum(psi)

    out_list <- list()
    for(lag in 0:lag_k){
      out_list[[paste0("w",lag)]] <- w[lag + 1]
    }

    out_list[["beta1"]] <- beta1

    return(out_list)

  }
  graph = function() {
    G = Matrix::Diagonal(n = length(x$lag0), x=1)
    return(G)
  }
  Q = function() {
    Q = prec.high * graph()
    return(Q)
  }
  mu = function() {
    par = interpret.theta()

    lag_cols <- grep("^lag[0-9]+$", names(x), value = TRUE)
    lag_k <- length(lag_cols) - 1L

    compile_lag_label <- c()
    for(lag in 0:lag_k){
      compile_lag_label <- c(compile_lag_label, paste0("lag",lag))
    }
    compile_w_label <- c()
    for(lag in 0:lag_k){
      compile_w_label <- c(compile_w_label, paste0("w",lag))
    }

    agg <- 0
    for(lag in 0:lag_k){
      agg <- agg + par[[compile_w_label[[lag+1]]]] * x[,which(names(x) == compile_lag_label[[lag+1]])]
    }

    return(par$beta1 * agg)
  }
  log.norm.const = function() {
    return(numeric(0))
  }
  log.prior = function() {
    val = (stats::dnorm(theta[1L], mean=0, sd=1, log=TRUE) +
             stats::dnorm(theta[2L], mean=0, sd=1, log=TRUE))
    return(val)
  }
  initial = function() {
    return(rep(1, 2))
  }
  quit = function() {
    return(invisible())
  }

  val = do.call(match.arg(cmd), args = list())
  return(val)
}





#' Rgeneric MIDAS spatially-varying coefficient model with hyperbolic
#' scheme constraint
#'
#' Defines a custom \code{rgeneric} model for use with the \code{INLA} framework,
#' implementing MIDAS-type lag weights using a Gaussian structure,
#' and where \code{beta[i]} varies for each location
#' following an iid model
#'

#' @param cmd Character string indicating the INLA command.
#'   One of \code{"graph"}, \code{"Q"}, \code{"mu"}, \code{"initial"},
#'   \code{"log.norm.const"}, \code{"log.prior"}, or \code{"quit"}.
#' @param theta Numeric vector of hyperparameters controlling the MIDAS weights.
#'
#' @return Depends on \code{cmd}:
#' \itemize{
#'   \item \code{graph}: Sparse precision structure
#'   \item \code{Q}: Precision matrix
#'   \item \code{mu}: Mean vector
#'   \item \code{initial}: Initial values for \code{theta}
#'   \item \code{log.norm.const}: Normalizing constant (numeric(0))
#'   \item \code{log.prior}: Log prior density
#' }
#'
#' @details
#' The MIDAS lag weights are constructed using a normalized Hyperbolic scheme
#' polynomial transformation of lag indices. The parameter \code{theta[1]}
#' control the shape of the lag weighting function.
#'
#' @importFrom Matrix Diagonal
#'
rgeneric.svc.Hyperbolic.midas.iid = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                                                     "log.prior", "quit"),
                                             theta = NULL){

  envir <- parent.env(environment())

  x <- envir$x

  regions <- sort(unique(x$region))
  n_regions <- length(regions)

  region_id <- match(x$region, regions)

  if (any(is.na(region_id))) {
    stop("Some observations have region labels not matched to 'regions'.")
  }
  if (!all(region_id %in% seq_len(n_regions))) {
    stop("region_id values are outside 1:n_regions.")
  }

  interpret.theta <- function() {

    lag_cols <- grep("^lag[0-9]+$", names(x), value = TRUE)
    lag_k <- length(lag_cols) - 1L

    gamma_val <- stats::plogis(theta[1L])
    tau_beta  <- exp(theta[2L])

    psi <- numeric(lag_k + 1L)
    for (lag in 0:lag_k) {
      psi[lag + 1L] <- exp(lgamma(lag + gamma_val) -
                             lgamma(lag + 1L) -
                             lgamma(gamma_val))
    }

    w <- psi / sum(psi)

    out_list <- list()
    out_list[["w"]] <- w
    out_list[["gamma"]] <- gamma_val
    out_list[["tau_beta"]] <- tau_beta

    return(out_list)
  }

  graph <- function() {

    n_obs <- nrow(x)
    N <- n_regions

    # region-specific beta_i block (iid => diagonal graph)
    Sbb <- Matrix::Diagonal(n = N, x = 1)

    # beta_i to observation links
    Sbe <- Matrix::sparseMatrix(
      i = region_id,
      j = seq_len(n_obs),
      x = 1,
      dims = c(N, n_obs)
    )

    # observation block
    See <- Matrix::Diagonal(n = n_obs, x = 1)

    G <- rbind(
      cbind(Sbb, Sbe),
      cbind(Matrix::t(Sbe), See)
    )

    G@x[] <- 1
    return(G)
  }

  Q <- function() {

    par <- interpret.theta()

    n_obs <- nrow(x)
    N <- n_regions
    prec.high <- exp(15)

    lag_cols <- grep("^lag[0-9]+$", names(x), value = TRUE)
    Zmat <- as.matrix(x[, lag_cols, drop = FALSE])
    z <- as.vector(Zmat %*% par$w)

    # Q_{beta beta} = tau_beta I + tau_eta A' D_z^2 A
    beta_diag <- rowsum(z^2, group = region_id, reorder = FALSE)
    beta_diag <- as.numeric(beta_diag)

    if (length(beta_diag) < N) {
      tmp <- numeric(N)
      tmp[sort(unique(region_id))] <- beta_diag
      beta_diag <- tmp
    }

    Qbb <- Matrix::Diagonal(n = N, x = par$tau_beta + prec.high * beta_diag)

    # Q_{beta eta} = -tau_eta A' D_z
    Qbe <- Matrix::sparseMatrix(
      i = region_id,
      j = seq_len(n_obs),
      x = -prec.high * z,
      dims = c(N, n_obs)
    )

    # Q_{eta eta} = tau_eta I
    Qee <- Matrix::Diagonal(n = n_obs, x = prec.high)

    Qfull <- rbind(
      cbind(Qbb, Qbe),
      cbind(Matrix::t(Qbe), Qee)
    )

    return(Qfull)
  }

  mu <- function() {
    n_obs <- nrow(x)
    N <- n_regions
    return(numeric(N + n_obs))
  }

  log.norm.const <- function() {
    return(numeric(0))
  }

  log.prior <- function() {
    val <- (stats::dnorm(theta[1L], mean = 0, sd = 1, log = TRUE) +
              stats::dnorm(theta[2L], mean = 0, sd = 1, log = TRUE))
    return(val)
  }

  initial <- function() {
    return(c(0, 0))
  }

  quit <- function() {
    return(invisible())
  }

  val <- do.call(match.arg(cmd), args = list())

  return(val)
}




#' Rgeneric MIDAS spatially-varying coefficient model with gaussian constraint
#'
#' Defines a custom \code{rgeneric} model for use with the \code{INLA} framework,
#' implementing MIDAS-type lag weights using a hyperbolic scheme structure,
#' and where \code{beta} is constant for all locations
#'

#' @param cmd Character string indicating the INLA command.
#'   One of \code{"graph"}, \code{"Q"}, \code{"mu"}, \code{"initial"},
#'   \code{"log.norm.const"}, \code{"log.prior"}, or \code{"quit"}.
#' @param theta Numeric vector of hyperparameters controlling the MIDAS weights.
#'
#' @return Depends on \code{cmd}:
#' \itemize{
#'   \item \code{graph}: Sparse precision structure
#'   \item \code{Q}: Precision matrix
#'   \item \code{mu}: Mean vector
#'   \item \code{initial}: Initial values for \code{theta}
#'   \item \code{log.norm.const}: Normalizing constant (numeric(0))
#'   \item \code{log.prior}: Log prior density
#' }
#'
#' @details
#' The MIDAS lag weights are constructed using a normalized Gaussian constraint
#' The parameters \code{theta[1]} and \code{theta[2]} control the shape of the
#' lag weighting function.
#'
#' @importFrom Matrix Diagonal
#'
rgeneric.globalbeta.Gaussian.midas = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                                                        "log.prior", "quit"),
                                                theta = NULL){

  envir = parent.env(environment())
  x <- envir$x

  ## artificial high precision to be added to the mean-model
  prec.high = exp(15)

  interpret.theta = function() {

    lag_cols <- grep("^lag[0-9]+$", names(x), value = TRUE)
    lag_k <- length(lag_cols) - 1L

    mu_val <- lag_k * (1 / (1 + exp(-theta[1L])))
    sigma_val = exp(theta[2L])
    beta1 <- theta[3L]

    psi <- numeric(lag_k + 1)
    for(lag in 0:lag_k){
      psi[lag + 1] <- exp( -(lag-mu_val)^2 / (2 * (sigma_val ^ 2)) )
    }

    w <- psi / sum(psi)

    out_list <- list()
    for(lag in 0:lag_k){
      out_list[[paste0("w",lag)]] <- w[lag + 1]
    }

    out_list[["beta1"]] <- beta1

    return(out_list)

  }
  graph = function() {
    G = Matrix::Diagonal(n = length(x$lag0), x=1)
    return(G)
  }
  Q = function() {
    Q = prec.high * graph()
    return(Q)
  }
  mu = function() {
    par = interpret.theta()

    lag_cols <- grep("^lag[0-9]+$", names(x), value = TRUE)
    lag_k <- length(lag_cols) - 1L

    compile_lag_label <- c()
    for(lag in 0:lag_k){
      compile_lag_label <- c(compile_lag_label, paste0("lag",lag))
    }
    compile_w_label <- c()
    for(lag in 0:lag_k){
      compile_w_label <- c(compile_w_label, paste0("w",lag))
    }

    agg <- 0
    for(lag in 0:lag_k){
      agg <- agg + par[[compile_w_label[[lag+1]]]] * x[,which(names(x) == compile_lag_label[[lag+1]])]
    }

    return(par$beta1 * agg)
  }
  log.norm.const = function() {
    return(numeric(0))
  }
  log.prior = function() {
    val = (stats::dnorm(theta[1L], mean=0, sd=1, log=TRUE) +
             stats::dnorm(theta[2L], mean=0, sd=1, log=TRUE) +
             stats::dnorm(theta[3L], mean=0, sd=1, log=TRUE))
    return(val)
  }
  initial = function() {
    return(rep(1, 3))
  }
  quit = function() {
    return(invisible())
  }

  val = do.call(match.arg(cmd), args = list())
  return(val)
}






#' Rgeneric MIDAS spatially-varying coefficient model with Gaussian constraint
#'
#' Defines a custom \code{rgeneric} model for use with the \code{INLA} framework,
#' implementing MIDAS-type lag weights using a Gaussian structure,
#' and where \code{beta[i]} varies for each location
#' following an iid model
#'

#' @param cmd Character string indicating the INLA command.
#'   One of \code{"graph"}, \code{"Q"}, \code{"mu"}, \code{"initial"},
#'   \code{"log.norm.const"}, \code{"log.prior"}, or \code{"quit"}.
#' @param theta Numeric vector of hyperparameters controlling the MIDAS weights.
#'
#' @return Depends on \code{cmd}:
#' \itemize{
#'   \item \code{graph}: Sparse precision structure
#'   \item \code{Q}: Precision matrix
#'   \item \code{mu}: Mean vector
#'   \item \code{initial}: Initial values for \code{theta}
#'   \item \code{log.norm.const}: Normalizing constant (numeric(0))
#'   \item \code{log.prior}: Log prior density
#' }
#'
#' @details
#' The MIDAS lag weights are constructed using a normalized Gaussian
#' polynomial transformation of lag indices. The parameters \code{theta[1]} and
#' \code{theta[2]} control the shape of the lag weighting function.
#'
#' @importFrom Matrix Diagonal
#'
rgeneric.svc.Gaussian.midas.iid = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                                                     "log.prior", "quit"),
                                             theta = NULL){

  envir <- parent.env(environment())

  x <- envir$x

  regions <- sort(unique(x$region))
  n_regions <- length(regions)

  region_id <- match(x$region, regions)

  if (any(is.na(region_id))) {
    stop("Some observations have region labels not matched to 'regions'.")
  }
  if (!all(region_id %in% seq_len(n_regions))) {
    stop("region_id values are outside 1:n_regions.")
  }

  interpret.theta <- function() {

    lag_cols <- grep("^lag[0-9]+$", names(x), value = TRUE)
    lag_k <- length(lag_cols) - 1L

    mu_val <- lag_k * (1 / (1 + exp(-theta[1L])))
    sigma_val = exp(theta[2L])

    psi <- numeric(lag_k + 1)
    for(lag in 0:lag_k){
      psi[lag + 1] <- exp( -(lag-mu_val)^2 / (2 * (sigma_val ^ 2)) )
    }

    tau_beta  <- exp(theta[3L])

    w <- psi / sum(psi)

    out_list <- list()
    out_list[["w"]] <- w
    out_list[["mu_val"]] <- mu_val
    out_list[["sigma_val"]] <- sigma_val
    out_list[["tau_beta"]] <- tau_beta

    return(out_list)
  }

  graph <- function() {

    n_obs <- nrow(x)
    N <- n_regions

    # region-specific beta_i block (iid => diagonal graph)
    Sbb <- Matrix::Diagonal(n = N, x = 1)

    # beta_i to observation links
    Sbe <- Matrix::sparseMatrix(
      i = region_id,
      j = seq_len(n_obs),
      x = 1,
      dims = c(N, n_obs)
    )

    # observation block
    See <- Matrix::Diagonal(n = n_obs, x = 1)

    G <- rbind(
      cbind(Sbb, Sbe),
      cbind(Matrix::t(Sbe), See)
    )

    G@x[] <- 1
    return(G)
  }

  Q <- function() {

    par <- interpret.theta()

    n_obs <- nrow(x)
    N <- n_regions
    prec.high <- exp(15)

    lag_cols <- grep("^lag[0-9]+$", names(x), value = TRUE)
    Zmat <- as.matrix(x[, lag_cols, drop = FALSE])
    z <- as.vector(Zmat %*% par$w)

    # Q_{beta beta} = tau_beta I + tau_eta A' D_z^2 A
    beta_diag <- rowsum(z^2, group = region_id, reorder = FALSE)
    beta_diag <- as.numeric(beta_diag)

    if (length(beta_diag) < N) {
      tmp <- numeric(N)
      tmp[sort(unique(region_id))] <- beta_diag
      beta_diag <- tmp
    }

    Qbb <- Matrix::Diagonal(n = N, x = par$tau_beta + prec.high * beta_diag)

    # Q_{beta eta} = -tau_eta A' D_z
    Qbe <- Matrix::sparseMatrix(
      i = region_id,
      j = seq_len(n_obs),
      x = -prec.high * z,
      dims = c(N, n_obs)
    )

    # Q_{eta eta} = tau_eta I
    Qee <- Matrix::Diagonal(n = n_obs, x = prec.high)

    Qfull <- rbind(
      cbind(Qbb, Qbe),
      cbind(Matrix::t(Qbe), Qee)
    )

    return(Qfull)
  }

  mu <- function() {
    n_obs <- nrow(x)
    N <- n_regions
    return(numeric(N + n_obs))
  }

  log.norm.const <- function() {
    return(numeric(0))
  }

  log.prior <- function() {
    val <- (stats::dnorm(theta[1L], mean = 0, sd = 1, log = TRUE) +
              stats::dnorm(theta[2L], mean = 0, sd = 1, log = TRUE) +
              stats::dnorm(theta[3L], mean = 0, sd = 1, log = TRUE))
    return(val)
  }

  initial <- function() {
    return(c(0, 0, 0))
  }

  quit <- function() {
    return(invisible())
  }

  val <- do.call(match.arg(cmd), args = list())

  return(val)
}





#' Rgeneric MIDAS spatially-varying coefficient model with Gaussian constraint
#'
#' Defines a custom \code{rgeneric} model for use with the \code{INLA} framework,
#' implementing MIDAS-type lag weights using a Gaussian structure,
#' and where \code{beta[i]} varies for each location following an iCAR model
#'

#' @param cmd Character string indicating the INLA command.
#'   One of \code{"graph"}, \code{"Q"}, \code{"mu"}, \code{"initial"},
#'   \code{"log.norm.const"}, \code{"log.prior"}, or \code{"quit"}.
#' @param theta Numeric vector of hyperparameters controlling the MIDAS weights.
#'
#' @return Depends on \code{cmd}:
#' \itemize{
#'   \item \code{graph}: Sparse precision structure
#'   \item \code{Q}: Precision matrix
#'   \item \code{mu}: Mean vector
#'   \item \code{initial}: Initial values for \code{theta}
#'   \item \code{log.norm.const}: Normalizing constant (numeric(0))
#'   \item \code{log.prior}: Log prior density
#' }
#'
#' @details
#' The MIDAS lag weights are constructed using a normalized Gaussian scheme
#' polynomial transformation of lag indices. The parameters \code{theta[1]} and
#' \code{theta[2]} control the shape of the lag weighting function.
#'
#' @importFrom Matrix Diagonal
#'
rgeneric.svc.Gaussian.midas.icar = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                                                      "log.prior", "quit"),
                                              theta = NULL){

  envir <- parent.env(environment())

  x <- envir$x
  g <- envir$g

  regions <- sort(unique(x$region))
  n_regions <- length(regions)

  region_id <- match(x$region, regions)

  if (any(is.na(region_id))) {
    stop("Some observations have region labels not matched to 'regions'.")
  }
  if (!all(region_id %in% seq_len(n_regions))) {
    stop("region_id values are outside 1:n_regions.")
  }

  build_Q_icar_from_nbs <- function(g) {

    nbs <- g$nbs
    N <- length(nbs)

    ii <- integer(0)
    jj <- integer(0)

    for (i in seq_len(N)) {
      nei <- nbs[[i]]
      if (length(nei) > 0L) {
        ii <- c(ii, rep(i, length(nei)))
        jj <- c(jj, nei)
      }
    }

    W <- Matrix::sparseMatrix(
      i = ii,
      j = jj,
      x = 1,
      dims = c(N, N)
    )

    W <- (W + Matrix::t(W)) > 0
    W <- Matrix::Matrix(W * 1, sparse = TRUE)

    D <- Matrix::Diagonal(x = Matrix::rowSums(W))
    Q_icar <- D - W

    return(Q_icar)
  }

  Q_icar <- build_Q_icar_from_nbs(g)
  stopifnot(n_regions == nrow(Q_icar))
  stopifnot(all(region_id %in% seq_len(n_regions)))

  prec.high <- exp(15)
  prec_beta <- 1 # this is the precision for the global beta parameter

  interpret.theta <- function() {

    lag_cols <- grep("^lag[0-9]+$", names(x), value = TRUE)
    lag_k <- length(lag_cols) - 1L

    mu_val <- lag_k * (1 / (1 + exp(-theta[1L])))
    sigma_val = exp(theta[2L])

    psi <- numeric(lag_k + 1)
    for(lag in 0:lag_k){
      psi[lag + 1] <- exp( -(lag-mu_val)^2 / (2 * (sigma_val ^ 2)) )
    }

    tau_icar  <- exp(theta[3L])

    w <- psi / sum(psi)

    out_list <- list()

    out_list[["w"]] <- w
    out_list[["mu_val"]] <- mu_val
    out_list[["sigma_val"]] <- sigma_val
    out_list[["tau_icar"]] <- tau_icar

    return(out_list)
  }

  graph <- function() {

    n_obs <- nrow(x)
    N <- nrow(Q_icar)

    # global beta node
    Sgg <- Matrix::Matrix(1, nrow = 1, ncol = 1, sparse = TRUE)

    # global-to-region links: must be present because Qgb is nonzero
    Sgb <- Matrix::Matrix(1, nrow = 1, ncol = N, sparse = TRUE)

    # region random effects block
    Sbb <- Q_icar
    Sbb@x[] <- 1
    Sbb <- Sbb + Matrix::Diagonal(n = N, x = 1)
    Sbb@x[] <- 1

    # global beta to obs links
    Sge <- Matrix::Matrix(1, nrow = 1, ncol = n_obs, sparse = TRUE)

    # region beta_i to obs links
    Sbe <- Matrix::sparseMatrix(
      i = region_id,
      j = seq_len(n_obs),
      x = 1,
      dims = c(N, n_obs)
    )

    # observation block
    See <- Matrix::Diagonal(n = n_obs, x = 1)

    G <- rbind(
      cbind(Sgg, Sgb, Sge),
      cbind(Matrix::t(Sgb), Sbb, Sbe),
      cbind(Matrix::t(Sge), Matrix::t(Sbe), See)
    )

    G@x[] <- 1
    return(G)
  }


  Q <- function() {

    par <- interpret.theta()

    n_obs <- nrow(x)
    N <- nrow(Q_icar)

    lag_cols <- grep("^lag[0-9]+$", names(x), value = TRUE)

    Zmat <- as.matrix(x[, lag_cols, drop = FALSE])
    z <- as.vector(Zmat %*% par$w)

    prec_beta <- 1

    # global beta prior
    Qgg <- Matrix::Matrix(prec_beta + prec.high * sum(z^2),
                          nrow = 1, ncol = 1, sparse = TRUE)

    # global-region block
    zb_sum_by_region <- rowsum(z^2, group = region_id, reorder = FALSE)
    zb_sum_by_region <- as.numeric(zb_sum_by_region)

    if (length(zb_sum_by_region) < N) {
      tmp <- numeric(N)
      tmp[sort(unique(region_id))] <- zb_sum_by_region
      zb_sum_by_region <- tmp
    }

    Qgb <- Matrix::Matrix(prec.high * zb_sum_by_region,
                          nrow = 1, ncol = N, sparse = TRUE)

    # global-observation block
    Qge <- Matrix::sparseMatrix(
      i = rep(1L, n_obs),
      j = seq_len(n_obs),
      x = -prec.high * z,
      dims = c(1, n_obs)
    )

    # region-region block
    beta_diag <- rowsum(z^2, group = region_id, reorder = FALSE)
    beta_diag <- as.numeric(beta_diag)

    if (length(beta_diag) < N) {
      tmp <- numeric(N)
      tmp[sort(unique(region_id))] <- beta_diag
      beta_diag <- tmp
    }

    Qbb <- par$tau_icar * Q_icar + Matrix::Diagonal(n = N, x = prec.high * beta_diag)

    # region-observation block
    Qbe <- Matrix::sparseMatrix(
      i = region_id,
      j = seq_len(n_obs),
      x = -prec.high * z,
      dims = c(N, n_obs)
    )

    # observation block
    Qee <- Matrix::Diagonal(n = n_obs, x = prec.high)

    # assemble
    Qfull <- rbind(
      cbind(Qgg, Qgb, Qge),
      cbind(Matrix::t(Qgb), Qbb, Qbe),
      cbind(Matrix::t(Qge), Matrix::t(Qbe), Qee)
    )

    return(Qfull)
  }

  mu <- function() {
    n_obs <- nrow(x)
    N <- nrow(Q_icar)
    return(numeric(1 + N + n_obs))
  }

  log.norm.const <- function() {
    return(numeric(0))
  }

  log.prior <- function() {
    val <- (stats::dnorm(theta[1L], mean=0, sd=1, log=TRUE) +
              stats::dnorm(theta[2L], mean=0, sd=1, log=TRUE) +
              stats::dnorm(theta[3L], mean=0, sd=1, log=TRUE))
    return(val)
  }

  initial <- function() {
    return(c(0, 0, 0))
  }

  quit <- function() {
    return(invisible())
  }

  val <- do.call(match.arg(cmd), args = list())

  return(val)

}


#' Rgeneric MIDAS spatially-varying coefficient model with beta constraint.
#' Defines a custom \code{rgeneric} model for use with the \code{INLA} framework,
#' implementing MIDAS-type lag weights using a beta structure of dimension 1,
#' and where \code{beta[i]} varies for each location following an iid model.
#'

#' @param cmd Character string indicating the INLA command.
#'   One of \code{"graph"}, \code{"Q"}, \code{"mu"}, \code{"initial"},
#'   \code{"log.norm.const"}, \code{"log.prior"}, or \code{"quit"}.
#' @param theta Numeric vector of hyperparameters controlling the MIDAS weights.
#'
#' @return Depends on \code{cmd}:
#' \itemize{
#'   \item \code{graph}: Sparse precision structure
#'   \item \code{Q}: Precision matrix
#'   \item \code{mu}: Mean vector
#'   \item \code{initial}: Initial values for \code{theta}
#'   \item \code{log.norm.const}: Normalizing constant (numeric(0))
#'   \item \code{log.prior}: Log prior density
#' }
#'
#' @details
#' The MIDAS lag weights are constructed using a beta constraint with dimension
#' 1. The parameter \code{theta[1]}
#' control the shape of the lag weighting function.
#'
#' @importFrom Matrix Diagonal
#'
rgeneric.svc.Beta1.midas.iid = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                                                "log.prior", "quit"),
                                        theta = NULL){

  envir <- parent.env(environment())

  x <- envir$x

  regions <- sort(unique(x$region))
  n_regions <- length(regions)

  region_id <- match(x$region, regions)

  if (any(is.na(region_id))) {
    stop("Some observations have region labels not matched to 'regions'.")
  }
  if (!all(region_id %in% seq_len(n_regions))) {
    stop("region_id values are outside 1:n_regions.")
  }

  interpret.theta <- function() {

    lag_cols <- grep("^lag[0-9]+$", names(x), value = TRUE)
    lag_k <- length(lag_cols) - 1L

    gamma2 <- exp(theta[1L]) + 1
    tau_beta  <- exp(theta[2L])

    psi <- numeric(lag_k + 1L)
    xi <- 1e-4
    for (lag in 0:lag_k) {
      #x_temp <- 0.0001 + (1-0.0001)*((lag-1)/(lag_k-1))
      x_temp <- xi + (1 - 2 * xi) * (lag / lag_k)
      psi[lag + 1L] <- 1 + (1-x_temp)^(gamma2-1)
    }

    w <- psi / sum(psi)

    out_list <- list()
    out_list[["w"]] <- w
    out_list[["gamma2"]] <- gamma2
    out_list[["tau_beta"]] <- tau_beta

    return(out_list)
  }

  graph <- function() {

    n_obs <- nrow(x)
    N <- n_regions

    # region-specific beta_i block (iid => diagonal graph)
    Sbb <- Matrix::Diagonal(n = N, x = 1)

    # beta_i to observation links
    Sbe <- Matrix::sparseMatrix(
      i = region_id,
      j = seq_len(n_obs),
      x = 1,
      dims = c(N, n_obs)
    )

    # observation block
    See <- Matrix::Diagonal(n = n_obs, x = 1)

    G <- rbind(
      cbind(Sbb, Sbe),
      cbind(Matrix::t(Sbe), See)
    )

    G@x[] <- 1
    return(G)
  }

  Q <- function() {

    par <- interpret.theta()

    n_obs <- nrow(x)
    N <- n_regions
    prec.high <- exp(15)

    lag_cols <- grep("^lag[0-9]+$", names(x), value = TRUE)
    Zmat <- as.matrix(x[, lag_cols, drop = FALSE])
    z <- as.vector(Zmat %*% par$w)

    # Q_{beta beta} = tau_beta I + tau_eta A' D_z^2 A
    beta_diag <- rowsum(z^2, group = region_id, reorder = FALSE)
    beta_diag <- as.numeric(beta_diag)

    if (length(beta_diag) < N) {
      tmp <- numeric(N)
      tmp[sort(unique(region_id))] <- beta_diag
      beta_diag <- tmp
    }

    Qbb <- Matrix::Diagonal(n = N, x = par$tau_beta + prec.high * beta_diag)

    # Q_{beta eta} = -tau_eta A' D_z
    Qbe <- Matrix::sparseMatrix(
      i = region_id,
      j = seq_len(n_obs),
      x = -prec.high * z,
      dims = c(N, n_obs)
    )

    # Q_{eta eta} = tau_eta I
    Qee <- Matrix::Diagonal(n = n_obs, x = prec.high)

    Qfull <- rbind(
      cbind(Qbb, Qbe),
      cbind(Matrix::t(Qbe), Qee)
    )

    return(Qfull)
  }

  mu <- function() {
    n_obs <- nrow(x)
    N <- n_regions
    return(numeric(N + n_obs))
  }

  log.norm.const <- function() {
    return(numeric(0))
  }

  log.prior <- function() {
    val <- (stats::dnorm(theta[1L], mean = 0, sd = 1, log = TRUE) +
              stats::dnorm(theta[2L], mean = 0, sd = 1, log = TRUE))
    return(val)
  }

  initial <- function() {
    return(c(0, 0))
  }

  quit <- function() {
    return(invisible())
  }

  val <- do.call(match.arg(cmd), args = list())

  return(val)
}





#' Rgeneric MIDAS spatially-varying coefficient model with beta constraint.
#' Defines a custom \code{rgeneric} model for use with the \code{INLA} framework,
#' implementing MIDAS-type lag weights using a beta structure of dimension 1,
#' and where \code{beta[i]} varies for each location following an iCAR model.
#'

#' @param cmd Character string indicating the INLA command.
#'   One of \code{"graph"}, \code{"Q"}, \code{"mu"}, \code{"initial"},
#'   \code{"log.norm.const"}, \code{"log.prior"}, or \code{"quit"}.
#' @param theta Numeric vector of hyperparameters controlling the MIDAS weights.
#'
#' @return Depends on \code{cmd}:
#' \itemize{
#'   \item \code{graph}: Sparse precision structure
#'   \item \code{Q}: Precision matrix
#'   \item \code{mu}: Mean vector
#'   \item \code{initial}: Initial values for \code{theta}
#'   \item \code{log.norm.const}: Normalizing constant (numeric(0))
#'   \item \code{log.prior}: Log prior density
#' }
#'
#' @details
#' The MIDAS lag weights are constructed using a beta constraint with dimension
#' 1. The parameter \code{theta[1]}
#' control the shape of the lag weighting function.
#'
#' @importFrom Matrix Diagonal
#'
rgeneric.svc.Beta1.midas.icar = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                                                 "log.prior", "quit"),
                                         theta = NULL){

  envir <- parent.env(environment())

  x <- envir$x
  g <- envir$g

  regions <- sort(unique(x$region))
  n_regions <- length(regions)

  region_id <- match(x$region, regions)

  if (any(is.na(region_id))) {
    stop("Some observations have region labels not matched to 'regions'.")
  }
  if (!all(region_id %in% seq_len(n_regions))) {
    stop("region_id values are outside 1:n_regions.")
  }

  build_Q_icar_from_nbs <- function(g) {

    nbs <- g$nbs
    N <- length(nbs)

    ii <- integer(0)
    jj <- integer(0)

    for (i in seq_len(N)) {
      nei <- nbs[[i]]
      if (length(nei) > 0L) {
        ii <- c(ii, rep(i, length(nei)))
        jj <- c(jj, nei)
      }
    }

    W <- Matrix::sparseMatrix(
      i = ii,
      j = jj,
      x = 1,
      dims = c(N, N)
    )

    W <- (W + Matrix::t(W)) > 0
    W <- Matrix::Matrix(W * 1, sparse = TRUE)

    D <- Matrix::Diagonal(x = Matrix::rowSums(W))
    Q_icar <- D - W

    return(Q_icar)
  }

  Q_icar <- build_Q_icar_from_nbs(g)
  stopifnot(n_regions == nrow(Q_icar))
  stopifnot(all(region_id %in% seq_len(n_regions)))

  prec.high <- exp(15)
  prec_beta <- 1 # this is the precision for the global beta parameter

  interpret.theta <- function() {

    lag_cols <- grep("^lag[0-9]+$", names(x), value = TRUE)
    lag_k <- length(lag_cols) - 1L

    gamma2 <- exp(theta[1L]) + 1
    tau_icar  <- exp(theta[2L])

    psi <- numeric(lag_k + 1L)
    xi <- 1e-4
    for (lag in 0:lag_k) {
      #x_temp <- 0.0001 + (1-0.0001)*((lag-1)/(lag_k-1))
      x_temp <- xi + (1 - 2 * xi) * (lag / lag_k)
      psi[lag + 1L] <- 1 + (1-x_temp)^(gamma2-1)
    }

    w <- psi / sum(psi)

    out_list <- list()

    out_list[["w"]] <- w
    out_list[["gamma2"]] <- gamma2
    out_list[["tau_icar"]] <- tau_icar

    return(out_list)
  }

  graph <- function() {

    n_obs <- nrow(x)
    N <- nrow(Q_icar)

    # global beta node
    Sgg <- Matrix::Matrix(1, nrow = 1, ncol = 1, sparse = TRUE)

    # global-to-region links: must be present because Qgb is nonzero
    Sgb <- Matrix::Matrix(1, nrow = 1, ncol = N, sparse = TRUE)

    # region random effects block
    Sbb <- Q_icar
    Sbb@x[] <- 1
    Sbb <- Sbb + Matrix::Diagonal(n = N, x = 1)
    Sbb@x[] <- 1

    # global beta to obs links
    Sge <- Matrix::Matrix(1, nrow = 1, ncol = n_obs, sparse = TRUE)

    # region beta_i to obs links
    Sbe <- Matrix::sparseMatrix(
      i = region_id,
      j = seq_len(n_obs),
      x = 1,
      dims = c(N, n_obs)
    )

    # observation block
    See <- Matrix::Diagonal(n = n_obs, x = 1)

    G <- rbind(
      cbind(Sgg, Sgb, Sge),
      cbind(Matrix::t(Sgb), Sbb, Sbe),
      cbind(Matrix::t(Sge), Matrix::t(Sbe), See)
    )

    G@x[] <- 1
    return(G)
  }


  Q <- function() {

    par <- interpret.theta()

    n_obs <- nrow(x)
    N <- nrow(Q_icar)

    lag_cols <- grep("^lag[0-9]+$", names(x), value = TRUE)

    Zmat <- as.matrix(x[, lag_cols, drop = FALSE])
    z <- as.vector(Zmat %*% par$w)

    prec_beta <- 1

    # global beta prior
    Qgg <- Matrix::Matrix(prec_beta + prec.high * sum(z^2),
                          nrow = 1, ncol = 1, sparse = TRUE)

    # global-region block
    zb_sum_by_region <- rowsum(z^2, group = region_id, reorder = FALSE)
    zb_sum_by_region <- as.numeric(zb_sum_by_region)

    if (length(zb_sum_by_region) < N) {
      tmp <- numeric(N)
      tmp[sort(unique(region_id))] <- zb_sum_by_region
      zb_sum_by_region <- tmp
    }

    Qgb <- Matrix::Matrix(prec.high * zb_sum_by_region,
                          nrow = 1, ncol = N, sparse = TRUE)

    # global-observation block
    Qge <- Matrix::sparseMatrix(
      i = rep(1L, n_obs),
      j = seq_len(n_obs),
      x = -prec.high * z,
      dims = c(1, n_obs)
    )

    # region-region block
    beta_diag <- rowsum(z^2, group = region_id, reorder = FALSE)
    beta_diag <- as.numeric(beta_diag)

    if (length(beta_diag) < N) {
      tmp <- numeric(N)
      tmp[sort(unique(region_id))] <- beta_diag
      beta_diag <- tmp
    }

    Qbb <- par$tau_icar * Q_icar + Matrix::Diagonal(n = N, x = prec.high * beta_diag)

    # region-observation block
    Qbe <- Matrix::sparseMatrix(
      i = region_id,
      j = seq_len(n_obs),
      x = -prec.high * z,
      dims = c(N, n_obs)
    )

    # observation block
    Qee <- Matrix::Diagonal(n = n_obs, x = prec.high)

    # assemble
    Qfull <- rbind(
      cbind(Qgg, Qgb, Qge),
      cbind(Matrix::t(Qgb), Qbb, Qbe),
      cbind(Matrix::t(Qge), Matrix::t(Qbe), Qee)
    )

    return(Qfull)
  }

  mu <- function() {
    n_obs <- nrow(x)
    N <- nrow(Q_icar)
    return(numeric(1 + N + n_obs))
  }

  log.norm.const <- function() {
    return(numeric(0))
  }

  log.prior <- function() {
    val <- (stats::dnorm(theta[1L], mean=0, sd=1, log=TRUE) +
              stats::dnorm(theta[2L], mean=0, sd=1, log=TRUE))
    return(val)
  }

  initial <- function() {
    return(c(0, 0))
  }

  quit <- function() {
    return(invisible())
  }

  val <- do.call(match.arg(cmd), args = list())

  return(val)

}


#' Rgeneric MIDAS spatially-varying coefficient model with beta constraint.
#' Defines a custom \code{rgeneric} model for use with the \code{INLA} framework,
#' implementing MIDAS-type lag weights using a beta structure of dimension 1,
#' and where \code{beta} is constant for all locations.
#'

#' @param cmd Character string indicating the INLA command.
#'   One of \code{"graph"}, \code{"Q"}, \code{"mu"}, \code{"initial"},
#'   \code{"log.norm.const"}, \code{"log.prior"}, or \code{"quit"}.
#' @param theta Numeric vector of hyperparameters controlling the MIDAS weights.
#'
#' @return Depends on \code{cmd}:
#' \itemize{
#'   \item \code{graph}: Sparse precision structure
#'   \item \code{Q}: Precision matrix
#'   \item \code{mu}: Mean vector
#'   \item \code{initial}: Initial values for \code{theta}
#'   \item \code{log.norm.const}: Normalizing constant (numeric(0))
#'   \item \code{log.prior}: Log prior density
#' }
#'
#' @details
#' The MIDAS lag weights are constructed using a beta constraint with dimension
#' 1. The parameter \code{theta[1]} control the shape of the lag weighting
#' function.
#'
#' @importFrom Matrix Diagonal
#'
rgeneric.globalbeta.Beta1.midas = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                                                   "log.prior", "quit"),
                                           theta = NULL){

  envir = parent.env(environment())
  x <- envir$x

  ## artificial high precision to be added to the mean-model
  prec.high = exp(15)

  interpret.theta = function() {

    lag_cols <- grep("^lag[0-9]+$", names(x), value = TRUE)
    lag_k <- length(lag_cols) - 1L

    gamma2 <- exp(theta[1L]) + 1
    beta1 <- theta[2L]

    psi <- numeric(lag_k + 1L)
    xi <- 1e-4
    for (lag in 0:lag_k) {
      #x_temp <- 0.0001 + (1-0.0001)*((lag-1)/(lag_k-1))
      x_temp <- xi + (1 - 2 * xi) * (lag / lag_k)
      psi[lag + 1L] <- 1 + (1-x_temp)^(gamma2-1)
    }

    w <- psi / sum(psi)

    out_list <- list()
    for(lag in 0:lag_k){
      out_list[[paste0("w",lag)]] <- w[lag + 1]
    }

    out_list[["beta1"]] <- beta1

    return(out_list)

  }
  graph = function() {
    G = Matrix::Diagonal(n = length(x$lag0), x=1)
    return(G)
  }
  Q = function() {
    Q = prec.high * graph()
    return(Q)
  }
  mu = function() {
    par = interpret.theta()

    lag_cols <- grep("^lag[0-9]+$", names(x), value = TRUE)
    lag_k <- length(lag_cols) - 1L

    compile_lag_label <- c()
    for(lag in 0:lag_k){
      compile_lag_label <- c(compile_lag_label, paste0("lag",lag))
    }
    compile_w_label <- c()
    for(lag in 0:lag_k){
      compile_w_label <- c(compile_w_label, paste0("w",lag))
    }

    agg <- 0
    for(lag in 0:lag_k){
      agg <- agg + par[[compile_w_label[[lag+1]]]] * x[,which(names(x) == compile_lag_label[[lag+1]])]
    }

    return(par$beta1 * agg)
  }
  log.norm.const = function() {
    return(numeric(0))
  }
  log.prior = function() {
    val = (stats::dnorm(theta[1L], mean=0, sd=1, log=TRUE) +
             stats::dnorm(theta[2L], mean=0, sd=1, log=TRUE))
    return(val)
  }
  initial = function() {
    return(rep(1, 2))
  }
  quit = function() {
    return(invisible())
  }

  val = do.call(match.arg(cmd), args = list())
  return(val)
}




#' Rgeneric MIDAS spatially-varying coefficient model with beta constraint.
#' Defines a custom \code{rgeneric} model for use with the \code{INLA} framework,
#' implementing MIDAS-type lag weights using a beta structure of dimension 2,
#' and where \code{beta[i]} varies for each location following an iid model.
#'

#' @param cmd Character string indicating the INLA command.
#'   One of \code{"graph"}, \code{"Q"}, \code{"mu"}, \code{"initial"},
#'   \code{"log.norm.const"}, \code{"log.prior"}, or \code{"quit"}.
#' @param theta Numeric vector of hyperparameters controlling the MIDAS weights.
#'
#' @return Depends on \code{cmd}:
#' \itemize{
#'   \item \code{graph}: Sparse precision structure
#'   \item \code{Q}: Precision matrix
#'   \item \code{mu}: Mean vector
#'   \item \code{initial}: Initial values for \code{theta}
#'   \item \code{log.norm.const}: Normalizing constant (numeric(0))
#'   \item \code{log.prior}: Log prior density
#' }
#'
#' @details
#' The MIDAS lag weights are constructed using a beta constraint with dimension
#' 1. The parameters \code{theta[1]} and \code{theta[2]}
#' control the shape of the lag weighting function.
#'
#' @importFrom Matrix Diagonal
#'
rgeneric.svc.Beta2.midas.iid = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                                                "log.prior", "quit"),
                                        theta = NULL){

  envir <- parent.env(environment())

  x <- envir$x

  regions <- sort(unique(x$region))
  n_regions <- length(regions)

  region_id <- match(x$region, regions)

  if (any(is.na(region_id))) {
    stop("Some observations have region labels not matched to 'regions'.")
  }
  if (!all(region_id %in% seq_len(n_regions))) {
    stop("region_id values are outside 1:n_regions.")
  }

  interpret.theta <- function() {

    lag_cols <- grep("^lag[0-9]+$", names(x), value = TRUE)
    lag_k <- length(lag_cols) - 1L

    gamma1 <- exp(theta[1L]) + 1
    gamma2 <- exp(theta[2L]) + 1
    tau_beta  <- exp(theta[3L])

    psi <- numeric(lag_k + 1L)
    xi <- 1e-4
    for(lag in 0:lag_k){
      #x_temp <- 0.0001 + (1-0.0001)*((lag-1)/(lag_k-1))
      x_temp <- xi + (1 - 2 * xi) * (lag / lag_k)
      psi[lag + 1L] <- (x_temp^(gamma1-1))*((1-x_temp)^(gamma2-1))
    }

    w <- psi / sum(psi)

    out_list <- list()
    out_list[["w"]] <- w
    out_list[["gamma1"]] <- gamma1
    out_list[["gamma2"]] <- gamma2
    out_list[["tau_beta"]] <- tau_beta

    return(out_list)
  }

  graph <- function() {

    n_obs <- nrow(x)
    N <- n_regions

    # region-specific beta_i block (iid => diagonal graph)
    Sbb <- Matrix::Diagonal(n = N, x = 1)

    # beta_i to observation links
    Sbe <- Matrix::sparseMatrix(
      i = region_id,
      j = seq_len(n_obs),
      x = 1,
      dims = c(N, n_obs)
    )

    # observation block
    See <- Matrix::Diagonal(n = n_obs, x = 1)

    G <- rbind(
      cbind(Sbb, Sbe),
      cbind(Matrix::t(Sbe), See)
    )

    G@x[] <- 1
    return(G)
  }

  Q <- function() {

    par <- interpret.theta()

    n_obs <- nrow(x)
    N <- n_regions
    prec.high <- exp(15)

    lag_cols <- grep("^lag[0-9]+$", names(x), value = TRUE)
    Zmat <- as.matrix(x[, lag_cols, drop = FALSE])
    z <- as.vector(Zmat %*% par$w)

    # Q_{beta beta} = tau_beta I + tau_eta A' D_z^2 A
    beta_diag <- rowsum(z^2, group = region_id, reorder = FALSE)
    beta_diag <- as.numeric(beta_diag)

    if (length(beta_diag) < N) {
      tmp <- numeric(N)
      tmp[sort(unique(region_id))] <- beta_diag
      beta_diag <- tmp
    }

    Qbb <- Matrix::Diagonal(n = N, x = par$tau_beta + prec.high * beta_diag)

    # Q_{beta eta} = -tau_eta A' D_z
    Qbe <- Matrix::sparseMatrix(
      i = region_id,
      j = seq_len(n_obs),
      x = -prec.high * z,
      dims = c(N, n_obs)
    )

    # Q_{eta eta} = tau_eta I
    Qee <- Matrix::Diagonal(n = n_obs, x = prec.high)

    Qfull <- rbind(
      cbind(Qbb, Qbe),
      cbind(Matrix::t(Qbe), Qee)
    )

    return(Qfull)
  }

  mu <- function() {
    n_obs <- nrow(x)
    N <- n_regions
    return(numeric(N + n_obs))
  }

  log.norm.const <- function() {
    return(numeric(0))
  }

  log.prior <- function() {
    val <- (stats::dnorm(theta[1L], mean = 0, sd = 1, log = TRUE) +
              stats::dnorm(theta[2L], mean = 0, sd = 1, log = TRUE) +
              stats::dnorm(theta[3L], mean = 0, sd = 1, log = TRUE))
    return(val)
  }

  initial <- function() {
    return(c(0, 0, 0))
  }

  quit <- function() {
    return(invisible())
  }

  val <- do.call(match.arg(cmd), args = list())

  return(val)
}



#' Rgeneric MIDAS spatially-varying coefficient model with beta constraint.
#' Defines a custom \code{rgeneric} model for use with the \code{INLA} framework,
#' implementing MIDAS-type lag weights using a beta structure of dimension 2,
#' and where \code{beta[i]} varies for each location following an iCAR model.
#'

#' @param cmd Character string indicating the INLA command.
#'   One of \code{"graph"}, \code{"Q"}, \code{"mu"}, \code{"initial"},
#'   \code{"log.norm.const"}, \code{"log.prior"}, or \code{"quit"}.
#' @param theta Numeric vector of hyperparameters controlling the MIDAS weights.
#'
#' @return Depends on \code{cmd}:
#' \itemize{
#'   \item \code{graph}: Sparse precision structure
#'   \item \code{Q}: Precision matrix
#'   \item \code{mu}: Mean vector
#'   \item \code{initial}: Initial values for \code{theta}
#'   \item \code{log.norm.const}: Normalizing constant (numeric(0))
#'   \item \code{log.prior}: Log prior density
#' }
#'
#' @details
#' The MIDAS lag weights are constructed using a beta constraint with dimension
#' 1. The parameters \code{theta[1]} and \code{theta[]}
#' control the shape of the lag weighting function.
#'
#' @importFrom Matrix Diagonal
#'
rgeneric.svc.Beta2.midas.icar = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                                                 "log.prior", "quit"),
                                         theta = NULL){

  envir <- parent.env(environment())

  x <- envir$x
  g <- envir$g

  regions <- sort(unique(x$region))
  n_regions <- length(regions)

  region_id <- match(x$region, regions)

  if (any(is.na(region_id))) {
    stop("Some observations have region labels not matched to 'regions'.")
  }
  if (!all(region_id %in% seq_len(n_regions))) {
    stop("region_id values are outside 1:n_regions.")
  }

  build_Q_icar_from_nbs <- function(g) {

    nbs <- g$nbs
    N <- length(nbs)

    ii <- integer(0)
    jj <- integer(0)

    for (i in seq_len(N)) {
      nei <- nbs[[i]]
      if (length(nei) > 0L) {
        ii <- c(ii, rep(i, length(nei)))
        jj <- c(jj, nei)
      }
    }

    W <- Matrix::sparseMatrix(
      i = ii,
      j = jj,
      x = 1,
      dims = c(N, N)
    )

    W <- (W + Matrix::t(W)) > 0
    W <- Matrix::Matrix(W * 1, sparse = TRUE)

    D <- Matrix::Diagonal(x = Matrix::rowSums(W))
    Q_icar <- D - W

    return(Q_icar)
  }

  Q_icar <- build_Q_icar_from_nbs(g)
  stopifnot(n_regions == nrow(Q_icar))
  stopifnot(all(region_id %in% seq_len(n_regions)))

  prec.high <- exp(15)
  prec_beta <- 1 # this is the precision for the global beta parameter

  interpret.theta <- function() {

    lag_cols <- grep("^lag[0-9]+$", names(x), value = TRUE)
    lag_k <- length(lag_cols) - 1L

    gamma1 <- exp(theta[1L]) + 1
    gamma2 <- exp(theta[2L]) + 1
    tau_icar  <- exp(theta[3L])


    psi <- numeric(lag_k + 1L)
    xi <- 1e-4
    for(lag in 0:lag_k){
      #x_temp <- 0.0001 + (1-0.0001)*((lag-1)/(lag_k-1))
      x_temp <- xi + (1 - 2 * xi) * (lag / lag_k)
      psi[lag + 1L] <- (x_temp^(gamma1-1))*((1-x_temp)^(gamma2-1))
    }

    w <- psi / sum(psi)

    out_list <- list()

    out_list[["w"]] <- w
    out_list[["gamma1"]] <- gamma1
    out_list[["gamma2"]] <- gamma2
    out_list[["tau_icar"]] <- tau_icar

    return(out_list)
  }

  graph <- function() {

    n_obs <- nrow(x)
    N <- nrow(Q_icar)

    # global beta node
    Sgg <- Matrix::Matrix(1, nrow = 1, ncol = 1, sparse = TRUE)

    # global-to-region links: must be present because Qgb is nonzero
    Sgb <- Matrix::Matrix(1, nrow = 1, ncol = N, sparse = TRUE)

    # region random effects block
    Sbb <- Q_icar
    Sbb@x[] <- 1
    Sbb <- Sbb + Matrix::Diagonal(n = N, x = 1)
    Sbb@x[] <- 1

    # global beta to obs links
    Sge <- Matrix::Matrix(1, nrow = 1, ncol = n_obs, sparse = TRUE)

    # region beta_i to obs links
    Sbe <- Matrix::sparseMatrix(
      i = region_id,
      j = seq_len(n_obs),
      x = 1,
      dims = c(N, n_obs)
    )

    # observation block
    See <- Matrix::Diagonal(n = n_obs, x = 1)

    G <- rbind(
      cbind(Sgg, Sgb, Sge),
      cbind(Matrix::t(Sgb), Sbb, Sbe),
      cbind(Matrix::t(Sge), Matrix::t(Sbe), See)
    )

    G@x[] <- 1
    return(G)
  }


  Q <- function() {

    par <- interpret.theta()

    n_obs <- nrow(x)
    N <- nrow(Q_icar)

    lag_cols <- grep("^lag[0-9]+$", names(x), value = TRUE)

    Zmat <- as.matrix(x[, lag_cols, drop = FALSE])
    z <- as.vector(Zmat %*% par$w)

    prec_beta <- 1

    # global beta prior
    Qgg <- Matrix::Matrix(prec_beta + prec.high * sum(z^2),
                          nrow = 1, ncol = 1, sparse = TRUE)

    # global-region block
    zb_sum_by_region <- rowsum(z^2, group = region_id, reorder = FALSE)
    zb_sum_by_region <- as.numeric(zb_sum_by_region)

    if (length(zb_sum_by_region) < N) {
      tmp <- numeric(N)
      tmp[sort(unique(region_id))] <- zb_sum_by_region
      zb_sum_by_region <- tmp
    }

    Qgb <- Matrix::Matrix(prec.high * zb_sum_by_region,
                          nrow = 1, ncol = N, sparse = TRUE)

    # global-observation block
    Qge <- Matrix::sparseMatrix(
      i = rep(1L, n_obs),
      j = seq_len(n_obs),
      x = -prec.high * z,
      dims = c(1, n_obs)
    )

    # region-region block
    beta_diag <- rowsum(z^2, group = region_id, reorder = FALSE)
    beta_diag <- as.numeric(beta_diag)

    if (length(beta_diag) < N) {
      tmp <- numeric(N)
      tmp[sort(unique(region_id))] <- beta_diag
      beta_diag <- tmp
    }

    Qbb <- par$tau_icar * Q_icar + Matrix::Diagonal(n = N, x = prec.high * beta_diag)

    # region-observation block
    Qbe <- Matrix::sparseMatrix(
      i = region_id,
      j = seq_len(n_obs),
      x = -prec.high * z,
      dims = c(N, n_obs)
    )

    # observation block
    Qee <- Matrix::Diagonal(n = n_obs, x = prec.high)

    # assemble
    Qfull <- rbind(
      cbind(Qgg, Qgb, Qge),
      cbind(Matrix::t(Qgb), Qbb, Qbe),
      cbind(Matrix::t(Qge), Matrix::t(Qbe), Qee)
    )

    return(Qfull)
  }

  mu <- function() {
    n_obs <- nrow(x)
    N <- nrow(Q_icar)
    return(numeric(1 + N + n_obs))
  }

  log.norm.const <- function() {
    return(numeric(0))
  }

  log.prior <- function() {
    val <- (stats::dnorm(theta[1L], mean=0, sd=1, log=TRUE) +
              stats::dnorm(theta[2L], mean=0, sd=1, log=TRUE) +
              stats::dnorm(theta[3L], mean=0, sd=1, log=TRUE))
    return(val)
  }

  initial <- function() {
    return(c(0, 0, 0))
  }

  quit <- function() {
    return(invisible())
  }

  val <- do.call(match.arg(cmd), args = list())

  return(val)

}



#' Rgeneric MIDAS spatially-varying coefficient model with beta constraint.
#' Defines a custom \code{rgeneric} model for use with the \code{INLA} framework,
#' implementing MIDAS-type lag weights using a beta structure of dimension 2,
#' and where \code{beta} is constant for all locations.
#'

#' @param cmd Character string indicating the INLA command.
#'   One of \code{"graph"}, \code{"Q"}, \code{"mu"}, \code{"initial"},
#'   \code{"log.norm.const"}, \code{"log.prior"}, or \code{"quit"}.
#' @param theta Numeric vector of hyperparameters controlling the MIDAS weights.
#'
#' @return Depends on \code{cmd}:
#' \itemize{
#'   \item \code{graph}: Sparse precision structure
#'   \item \code{Q}: Precision matrix
#'   \item \code{mu}: Mean vector
#'   \item \code{initial}: Initial values for \code{theta}
#'   \item \code{log.norm.const}: Normalizing constant (numeric(0))
#'   \item \code{log.prior}: Log prior density
#' }
#'
#' @details
#' The MIDAS lag weights are constructed using a beta constraint with dimension
#' 1. The parameters \code{theta[1]} and \code{theta[2]} control the shape of
#' the lag weighting function.
#'
#' @importFrom Matrix Diagonal
#'
rgeneric.globalbeta.Beta2.midas = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                                                   "log.prior", "quit"),
                                           theta = NULL){

  envir = parent.env(environment())
  x <- envir$x

  ## artificial high precision to be added to the mean-model
  prec.high = exp(15)

  interpret.theta = function() {

    lag_cols <- grep("^lag[0-9]+$", names(x), value = TRUE)
    lag_k <- length(lag_cols) - 1L

    gamma1 <- exp(theta[1L]) + 1
    gamma2 <- exp(theta[2L]) + 1
    beta1 <- theta[3L]

    psi <- numeric(lag_k + 1L)
    xi <- 1e-4
    for(lag in 0:lag_k){
      #x_temp <- 0.0001 + (1-0.0001)*((lag-1)/(lag_k-1))
      x_temp <- xi + (1 - 2 * xi) * (lag / lag_k)
      psi[lag + 1L] <- (x_temp^(gamma1-1))*((1-x_temp)^(gamma2-1))
    }

    w <- psi / sum(psi)

    out_list <- list()
    for(lag in 0:lag_k){
      out_list[[paste0("w",lag)]] <- w[lag + 1]
    }
    out_list[["beta1"]] <- beta1

    return(out_list)

  }
  graph = function() {
    G = Matrix::Diagonal(n = length(x$lag0), x=1)
    return(G)
  }
  Q = function() {
    Q = prec.high * graph()
    return(Q)
  }
  mu = function() {
    par = interpret.theta()

    lag_cols <- grep("^lag[0-9]+$", names(x), value = TRUE)
    lag_k <- length(lag_cols) - 1L

    compile_lag_label <- c()
    for(lag in 0:lag_k){
      compile_lag_label <- c(compile_lag_label, paste0("lag",lag))
    }
    compile_w_label <- c()
    for(lag in 0:lag_k){
      compile_w_label <- c(compile_w_label, paste0("w",lag))
    }

    agg <- 0
    for(lag in 0:lag_k){
      agg <- agg + par[[compile_w_label[[lag+1]]]] * x[,which(names(x) == compile_lag_label[[lag+1]])]
    }

    return(par$beta1 * agg)
  }
  log.norm.const = function() {
    return(numeric(0))
  }
  log.prior = function() {
    val = (stats::dnorm(theta[1L], mean=0, sd=1, log=TRUE) +
             stats::dnorm(theta[2L], mean=0, sd=1, log=TRUE) +
             stats::dnorm(theta[3L], mean=0, sd=1, log=TRUE))
    return(val)
  }
  initial = function() {
    return(rep(1, 3))
  }
  quit = function() {
    return(invisible())
  }

  val = do.call(match.arg(cmd), args = list())
  return(val)
}



#' Rgeneric MIDAS spatially-varying coefficient model with Almon polynomical
#' constraint, d = 2
#'
#' Defines a custom \code{rgeneric} model for use with the \code{INLA} framework,
#' implementing MIDAS-type lag weights using Almon lag structure,
#' and where \code{beta[i]} varies for each location following an iid model.
#'

#' @param cmd Character string indicating the INLA command.
#'   One of \code{"graph"}, \code{"Q"}, \code{"mu"}, \code{"initial"},
#'   \code{"log.norm.const"}, \code{"log.prior"}, or \code{"quit"}.
#' @param theta Numeric vector of hyperparameters controlling the MIDAS weights.
#'
#' @return Depends on \code{cmd}:
#' \itemize{
#'   \item \code{graph}: Sparse precision structure
#'   \item \code{Q}: Precision matrix
#'   \item \code{mu}: Mean vector
#'   \item \code{initial}: Initial values for \code{theta}
#'   \item \code{log.norm.const}: Normalizing constant (numeric(0))
#'   \item \code{log.prior}: Log prior density
#' }
#'
#' @details
#' The MIDAS lag weights are constructed using a normalized Almon constraint.
#' The parameters \code{theta[1]} and and \code{theta[2]}
#' control the shape of the lag weighting function.
#'
#' @importFrom Matrix Diagonal
#'
rgeneric.svc.Almon2.midas.iid = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                                                 "log.prior", "quit"),
                                         theta = NULL){

  envir <- parent.env(environment())

  x <- envir$x

  regions <- sort(unique(x$region))
  n_regions <- length(regions)

  region_id <- match(x$region, regions)

  if (any(is.na(region_id))) {
    stop("Some observations have region labels not matched to 'regions'.")
  }
  if (!all(region_id %in% seq_len(n_regions))) {
    stop("region_id values are outside 1:n_regions.")
  }

  interpret.theta <- function() {

    lag_cols <- grep("^lag[0-9]+$", names(x), value = TRUE)
    lag_k <- length(lag_cols) - 1L

    #gamma1 <- 0.01*sin(theta[1L])
    #gamma2 <- 0.01*sin(theta[2L])

    gamma1 <- 0.01 * tanh(theta[1])
    gamma2 <- 0.01 * tanh(theta[2])
    tau_beta  <- exp(theta[3L])

    psi <- numeric(lag_k + 1L)
    for(lag in 0:lag_k){
      psi[lag + 1L] <- exp(gamma1*(lag^1) + gamma2*(lag^2))
    }

    w <- psi / sum(psi)

    out_list <- list()
    out_list[["w"]] <- w
    out_list[["gamma1"]] <- gamma1
    out_list[["gamma2"]] <- gamma2
    out_list[["tau_beta"]] <- tau_beta

    return(out_list)
  }

  graph <- function() {

    n_obs <- nrow(x)
    N <- n_regions

    # region-specific beta_i block (iid => diagonal graph)
    Sbb <- Matrix::Diagonal(n = N, x = 1)

    # beta_i to observation links
    Sbe <- Matrix::sparseMatrix(
      i = region_id,
      j = seq_len(n_obs),
      x = 1,
      dims = c(N, n_obs)
    )

    # observation block
    See <- Matrix::Diagonal(n = n_obs, x = 1)

    G <- rbind(
      cbind(Sbb, Sbe),
      cbind(Matrix::t(Sbe), See)
    )

    G@x[] <- 1
    return(G)
  }

  Q <- function() {

    par <- interpret.theta()

    n_obs <- nrow(x)
    N <- n_regions
    prec.high <- exp(15)

    lag_cols <- grep("^lag[0-9]+$", names(x), value = TRUE)
    Zmat <- as.matrix(x[, lag_cols, drop = FALSE])
    z <- as.vector(Zmat %*% par$w)

    # Q_{beta beta} = tau_beta I + tau_eta A' D_z^2 A
    beta_diag <- rowsum(z^2, group = region_id, reorder = FALSE)
    beta_diag <- as.numeric(beta_diag)

    if (length(beta_diag) < N) {
      tmp <- numeric(N)
      tmp[sort(unique(region_id))] <- beta_diag
      beta_diag <- tmp
    }

    Qbb <- Matrix::Diagonal(n = N, x = par$tau_beta + prec.high * beta_diag)

    # Q_{beta eta} = -tau_eta A' D_z
    Qbe <- Matrix::sparseMatrix(
      i = region_id,
      j = seq_len(n_obs),
      x = -prec.high * z,
      dims = c(N, n_obs)
    )

    # Q_{eta eta} = tau_eta I
    Qee <- Matrix::Diagonal(n = n_obs, x = prec.high)

    Qfull <- rbind(
      cbind(Qbb, Qbe),
      cbind(Matrix::t(Qbe), Qee)
    )

    return(Qfull)
  }

  mu <- function() {
    n_obs <- nrow(x)
    N <- n_regions
    return(numeric(N + n_obs))
  }

  log.norm.const <- function() {
    return(numeric(0))
  }

  log.prior <- function() {
    val <- (stats::dnorm(theta[1L], mean = 0, sd = 1, log = TRUE) +
              stats::dnorm(theta[2L], mean = 0, sd = 1, log = TRUE) +
              stats::dnorm(theta[3L], mean = 0, sd = 1, log = TRUE))
    return(val)
  }

  initial <- function() {
    return(c(0, 0, 0))
  }

  quit <- function() {
    return(invisible())
  }

  val <- do.call(match.arg(cmd), args = list())

  return(val)
}



#' Rgeneric MIDAS spatially-varying coefficient model with Almon polynomical
#' constraint, d = 2
#'
#' Defines a custom \code{rgeneric} model for use with the \code{INLA} framework,
#' implementing MIDAS-type lag weights using Almon lag structure,
#' and where \code{beta[i]} varies for each location following an iCAR model.
#'

#' @param cmd Character string indicating the INLA command.
#'   One of \code{"graph"}, \code{"Q"}, \code{"mu"}, \code{"initial"},
#'   \code{"log.norm.const"}, \code{"log.prior"}, or \code{"quit"}.
#' @param theta Numeric vector of hyperparameters controlling the MIDAS weights.
#'
#' @return Depends on \code{cmd}:
#' \itemize{
#'   \item \code{graph}: Sparse precision structure
#'   \item \code{Q}: Precision matrix
#'   \item \code{mu}: Mean vector
#'   \item \code{initial}: Initial values for \code{theta}
#'   \item \code{log.norm.const}: Normalizing constant (numeric(0))
#'   \item \code{log.prior}: Log prior density
#' }
#'
#' @details
#' The MIDAS lag weights are constructed using a normalized Almon constraint.
#' The parameters \code{theta[1]} and and \code{theta[2]}
#' control the shape of the lag weighting function.
#'
#' @importFrom Matrix Diagonal
#'
rgeneric.svc.Almon2.midas.icar = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                                                  "log.prior", "quit"),
                                          theta = NULL){

  envir <- parent.env(environment())

  x <- envir$x
  g <- envir$g

  regions <- sort(unique(x$region))
  n_regions <- length(regions)

  region_id <- match(x$region, regions)

  if (any(is.na(region_id))) {
    stop("Some observations have region labels not matched to 'regions'.")
  }
  if (!all(region_id %in% seq_len(n_regions))) {
    stop("region_id values are outside 1:n_regions.")
  }

  build_Q_icar_from_nbs <- function(g) {

    nbs <- g$nbs
    N <- length(nbs)

    ii <- integer(0)
    jj <- integer(0)

    for (i in seq_len(N)) {
      nei <- nbs[[i]]
      if (length(nei) > 0L) {
        ii <- c(ii, rep(i, length(nei)))
        jj <- c(jj, nei)
      }
    }

    W <- Matrix::sparseMatrix(
      i = ii,
      j = jj,
      x = 1,
      dims = c(N, N)
    )

    W <- (W + Matrix::t(W)) > 0
    W <- Matrix::Matrix(W * 1, sparse = TRUE)

    D <- Matrix::Diagonal(x = Matrix::rowSums(W))
    Q_icar <- D - W

    return(Q_icar)
  }

  Q_icar <- build_Q_icar_from_nbs(g)
  stopifnot(n_regions == nrow(Q_icar))
  stopifnot(all(region_id %in% seq_len(n_regions)))

  prec.high <- exp(15)
  prec_beta <- 1 # this is the precision for the global beta parameter

  interpret.theta <- function() {

    lag_cols <- grep("^lag[0-9]+$", names(x), value = TRUE)
    lag_k <- length(lag_cols) - 1L

    gamma1 <- 0.01 * tanh(theta[1])
    gamma2 <- 0.01 * tanh(theta[2])
    tau_icar  <- exp(theta[3L])

    psi <- numeric(lag_k + 1L)
    for(lag in 0:lag_k){
      psi[lag + 1L] <- exp(gamma1*(lag^1) + gamma2*(lag^2))
    }

    w <- psi / sum(psi)

    out_list <- list()

    out_list[["w"]] <- w
    out_list[["gamma1"]] <- gamma1
    out_list[["gamma2"]] <- gamma2
    out_list[["tau_icar"]] <- tau_icar

    return(out_list)
  }

  graph <- function() {

    n_obs <- nrow(x)
    N <- nrow(Q_icar)

    # global beta node
    Sgg <- Matrix::Matrix(1, nrow = 1, ncol = 1, sparse = TRUE)

    # global-to-region links: must be present because Qgb is nonzero
    Sgb <- Matrix::Matrix(1, nrow = 1, ncol = N, sparse = TRUE)

    # region random effects block
    Sbb <- Q_icar
    Sbb@x[] <- 1
    Sbb <- Sbb + Matrix::Diagonal(n = N, x = 1)
    Sbb@x[] <- 1

    # global beta to obs links
    Sge <- Matrix::Matrix(1, nrow = 1, ncol = n_obs, sparse = TRUE)

    # region beta_i to obs links
    Sbe <- Matrix::sparseMatrix(
      i = region_id,
      j = seq_len(n_obs),
      x = 1,
      dims = c(N, n_obs)
    )

    # observation block
    See <- Matrix::Diagonal(n = n_obs, x = 1)

    G <- rbind(
      cbind(Sgg, Sgb, Sge),
      cbind(Matrix::t(Sgb), Sbb, Sbe),
      cbind(Matrix::t(Sge), Matrix::t(Sbe), See)
    )

    G@x[] <- 1
    return(G)
  }


  Q <- function() {

    par <- interpret.theta()

    n_obs <- nrow(x)
    N <- nrow(Q_icar)

    lag_cols <- grep("^lag[0-9]+$", names(x), value = TRUE)

    Zmat <- as.matrix(x[, lag_cols, drop = FALSE])
    z <- as.vector(Zmat %*% par$w)

    prec_beta <- 1

    # global beta prior
    Qgg <- Matrix::Matrix(prec_beta + prec.high * sum(z^2),
                          nrow = 1, ncol = 1, sparse = TRUE)

    # global-region block
    zb_sum_by_region <- rowsum(z^2, group = region_id, reorder = FALSE)
    zb_sum_by_region <- as.numeric(zb_sum_by_region)

    if (length(zb_sum_by_region) < N) {
      tmp <- numeric(N)
      tmp[sort(unique(region_id))] <- zb_sum_by_region
      zb_sum_by_region <- tmp
    }

    Qgb <- Matrix::Matrix(prec.high * zb_sum_by_region,
                          nrow = 1, ncol = N, sparse = TRUE)

    # global-observation block
    Qge <- Matrix::sparseMatrix(
      i = rep(1L, n_obs),
      j = seq_len(n_obs),
      x = -prec.high * z,
      dims = c(1, n_obs)
    )

    # region-region block
    beta_diag <- rowsum(z^2, group = region_id, reorder = FALSE)
    beta_diag <- as.numeric(beta_diag)

    if (length(beta_diag) < N) {
      tmp <- numeric(N)
      tmp[sort(unique(region_id))] <- beta_diag
      beta_diag <- tmp
    }

    Qbb <- par$tau_icar * Q_icar + Matrix::Diagonal(n = N, x = prec.high * beta_diag)

    # region-observation block
    Qbe <- Matrix::sparseMatrix(
      i = region_id,
      j = seq_len(n_obs),
      x = -prec.high * z,
      dims = c(N, n_obs)
    )

    # observation block
    Qee <- Matrix::Diagonal(n = n_obs, x = prec.high)

    # assemble
    Qfull <- rbind(
      cbind(Qgg, Qgb, Qge),
      cbind(Matrix::t(Qgb), Qbb, Qbe),
      cbind(Matrix::t(Qge), Matrix::t(Qbe), Qee)
    )

    return(Qfull)
  }

  mu <- function() {
    n_obs <- nrow(x)
    N <- nrow(Q_icar)
    return(numeric(1 + N + n_obs))
  }

  log.norm.const <- function() {
    return(numeric(0))
  }

  log.prior <- function() {
    val <- (stats::dnorm(theta[1L], mean=0, sd=1, log=TRUE) +
              stats::dnorm(theta[2L], mean=0, sd=1, log=TRUE) +
              stats::dnorm(theta[3L], mean=0, sd=1, log=TRUE))
    return(val)
  }

  initial <- function() {
    return(c(0, 0, 0))
  }

  quit <- function() {
    return(invisible())
  }

  val <- do.call(match.arg(cmd), args = list())

  return(val)

}


#' Rgeneric MIDAS spatially-varying coefficient model with Almon polynomical
#' constraint, d = 2
#'
#' Defines a custom \code{rgeneric} model for use with the \code{INLA} framework,
#' implementing MIDAS-type lag weights using Almon lag structure,
#' and where \code{beta} is constant for all locations.
#'

#' @param cmd Character string indicating the INLA command.
#'   One of \code{"graph"}, \code{"Q"}, \code{"mu"}, \code{"initial"},
#'   \code{"log.norm.const"}, \code{"log.prior"}, or \code{"quit"}.
#' @param theta Numeric vector of hyperparameters controlling the MIDAS weights.
#'
#' @return Depends on \code{cmd}:
#' \itemize{
#'   \item \code{graph}: Sparse precision structure
#'   \item \code{Q}: Precision matrix
#'   \item \code{mu}: Mean vector
#'   \item \code{initial}: Initial values for \code{theta}
#'   \item \code{log.norm.const}: Normalizing constant (numeric(0))
#'   \item \code{log.prior}: Log prior density
#' }
#'
#' @details
#' The MIDAS lag weights are constructed using a normalized Almon constraint.
#' The parameters \code{theta[1]} and and \code{theta[2]}
#' control the shape of the lag weighting function.
#'
#' @importFrom Matrix Diagonal
#'
rgeneric.globalbeta.Almon2.midas = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                                                    "log.prior", "quit"),
                                            theta = NULL){

  envir = parent.env(environment())
  x <- envir$x

  ## artificial high precision to be added to the mean-model
  prec.high = exp(15)

  interpret.theta = function() {

    lag_cols <- grep("^lag[0-9]+$", names(x), value = TRUE)
    lag_k <- length(lag_cols) - 1L

    gamma1 <- 0.01 * tanh(theta[1])
    gamma2 <- 0.01 * tanh(theta[2])
    beta1 <- theta[3L]

    psi <- numeric(lag_k + 1L)
    for(lag in 0:lag_k){
      psi[lag + 1L] <- exp(gamma1*(lag^1) + gamma2*(lag^2))
    }

    w <- psi / sum(psi)

    out_list <- list()
    for(lag in 0:lag_k){
      out_list[[paste0("w",lag)]] <- w[lag + 1]
    }

    out_list[["beta1"]] <- beta1

    return(out_list)

  }
  graph = function() {
    G = Matrix::Diagonal(n = length(x$lag0), x=1)
    return(G)
  }
  Q = function() {
    Q = prec.high * graph()
    return(Q)
  }
  mu = function() {
    par = interpret.theta()

    lag_cols <- grep("^lag[0-9]+$", names(x), value = TRUE)
    lag_k <- length(lag_cols) - 1L

    compile_lag_label <- c()
    for(lag in 0:lag_k){
      compile_lag_label <- c(compile_lag_label, paste0("lag",lag))
    }
    compile_w_label <- c()
    for(lag in 0:lag_k){
      compile_w_label <- c(compile_w_label, paste0("w",lag))
    }

    agg <- 0
    for(lag in 0:lag_k){
      agg <- agg + par[[compile_w_label[[lag+1]]]] * x[,which(names(x) == compile_lag_label[[lag+1]])]
    }

    return(par$beta1 * agg)
  }
  log.norm.const = function() {
    return(numeric(0))
  }
  log.prior = function() {
    val = (stats::dnorm(theta[1L], mean=0, sd=1, log=TRUE) +
             stats::dnorm(theta[2L], mean=0, sd=1, log=TRUE) +
             stats::dnorm(theta[3L], mean=0, sd=1, log=TRUE))
    return(val)
  }
  initial = function() {
    return(rep(1, 3))
  }
  quit = function() {
    return(invisible())
  }

  val = do.call(match.arg(cmd), args = list())
  return(val)
}




